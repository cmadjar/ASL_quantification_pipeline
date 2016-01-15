#! /usr/bin/env perl

require 5.001;
use strict;
use warnings;
use File::Copy;
use File::Path 'make_path';
use Getopt::Tabular;
use File::Basename;
use FindBin;
use ASL;

my $Usage = <<USAGE;

This pipeline create CBF quantification maps from row ASL data using a .xml file as input
for Neurolens analysis. This takes DICOM images as input.

Usage $0 [options]

-help for options

USAGE

my $data_dir 	 = '/Users/deback/myfiles/postdoc/experiments/PreventAD/Transfers/ASL/'
					. 'From_dcm_files/DATA_RELEASE_1.0_candidates';
my $script_dir   = '/Users/deback/myfiles/postdoc/experiments/PreventAD/Scripts/ASLquantification';
my $log_dir     =   '/Users/deback/myfiles/postdoc/experiments/PreventAD/analysis/ASL_quantification/logs/ASLquantFromDICOM_v2-094_Elcapitan/ASLquant_logs';
my $out         =   '/Users/deback/myfiles/postdoc/experiments/PreventAD/output/ASL/From_dcm_files/31LP_candidates/ANALYSES/ASLquant_v2-094_Elcapitan';
my $xml_template=   '/Users/deback/myfiles/postdoc/experiments/PreventAD/analysis/ASL_quantification/ASL_quantification_pipeline/XML_analyses_parameters/ASLparameters_v2.0_2013-10-28.xml';
my ($list,@args);
my $unwarp  	 = 0; # set fieldmap correction to none

my @args_table = (
 ["-list",    "string",  1, \$list,         "list of directories with nlvolume files"],
 ["-log_dir", "string",  1, \$log_dir,      "directory for log files"],
 ["-out",     "string",  1, \$out,          "directory where the ANALYSES should be ran"],
 ["-xml",     "string",  1, \$xml_template, "xml file with ASL analysis parameters"],
 ["-unwarp",  "boolean", 1, \$unwarp, 		"whether fieldmap correction should be run"]
);

Getopt::Tabular::SetHelp ($Usage, '');
GetOptions(\@args_table, \ @ARGV, \@args) || exit 1;

# needed for log file
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
my $date = sprintf("%4d-%02d-%02d_%02d:%02d:%02d",
					$year+1900,
					$mon+1,
					$mday,
					$hour,
					$min,
					$sec
				   );
my $log = "$log_dir/ASLQuantification$date.log";
#open (LOG, ">>$log");
print "Log file, $date\n\n";
#close (LOG);

if(!$list) {
	print "You need to specify a file with the list of directory to analyze.\n";
	exit 1;
}

# read list of candID and visit lable
open(DIRS,"<$list");
my @dirs = <DIRS>;
close(DIRS);

foreach my $natdir (@dirs){
	chomp ($natdir);
	my ($site, $candID, $visit) = &ASL::getSiteSubjectVisitIDs($natdir);
	next if(!$site);
	print "Site \t\tCandID \t\tVisit\n$site \t$candID \t\t$visit\n";

    my $outdir   = $out . "/" . $candID . "/" . $visit;
    my $xml_file = $outdir . "/" . basename($xml_template);
    make_path($outdir,0,0755) unless -e $outdir;
    `cp $xml_template $xml_file` unless -e $xml_file;

    # read the XML file with analysis' options and get the list of plugins to run
    my $nldo_opt = &ASL::readNeurolensXMLfile($xml_file);
    my @plugin_list;
    foreach my $plug (@{$nldo_opt->{plugin}}){
        push(@plugin_list, $plug->{name});
    }

    # parse the native directory for ASL files
    opendir(NATDIR, $natdir);
    my @files = readdir(NATDIR);
	closedir(NATDIR);

	# loop through found files in native directory
	my @natfiles = grep(/nlvolume$/, @files);
    foreach my $filename (@natfiles){
        if ($filename!~/_ASL_/i){
            print "$filename is not an ASL file...\n";
            next;
        }

        # Get the names of the output files
        my ($MC_nlvolume,         $MC_minc,
    		$preprocessed_flow,   $preprocessed_even,
    		$flow_eff, 			  $even_eff,
    		$flow_se_eff,         $even_se_eff,
    		$cbf_map_nlvol,       $cbf_map_mnc,
    		$fmap_rads,           $phase_map,
    		$mag_map,             $MC_unwarp,
    		$flow_snr,            $even_snr
    	   ) = &ASL::getOutputNames($filename, $outdir, $nldo_opt, $natdir, $unwarp);

        if ((-e $preprocessed_flow) && (-e $preprocessed_even)
            && (-e $flow_eff)       && (-e $even_eff)
            && (-e $cbf_map_nlvol)  && (-e $cbf_map_mnc)) {
            print "$filename has already been processed\n";
            next;
        }

        if (-e $preprocessed_flow
        	|| -e $preprocessed_even || -e $flow_eff
        	|| -e $even_eff 		 || -e $cbf_map_nlvol
        	|| -e $cbf_map_mnc){
            print "Only part of the outputs are already present for $candID $visit \n"
            		. " Deleting all partial outputs...\n"
            		. "rm $outdir/*";
        }

        # Run preprocessing and GLM on the flow series
        my $command_open     = "nldo open " . $natdir . "/" . $filename;
        my $command_close    = "nldo close LAST";
        my $command_closeALL = "nldo close ALL";
        system($command_open);

        # run motion correction
        foreach my $plug (@plugin_list){
            next unless ($plug eq "Motion Correction");
            print "Running $plug on $candID $visit ...\n";
            my ($options) = &ASL::getParameters($nldo_opt, $plug);
            my $command   = "nldo run '$plug' -inputDataset LAST $options";
            system ($command);
        }

        # save MC file as nlvolume and minc
        my $command_saveMCnlvol = "nldo save LAST " . $MC_nlvolume;
        my $command_saveMCminc  = "nldo save LAST " . $MC_minc;
        system ($command_saveMCnlvol);
        system ($command_saveMCminc);
        system ($command_closeALL);

        # Do fieldmap correction if unwarp is set
        if ($unwarp) {
        	my $success = &ASL::runFieldmapCorrection($fmap_rads, $phase_map,
        	                                          $mag_map,   $MC_minc,
        	                                          $MC_unwarp
        			   					             );
        }

        # open unwarped image if fieldmap correction, or MC file otherwise
        $command_open = "nldo open ";
        if ($unwarp) {
        	$command_open .= $MC_unwarp;
        } else {
        	$command_open .= $MC_nlvolume;
        }
				system($command_open);

        # run CBF analyses
        foreach my $plug (@plugin_list){
            next if ($plug eq "ASL Quantification" || $plug eq "Motion Correction");
            print "Running $plug on $candID $visit ...\n";
            my ($options)   = &ASL::getParameters($nldo_opt,$plug);
            my $command     = "nldo run '$plug' -inputDataset LAST $options";
            system($command);
        }

        # Save and close the preprocessed images and the effect map for the flow series
        my $command_saveSNRFlow     = "nldo save LAST $flow_snr"  if ($flow_snr);
        my $command_saveEffFlow     = "nldo save LAST $flow_eff";
        my $command_saveSeEffFlow   = "nldo save LAST $flow_se_eff";
        my $command_savePreprocFlow = "nldo save LAST $preprocessed_flow";
        system($command_saveSNRFlow) if ($flow_snr);
        system($command_close)       if ($flow_snr);
        system($command_saveEffFlow);
        system($command_close); # close the effet flow map
        system($command_saveSeEffFlow);
        system($command_close); # close the Standard Error flow map
        system($command_savePreprocFlow);
        system($command_close); # close the preprocessed flow map
        system($command_close); # close the flow series

        # run Spatial Filtering and GLM Time Series plugins on the even series
        foreach my $plug (@plugin_list){
            next unless ($plug eq "GLM Time Series" || $plug eq "Spatial Filtering");
            print "Running $plug on $candID $visit even series ...\n";
            my ($options)   =   &ASL::getParameters($nldo_opt,$plug);
            my $command     =   "nldo run '$plug' -inputDataset LAST $options";
            system($command);
        }

        # Save and close the preprocessed images and the effect map for the even series
        my $command_saveSNREven     = "nldo save LAST $even_snr" if ($even_snr);
        my $command_saveEffEven     = "nldo save LAST $even_eff";
        my $command_saveSeEffEven   = "nldo save LAST $even_se_eff";
        my $command_savePreprocEven = "nldo save LAST $preprocessed_even";
        system($command_saveSNREven) if ($even_snr);
        system($command_close)       if ($even_snr);
        system($command_saveEffEven);
        system($command_close);
        system($command_saveSeEffEven);
        system($command_close);
        system($command_savePreprocEven);
        system($command_close);
        system($command_close);
        system($command_closeALL);

        # Open the even and flow maps and run ASL quantification on them
        my $command_open_flow = "nldo open ".$flow_eff;
        my $command_open_even = "nldo open ".$even_eff;
        system($command_open_flow);
        system($command_open_even);
        foreach my $plug (@plugin_list){
            next unless ($plug eq "ASL Quantification");
            print "Running $plug on $candID $visit even series ...\n";
            my ($options) = &ASL::getParameters($nldo_opt,$plug);
            my $cmd = "nldo run '$plug' -m0Estimate $even_eff -deltaM $flow_eff $options";
            system($cmd);
        }
        my $command_saveCBF_nlvol = "nldo save LAST $cbf_map_nlvol";
        my $command_saveCBF_mnc   = "nldo save LAST $cbf_map_mnc";
        system($command_saveCBF_nlvol);
        system($command_saveCBF_mnc);
        system($command_closeALL);
    }
}

exit 0;
