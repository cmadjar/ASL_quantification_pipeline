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

This script computes ROI averaging with Neurolens parameters given in an xml file.

Usage $0 [options]

-help for options

USAGE
# define root directory as well as the directory where the scripts and the data are
# located
my $root_dir    = '/Users/deback/myfiles/postdoc/experiments/PreventAD';
my $scripts_dir = $root_dir . '/analysis/ASL_quantification';
my $data_dir  	= $root_dir . '/output/ASL/From_dcm_files/DATA_RELEASE_1.0_candidates';
# defaults
my $log_dir   = $scripts_dir . '/logs';
my $cbf_dir   = $data_dir    . '/ANALYSES/ASLquantFromDICOM_v2-094_Elcapitan_noFieldmapCorrection';
my $masks_dir = $data_dir    . '/MASK_DKT_mindboggle/partial_volume_correction';
my $avg_dir   = $data_dir    . '/EXTRACTED_CBF';
my $DKT_dir	  = $scripts_dir . '/Templates/DKT_mindboggle_101/parcel_files';
my $xml_file  = $scripts_dir . '/XML_analyses_parameters/ROIaveragingParameters_v2.0_2013-10-28.xml';
my $nlavg	  = undef;  # set compute neurolens average to none
my $mincstats = undef;  # set compute minc stats to none
my ($list, @args);

# define descriptions
my $list_desc    = "list of subject IDs and visit label";
my $log_desc     = "directory for the log files";
my $masks_desc   = "directory containing masks of grey matter and ROIs";
my $cbf_desc     = "directory containing CBF maps";
my $avg_desc     = "directory where CSV files with average CBF will be saved";
my $DKT_desc     = "directory with the DKT parcels in ICBM";
my $xml_desc     = "xml file summarizing ASL analysis parameters";
my $nlavg_desc   = "if set, compute neurolens average";
my $mncstat_desc = "if set, compute mincstats";

# set the args table
my @args_table = (
	["-list",      "string",  1, \$list,      $list_desc  	],
	["-log_dir",   "string",  1, \$log_dir,   $log_desc   	],
	["-masksdir",  "string",  1, \$masks_dir, $masks_desc 	],
	["-cbfdir",    "string",  1, \$cbf_dir,   $cbf_desc   	],
	["-avgdir",    "string",  1, \$avg_dir,   $avg_desc   	],
	["-DKTdir",    "string",  1, \$DKT_dir,   $DKT_desc   	],
	["-xml",	   "string",  1, \$xml_file,  $xml_desc   	],
	["-nlavg",	   "boolean", 1, \$nlavg,     $nlavg_desc	],
	["-mincstats", "boolean", 1, \$mincstats, $mncstat_desc  ]
);

Getopt::Tabular::SetHelp ($Usage,'');
GetOptions(\@args_table,\ @ARGV,\@args) || exit 1;

# needed for log file
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
my $date    =   sprintf("%4d-%02d-%02d",$year+1900,$mon+1,$mday);
my $log     =   "$log_dir/ASLQuantification$date.log";
open(LOG,">>$log");
print LOG "Log file, $date\n\n";

my $message;
if(!$list) {
	$message = "You need to specify a file with the list of directory to compute "
				. "CBF ROI averaging.\n";
	print $message;
	exit 1;
}

# Exits if neither nlavg or mincstats have been set
if (!$nlavg && !$mincstats) {
	$message = "You need to use either -nlavg or mincstats option.\n";
	print $message;
	exit 1;
}

# read the list of subjects into @rows
open(FILE,"<$list");
my  @rows = <FILE>;
close(FILE);

my $count = 0;
# loop through each row (aka subjID/visit)
foreach my $row (@rows) {

    # remove last \n
    chomp($row);

    # determine the site, candidate ID and visit label base on directory path
    my ($site, $candID, $visit) = &ASL::getSiteSubjectVisitIDs($row);
    next if (!$site);
    $message = "Site \t\tCandID \t\tVisit\n$site \t$candID \t\t$visit\n";
    print LOG $message;
    print $message;

    # go to the next one if no folder for the candID and visit label exists in the GM and
    # ASL Analyses directories
    my $cand_masks_dir = $masks_dir . "/" . $candID . "/" . $visit;
    my $cand_cbf_dir   = $cbf_dir   . "/" . $candID . "/" . $visit;
	next if((!$cand_masks_dir) || (!$cand_cbf_dir));

    # get list of CBF maps in $cand_cbf_dir
    my ($cbf_maps) = &ASL::getMincs($cand_cbf_dir,"cbf.mnc");
    next if (@$cbf_maps == 0); # go to the next candidate if no CBF maps in $cand_cbf_dir

    # read XML neurolens file if compute nlavg is set
    my ($plugin, $nloptions);
    if ($nlavg) {
    	$plugin = "ROI Averaging";
		# read neurolens XML file
		my ($nldo_opt) = &ASL::readNeurolensXMLfile($xml_file) if ($nlavg);
		# determine options for ROI Averaging
    	($nloptions)= &ASL::getParameters($nldo_opt, $plugin);
    }

    # determine mincstats options if $mincstats is set
    my $mncoptions = "-stddev -min -max -count -volume";

	# determine list of ROI to average (including overall GM mask)
	my ($roi_list) = &ASL::getROIsList($cand_masks_dir,
										$DKT_dir,
										"pve_exactgm_brain_asl.mnc"
									  );

	# create the CSV file in which avg CBF will be saved
	my $csv_gp = $avg_dir . "/GM_avg_and_DKT_CBF_" . $date . ".csv";

	# if first line, create CSV file and write title
	if ($count < 1) {
		# Determine title row and print it into the CSV file
		my ($title_row) = &ASL::createTitleRow($roi_list, $nlavg, $mincstats);
		open (CSV, ">$csv_gp") or die "$!";
		print CSV $title_row;
		close (CSV);
	}

	# loop through list of CBF maps (each CBF map will be a row of the spreadsheet)
	foreach my $cbf (@$cbf_maps) {
		my ($row) = &ASL::createSpreadsheetRow($roi_list,   $candID,   $visit,	   $cbf,
											   $masks_dir,  $plugin,   $nloptions, $nlavg,
											   $mncoptions, $mincstats
											  );

		# print row into the spreadsheet
		open (CSV, ">>$csv_gp");
		print CSV $row;
		close (CSV);
	}

	$count ++;
}
