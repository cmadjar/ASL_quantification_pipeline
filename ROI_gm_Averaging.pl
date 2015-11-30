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
my $log_dir     =   '/Users/cmadjar/Documents/McGill/PreventAD/Scripts/ASLquantification/logs';
my $cbf_dir     =   '/Users/cmadjar/Documents/McGill/PreventAD/Transfers/ASL/From_dcm_files/31LP_candidates/ANALYSES/ASLquant_v2-054_Maverick';
my $gm_dir      =   '/Users/cmadjar/Documents/McGill/PreventAD/Transfers/ASL/From_dcm_files/31LP_candidates/MASKS';
my $avg_dir     =   '/Users/cmadjar/Documents/McGill/PreventAD/Transfers/ASL/From_dcm_files/31LP_candidates/EXTRACTED_CBF';
my $xml_file    =   '/Users/cmadjar/Documents/McGill/PreventAD/Scripts/ASLquantification/ROIaveragingParameters_v2.0_2013-10-28.xml';
my ($list,@args);

my @args_table = (["-list", "string",   1,  \$list,     "list of directories to look in for dicom files"            ],
                  ["-log_dir","string", 1,  \$log_dir,  "directory for log files"                                   ],
                  ["-gmdir","string",   1,  \$gm_dir,   "directory containing GM masks"                             ],
                  ["-cbfdir","string",  1,  \$cbf_dir,  "directory containing CBF maps"                             ],
                  ["-avgdir","string",  1,  \$avg_dir,  "directory where CSV files with average CBF will be saved"  ],
                  ["-xml","string",     1,  \$xml_file, "xml file summarizing ASL analysis parameters"              ]
                 );

Getopt::Tabular::SetHelp ($Usage,'');
GetOptions(\@args_table,\ @ARGV,\@args) || exit 1;

# needed for log file
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
my $date    =   sprintf("%4d-%02d-%02d",$year+1900,$mon+1,$mday);
my $log     =   "$log_dir/ASLQuantification$date.log";
open(LOG,">>$log");
print LOG "Log file, $date\n\n";

if(!$list) { print "You need to specify a file with the list of directory to compute CBF ROI averaging.\n"; exit 1;}

# create the group CSV file in which avg CBF will be saved 
my  $csv_gp =   $avg_dir."/GM_avg_CBF_".$date.".csv";
open(CSV,">$csv_gp") or die "$!";
print CSV "CandID,VisitLabel,CBF_map,GM_mask,AveragedCBF\n";
close(CSV);
    
open(DIRS,"<$list");
my  @dirs   =   <DIRS>;
close(DIRS);

foreach my $dir (@dirs) {
    chomp($dir);

    # determine the site, candidate ID and visit label base on directory path
    my  ($site,$candID,$visit)   =   &ASL::getSiteSubjectVisitIDs($dir);
    next    if (!$site);
    print LOG "Site \t\tCandID \t\tVisit\n$site \t$candID \t\t$visit\n";
    print "Site \t\tCandID \t\tVisit\n$site \t$candID \t\t$visit\n";

    # go to the next one if no folder for the candID and visit label exists in the GM and ASL Analyses directories
    my  $cand_gm_dir    =   $gm_dir."/".$candID."/".$visit;
    print "Cand: $cand_gm_dir\n";
    my  $cand_cbf_dir   =   $cbf_dir."/".$candID."/".$visit;
    next    if((!$cand_gm_dir) || (!$cand_cbf_dir));

    # get list of GM masks in $cand_gm_dir
    my  ($gm_masks)     =   &ASL::getMincs($cand_gm_dir, "pve_exactgm_brain_asl.mnc");
    my  ($cbf_maps)     =   &ASL::getMincs($cand_cbf_dir,"cbf.mnc");
    next    if ((@$gm_masks  ==  0) || (@$cbf_maps == 0)); # go to the next candidate if no GM mask in $cand_gm_dir or no CBF maps in $cand_cbf_dir
    
    # resample GM mask to ASL resolution and get the name of resampled GM masks
#    my  ($gm_ASLrspld)  =   ASL::resampleGMtoASL($gm_masks,$cbf_maps);
    my  ($gm_ASLrspld)     =   $gm_masks;
   
    # read the XML file with Neurolens' options and get the list of pluggins to run (should only be ROI Averaging)
    my  $xml        =   new XML::Simple (KeyAttr=>[]);
    my  $nldo_opt   =   $xml->XMLin($xml_file);

    my  ($plugin,$options);
    foreach my $plug ($nldo_opt->{plugin}){
        if  ($plug->{name} eq "ROI Averaging") {
            $plugin  =   $plug->{name};
            my @parameters_list =   @{$plug->{parameter}};
            foreach my $p (@parameters_list) {
                if($p->{name} eq "-maskOperation") { 
                    $options = $options." ".$p->{name}." \"".$p->{value}."\"";
                    next;
                }
                if ($options) {
                    $options = $options." ".$p->{name}." ".$p->{value};
                } else {
                    $options = $p->{name}." ".$p->{value};
                }
            }
        }
    }
    &ASL::computeROIs($candID,$visit,$gm_ASLrspld,$cbf_maps,$csv_gp,$plugin,$options);
}    


