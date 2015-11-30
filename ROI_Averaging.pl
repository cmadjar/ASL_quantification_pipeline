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
my $avg_dir     =   '/Users/cmadjar/Documents/McGill/PreventAD/Transfers/ASL/From_dcm_files/31LP_candidates/EXTRACTED_CBF';
my $xml_file    =   '/Users/cmadjar/Documents/McGill/PreventAD/Scripts/ASLquantification/ROIaveragingParameters_v2.0_2013-10-28.xml';
my ($cbf_dir, $roi_dir, $roi_basename);
my ($list,@args);

my @args_table = (["-list",         "string", 1,  \$list,         "list of directories to look in for dicom files"            ],
                  ["-log_dir",      "string", 1,  \$log_dir,      "directory for log files"                                   ],
                  ["-roi_dir",      "string", 1,  \$roi_dir,      "directory containing GM masks"                             ],
                  ["-cbf_dir",      "string", 1,  \$cbf_dir,      "directory containing CBF maps"                             ],
                  ["-avg_dir",      "string", 1,  \$avg_dir,      "directory where CSV files with average CBF will be saved"  ],
                  ["-xml",          "string", 1,  \$xml_file,     "xml file summarizing ASL analysis parameters"              ],
                  ["-roi_basename", "string", 1,  \$roi_basename, "Basename of the ROI to use for averaging"                  ]
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
my  $csv_gp =   $avg_dir."/GM_avg_CBF_".$roi_basename."_".$date.".csv";
open(CSV,">$csv_gp") or die "$!";
print CSV "CandID,"  
            . "VisitLabel," 
            . "CBF_map,"  
            . $roi_basename . "_mask,"  
            . $roi_basename . "_Weighted_Average_CBF,"
            . $roi_basename . "_Standard_Deviation_CBF,"
            . $roi_basename . "_Minimum_CBF,"
            . $roi_basename . "_Maximum_CBF,"
            . $roi_basename . "_Voxel_Count,"
            . $roi_basename . "_Volume_Size"
            . "\n";
close(CSV);
    
open(DIRS,"<$list");
my  @dirs   =   <DIRS>;
close(DIRS);

foreach my $dir (@dirs) {
    chomp($dir);

    # determine the site, candidate ID and visit label base on directory path
    my  ($site,$candID,$visit)   =   ASL::getSiteSubjectVisitIDs($dir);
    next    if (!$site);
    print LOG "Site \t\tCandID \t\tVisit\n$site \t$candID \t\t$visit\n";
    print "Site \t\tCandID \t\tVisit\n$site \t$candID \t\t$visit\n";

    # go to the next one if no folder for the candID and visit label exists in the GM and ASL Analyses directories
    my  $cand_roi_dir = $roi_dir."/".$candID."/".$visit;
    print "Cand: $cand_roi_dir\n";
    my  $cand_cbf_dir = $cbf_dir."/".$candID."/".$visit;
    next    if((!$cand_roi_dir) || (!$cand_cbf_dir));

    # get list of GM masks in $cand_roi_dir
    my  ($roi_masks)  = &ASL::getMincs($cand_roi_dir, $roi_basename);
    my  ($cbf_maps)   = &ASL::getMincs($cand_cbf_dir,"cbf.mnc");
    next    if ((@$roi_masks  ==  0) || (@$cbf_maps == 0)); # go to the next candidate if no GM mask in $cand_roi_dir or no CBF maps in $cand_cbf_dir
    
    # resample GM mask to ASL resolution and get the name of resampled GM masks
    my  ($roi_ASLrspld) = &ASL::resampleGMtoASL($roi_masks,$cbf_maps);
    my  ($roi_ASLrspld) = $roi_masks;
   
    # read the XML file with Neurolens' options and get the list of pluggins to run (should only be ROI Averaging)
    my  $xml        = new XML::Simple (KeyAttr=>[]);
    my  $nldo_opt   = $xml->XMLin($xml_file);

    my  ($plugin); 
    my  $options    =   "";
    foreach my $plug ($nldo_opt->{plugin}){
        if  ($plug->{name} eq "ROI Averaging") {
            $plugin  =   $plug->{name};
            my @parameters_list =   @{$plug->{parameter}};
            foreach my $p (@parameters_list) {
                if($p->{name} eq "-maskOperation") { 
                    $options   =   $options." ".$p->{name}." \"".$p->{value}."\"";
                    next;
                }
                $options    =   $options." ".$p->{name}." ".$p->{value};
            }
        }
    }
    &ASL::computeROIs($candID,$visit,$roi_ASLrspld,$cbf_maps,$csv_gp,$plugin,$options);
}    


