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

This pipeline 
    1. creates XFM transforms from the stereotaxic (stx) space to the ASL space using CIVET outputs and motion corrected ASL series.
    2. convert GM maps into the ASL space.
    3. convert specified regions of interest back into the ASL space.

Usage $0 [options]

-help for options

USAGE
my $log_dir   = '/Users/cmadjar/Documents/McGill/PreventAD/Scripts/ASLquantification/logs/ASLquantFromDICOM_v2-054_Maverick/Create_Masks_logs';
my $root_dir  = '/Users/cmadjar/Documents/McGill/PreventAD/Transfers/ASL/From_dcm_files/DATA_RELEASE_1.0_candidates';
my $civet_root_dir = '/Users/cmadjar/Documents/McGill/PreventAD/Transfers/ASL/From_dcm_files/DATA_RELEASE_1.0_candidates/CIVET_1.12_icbm152nl_AAL_N3_125';
my $XFMcreate = 0; 
my ($mask_root_dir, $list, $roi, @args);

my @args_table = (
    ["-list",      "string", 1, \$list,           "list of 'candIDs/VisitLabel' (aka: 123456/NAPBL00)."],
    ["-log_dir",   "string", 1, \$log_dir,        "directory for log files."],
    ["-root_dir",  "string", 1, \$root_dir,       "directory where the ANALYSES, TRANSFORMS folders are stored."],
    ["-mask_dir",  "string", 1, \$mask_root_dir,  "directory where the masks will be created"],
    ["-civet_dir", "string", 1, \$civet_root_dir, "directory where the CIVET folders are stored."],
    ["-createXFM", "boolean", 1, \$XFMcreate,     "if set, create XFM transform from the stx space to the ASL space."],
    ["-roi",       "string", 1, \$roi,            "roi to be used to create specific brain roi masks"]
);

Getopt::Tabular::SetHelp ($Usage,'');
GetOptions(\@args_table,\ @ARGV,\@args) || exit 1;

# needed for log file
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date    = sprintf("%4d-%02d-%02d_%02d:%02d:%02d",$year+1900,$mon+1,$mday,$hour,$min,$sec);
my $roi_base  = basename($roi, '.mnc');
my $log     = "$log_dir/Create_$roi_base\_$date.log";
open(LOG,">>$log") or die "Can't write to file '$log' [$!]\n";
print LOG "Log file, $date\n\n";

# quit if $list not set or root directory mentionned not found on the file system
if (!$list) { 
    print "You need to specify a file with the list of candID/VisitLabel to use.\n"; 
    exit 1;
}
unless (($mask_root_dir) && (-e $mask_root_dir)) { 
    print "You need to specify a valid directory where the masks will be created.\n"; 
    exit 1;
}
unless (-e $root_dir) {
    print "Specified root_dir $root_dir not found. Please verify path";
    exit 1;
}

# read list of candidates and visit labels
open(SUBJIDs,"<$list");
my @subjIDs    =   <SUBJIDs>;
close(SUBJIDs);

# loop through candidates and visit labels
foreach my $subjID (@subjIDs) {
    chomp ($subjID); # remove return character from the line

    # grep dataset's candID and visit label
    my ($site, $candID, $visit) = &ASL::getSiteSubjectVisitIDs($subjID);
    next    if(!$site);
    print LOG "Site \t\tCandID \t\tVisit\n$site \t$candID \t\t$visit\n";

    # set paths
    $root_dir          =~ s/\/$//;   # remove last '/' from root_dir
    $mask_root_dir     =~ s/\/$//;   # remove last '/' from mask_root_dir
    my $nls_path       = $root_dir . '/ANALYSES/ASLquant_v2-054_Maverick';
    my $analyses_dir   = $nls_path . '/'            . $candID . '/' . $visit;
    my $mask_dir       = $mask_root_dir . '/'       . $candID . '/' . $visit;
    my $transform_dir  = $root_dir . '/TRANSFORMS/' . $candID . '/' . $visit;
    # determine civet complete directory (aka. DCCID_Visit_adniT1_001)
    my $civet_pattern  = $candID . "_" . $visit . "_adniT1";
    my $civet_dir      = &grepFiles($civet_root_dir, $civet_pattern);

    # create mask_dir and transform_dir if do not exist
    make_path($mask_dir,      0, 0755) unless (-e $mask_dir);
    make_path($transform_dir, 0, 0755) unless (-e $transform_dir);

    # grep files
    my $MC_series   = &grepFiles($analyses_dir, '-MC.mnc$');
    my $cbf_map     = &grepFiles($analyses_dir, '-cbf.mnc$');
#    my $gm_stx      = &grepFiles($civet_dir,    '_pve_exactgm.mnc$', 'classify');
    my $t1_2_stx    = &grepFiles($civet_dir,    '_t1_tal.xfm$',      'transforms/linear');
    my $nat_t1      = &grepFiles($civet_dir,    '_t1.mnc$',          'native');
    my $civet_brain = &grepFiles($civet_dir,    '_brain_mask.mnc$',  'mask');
    my $t1_2_asl    = &grepFiles($transform_dir,'_t1_asl.xfm$');

    # check that stx_2_func exists if flag -createXFM is not set
    if ((!$XFMcreate) && (!$t1_2_asl)) {
        my $message = "ERROR: XFM file *t1_asl.xfm does not exist. Consider setting -createXFM option.\n";
        print $message;
        print LOG $message;
        next;
    }

    # check that all necessary files were found on the filesystem
    my $message_intro = "Could not find file ";
    my $message_end   = " on the filesystem. Check that the path are set correctly.\n";
    print LOG $message_intro . "MC series"       . $message_end  unless (-e $MC_series);
    print LOG $message_intro . "CBF map"         . $message_end  unless (-e $cbf_map);
#    print LOG $message_intro . "GM mask in stx"  . $message_end  unless (-e $gm_stx);
    print LOG $message_intro . "t1 2 tal XFM"    . $message_end  unless (-e $t1_2_stx);
    print LOG $message_intro . "native T1"       . $message_end  unless (-e $nat_t1);
    print LOG $message_intro . "CIVET brain mask". $message_end  unless (-e $civet_brain);
#    next unless ((-e $MC_series) && (-e $cbf_map) && (-e $gm_stx) && (-e $t1_2_stx) && (-e $nat_t1) && ($civet_brain));
    next unless ((-e $MC_series) && (-e $cbf_map) && (-e $t1_2_stx) && (-e $nat_t1) && ($civet_brain));
    
    # set new file names
    my $t1_base   = basename($nat_t1, '.mnc');
    $t1_base      =~ s/_t1$//;
    my $mask_base = $mask_dir  . '/' . $t1_base;
#    my $gm_t1     = $mask_base . '_pve_exactgm_brain_t1.mnc';
#    my $gm_asl    = $mask_base . '_pve_exactgm_brain_asl.mnc';
#    my $gm_brain  = $mask_base . '_pve_exactgm_brain.mnc';
    my ($roi_t1, $roi_asl, $roi_gm_stx, $roi_stx);
    if ($roi) {
        $roi_stx   = $mask_dir  . '/' . $roi_base . '_stx.mnc';
#        $roi_t1    = $mask_base . '_' . $roi_base . '_pve_exactgm_t1.mnc';
#        $roi_asl   = $mask_base . '_' . $roi_base . '_pve_exactgm_asl.mnc';
#        $roi_gm_stx= $mask_base . '_' . $roi_base . '_pve_exactgm_stx.mnc';
        $roi_t1    = $mask_base . '_' . $roi_base . '_t1.mnc';
        $roi_asl   = $mask_base . '_' . $roi_base . '_asl.mnc';
        $roi_gm_stx= $mask_base . '_' . $roi_base . '_stx.mnc';
    }
    $t1_2_asl     = $transform_dir . '/' . $t1_base . '_t1_asl.xfm'  if (!$t1_2_asl);

    # create t1_2_asl.xfm transform if createXFM set
    if ($XFMcreate) {
        my ($t1_2_asl_success) = &createT12aslXFM( $nat_t1, 
                                                   $MC_series,
                                                   $t1_2_asl                              
                                                 );
        next if (!$t1_2_asl_success);
        print LOG "T1 to ASL XFM successfully created.\n";
    }
    
    # resample masks into T1 native space
    my $stx_to_t1_options = "-use_input_sampling " .
                                "-float " .
                                "-invert_transformation " .
                                "-transformation " . $t1_2_stx;
#    my ($gm_t1_success)   = &resampleMask( $stx_to_t1_options,
#                                           $gm_stx,
#                                           $gm_t1,
#                                           $civet_brain,
#                                           $gm_brain
#                                         );
    
#    # resample masks into ASL native space
    my $t1_to_asl_options = "-transformation " . $t1_2_asl .
                                " -like "      . $cbf_map  .
                                " -float ";
#    my ($gm_asl_success)  = &resampleMask( $t1_to_asl_options,
#                                           $gm_t1,
#                                           $gm_asl
#                                         );

    # If roi set
    # 1. resample roi mask into Stereotaxic space (to be able to multiply GM and ROI in stx space)
    # 2. resample roi mask into T1 native space
    # 3. resample roi mask into ASL native space
    if ($roi) {
#        my $roi_to_stx_options = "-like " . $gm_brain . " -float ";
#        my ($roi_stx_success)= &resampleMask( $roi_to_stx_options,
#                                              $roi,
#                                              $roi_stx
#                                            );
        my ($roi_t1_success) = &resampleMask( $stx_to_t1_options,
#                                              $roi_stx,
                                              $roi,
#                                              $roi_t1,
                                              $roi_t1
#                                              $gm_brain,
#                                              $roi_gm_stx
                                            );
        my ($roi_asl_success)= &resampleMask( $t1_to_asl_options,
                                              $roi_t1,
                                              $roi_asl
                                            );
    }
}

exit 0;

# Functions

=pod
Grep file matching pattern $pattern in directory $dir
and return it. Return undef if could not find file
matching pattern $pattern.
=cut
sub grepFiles {
    my ($dir, $pattern, $civet_subdir) = @_;

    # if $civet_subdir is defined, add it to the civet $dir
    $dir = $dir . '/' . $civet_subdir if ($civet_subdir);

    # read directory
    opendir (DIR, $dir);
    my @files = readdir (DIR);
    closedir (DIR);

    # grep only the file matching pattern $pattern
    my @matches = grep(/$pattern/, @files);
    my $file    = $dir . "/" . $matches[0] unless (!@matches);

    return undef unless (($file) && (-e $file));
    return ($file);
}    


=pod
Create the t1_2_asl XFM file using mritoself command.
Returns undef if XFM was not created 
or 1 if creation is a success.
=cut
sub createT12aslXFM {
    my ($nat_t1, $MC_series, $t1_2_asl) = @_;

    # run mritoself onto T1 and MC series to create t1_2_asl XFM file
    my $cmd = "mritoself " . $nat_t1 . " " . $MC_series . " " . $t1_2_asl;
    system ($cmd) unless (-e $t1_2_asl);

    return undef unless (-e $t1_2_asl);
    return 1;
}


=pod
Resample mask to T1 space.
Returns undef if mask was not created
or 1 if creation was a success.
=cut
sub resampleMask {
    my ($options, $input_mask, $output_mask, $civet_brain_mask, $gm_brain_mask) = @_;

    # If civet_brain_mask (brain mask produced by CIVET)
    # and gm_brain_mask (gm mask multiplied by CIVET brain mask) name provided
    # create the extracted gm mask and rename input_mask to this gm_brain_mask
    if (($civet_brain_mask) && ($gm_brain_mask)) {
        my $bet_cmd = "mincmath -mult " .
                      $input_mask       . " " .
                      $civet_brain_mask . " " .
                      $gm_brain_mask;
        # run command unless file already exists
        system($bet_cmd) unless (-e $gm_brain_mask);
        # if file was not create, then return undef
        return undef unless (-e $gm_brain_mask);
        # overwrite input_mask to $gm_brain_mask that should be used for mincresample
        $input_mask = $gm_brain_mask;
    }

    # run mincresample on stx mask to put it in the T1 native space
    my $cmd = "mincresample "     .
                $options    . " " . 
                $input_mask . " " . 
                $output_mask;
    system ($cmd) unless (-e $output_mask);   
             
    return undef unless (-e $output_mask);
    return 1;
}
