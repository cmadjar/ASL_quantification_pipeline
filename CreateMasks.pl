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

This pipeline:
    1. creates XFM transforms from the stereotaxic (stx) space to the ASL space using
    	CIVET outputs and motion corrected ASL series.
    2. convert GM maps into the ASL space.
    3. convert specified regions of interest back into the ASL space.

Usage $0 [options]

-help for options

USAGE
my $log_dir   = '/Users/deback/myfiles/postdoc/experiments/PreventAD/analysis/ASL_quantification/'
					. 'logs/ASLquantFromDICOM_v2-094_Elcapitan/Create_Masks_logs';
my $root_dir  = '/Users/deback/myfiles/postdoc/experiments/PreventAD/output/ASL/From_dcm_files/'
     				. 'DATA_RELEASE_1.0_candidates';
my $civet_root_dir = '/Users/deback/myfiles/postdoc/experiments/PreventAD/output/ASL/'
					. 'From_dcm_files/DATA_RELEASE_1.0_candidates/'
					. 'CIVET_1.12_icbm152nl_AAL_N3_125';
my $XFMcreate = 0;
my ($mask_root_dir, $list, $roi_dir, @args);

# description of the table
my $list_desc       = "list of 'candIDs/VisitLabel' (aka: 123456/NAPBL00).";
my $log_dir_desc    = "directory for log files.";
my $root_dir_desc   = "directory where the ANALYSES, TRANSFORMS folders are stored.";
my $mask_dir_desc   = "directory where the masks will be created";
my $civet_dir_desc  = "directory where the CIVET folders are stored.";
my $createXMFM_desc = "if set, create XFM transform from the stx space to the ASL space.";
my $roi_dir_desc    = "directory with a list of ROIs from a template to create specific "
						. "brain roi masks. (e.g. directory to parcel files from DKT "
						. "template (~/Documents/McGill/PreventAD/Scripts/"
						. "ASLquantification/Templates/DKT_mindboggle_101/parcel_files";
my @args_table = (
    [ "-list",      "string",  1, \$list,           $list_desc       ],
    [ "-log_dir",   "string",  1, \$log_dir,        $log_dir_desc    ],
    [ "-root_dir",  "string",  1, \$root_dir,       $root_dir_desc   ],
    [ "-mask_dir",  "string",  1, \$mask_root_dir,  $mask_dir_desc   ],
    [ "-civet_dir", "string",  1, \$civet_root_dir, $civet_dir_desc  ],
    [ "-createXFM", "boolean", 1, \$XFMcreate,      $createXMFM_desc ],
    [ "-roi_dir",   "string",  1, \$roi_dir,        $roi_dir_desc    ]
);

Getopt::Tabular::SetHelp ($Usage,'');
GetOptions(\@args_table,\ @ARGV,\@args) || exit 1;

# needed for log file
my ($log, $date) = &createLogFile($log_dir);
open(LOG,">>$log") or die "Can't write to file '$log' [$!]\n";
my $message = "Log file, $date\n\n";
print LOG $message;

# quit if $list not set or root directory mentionned not found on the file system
if (!$list) {
	$message = "You need to specify a file with the list of candID/VisitLabel to use.\n";
    print $message; print LOG $message;
    exit 1;
}
unless (($mask_root_dir) && (-e $mask_root_dir)) {
	$message = "You need to specify a valid directory where the masks will be created.\n";
    print $message; print LOG $message;
    exit 1;
}
unless (-e $root_dir) {
	$message = "Specified root_dir $root_dir not found. Please verify path";
    print $message; print LOG $message;
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
    next if(!$site);
    print LOG "Site \t\tCandID \t\tVisit\n$site \t$candID \t\t$visit\n";

    # set paths
    $root_dir          =~ s/\/$//;   # remove last '/' from root_dir
    $mask_root_dir     =~ s/\/$//;   # remove last '/' from mask_root_dir
    my ($nls_path, $analyses_dir, $mask_dir, $transform_dir, $civet_dir)
    	= &getPipelinePaths($root_dir, $mask_root_dir, $civet_root_dir, $candID, $visit);

    # grep input files
    my ($MC_series, $cbf_map, $gm_stx, $t1_2_stx, $nat_t1, $civet_brain, $t1_2_asl)
    	= &grepInputFiles($analyses_dir, $civet_dir, $transform_dir);

	# check that all necessary files were found on the filesystem
    next unless ((-e $MC_series) && (-e $cbf_map) && (-e $gm_stx)
    				&& (-e $t1_2_stx) && (-e $nat_t1) && ($civet_brain));

    # check that stx_2_func exists if flag -createXFM is not set
    if ((!$XFMcreate) && (!$t1_2_asl)) {
        $message = "ERROR: XFM file *t1_asl.xfm does not exist. "
        			. "Consider setting -createXFM option.\n";
        print $message;
        print LOG $message;
        next;
    }

    # set new file names
    my ($mask_base, $gm_t1, $gm_asl, $gm_brain);
    ($mask_base, $gm_t1, $gm_asl, $gm_brain, $t1_2_asl)
    		= &setNewFileName($nat_t1, $mask_dir, $transform_dir, $t1_2_asl);

    # create t1_2_asl.xfm transform if createXFM is set
    if ($XFMcreate) {
        my ($t1_2_asl_success) = &createT12aslXFM( $nat_t1,
                                                   $MC_series,
                                                   $t1_2_asl
                                                 );
        next if (!$t1_2_asl_success);
        print LOG "T1 to ASL XFM successfully created.\n";
    }

    # resample GM mask to ASL space
    unless (-e $gm_asl) {
    	my ($resample_success) = &resampleMasksToASLspace($t1_2_stx, $gm_stx, $gm_t1,
    								$civet_brain, $gm_brain, $t1_2_asl, $cbf_map, $gm_asl);
    	# go to next subject unless resample was a success
    	next unless ($resample_success);
    }

	# next unless $roi_dir was set and exists
	next unless ($roi_dir && (-e $roi_dir));

    # Grep all the ROI files from the parcel template directory
    my @roi_list = ();
   	# read ROI dir into an array
   	@roi_list = &readDirectory($roi_dir, 'mnc$');

    # For each ROI found in $roi_dir, do:
    	# 1. resample roi mask into Stereotaxic space to use same resolution between
          # ROI and gm mask before multiplication
    	# 2. resample roi mask into T1 native space
    	# 3. resample roi mask into ASL native space
    # Loop through ROI files
    foreach my $roi (@roi_list) {

    	# define the name of the ROI masks to be created
  	 	my ($roi_stx, $roi_t1, $roi_asl, $roi_gm_stx)
  	 		= &setROInames($roi, $mask_dir, $mask_base);

    	# Multiply the ROI by the GM mask to get only the GM within that ROI
		my $roi_to_stx_options = "-like " . $gm_brain . " -float ";
		my ($roi_stx_success)= &resampleMask( $roi_to_stx_options,
											  $roi,
											  $roi_stx
											);

		# resample the ROI into the ASL space
		my ($resample_success) = &resampleMasksToASLspace($t1_2_stx, $roi_stx,
									$roi_t1, $gm_brain, $roi_gm_stx, $t1_2_asl,
									$cbf_map, $roi_asl);

		# if resample was a success, remove temporary masks that are in the stx space
		if ($resample_success) {
			unlink ($roi_stx);
			unlink ($roi_gm_stx);
		}
    }

}

exit 0;

# Functions

=pod
Will find out which path to use to create the LOG file.
Input:  - $log_dir: directory where the logs should be created
Return: - $log: path to the log file to use and store information about this pipeline
=cut
sub createLogFile {

	# grep the log_dir path
	my ($log_dir) = @_;

	# get the local time information
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	my $date    = sprintf("%4d-%02d-%02d_%02d:%02d:%02d",
						  $year+1900,
						  $mon+1,
						  $mday,
						  $hour,
						  $min,
					  	  $sec);

	# determine the path of the log file
	my $log     = "$log_dir/Create_$date.log";

	# return the log path
	return ($log, $date);
}

=pod
Find and create pipeline paths if they do not exist.
Inputs: - $root_dir:       root directory of the analyses
		- $mask_root_dir:  root directory of the masks to be created
		- $civet_root_dir: root directory of the CIVET outputs
		- $candID: candidate ID
		- $visit:  visit label
Return: - $nls_path:      path to Neurolens root directory
		- $analyses_dir:  path to visit level Neurolens directory
		- $mask_dir:      path to the mask directory of the visit
		- $transform_dir: path to the transform directory of the visit
		- $civet_dir: 	  path to the CIVET directory of the visit
=cut
sub getPipelinePaths {

	# Initialize argument variables
	my ($root_dir, $mask_root_dir, $civet_root_dir, $candID, $visit) = @_;

	# determine CBF analyses paths
    my $nls_path       = $root_dir . '/ANALYSES/ASLquant_v2-094_Elcapitan';
    my $analyses_dir   = $nls_path . '/'            . $candID . '/' . $visit;
    my $mask_dir       = $mask_root_dir . '/'       . $candID . '/' . $visit;
    my $transform_dir  = $root_dir . '/TRANSFORMS/' . $candID . '/' . $visit;

    # determine civet complete directory (aka. DCCID_Visit_adniT1_001)
    my $civet_pattern  = $candID . "_" . $visit . "_adniT1";
    my $civet_dir      = &grepFile($civet_root_dir, $civet_pattern);

    # create mask_dir and transform_dir if do not exist
    make_path($mask_dir,      0, 0755) unless (-e $mask_dir);
    make_path($transform_dir, 0, 0755) unless (-e $transform_dir);

	# Return paths
	return ($nls_path, $analyses_dir, $mask_dir, $transform_dir, $civet_dir);

}

=pod
Function that greps all the input files needed for the creation of masks
Inputs: - $analyses_dir:  directory with the CBF outputs
		- $civet_dir:     directory with the CIVET outputs
		- $transform_dir: directory with the XFM outputs
Return: - $MC_series:   motion corrected ASL series
   		- $cbf_map:     CBF map
   		- $gm_stx:      grey matter mask in stereotaxic (stx) space
   		- $t1_2_stx:    transformation from T1 space to stx space
   		- $nat_t1:      native T1 file
   		- $civet_brain: CIVET mask of the brain
   		- $t1_2_asl:    transformation from the T1 space to the ASL/CBF space
=cut
sub grepInputFiles {

	# Initialize need variables
	my ($analyses_dir, $civet_dir, $transform_dir) = @_;

	# Grep files using greFile function
	my $MC_series   = &grepFile($analyses_dir,  '-MC.mnc$');
    my $cbf_map     = &grepFile($analyses_dir,  '-cbf.mnc$');
    my $gm_stx      = &grepFile($civet_dir,     '_pve_exactgm.mnc$', 'classify');
    my $t1_2_stx    = &grepFile($civet_dir,     '_t1_tal.xfm$',      'transforms/linear');
    my $nat_t1      = &grepFile($civet_dir,     '_t1.mnc$',          'native');
    my $civet_brain = &grepFile($civet_dir,     '_brain_mask.mnc$',  'mask');
    my $t1_2_asl    = &grepFile($transform_dir, '_t1_asl.xfm$');

	# Returns found files
	if ($MC_series && $cbf_map && $gm_stx && $t1_2_stx
	    && $nat_t1 && $civet_brain) {
		return ($MC_series, $cbf_map, $gm_stx, $t1_2_stx,
		 		$nat_t1, $civet_brain, $t1_2_asl);
	} else {
		return undef;
	}

}

=pod
Grep the first file matching pattern $pattern in directory $dir
and return it. Return undef if could not find file
matching pattern $pattern.
TODO: improve this to grep multiple files?
=cut
sub grepFile {

	# Initialize arguments
    my ($dir, $pattern, $civet_subdir) = @_;

    # If $civet_subdir is defined, add it to the civet $dir
    $dir = $dir . '/' . $civet_subdir if ($civet_subdir);

	# Read directory and return it into $matches dereferenced array
	my @matches = &readDirectory($dir, $pattern);

    # Grep only the first file matching pattern $pattern
    my $file    = $matches[0] unless (!@matches);

	# Print in the log if could not find the file
	my $message = "Could not find file "
					. $file
					. " on the filesystem. Check that the path are set correctly.\n";
   	print LOG $message unless (-e $file);

	# returns undef if file could not be found on the file system
    return undef unless (($file) && (-e $file));
    return ($file);
}

=pod
Set the name of the new files to be created.
Inputs: - $nat_t1:        native T1 file
        - $mask_dir:      directory where the masks should be created
        - $transform_dir: directory where the transforms will be created
        - $t1_2_asl:      transform from the T1 to the ASL space
Return: - $mask_base: base name of the masks to be created
		- $gm_t1:     name of the grey matter mask in the T1 space
		- $gm_asl:    name of the grey matter mask in the ASL space
		- $gm_brain:  name of the grey matter brain mask (in stx space?)
		- $t1_2_asl:  transform from the t1 to the ASL space
=cut
sub setNewFileName {

	# grep the arguments
	my ($nat_t1, $mask_dir, $transform_dir, $t1_2_asl) = @_;

	# grep the name of the T1
    my $t1_base   = basename($nat_t1, '.mnc');
    $t1_base      =~ s/_t1$//;

    # determine the base name of the masks that will be created
    my $mask_base = $mask_dir  . '/' . $t1_base;

    # determine the names of the different masks that will be created
    my $gm_t1     = $mask_base . '_pve_exactgm_brain_t1.mnc';
    my $gm_asl    = $mask_base . '_pve_exactgm_brain_asl.mnc';
    my $gm_brain  = $mask_base . '_pve_exactgm_brain.mnc';

    # determine the name of the t1 to asl space XFM file if it is not defined
    $t1_2_asl     = $transform_dir . '/' . $t1_base . '_t1_asl.xfm'  if (!$t1_2_asl);

	# return statement
	return ($mask_base, $gm_t1, $gm_asl, $gm_brain, $t1_2_asl);

}

=pod
Resample a mask from stx space to ASL space.
Inputs: - $t1_2_stx:    transform from t1 to stx space
		- $gm_stx:      mask in stx space (CIVET GM mask)
		- $gm_t1:       mask in t1 space
		- $civet_brain: brain mask from CIVET (or GM mask in stx space if roi respamling)
		- $gm_brain:    mask in the stx space (result of multiplication between grey
							matter and ROI mask or brain_mask)
		- $t1_2_asl:    transform from t1 to ASL space
		- $cbf_map:     CBF map
		- $gm_asl:      mask in ASL space
Return: - True if resampled mask in the ASL space exists in the file system
=cut
sub resampleMasksToASLspace {

	# get the arguments of the function
	my ($t1_2_stx, $gm_stx, $gm_t1, $civet_brain,
		$gm_brain, $t1_2_asl, $cbf_map, $gm_asl) = @_;

    # resample masks into T1 native space
    my $stx_to_t1_options = "-use_input_sampling " .
                                "-float " .
                                "-invert_transformation " .
                                "-transformation " . $t1_2_stx;
    my ($gm_t1_success)   = &resampleMask( $stx_to_t1_options,
                                           $gm_stx,
                                           $gm_t1,
                                           $civet_brain,
                                           $gm_brain
                                         );

    # resample masks into ASL native space
    my $t1_to_asl_options = "-transformation " . $t1_2_asl .
                                " -like "      . $cbf_map  .
                                " -float ";
    my ($gm_asl_success)  = &resampleMask( $t1_to_asl_options,
                                           $gm_t1,
                                           $gm_asl
                                         );

	# if mask in ASL space was created successfully, then return True, otherwise, False
	if (-e $gm_asl) {
		print LOG "Successfully created $gm_asl mask.\n";
		unlink ($gm_t1); # remove the GM mask in the t1 space as it won't be used anymore
		return 1;
	} else {
		print LOG "Could not create $gm_asl mask.\n";
		return undef;
	}

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

=pod
Function to read the content of a directory and grep only files that match the criteria
$match.
Inputs:  - $dir:    directory containing files to be fetch
         - $match:  string used to look for files to be returned
Outputs: - @files_list: list of files found in $dir matching $match
=cut
sub readDirectory {

	# get the directory to read from the function's arguments
	my ($dir, $match) = @_;

	# Initialize the array of the list of files to be returned
	my (@files_list) = ();

	# Read directory $dir and store its content in @entries
	opendir (DIR, "$dir") || die "Cannot open $dir\n";
	my @entries = readdir(DIR);
	closedir(DIR);

	# Keep only files that matches the regex stored in $match
	@files_list = grep (/$match/i, @entries);
	# Add directory path to each element (file) of the array
	@files_list = map {"$dir/" . $_} @files_list;

	# return the list of files
	return (@files_list);

}

=pod
Define the name of the ROI masks to be created.
Inputs: - $roi:       ROI file (with full path)
		- $mask_dir:  path to the masks root directory
		- $mask_base: base name of the mask to be created
Return: - $roi_stx:    ROI in the stereotaxic space
		- $roi_t1:	   ROI in the T1 space
		- $roi_asl:    ROI in the ASL space
		- $roi_gm_stx: ROI in the stereotaxic space multiplied by the grey matter
=cut
sub setROInames {

		# Get the arguments of the function
		my ($roi, $mask_dir, $mask_base) = @_;

    	# Get the base name of the ROI
    	my $roi_base = basename($roi, '.mnc');

    	# Set output filenames based on the ROI name
        my $roi_stx   = $mask_dir  . '/' . $roi_base . '_stx.mnc';
        my $roi_t1    = $mask_base . '_' . $roi_base . '_pve_exactgm_t1.mnc';
        my $roi_asl   = $mask_base . '_' . $roi_base . '_pve_exactgm_asl.mnc';
        my $roi_gm_stx= $mask_base . '_' . $roi_base . '_pve_exactgm_stx.mnc';

		# return names
		return ($roi_stx, $roi_t1, $roi_asl, $roi_gm_stx);

}
