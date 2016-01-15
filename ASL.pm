=pod

=head1 NAME

ASL -- A set of utility functions to perform ASL quantification.

=head1 SYNOPSIS

 use ASL;

=head1 DESCRIPTION

This is a mismatch of functions that are used to run ASL quantification.
This is called by ASLpreprocessing.pl.

=head1 METHODS

=cut

package ASL;

use Exporter();
use File::Temp qw(tempdir);
use XML::Simple;
use File::Basename;

$VERSION    =   0.0;
@ISA        =   qw(Exporter);

@EXPORT     =   qw();
@EXPORT_OK  =   qw(getSiteSubjectVisitIDs getOutputNames getParameter);


=pod

This function extracts the site, candid and visit from the path to the ASL file given to the script.

=cut

sub getSiteSubjectVisitIDs {
    my  ($d)=   @_;

    if (($d  =~  m/(\d+)\/([N,P][A,R][P,E][B,F][L,U]\d+)/i) || ($d =~ m/(\d+)\/(living_phantom_MTL_MS_\d+)/i)){
        my  $site   ="PreventAD";
        my  $candID =$1;
        my  $visit  =$2;
        return  ($site,$candID,$visit);
    } else {
        return  undef;
    }
}

=pod

This function determines the output names based on which plugin will be run.

=cut

sub getOutputNames {
    my  ($filename, $outdir, $nldo_opt, $natdir, $unwarp) = @_;

    my  ($pre_flow_suffix, $pre_even_suffix, $create_snr);

    foreach my $plug (@{$nldo_opt->{plugin}}) {

        if ($plug->{name} eq 'Motion Correction'){
            $pre_flow_suffix = "-MC";
            $pre_even_suffix = "-MC";
            next unless ($unwarp);
			$pre_flow_suffix .= "-unwarp";
			$pre_even_suffix .= "-unwarp";
			next;
        }



        if ($plug->{name} eq 'ASL Subtraction'){
            $pre_flow_suffix .= "-flow";
            $pre_even_suffix .= "-even";
            next;
        }

        if ($plug->{name} eq 'Spatial Filtering'){
            $pre_flow_suffix .= "-SM";
            $pre_even_suffix .= "-SM";
            next;
        }

        if ($plug->{name} eq 'GLM Time Series'){
            my $options = &ASL::getParameters($nldo_opt,$plug->{name});
            if ($options =~ /snr 1/i) {
                $create_snr = 1;
            }
        }

    }

    # preprocessed file names
    my $file_basename = $outdir . "/" . substr(basename($filename), 0, -9);
    my $MC_nlvolume   = $file_basename . "-MC.nlvolume";
    my $MC_minc       = $file_basename . "-MC.mnc";
    my $preprocessed_flow = $file_basename . $pre_flow_suffix . ".nlvolume";
    my $preprocessed_even = $file_basename . $pre_even_suffix . ".nlvolume";

    # for fieldmap correction
    my ($MC_nii, $MC_unwarp, $phase_map, $mag_map, $fmap_rad) = undef;
    if ($unwarp) {
    	($phase_map, $mag_map) = &findFieldmaps($natdir);
        $MC_unwarp = $file_basename . "-MC_unwarp.nlvolume";
    	$fmap_rads = substr(basename($phase_map), 0, -4) . "_fmap_rad.nii.gz";
    }

    # flow file names
    my $flow_basename = substr($preprocessed_flow, 0, -9);
    my $flow_snr      = $flow_basename . "-snr.nlvolume" if ($create_snr);
    my $flow_eff      = $flow_basename . "-eff.nlvolume";
    my $flow_se_eff   = $flow_basename . "-se_eff.nlvolume";

    # even file names
    my $even_basename    = substr($preprocessed_even, 0, -9);
	my $even_snr         = $even_basename . "-snr.nlvolume" if ($create_snr);
    my $even_eff         = $even_basename . "-eff.nlvolume";
    my $even_se_eff      = $even_basename . "-se_eff.nlvolume";
    my $cbf_map_nlvolume = substr($flow_eff,0,-9) . "-cbf.nlvolume";
    my $cbf_map_minc     = substr($flow_eff,0,-9) . "-cbf.mnc";

    print "\n\n$cbf_map\n\n";

    return ($MC_nlvolume,         $MC_minc,
    		$preprocessed_flow,   $preprocessed_even,
    		$flow_eff, 			  $even_eff,
    		$flow_se_eff,         $even_se_eff,
    		$cbf_map_nlvolume,    $cbf_map_minc,
    		$fmap_rads,           $phase_map,
    		$mag_map,             $MC_unwarp,
    		$flow_snr,            $even_snr
    	   );
}

=pod
Grep the fieldmaps from the native directory
=cut
sub findFieldmaps {
	my ($natdir) = @_;

	# read native directory
	opendir (NATDIR, $natdir);
	my @files = readdir (NATDIR);
	closedir (NATDIR);

	# grep first phase
	my @phase_files = grep (/Fieldmap_001\.mnc|Fieldmap_019\.mnc/, @files);
	my @mag_files   = grep (/Fieldmap_002\.mnc|Fieldmap_020\.mnc/, @files);

	# grep first phase
	my $phase_map = $natdir . "/" . $phase_files[0];
	my $mag_map   = $natdir . "/" . $mag_files[0];

	return ($phase_map, $mag_map);
}

=pod
=cut
sub convert2nifti {
	my ($minc_file) = @_;

	my $nii_file = substr($minc_file, 0, -4) . "-nii.nii";

	my $mnc2nii_cmd = "mnc2nii $minc_file $nii_file";
	system ($mnc2nii_cmd) unless (-e $nii_file);

	return $nii_file;
}

=pod
Run fieldmap correction.
=cut
sub runFieldmapCorrection {
	my ($fmap_rads, $phase_map, $mag_map, $MC_minc, $MC_unwarp) = @_;

	# return success if fmap_rad exists already
	return 1 if (-e $MC_unwarp);

	# convert to fieldmaps and MC series to nifti format
	my ($MC_nii)    = &convert2nifti($MC_minc);
	my ($phase_nii) = &convert2nifti($phase_map);
	my ($mag_nii)   = &convert2nifti($mag_map);

	# define out directory
	my $outdir = dirname($MC_nii);

	# Brain extraction of the magnitude file
	my $mag_name  = substr(basename($mag_nii), 0, -8);
	my $mag_brain = $outdir . "/" . $mag_name . "_brain.nii.gz";
	my $bet_cmd   = "bet $mag_nii $mag_brain -B -f 0.5 -g 0";
	system ($bet_cmd) unless (-e $mag_brain);

	# transform the phase map in radians
	my $phase_name  = substr(basename($phase_nii), 0, -8);
	my $phase_rad   = $outdir . "/" . $phase_name . "_phase_rad.nii.gz";
	my $convert_rad = "fslmaths $phase_nii -div 8190 -mul 3.14159 $phase_rad -odt float";
	system ($convert_rad) unless (-e $phase_rad);

	# unwarp the phase map
	my $uphase_rad =  $outdir . "/" . $phase_name . "_uphase_rad.nii.gz";
	my $unwarp_cmd = "prelude -p $phase_rad -a $mag_brain -o $uphase_rad";
	system ($unwarp_cmd) unless (-e $uphase_rad);

	# transform the fieldmap in rad/s
	my $transf_cmd = "fslmaths $uphase_rad -div 0.00246 $fmap_rads";
	system ($transf_cmd) unless (-e $fmap_rads);

	# unwarp MC ASL file
	my $MC_unwarp_name  = substr($MC_unwarp, 0, -9);
	my $MC_unwarp_nii   = $MC_unwarp_name . ".nii";
  my $MC_unwarp_nlvol = $MC_unwarp_name . ".nlvolume";

	my $unwarp_cmd = "fugue -i $MC_nii --loadfmap=$fmap_rads --dwell=0.000209999328 "
						. "--asym=0.00246 --unwarpdir=y- -s 0.5 -u $MC_unwarp_nii";
	system ($unwarp_cmd);

	# gunzip the $MC_unwarp output
	my $gunzip_cmd = "gunzip $MC_unwarp_nii.gz";
	system ($gunzip_cmd);

	# convert MC_unwarp file back to minc and reinsert minc header information
	my ($MC_unwarp_mnc)	= &convert2minc($MC_unwarp_nii, $MC_minc);

	# open with neurolens, save it as nlvolume
  my ($open_mnc) = "nldo open $MC_unwarp_mnc";
  my ($nldo_save) = "nldo save LAST $MC_unwarp_nlvol";
  my ($nldo_close) = "nldo close ALL";
  system($open_mnc);
  system($nldo_save);
  system($nldo_close);

	# return undef if $fmap_rad does not exist, 1 otherwise
	return undef unless (-e $MC_unwarp);

	# remove temporary files
	my $rm_cmd = "rm $mag_brain $phase_rad $uphase_rad $fmap_rads $MC_nii *_brain_mask.nii.gz";
	system ($rm_cmd);

	return 1;

}


=pod

=cut
sub convert2minc {
	my ($MC_unwarp_nii, $MC_mnc) = @_;

	# determine minc name
	my $MC_unwarp_name = substr($MC_unwarp_nii, 0, -4);
	my $MC_unwarp_mnc  = $MC_unwarp_name . ".mnc";

	# nii2mnc conversion
	my $mnc2nii_cmd = "nii2mnc $MC_unwarp_nii $MC_unwarp_mnc";
	system($mnc2nii_cmd);

	# reinsert mincheaders in the file
	my ($success) = &insertMincHeader($MC_mnc, $MC_unwarp_mnc);

	return ($MC_mnc) if ($success && (-e $MC_mnc));
	return undef;
}


=pod
Read neurolens XML file
=cut
sub readNeurolensXMLfile {

	my ($xml_file) = @_;

	# read the XML file with Neurolens' options and get the list of plugins to run
    # (should only be ROI Averaging)
    my $xml      = new XML::Simple (KeyAttr=>[]);
    my $nldo_opt = $xml->XMLin($xml_file);

    return ($nldo_opt);

}


=pod

This function reads the xml file filled with analysis' options to use and returns the options in a string.

=cut

sub getParameters{
    my  ($nldo_opt, $plug) = @_;

    my  (%parameters_list, $outputs_list, $options);

    foreach my $plugin (@{$nldo_opt->{plugin}}){
        next unless ($plugin->{name} eq $plug);

        my @parameters_list = @{$plugin->{parameter}};
        foreach my $p (@parameters_list){

            if($p->{name} eq "-subtractionOrder" ||
               $p->{name} eq "-kernelType"       ||
               $p->{name} eq "-contrastList"     ||
               $p->{name} eq "-interpolationType"||
               $p->{name} eq "-aslType"			 ||
               $p->{name} eq "-maskOperation"
              ) {
                $options = $options . " " . $p->{name} . " \"" . $p->{value} . "\"";
                next;
            }
            $options = $options . " " . $p->{name} . " " . $p->{value};

        }
    }

    return  $options;
}

=pod
Return the list of GM masks in GM_dir/candID/visitlabel.
=cut

sub getMincs {
    my ($dir,$pattern)    =   @_;

    opendir(DIR,$dir);
    my  @files   =  readdir(DIR);
    closedir(DIR);

    my  (@keep)   =   grep( /$pattern/i,   @files);
    my  @mincs;
    foreach my $k (@keep) {
        push(@mincs,"$dir/$k");
    }

    return  (\ @mincs);
}

=pod
Resample each GM mask to CBF resolution.
=cut

sub resampleGMtoASL {
    my  ($gm_masks,$cbf_maps)   =   @_;

    my  ($gm_rspled);
    foreach my $gm  (@$gm_masks) {
        ($gm_rspled)    =   &resample($gm, $cbf_maps);   # determine the array of resample GM name
    }

    return  ($gm_rspled);
}

=pod
Determine name of resampled GM masks based on ASL and GM native mask's names and execute mincresample.
=cut

sub resample {
    my  ($gm, $cbf_maps)    =   @_;

    my  @gm_rspled;
    foreach my $cbf (@$cbf_maps) {
        my  $respled_mask   =   $gm;
        if  (($cbf =~ m/_(ASL)_(\d\d\d)/i) && ($gm !~ m/_ASL\d\d\d/i)) {
            my  $ext        =   "_$1$2res.nlvolume";
            $respled_mask   =~  s/.nlvolume$/$ext/g;
        }

        my $command =   "mincresample -like $cbf $gm $respled_mask -clobber";
        system($command)    unless (-e $respled_mask);
        push(@gm_rspled,$respled_mask);
    }

    return  (\ @gm_rspled);
}


=pod
Determine list of ROIs to use for ROI averaging (including overall GM mask)
=cut
sub getROIsList {
	my ($cand_masks_dir, $DKT_dir, $gm_match) = @_;

	my ($gm_masks) = &getMincs($cand_masks_dir, $gm_match);
	my @roi_list = $$gm_masks[0];
	opendir (DIR, $DKT_dir) or die "$!";
	my @files = readdir(DIR);
	closedir(DIR);
	foreach my $file (@files) {
		next if ($file eq ".") || ($file eq "..") || ($file eq ".DS_Store");
		$file =~ s/\.mnc//i;
		push (@roi_list, $file);
	}

	return (\ @roi_list);

}


=pod
Determine title row of the spreadsheet to create with ROI averages
=cut
sub createTitleRow {
	my ($roi_list, $nlavg, $mincstats) = @_;

	my $title_row = "CandID, Visit, CBF_map";

	# Loop through ROI list to determine column names
	foreach my $roi (@$roi_list) {
		my $roi_name = $roi;
		my $dkt_basename = "corrected_rwOASIS-TRT-"
							. "20_DKT31_CMA_jointfusion_labels_in_MNI152_onlyGM_";

		# determine ROI name to use for spreadsheet
		$roi_name = "GM" if ($roi =~ m/pve_exactgm_brain_asl.mnc/i);
		$roi_name =~ s/$dkt_basename//i if ($roi =~ m/$dkt_basename/i);
		$roi_name =~ s/_pve_exactgm_asl//i;

		# Push name of the ROI in the title row with average if nlavg is set
		$title_row = $title_row . "," . $roi_name . "_average" if $nlavg;

		# Push name of the ROI in the title row with mincstats fields if mincstats is set
		if ($mincstats) {
			$title_row = $title_row  . ","
						 . $roi_name . "_stddev,"
						 . $roi_name . "_min,"
						 . $roi_name . "_max,"
						 . $roi_name . "_number_of_voxels,"
						 . $roi_name . "_volume_mm3";
		}
	}

	$title_row .= "\n";

	return $title_row;

}


=pod
Create spreadsheet row for the CBF map
=cut
sub createSpreadsheetRow {
	my ($roi_list,   $candID,   $visit,		$cbf,
		$masks_dir,  $plugin,   $nloptions, $nlavg,
		$mncoptions, $mincstats
	   ) = @_;

	my $row = $candID . ", " . $visit . ", " . basename($cbf);

	foreach my $roi (@$roi_list) {

		unless ($roi =~ m/pve_exactgm_brain_asl.mnc/i) {
			my $subject_dir = $masks_dir . "/" . $candID . "/" . $visit;
			($roi) = &getMincs($subject_dir, $roi);
			$roi = $$roi[0];
		}
		my ($roi_values) = &computeROIs($roi,     $cbf,    	   $plugin,
										$nloptions, $mncoptions, $nlavg,
										$mincstats
									   );

		# if $nlavg is set, add average CBF to the row of the spreadsheet
		$row .= "," . $$roi_values{'Average'} if ($nlavg);

		# if #mincstats is set, add minc stats of the ROI to the row of the spreadsheet
		if ($mincstats) {
			$row = $row . ","
					. $$roi_values{'Stddev'} . ","
					. $$roi_values{'Min'} . ","
					. $$roi_values{'Max'} . ","
					. $$roi_values{'# voxels'} . ","
					. $$roi_values{'Volume (mm3)'};
		}

	}

	$row .= "\n";

	return $row;

}

=pod
Compute ROI averaging in GM masks and write it on
=cut
sub computeROIs {
    my ($roi, $cbf_map, $plugin, $nl_options, $mnc_options, $nlavg, $mincstats) = @_;

    my ($values) = &executeMincstats($cbf_map, $roi, $mnc_options) if ($mincstats);

    if ($nlavg) {
    	my ($average) = &executeNL($cbf_map, $roi, $plugin, $nl_options);
    	# add average to the hash of values
    	$values->{'Average'} = $average;
    }

	return ($values);
}

=pod
Execute Neurolens command to get the average CBF in GM mask.
=cut
sub executeNL {
    my  ($cbf_map,$gm_mask,$plugin,$options)    =   @_;

    my  $open_cmd     = "nldo open ".$cbf_map." ".$gm_mask;
    my  $closeALL_cmd = "nldo close ALL";
    my  $ROIavg_cmd   = "nldo run '$plugin' -maskDataset $gm_mask -inputDataset $cbf_map $options";

    system($open_cmd);
    my  $average      = `$ROIavg_cmd`;
    chomp($average);
    system($closeALL_cmd);

    return  ($average);
}

=pod
Execute mincstats command to get standard deviation, min, max ...
=cut
sub executeMincstats {
    my ($cbf_map, $gm_mask, $options) = @_;

    my $mincstats_cmd = "mincstats "
                          . $cbf_map
                          . " -mask "
                          . $gm_mask
                          . " -mask_range 0.5,1.5 "
                          . $options;
    my $results = `$mincstats_cmd`;
    my @vals    = split("\n", $results);
    my %values;
    foreach my $val (@vals) {
        my ($key, $value) = split(':', $val);
        $value =~ s/ //g;  # remove spaces from value
        $values{$key} = $value;
    }

    return (\ %values);
}


=pod
Write in CSV file CBF mean average over GM mask.
=cut
sub writeCSV {
    my ($candID,$visit,$cbf_map,$gm_mask,$csv_gp,$average,$stddev,$min,$max,$count,$volume) = @_;

    my $cbf = basename($cbf_map);
    my $gm  = basename($gm_mask);

    open(CSV,">>$csv_gp") or die "$!";
    print CSV $candID   . ","
            . $visit    . ","
            . $cbf      . ","
            . $gm       . ","
            . $average  . ","
            . $stddev   . ","
            . $min      . ","
            . $max      . ","
            . $count    . ","
            . $volume   . ","
            . "\n";
    close(CSV)

}


=pod
Insert in the minc header all the acquisition arguments.
Takes the raw MC ASL minc file as input and modify
the unwarped ASL minc file based on the raw minc file's argument.
Inputs:  - $raw_file: raw DTI minc file to grep header information
         - $processed_minc: processed minc file in which header info will be inserted
Outputs: - 1 if if all mincheader information was inserted
         - undef otherwise
=cut
sub insertMincHeader {
    my  ($raw_file, $processed_minc) = @_;

    # insert old acquisition, patient and study arguments
    my  ($acqInsert) = &insertFieldList($raw_file, $processed_minc, 'acquisition:');

    # insert patient information from the raw dataset into the processed files
    my ($patientInsert) = &insertFieldList($raw_file, $processed_minc, 'patient:');

    # insert study information from the raw dataset into the processed files
    my ($studyInsert) = &insertFieldList($raw_file, $processed_minc, 'study:');

    if (($acqInsert) && ($patientInsert) && ($studyInsert)) {
        return 1;
    } else {
        return undef;
    }
}


=pod
Insert information extracted from MC ASL dataset and insert it in the unwarped file.
If one of the value to insert is not defined, return undef, otherwise return 1.
Inputs:  - $raw_mnc: raw MC ASL minc file to grep header information from
         - $processed_minc: processed file in which header information will be inserted
         - $minc_field: minc field to be inserted in processed minc file
Outputs: - 1 if information were inserted in unwarped ASL minc
         - undef otherwise
=cut
sub insertFieldList {
    my ($raw_mnc, $processed_minc, $minc_field) = @_;

    # fetches list of arguments starting with $minc_field (i.e. 'patient:'; 'study:' ...)
    my ($arguments) = &fetch_header_info($minc_field, $raw_mnc, '$1, $2');

    # fetches list of values with arguments starting with $minc_field. Don't remove
    # semi_colon (last option of fetch_header_info).
    my ($values) = &fetch_header_info($minc_field,
    								  $raw_mnc,
    								  '$3, $4, $5, $6, $7',
    								  1
    								 );

    my ($arguments_list, $arguments_list_size) = &get_header_list('=', $arguments);
    my ($values_list, $values_list_size)       = &get_header_list(';', $values);

    my @insert_failure;
    if ($arguments_list_size == $values_list_size)  {
        for (my $i=0; $i<$arguments_list_size; $i++)   {
            my  $argument = @$arguments_list[$i];
            my  $value    = @$values_list[$i];
            my ($insert)  = &modify_header($argument,
            							   $value,
            							   $processed_minc,
            							   '$3, $4, $5, $6'
            							  );
            # store in array @insert_failure the arguments that were not successfully
            # inserted in the mincheader
            push (@insert_failure, $argument) if (!$insert);
        }
        # if at least one insertion failed, will return undef, otherwise 1.
        if ($#insert_failure >= 0) {
            return  undef;
        } else {
            return 1;
        }
    # if arguments_list and values_list do not have the same size, will return undef
    } else {
        return undef;
    }
}

=pod
Function that runs minc_modify_header and insert
minc header information if not already inserted.
Inputs:  - $argument: argument to be inserted in minc header
         - $value: value of the argument to be inserted in minc header
         - $minc: minc file
         - $awk: awk information to check if argument not already inserted in minc header
Outputs: - 1 if argument was indeed inserted into the minc file
         - undef otherwise
=cut
sub modify_header {
    my  ($argument, $value, $minc, $awk) = @_;

    # check if header information not already in minc file
    my $hdr_val = &fetch_header_info($argument, $minc, $awk);

    # insert mincheader unless mincheader field already inserted ($hdr_val eq $value)
    my  $cmd = "minc_modify_header -sinsert $argument=$value $minc";
    system($cmd) unless (($hdr_val) && ($value eq $hdr_val));

    # check if header information was indeed inserted in minc file
    my $hdr_val2 = &fetch_header_info($argument, $minc, $awk);
    if ($hdr_val2) {
        return 1;
    } else {
        return undef;
    }
}

=pod
Function that fetch header information in minc file
Inputs:  - $field: field to look for in minc header
         - $minc: minc file
         - $awk: awk information to check if argument not already inserted in minc header
         - $keep_semicolon: if defined, keep semicolon at the end of the value extracted
Outputs: - $value: value of the field found in the minc header
=cut
sub fetch_header_info {
    my ($field, $minc, $awk, $keep_semicolon) = @_;

    my $val   = `mincheader $minc | grep $field | awk '{print $awk}' | tr '\n' ' '`;
    my $value = $val if $val !~ /^\s*"*\s*"*\s*$/;
    if ($value) {
        $value =~ s/^\s+//;                       # remove leading spaces
        $value =~ s/\s+$//;                       # remove trailing spaces
        $value =~ s/;// unless ($keep_semicolon); # remove ";" unless $keep_semicolon is defined
    } else {
        return undef;
    }

    return ($value);
}

=pod
Get the list of arguments and values to insert into the mincheader (acquisition:*,
patient:* and study:*).
Inputs:  - $splitter: splitter used to separate list of fields stored in $fields
         - $fields: list header arguments and values to insert in the minc header
Outputs: - $list: array of header arguments and values' list
         - $list_size: size of the array $list
=cut
sub get_header_list {
    my ($splitter, $fields) = @_;

    my @tmp = split ($splitter, $fields);
    pop (@tmp);
    my @items;
    foreach my $item (@tmp) {
        $item =~ s/^\s+//;   # remove leading spaces
        $item =~ s/\s+$//;   # remove trailing spaces
        push (@items, $item);
    }
    my  $list      = \@items;
    my  $list_size = @$list;

    return ($list, $list_size);
}
