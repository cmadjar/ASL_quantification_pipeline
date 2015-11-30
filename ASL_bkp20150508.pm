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
    my  ($filename,$outdir,$nldo_opt)   =   @_;
    
    my  ($pre_flow_suffix,$pre_even_suffix, $create_snr);
    foreach my $plug (@{$nldo_opt->{plugin}}){

        if ($plug->{name} eq 'Motion Correction'){
            $pre_flow_suffix    =   "-MC";
            $pre_even_suffix    =   "-MC";
            next;
        }

        if ($plug->{name} eq 'ASL Subtraction'){
            $pre_flow_suffix    =   $pre_flow_suffix."-flow";
            $pre_even_suffix    =   $pre_even_suffix."-even";
            next;
        }

        if ($plug->{name} eq 'Spatial Filtering'){
            $pre_flow_suffix    =   $pre_flow_suffix."-SM";
            $pre_even_suffix    =   $pre_even_suffix."-SM";
            next;
        }

        if ($plug->{name} eq 'GLM Time Series'){
            my $options    = &ASL::getParameters($nldo_opt,$plug->{name});
            if ($options =~ /snr 1/i) {
                $create_snr = 1;
            }
        }
        
    }
    
    my $MC_nlvolume         =   $outdir."/".substr(basename($filename),0,-9)."-MC.nlvolume";
    my $MC_minc             =   $outdir."/".substr(basename($filename),0,-9)."-MC.mnc";
    my $preprocessed_flow   =   $outdir."/".substr(basename($filename),0,-9).$pre_flow_suffix.".nlvolume";
    my $preprocessed_even   =   $outdir."/".substr(basename($filename),0,-9).$pre_even_suffix.".nlvolume";
    my $flow_snr            =   substr($preprocessed_flow,0,-9)."-snr.nlvolume"  if ($create_snr);
    my $even_snr            =   substr($preprocessed_even,0,-9)."-snr.nlvolume"  if ($create_snr);
    my $flow_eff            =   substr($preprocessed_flow,0,-9)."-eff.nlvolume";
    my $flow_se_eff         =   substr($preprocessed_flow,0,-9)."-se_eff.nlvolume";
    my $even_eff            =   substr($preprocessed_even,0,-9)."-eff.nlvolume";
    my $even_se_eff         =   substr($preprocessed_even,0,-9)."-se_eff.nlvolume";
    my $cbf_map_nlvolume    =   substr($flow_eff,0,-9)."-cbf.nlvolume";
    my $cbf_map_minc        =   substr($flow_eff,0,-9)."-cbf.mnc";
    print "\n\n$cbf_map\n\n";
    
    return  ($MC_nlvolume, $MC_minc, $preprocessed_flow,$preprocessed_even,$flow_eff,$even_eff,$flow_se_eff,$even_se_eff,$cbf_map_nlvolume, $cbf_map_minc, $flow_snr, $even_snr);
}

=pod

This function reads the xml file filled with analysis' options to use and returns the options in a string.

=cut

sub getParameters{
    my  ($nldo_opt, $plug)  =   @_;
    
    # read the XML file with analyses parameters
#    my $xml = new XML::Simple (KeyAttr=>[]);
#    my $data = $xml->XMLin($xml_file);

    # dereference hash ref
    # access <plugin> array

    my  (%parameters_list,$outputs_list,$options);

    foreach my $plugin (@{$nldo_opt->{plugin}}){
        next    unless ($plugin->{name} eq $plug); 
        
        my @parameters_list =   @{$plugin->{parameter}};
        foreach my $p (@parameters_list){
            
            if($p->{name} eq "-subtractionOrder" ||
               $p->{name} eq "-kernelType"       ||
               $p->{name} eq "-contrastList"     ||
               $p->{name} eq "-interpolationType"||
               $p->{name} eq "-aslType"){
                $options   =   $options." ".$p->{name}." \"".$p->{value}."\"";
                next;
            }
            $options    =   $options." ".$p->{name}." ".$p->{value};
        
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
Compute ROI averaging in GM masks and write it on
=cut
sub computeROIs {
    my  ($candID,$visit,$gm_ASLrspld,$cbf_maps,$csv_gp,$plugin,$options)  =   @_;

    foreach my $cbf_map (@$cbf_maps) {
        foreach my $gm_mask (@$gm_ASLrspld){
            my  ($average)  =   &executeNL($cbf_map,$gm_mask,$plugin,$options);
            &writeCSV($candID,$visit,$cbf_map,$gm_mask,$csv_gp,$average);
        }
    }
}

=pod
Execute Neurolens commands to get the average CBF in GM mask.
=cut
sub executeNL {
    my  ($cbf_map,$gm_mask,$plugin,$options)    =   @_;

    my  $open_cmd       =   "nldo open ".$cbf_map." ".$gm_mask;  
    my  $closeALL_cmd   =   "nldo close ALL";
    my  $ROIavg_cmd     =   "nldo run '$plugin' -maskDataset $gm_mask -targetDataset $cbf_map $options";

    system($open_cmd);
    my  $average        =   `$ROIavg_cmd`;
    chomp($average);
    system($closeALL_cmd);

    return  ($average);
}

=pod
Write in CSV file CBF mean average over GM mask.
=cut
sub writeCSV {
    my  ($candID,$visit,$cbf_map,$gm_mask,$csv_gp,$average)  =   @_;

    my  $cbf    =   basename($cbf_map);
    my  $gm     =   basename($gm_mask);

    open(CSV,">>$csv_gp") or die "$!";
    print CSV "$candID,$visit,$cbf,$gm,$average\n";
    close(CSV)

}

