#! /usr/bin/env perl

require 5.001;
use strict;
use warnings;
use File::Copy;
use File::Path 'make_path';
use File::Temp qw/ tempdir /;
use Getopt::Tabular;
use File::Basename;
use FindBin;
use ASL;

my $Usage = <<USAGE;

This pipeline create CBF quantification maps from row ASL data using a .xml file as input for Neurolens analysis.

Usage $0 [options]

-help for options

USAGE
my $log_dir     =   '/Users/deback/myfiles/postdoc/experiments/PreventAD/analysis/ASL_quantification/logs/ASLquantFromDICOM_v2-094_Elcapitan/Convert_logs';
my ($list,@args);

my @args_table = (["-list",     "string",   1,  \$list,         "list of archives to look in for ASL dicom files"],
                  ["-log_dir",  "string",   1,  \$log_dir,      "directory for log files"]
                 );

Getopt::Tabular::SetHelp ($Usage,'');
GetOptions(\@args_table,\ @ARGV,\@args) || exit 1;

# needed for log file
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
my $date    =   sprintf("%4d-%02d-%02d_%02d:%02d:%02d",$year+1900,$mon+1,$mday,$hour,$min,$sec);
my $log     =   "$log_dir/ConvertDICOMs_$date.log";
open(LOG,">>$log");
print LOG "Log file, $date\n\n";

if(!$list) { print "You need to specify a file with the list of directory to analyze.\n"; exit 1;}

my $template 	= "ExtractDICOM-$hour-$min-XXXXXX"; # for tempdir
my $TmpDir 		= tempdir($template, TMPDIR => 1, CLEANUP => 1 );

open(TARS,"<$list");
my @tars    =   <TARS>;
close(TARS);
foreach my $tar (@tars){
	chomp($tar);
	my $natdir	= dirname($tar);
	my ($site,$candID,$visit)   =   ASL::getSiteSubjectVisitIDs($natdir);
	next    if(!$site);
	print LOG "Site \t\tCandID \t\tVisit\n$site \t$candID \t\t$visit\n";

    # Check if archive is valid
    unless ($tar =~ /tgz$/i){
    	print LOG "\nERROR: $tar is not an archive...\n";
        next;
    }

    # Convert
    my $filename;
    if ($tar =~ /pCASL_(\d+)/) {
		# Determine candidate information based on archive name
		my $run		= $1;
		$filename	= $natdir . "/" . $site . "_" . $candID . "_" . $visit . "_ASL_" . $run . ".nlvolume" ;
		# Extract archive
        `cd $TmpDir; tar -xzf $tar`;
        # Open DICOM and save run as nlvolume
        my $dirname			= basename($tar);
        $dirname			=~ s/_\d\.tgz$//i;
		my $command_open	= "nldo open " . $TmpDir . "/" . $dirname;
		my $command_save	= "nldo save LAST " . $filename;
		my $command_close	= "nldo close ALL";
		system($command_open);
		system($command_save) unless (-e $filename);
		system($command_close);
	} else {
		print LOG "\nERROR: $tar is not a valid archive containing ASL DICOM\n";
        next;
	}

	# Check that the file was successfully converted
	if (-e $filename) {
		print LOG "Successfully create $filename.\n";
	} else {
		print LOG "\nERROR: $filename was not created from $tar.\n";
	}
}

exit 0;
