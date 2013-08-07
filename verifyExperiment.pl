#!/usr/bin/perl -w 

use strict;
use diagnostics;
use Getopt::Long;
use File::Glob;
use IO::Handle;
use Data::Dumper;
use b2b::tools;


my $exp;
my 	 $samplesheet = "/home/laurentt/sampleSheet/gnomex_sampleannotations_SEIDMAN_030713.txt";
my	 $analysisDir = "/Data01/gnomex/Analysis/experiment/";

GetOptions(

	"exp=s"		=>	\$exp

	);

if(!defined($exp)) {die "You must specify an experiment to check with the --exp parameter";}

if ( $exp =~ m/R$/){
	print "$exp name is correctly formed";
	} elsif ($exp =~ m/\d+/) {
		$exp .= "R";
	} else {
		die "Invalid experiment identifier\t$exp";
	}

my $expDir = $analysisDir.$exp;

if( chdir "${expDir}" == 0 ) {die "${expDir} does not exist";}

my $expSampleHash = b2b::tools::makeExpSampleHash(
		sampleHash => b2b::tools::parseSampleSheet($samplesheet),
		exp        => $exp
	);

open my $REPORT, ">${exp}_RunReport.txt";


checkFastQC(
	fh => $REPORT,
	sampHash => $expSampleHash
	);


## parameter takes the sample hash and the report file handle and
## checks for the appropriate files in the fastqc directory
sub checkFastQC{
	my %args = @_;
	my $sampHash = $args{sampHash};
	my $fh = $args{fh};
	my @fastqFiles;
	for my $sample ( keys( %$sampHash )){
		my @files = split( " " , $sampHash->{$sample}->{"Associated Files"});
		for my $file ( @files ){
			chomp($file);
			print "$file\n";
			push(@fastqFiles, $file);
		}
	}
	print "fastq files\n";
	for my $file (@fastqFiles){
		print $file."\n";
	}

}