#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use diagnostics;
use File::Glob;
use b2b::qc;
use b2b::tools;


my $exp = "52R";
my $dry;
my	 $samplesheet = "/home/laurentt/sampleSheet/BruneauExperimentsGNOMex20130325.txt";
#my $samplesheet = "/home/laurentt/sampleSheet/gnomex_sampleannotations_SEIDMAN_030713.txt";

my $sampleHash = b2b::tools::parseSampleSheet($samplesheet);
my $expSampleHash = b2b::tools::makeExpSampleHash(
	exp => $exp,
	sampleHash => $sampleHash,
	);

my $species = b2b::tools::getSpecies(
	exp => $exp,
	sampleHash => $sampleHash,
	); 

my $analysisDir = "/Data01/gnomex/Analysis/experiment/";


my $expDir= $analysisDir.$exp."/";

my $inpattern = "$expDir/Tophat*/*marked_dup.bam";
print "inpattern\t$inpattern\n";
my @bams = glob $inpattern;

# for my $bam (@bams){
# 	chomp $bam;
# 	print $bam."\n";
# 	b2b::qc::collectDups( 
# 		file=>$bam, 
# 		);
# }

b2b::tools::runDRDS(
	bamID=> 	"dupes",
	sampleHash => $expSampleHash,
	analysisDir => $expDir,
	dry 	=>	$dry, 
	species => $species,
	);
