#!/usr/bin/perl -w
# use b2b::tools;
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename;
use File::Glob;
use IO::Handle;
use b2b::tools;
use diagnostics -verbose;
my $debug = 1;
my $dry = 0;

# my	 $samplesheet = "/home/laurentt/sampleSheet/BruneauExperimentsGNOMex20130325.txt";
my 	 $samplesheet = "/home/laurentt/sampleSheet/gnomex_sampleannotations_SEIDMAN_030713.txt";
my	 $dataPath = "/Data01/gnomex/ExperimentData/";
my 	 $exp; 
my	 $analysisDir = "/Data01/gnomex/Analysis/experiment/";
my 	 $fastqDir;
my 	 $nofastqc;
my 	 $USEQonly;
my	 $TRACKonly;
my   $noTophat;

GetOptions(
	"exp=s"			=>	\$exp,
	"nofastqc"		=>  \$nofastqc,
	"dry"			=>	\$dry,
	"USEQonly"		=>	\$USEQonly,
	"TRACKonly"		=>	\$TRACKonly,
	"noTophat"		=>	\$noTophat,
	);

if (!defined($exp)){
	print "enter experiment id:\n";
	$exp = <stdin>;
}
chomp $exp;




$analysisDir .= $exp."/";

system("mkdir -p $analysisDir");
chdir($analysisDir) or die "$!";

my $now_string = localtime;
$now_string =~ s/\s/_/g;
my $logfile = $analysisDir."runExp_".$now_string.".log";
open ERROR,  '>', $logfile  or die $!;
STDERR->fdopen( \*ERROR,  'w' ) or die $!;


my $fastqcDirCommand = "mkdir -p ${analysisDir}fastqc";
print "$fastqcDirCommand";
b2b::tools::runAndLog($fastqcDirCommand);

my $sampleHash = b2b::tools::parseSampleSheet($samplesheet);
my $expSampleHash = b2b::tools::makeExpSampleHash(
	exp => $exp,
	sampleHash => $sampleHash,
	);

my $species = b2b::tools::getSpecies(
	exp => $exp,
	sampleHash => $sampleHash,
	); 

## check if USEQ only
if ($USEQonly) {
	print "Running USEQ only\n";
	goto USEQ;
}

## chech for track only
if ($TRACKonly) {
	print "Making browser tracks only\n";
	goto TRACK;
}

$fastqDir = b2b::tools::findFastQPath (
	exp 		=> $exp,
	path 		=> $dataPath,
	sampleHash 	=> $sampleHash,
	) unless ($nofastqc && $noTophat);

## run FastQC on the fastqfiles
my $fastQCcommand = "fastqc --outdir=${analysisDir}fastqc ${fastqDir}*" unless ($nofastqc);
print "$fastQCcommand\n";
b2b::tools::runAndLog("$fastQCcommand") unless($dry || $nofastqc);


b2b::tools::runTophat(
	exp 	=> $exp,
	sampleHash => $expSampleHash,
	fastqDir 	=> $fastqDir,
	analysisDir => $analysisDir,
	dry 		=> $dry,
	) unless ($noTophat);

## run runQC_TL.pl on the samples

chdir($analysisDir);

print("running QC\n");
my $sample1 = $exp;
$sample1 =~ s/R//;
$sample1 .= "X1";


my $runQCcommand = "runQC_TL.pl --species $species";
print $runQCcommand."\n";
b2b::tools::runAndLog($runQCcommand) unless ($dry);

TRACK:
## make Genome tracks

my @files = glob("${analysisDir}Tophat*/*processed.bam");
for my $file (@files){
	my $trackGenCommand = "convert_SAM_or_BAM_for_Genome_Browser.pl --nosort $file";
	print "$trackGenCommand\n";
	b2b::tools::runAndLog($trackGenCommand) unless ($dry);
}

print "Finished making browser tracks\n";

## chech for track only
if ($TRACKonly) {
	
	goto FINISH;
}


USEQ:
print "Starting USEQ DefinedRegionDifferentialSeq\n";
## differential expression 
b2b::tools::runDRDS(
	sampleHash => $expSampleHash,
	analysisDir => $analysisDir,
	dry 	=>	$dry, 
	species => $species,
	);
if ($USEQonly){
	goto FINISH;
}



FINISH:
print "analysis complete\n";