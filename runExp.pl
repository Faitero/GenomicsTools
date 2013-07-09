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
use b2b::qc;
use diagnostics -verbose;
my $debug = 1;
my $dry = 0;

my	 $samplesheet = "/home/laurentt/sampleSheet/BruneauExperimentsGNOMex20130325.txt";
# my 	 $samplesheet = "/home/laurentt/sampleSheet/gnomex_sampleannotations_SEIDMAN_030713.txt";
my	 $dataPath = "/Data01/gnomex/ExperimentData/";
my 	 $exp; 
my	 $analysisDir = "/Data01/gnomex/Analysis/experiment/";
my 	 $fastqDir;
my 	 $nofastqc;
my 	 $USEQonly;
my	 $TRACKonly;
my   $noTophat;
my $inPattern = "Tophat*/accepted_hits.bam";
my $startingDir = "./";
my $outDir		= "output/";
my $noMarkDup; 						## this flag will skip the dedup step.
my $noClipping;						## this flag will skip the clipping profile
my $markedDupFile;
my $qcTabOnly;
my $mouserefgene = "/work/Common/Data/Annotation/mouse/mm9/Mus_musculus.NCBIM37.67.fixed.bed";
my $humanrefgene = "/work/Common/Data/Annotation/human/Homo_sapiens.GRCh37.71.fixed.bed";
my $refgene;
my $now_string = localtime;
$now_string =~ s/\s/_/g;
my $skipFlag = 0;
my $refgene;

GetOptions(
	"exp=s"			=>	\$exp,
	"nofastqc"		=>  \$nofastqc,
	"dry"			=>	\$dry,
	"USEQonly"		=>	\$USEQonly,
	"TRACKonly"		=>	\$TRACKonly,
	"noTophat"		=>	\$noTophat,
	"inPattern=s"	=>		\$inPattern,
	"startDir=s"	=>		\$startingDir,
	"outDir=s"		=>		\$outDir,
	"noMarkDup"		=>		\$noMarkDup,
	"noClipping"	=>		\$noClipping,
	"noGeneCoverage" =>		\$noGeneCoverage,
	"noJuncAnnot"	=>		\$noJuncAnnot,
	"noJuncSat"		=>		\$noJuncSat,
	"noReadDist"	=>		\$noReadDist,
	"noGC"			=>		\$noGC,
	"noNVC"			=>		\$noNVC,
	"noQuality"		=>		\$noQuality,
	"noFilter"		=>		\$noFilter,
	"noBamStat"		=>		\$noBamStat,
	"noQC"			=>		\$noQC,	
	"noCatTab"		=>		\$noCatTab,
	"species=s"			=>	\$species,
	"qcTabOnly"		=>		\$qcTabOnly,
	);

if (!defined($exp)){
	print "enter experiment id:\n";
	$exp = <stdin>;
}
chomp $exp;

if ( $exp =~ m/^\d+$/){
	$exp .= "R";
} elsif ( $exp =~ m/^\d+R$/){
} else{
	die "Invalid Experiment: must be of the form <experiment #>R";
}

open my $RUNLOG, '>', $analysisDir.$now_string."_".$exp."runlog.txt";
print $RUNLOG "$now_string\nrunlog for experiment $exp\n\n";

print "Experiment $exp\n";
$analysisDir .= $exp."/";
chdir($analysisDir);
system("mkdir -p $analysisDir");
chdir($analysisDir) or die "$!";

print $RUNLOG "$analysisDir created successfully";

my $logfile = $analysisDir."runExp_".$now_string.".log";

open (LOG, ">$logfile");

*STDERR = *LOG;
*STDOUT = *LOG;

# open ERROR,  '>', $logfile  or die $!;
# STDERR->fdopen( \*ERROR,  'w' ) or die $!;

my $fastqcDirCommand = "mkdir -p ${analysisDir}fastqc";
print "$fastqcDirCommand\n";
b2b::tools::runAndLog($fastqcDirCommand);

if(-d "${analysisDir}fastqc"){
	print $RUNLOG "${analysisDir}fastqc created successfully";
} else {
	die "unable to create ${analysisDir}fastqc";
}

unless ( -e $samplesheet ){
	die "Invalid Sample Sheet";
}

my $sampleHash = b2b::tools::parseSampleSheet($samplesheet);
my $expSampleHash = b2b::tools::makeExpSampleHash(
	exp => $exp,
	sampleHash => $sampleHash,
	);
my $species = b2b::tools::getSpecies(
	exp => $exp,
	sampleHash => $sampleHash,
	); 
if( defined($species) ){
	print $RUNLOG "Species:\t$species\n"; 
} else{
	print $RUNLOG "Species could not be determined from the sample sheet\n";
	die "No Specied defined";
}

## check if USEQ only
if ($USEQonly) {
	print "Running USEQ only\n";
	goto USEQ;
}

## check for track only
if ($TRACKonly) {
	print "Making browser tracks only\n";
	goto TRACK;
}

$fastqDir = b2b::tools::findFastQPath (
	exp 		=> $exp,
	path 		=> $dataPath,
	sampleHash 	=> $sampleHash,
	) unless ($nofastqc && $noTophat);

my $fastQCcommand = "fastqc --outdir=${analysisDir}fastqc ${fastqDir}*" unless ($nofastqc);
print "$fastQCcommand\n";
b2b::tools::runAndLog("$fastQCcommand") unless($dry || $nofastqc);

#check that all fastqc dirs were created
b2b::report::checkFastQC(
	RUNLOG =>	$RUNLOG,
	fastqDir => $fastqDir,
	analysisDir => $analysisDir
	);


b2b::tools::runTophat(
	exp 	=> $exp,
	sampleHash => $expSampleHash,
	fastqDir 	=> $fastqDir,
	analysisDir => $analysisDir,
	dry 		=> $dry,
	) unless ($noTophat);

## run runQC_TL.pl on the samples
b2b::report::checkTophat(
	RUNLOG => $RUNLOG,
	sampleHash => $expSampleHash,
	analysisDir => $analysisDir,
	);

print("running QC\n");

# ## this part needs refactoring:
# my $runQCcommand = "runQC_TL.pl --species $species";
# print $runQCcommand."\n";
# b2b::tools::runAndLog($runQCcommand) unless ($dry);

if (lc($species) eq "human"){
	print "species: human\n";
	$refgene = $humanrefgene;
} elsif (lc($species) eq "mouse"){
	print "species: mouse\n";
	$refgene = $mouserefgene;
} else {
	die "species not recognized -- contact bioinformatician\n";
}

print "refgene = $refgene\n";
print $RUNLOG "refgene = $refgene\n";

print "searching for inpattern = \t$inPattern";
my @inputFiles = glob $inPattern;

if (@inputFiles == 0 ){
	die "no Input Files matching the pattern\n";
}


for my $file (@inputFiles){
	print "Processing $file";
	my $outfolder = $file;
	if($outfolder =~ m/(Tophat.+\/)/){
		print("outfolder\t$1\n");
		$outfolder = $1;
	}
	my $markedDupFile = b2b::qc::markDup($file) unless $noMarkDup; 
    $file = $markedDupFile unless !$markedDupFile;
	b2b::report::checkIfExists( file=>$markedDupFile , RUNLOG=>$RUNLOG);
    my $indexFile = $markedDupFile;
    $indexFile =~ s/.bam$/.bai/;
	b2b::report::checkIfExists( file=>$indexFile , RUNLOG=>$RUNLOG);
    my $bamstattab = $file."bam_stat.tab";
    runAndLog("bam_stat.py -i $file 2> $file.bam_stat.tab") unless $noBamStat;   

	b2b::report::checkIfExists( file=>$bamstattab , RUNLOG=>$RUNLOG);
	
	b2b::qc::runQC({
		file => $file,
		refgene => $refgene,
		}) unless $noQC;	
	b2b::report::checkRunQC(
		outfolder => $file."_QC/",
		RUNLOG => $RUNLOG,
		);
}	


TRACK:
## make Genome tracks

my @files = glob("${analysisDir}Tophat*/*processed.bam");
for my $file (@files){
	my $trackGenCommand = "convert_SAM_or_BAM_for_Genome_Browser.pl --nosort $file";
	print "$trackGenCommand\n";
	b2b::tools::runAndLog($trackGenCommand) unless ($dry);
	b2b::report::checkBrowserTracks( file => $file , RUNLOG => $RUNLOG );
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

sub runAndLog{
	my $command = shift;
	my $time = localtime;
	print "$time\t$command";
	system($command);
}