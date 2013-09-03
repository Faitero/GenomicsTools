#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename;
use File::Glob;
use IO::Handle;
use tools::fq;
use tools::bam;
use tools::meta;
use tools::db;
use tools::report;
# use diagnostics -verbose;
our $debug = 1;
our $dry = 0;

my $log = 0;
# my	 $samplesheet = "/home/laurentt/sampleSheet/BruneauExperimentsGNOMex20130325.txt";
# my 	 $samplesheet = "/home/laurentt/sampleSheet/gnomex_sampleannotations_SEIDMAN_030713.txt";
# my 	 $samplesheet = "/home/laurentt/sampleSheet/Alisha144R_2013-08-14.txt";
my 	 $samplesheet= "/home/laurentt/sampleSheet/166R_2013-08-20.txt";
my 	 $alignGTF;						
my	 $dataPath = "/Data01/gnomex/ExperimentData/";
my 	 $exp; 
my	 $analysisDir = "/Data01/gnomex/Analysis/experiment/";
my 	 $fastqDir;
my 	 $nofastqc;
my 	 $fastqconly;
my 	 $useqonly;
my	 $trackonly;
my   $notophat;
my $inPattern = "Tophat*/accepted_hits.bam";
my $startingDir = "./";
my $outDir		= "output/";
my $noMarkDup; 						## this flag will skip the dedup step.
my $noClipping;						## this flag will skip the clipping profile
my $markedDupFile;
my $qcTabOnly;
my $noGeneCoverage;
my $noInnerDist;  					## this prevents the inner_distance.py from being executed
my $noJuncAnnot;					## this prevents the junciton_annotation.py from being executed				
my $noJuncSat;						## prevent junction_saturation.py from being executed
my $noReadDist;						## prevents read_distribution.py from being executed
my $noReadDup;						## prevents read_duplication.py from being executed
my $noGC;							## prevents read_GC.py from being executed
my $noNVC;							## prevent read_NVC.py from being called
my $noQuality;						## prevent read_quality from being called
my $noRPKMSat;						## prevent RPKM_saturation.py from being run
my $noFilter;						## prevent script from deDuping
my $noBamStat;						## prenent bam_stat.py from being run
my $noqc;							## prevent the runQC sub from being run
my $tophatonly;			
my $species;
my $mouserefgene = "/work/Common/Data/Annotation/mouse/mm9/Mus_musculus.NCBIM37.67.fixed.bed";
my $humanrefgene = "/work/Common/Data/Annotation/human/Homo_sapiens.GRCh37.71.fixed.bed";
my $refgene;
my $now_string = getCurrentTime();
my $skipFlag = 0;
my $reponly;
my $qconly;
my $qcon; 							## this goes from after tophat to the end of the script
my $trackon; 						## run from the making tracks onward
my $processonly;  
my $processon; 


GetOptions(
	"alignGTF=s"	=>		\$alignGTF,						## this overrides the default GTF options
	"exp=s"			=>		\$exp,
	"nofastqc"		=>  	\$nofastqc,
	"dry"			=>		\$dry,
	"useqonly"		=>		\$useqonly,
	"trackonly"		=>		\$trackonly,
	"notophat"		=>		\$notophat,
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
	"noqc"			=>		\$noqc,	
	"reponly"		=> 		\$reponly,
	"fastqconly"	=>		\$fastqconly,
	"tophatonly"    =>      \$tophatonly,
	"qconly"		=>		\$qconly,
	"qcon"			=>		\$qcon,
	"trackon"		=> 		\$trackon,
	"processonly"	=>		\$processonly,
	"processon" 	=>		\$processon,
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

open my $RUNLOG, '>', $analysisDir.$now_string."_".$exp."runlog.txt" or die "error opening run log: $!\n";
print $RUNLOG "$now_string\nrunlog for experiment $exp\n\n";

print "Experiment $exp\n";
$analysisDir .= $exp."/";
system("mkdir -p $analysisDir");
chdir($analysisDir) or die "$!";
print $RUNLOG "\n$analysisDir created successfully";

if ($log){
	my $logfile = $analysisDir."runExp_".$now_string.".log";
	## this code allows all STDOUT and STDERR to print to the log file. 
	open (LOG, ">$logfile");
	*STDERR = *LOG;
	*STDOUT = *LOG;
}
# open ERROR,  '>', $logfile  or die $!;
# STDERR->fdopen( \*ERROR,  'w' ) or die $!;

unless ( -e $samplesheet ){
	die "Invalid Sample Sheet";
}

my $sampleHash = tools::meta::makeExpSampleHash(
	exp => $exp,
	samplesheet => $samplesheet,
	);

$species = tools::meta::getSpecies(
	exp => $exp,
	sampleHash => $sampleHash,
	); 
print "species\t$species\n";

$fastqDir = tools::fq::findFastQPath(
	exp 		=> $exp,
	path 		=> $dataPath,
	sampleHash 	=> $sampleHash,
	) unless ($nofastqc && $notophat);

if ($reponly){
	print "report only\n";
	goto FINISH;
}

if ( $qconly ){
	print "QC only\n";
	goto QC;
}

if ( $qcon ){
	print "QC on -- Skip alignment";
	goto QC;
}



if( defined($species) ){
	print $RUNLOG "Species:\t$species\n"; 
} else{
	print $RUNLOG "Species could not be determined from the sample sheet\n";
	die "No Species defined";
}

if ($fastqconly){
	print "fastQC only\n";
	goto FASTQC;
}

if ($nofastqc){
	print "nofastqc\n";
	goto ENDFASTQC;
}

if ( $tophatonly ){
	print "Tophat only\n";
	goto ENDFASTQC;
}

## check if USEQ only
if ($useqonly) {
	print "Running USEQ only\n";
	goto USEQ;
}

## check for track only
if ($trackonly) {
	print "Making browser tracks only\n";
	goto TRACK;
}

if ($trackon) {
	print "Going from making tracks on\n";
	goto TRACK;
}

if ( $processonly || $processon ){
	goto PROCESS;
}

FASTQC:
## Step --- fastQC 
tools::fq::fastQC(
	analysisDir => $analysisDir,
	fastqDir => $fastqDir
	) unless $nofastqc;

#check that all fastqc dirs were created
print $RUNLOG "\n\nFASTQC VALIDATION\n";


if ($fastqconly){
	goto FINISH;
}
ENDFASTQC:

## Step -- TOPHAT

tools::fq::runTophat(
	exp 	=> $exp,
	sampleHash => $sampleHash,
	fastqDir 	=> $fastqDir,
	analysisDir => $analysisDir,
	species => $species,
	dry 		=> $dry,
	alignGTF => $alignGTF,
	) unless ($notophat);

print $RUNLOG "\n\nTOPHAT VALIDATION\n";



if ($tophatonly){
	goto FINISH;
}

QC:
## Step -- QC

print("running QC\n");

tools::bam::processBam(
	species => $species,
	inPattern => $inPattern,
	RUNLOG => $RUNLOG
	);



if ( $qconly ){
	goto FINISH;
}

PROCESS:

if ($processonly || $processon){
	tools::bam::makeProcessedBam( RUNLOG => $RUNLOG );
}

if ($processonly) {
	goto FINISH;
}


TRACK:
## make Genome tracks

tools::bam::makeTracks(
	analysisDir => $analysisDir,
	RUNLOG => $RUNLOG
	);
## check for track only
if ($trackonly) {
	print "Qconly - going to FINSH\n";
	goto FINISH;
}

## does differenial expression analysis
USEQ:
print "Starting USEQ DefinedRegionDifferentialSeq\n";

## differential expression 
tools::bam::runDRDS(
	bamID=> 	"marked_dup",
	sampleHash => $sampleHash,
	analysisDir => $analysisDir,
	dry 	=>	$dry, 
	species => $species,
	);

tools::bam::runDRDS(
	bamID=> 	"processed",
	sampleHash => $sampleHash,
	analysisDir => $analysisDir,
	dry 	=>	$dry, 
	species => $species,
	);

if ($useqonly){
	goto FINISH;
}

FINISH:

# ## addRows to sample database
tools::db::addSamples( sh => $sampleHash, fastqdir => $fastqDir, species => $species );

# tools::db::addBams{ 'sh' , $sampleHash, 'type' , "accepted_hits_marked_dup" };
tools::db::addBams( sh=>$sampleHash, type=>"accepted_hits_marked_dup" );
tools::db::addMarkDup(sh =>$sampleHash);
tools::db::addBamStat(sh => $sampleHash);
tools::db::addReadDistribution(sh => $sampleHash);

## not essential feature right now
# tools::report::checkFastQC(
# 	RUNLOG =>	$RUNLOG,
# 	analysisDir => $analysisDir,
# 	);
tools::report::checkTophat(
	RUNLOG => $RUNLOG,
	sampleHash => $sampleHash,
	analysisDir => $analysisDir,
	);

tools::report::checkRunQC(
	RUNLOG => $RUNLOG,
	sampleHash => $sampleHash,
	analysisDir => $analysisDir,
	);

tools::report::checkDRDS(
	bamID=>	"processed",
	sampleHash => $sampleHash,
	analysisDir => $analysisDir,
	);
tools::report::checkDRDS(
	bamID=>	"marked_dup",
	sampleHash => $sampleHash,
	analysisDir => $analysisDir,
	);
my $repDir = "QCreport/";
system("mkdir -p $repDir");
my $filename = "${repDir}/${exp}-QCreport.html";
print "filename: ".$filename."\n";
open my $QCreport, "> $filename" or die "could not open QC report\n";
tools::report::writeHeader(fh => $QCreport, exp => $exp, curtime => $now_string );
tools::report::writeLegend(fh => $QCreport, exp => $exp );
tools::report::writeBamStat(fh => $QCreport, exp => $exp  );
tools::report::writeMarkDup(fh => $QCreport, exp => $exp  );

tools::report::writeFigTables(fh => $QCreport, 
	sampleHash => $sampleHash,
	repdir => $repDir,
	analysisDir => $analysisDir,
	);

tools::report::writeFooter(fh => $QCreport);

print "analysis complete\n";
close( $RUNLOG );
close ( LOG );

sub runAndLog{
	my $command = shift;
	my $time = localtime;
	print "$time\t$command";
	system($command);
}

sub getCurrentTime{
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	$mon = $mon + 1;
	$mon = prependZero($mon);
	$mday = prependZero($mday);
	$hour = prependZero($hour);
	$min = prependZero($min);
	$sec = prependZero($sec);

	my $curTime = $year + 1900 ."-".( $mon)  ."-". $mday . "-" . $hour . "-" . $min . "-" . $sec ."-";
	$curTime =~ s/\s/_/g;
	return $curTime;
}

sub prependZero{
	my $num = shift;
	if (length ($num ) == 1 ){
		$num = "0" . $num;
	}
	return $num;
}
