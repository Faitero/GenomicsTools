#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename;
use File::Glob;
use IO::Handle;
use b2b::tools;
use b2b::qc;
use b2b::runval;
use b2b::db;
use b2b::report;
use TimUtil;
# use diagnostics -verbose;
my $debug = 1;
my $dry = 0;

my $log = 0;
my	 $samplesheet = "/home/laurentt/sampleSheet/BruneauExperimentsGNOMex20130325.txt";
# my 	 $samplesheet = "/home/laurentt/sampleSheet/gnomex_sampleannotations_SEIDMAN_030713.txt";
my	 $dataPath = "/Data01/gnomex/ExperimentData/";
my 	 $exp; 
my	 $analysisDir = "/Data01/gnomex/Analysis/experiment/";
my 	 $fastqDir;
my 	 $nofastqc;
my 	 $fastQConly;
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
my $noQC;							## prevent the runQC sub from being run
my $noCatTab;				
my $species;
my $mouserefgene = "/work/Common/Data/Annotation/mouse/mm9/Mus_musculus.NCBIM37.67.fixed.bed";
my $humanrefgene = "/work/Common/Data/Annotation/human/Homo_sapiens.GRCh37.71.fixed.bed";
my $refgene;
my $now_string = TimUtil::getCurrentTime();
$now_string =~ s/\s/_/g;
my $skipFlag = 0;
my $repOnly;

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
	"qcTabOnly"		=>		\$qcTabOnly,
	"repOnly"		=> 		\$repOnly,
	"fastQConly"	=>		\$fastQConly,
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
chdir($analysisDir);
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

my $sampleHash = b2b::tools::parseSampleSheet($samplesheet);
my $expSampleHash = b2b::tools::makeExpSampleHash(
	exp => $exp,
	sampleHash => $sampleHash,
	);

# print Dumper($expSampleHash);
# <STDIN>;

$species = b2b::tools::getSpecies(
	exp => $exp,
	sampleHash => $sampleHash,
	); 

$fastqDir = b2b::tools::findFastQPath (
	exp 		=> $exp,
	path 		=> $dataPath,
	sampleHash 	=> $sampleHash,
	) unless ($nofastqc && $noTophat);



# ## addRows to sample database
b2b::db::addSamples( sh => $expSampleHash, fastqdir => $fastqDir );


b2b::db::addBams(
	sh => $expSampleHash,
	type => "accepted_hits_marked_dup"
);

b2b::db::addMarkDup(sh =>$expSampleHash);
b2b::db::addBamStat(sh => $expSampleHash);
b2b::db::addReadDistribution(sh => $expSampleHash);


my $repDir = "QCreport";
system("mkdir -p $repDir");
my $filename = "${repDir}/${now_string}-${exp}-QCreport.html";
print "filename: ".$filename."\n";

open my $QCreport, "> $filename" or die "could not open QC report\n";
b2b::report::writeHeader(fh => $QCreport, exp => $exp );
b2b::report::writeLegend(fh => $QCreport, exp => $exp );
b2b::report::writeBamStat(fh => $QCreport, exp => $exp  );
b2b::report::writeMarkDup(fh => $QCreport, exp => $exp  );

b2b::report::writeFooter(fh => $QCreport);

die;

print "species\t$species\n";

if( defined($species) ){
	print $RUNLOG "Species:\t$species\n"; 
} else{
	print $RUNLOG "Species could not be determined from the sample sheet\n";
	die "No Species defined";
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

## Step --- fastQC 

my $fastqcDirCommand = "mkdir -p ${analysisDir}fastqc";
# print "$fastqcDirCommand\n";

b2b::tools::runAndLog($fastqcDirCommand);

if(-d "${analysisDir}fastqc"){
	print $RUNLOG "${analysisDir}fastqc created successfully";
} else {
	print $RUNLOG "unable to create ${analysisDir}fastqc -- Error\n";
}

my $fastQCcommand = "fastqc --outdir=${analysisDir}fastqc ${fastqDir}*" unless ($nofastqc);

b2b::tools::runAndLog("$fastQCcommand") unless($dry || $nofastqc);

#check that all fastqc dirs were created
print $RUNLOG "\n\nFASTQC VALIDATION\n";

b2b::runval::checkFastQC(
	RUNLOG =>	$RUNLOG,
	fastqDir => $fastqDir,
	analysisDir => $analysisDir
	);


if ($fastQConly){
	goto FINISH;
}

b2b::tools::runTophat(
	exp 	=> $exp,
	sampleHash => $expSampleHash,
	fastqDir 	=> $fastqDir,
	analysisDir => $analysisDir,
	dry 		=> $dry,
	) unless ($noTophat);

print $RUNLOG "\n\nTOPHAT VALIDATION\n";

b2b::runval::checkTophat(
	RUNLOG => $RUNLOG,
	sampleHash => $expSampleHash,
	analysisDir => $analysisDir,
	);

print("running QC\n");

if (lc($species) eq "human"){
	print "species: human\n";
	$refgene = $humanrefgene;
} elsif (lc($species) eq "mouse"){
	print "species: mouse\n";
	$refgene = $mouserefgene;
} else {
	print $RUNLOG "Invalid Species, dying\n";
	die "species not recognized -- contact bioinformatician\n";
}

print "refgene = $refgene\n";
print $RUNLOG "species = $species\trefgene = $refgene\n";

print "searching for inpattern = \t$inPattern";
my @inputFiles = glob $inPattern;

if (@inputFiles == 0 ){
	print $RUNLOG "Tried to run QC without valid BAM files that match the in Pattern, dying\n"; 
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
	
    print $RUNLOG "\n\nMARK DUP VALIDATION\n";
	b2b::runval::checkIfExists( file=>$markedDupFile , RUNLOG=>$RUNLOG);
    my $indexFile = $markedDupFile;
    $indexFile =~ s/.bam$/.bai/;
	b2b::runval::checkIfExists( file=>$indexFile , RUNLOG=>$RUNLOG);
    
    my $bamstattab = $file."bam_stat.tab";
    runAndLog("bam_stat.py -i $file 2> $file.bam_stat.tab") unless $noBamStat;   

    print $RUNLOG "\n\nBAMSTAT TAB VALIDATIONfor $file\n";
	b2b::runval::checkIfExists( file=>$bamstattab , RUNLOG=>$RUNLOG);
	
	b2b::qc::runQC({
		file => $file,
		refgene => $refgene,
		}) unless $noQC;	
	
	print $RUNLOG "\n\nrunQC VALIDATION for $file\n";
	b2b::runval::checkRunQC(
		outfolder => $file."_QC/",
		RUNLOG => $RUNLOG,
		);
}	

print "Beginning to process marked files\n";
my @markedFiles = glob "Tophat*/*marked_dup.bam";

b2b::db::addBams{
	sh => $expSampleHash,
	type => "accepted_hits_marked_dup"
};

b2b::db::addMarkDup(sh =>$expSampleHash);

for my $file (@markedFiles){	
	print "processing\t$file";
	my $baiFile = $file;
	$baiFile =~ s/.bam$/.bai/;	
	print $RUNLOG "\n\nmarkedDup BAM index file validation\n";
	b2b:runval::checkIfExists(file=>$baiFile, RUNLOG=>$RUNLOG);

	$file = filterBam($file);
	print "processed file:\t$file";
	
	buildBamIndex($file);

	$baiFile = $file;
	$baiFile =~ s/.bam$/.bai/;

	print $RUNLOG "\n\nprocessed BAM index file validation\n";
	b2b:runval::checkIfExists(file=>$baiFile, RUNLOG=>$RUNLOG);

# code below is unnecessary for the new analysis plan
	# runQC({
	# 	file => $file,
	# 	refgene => $refgene,
	# });

	# print $RUNLOG "\n\nrunQC VALIDATION for $file\n";
	# b2b::runval::checkRunQC(
	# 	outfolder => $file."_QC/",
	# 	RUNLOG => $RUNLOG,
	# 	);
}

TRACK:
## make Genome tracks

my @files = glob("${analysisDir}Tophat*/*marked_dup.bam");

print $RUNLOG "\n\nBROWSER TRACK VALIDATION\n";
for my $file (@files){
	
	my $trackGenCommand = "convert_SAM_or_BAM_for_Genome_Browser.pl --nosort $file";
	print "$trackGenCommand\n";
	b2b::tools::runAndLog($trackGenCommand) unless ($dry);
	b2b::runval::checkBrowserTracks( file => $file , RUNLOG => $RUNLOG );
}

print "Finished making browser tracks for markedDup.bam's\n";

@files = glob("${analysisDir}Tophat*/*processed.bam");

print $RUNLOG "\n\nBROWSER TRACK VALIDATION\n";
for my $file (@files){
	my $trackGenCommand = "convert_SAM_or_BAM_for_Genome_Browser.pl --nosort $file";
	print "$trackGenCommand\n";
	b2b::tools::runAndLog($trackGenCommand) unless ($dry);
	b2b::runval::checkBrowserTracks( file => $file , RUNLOG => $RUNLOG );
}

print "Finished making browser tracks for processed.bam's\n";

## check for track only
if ($TRACKonly) {
	goto FINISH;
}

## does differenial expression analysis
USEQ:
print "Starting USEQ DefinedRegionDifferentialSeq\n";

## differential expression 
b2b::tools::runDRDS(
	bamID=> 	"marked_dup",
	sampleHash => $expSampleHash,
	analysisDir => $analysisDir,
	dry 	=>	$dry, 
	species => $species,
	);

b2b::tools::runDRDS(
	bamID=> 	"processed",
	sampleHash => $expSampleHash,
	analysisDir => $analysisDir,
	dry 	=>	$dry, 
	species => $species,
	);

b2b::runval::checkDRDS(
	bamID=>	"processed",
	sampleHash => $expSampleHash,
	analysisDir => $analysisDir,
	);

b2b::runval::checkDRDS(
	bamID=>	"marked_dup",
	sampleHash => $expSampleHash,
	analysisDir => $analysisDir,
	);

if ($USEQonly){
	goto FINISH;
}


FINISH:
print "analysis complete\n";
close( $RUNLOG );
close ( LOG );

sub runAndLog{
	my $command = shift;
	my $time = localtime;
	print "$time\t$command";
	system($command);
}