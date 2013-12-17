package tools::bam;
# this library contains code required to run QC metrics and differential expression analysis on 
# Bam files 
## test to save
use tools::runval;
use strict;
use warnings;
use Getopt::Long;
use IO::Handle;
use diagnostics -verbose;
use Data::Dumper;
my $dry = 0;

my $mouserefgene = "/work/Common/Data/Annotation/mouse/mm9/Mus_musculus.NCBIM37.67.fixed.bed";
my $humanrefgene = "/work/Common/Data/Annotation/human/Homo_sapiens.GRCh37.71.fixed.bed";
my $chickenrefgene = "/work/Common/Data/Annotation/chicken/Gallus_gallus.Galgal4.72.bed";


## this begins the QC process. This sub assumes that you are in the directory that contain
## the bam files. If not, pass the optional parameter to the sub expdir => <somthing> to change to
## the working directory 
sub processBam{
	my $refgene;
	my %args = @_;
	my $species = $args{species};
	my $inPattern = $args{inPattern};
	my $noMarkDup = $args{noMarkDup};
	my $noQC = $args{noQC};
	my $RUNLOG = $args{RUNLOG};
	if (exists($args{expDir})){
		chdir($args{expDir}) or die "$!";
	}
	if (lc($species) eq "human"){
		print "species: human\n";
		$refgene = $humanrefgene;
	} elsif (lc($species) eq "mouse"){
		print "species: mouse\n";
		$refgene = $mouserefgene;
	} elsif( lc( $species) eq "chicken" ){
		print "species: chicken\n";
		$refgene = $chickenrefgene;
	}
	else {
	die "species not recognized -- contact bioinformatician\n";
	}
	print "refgene = $refgene\n";
	print "Bam Pattern = $inPattern\n";
	my @inputFiles = glob $inPattern;
	if ((scalar(@inputFiles) == 0 )){
		print $RUNLOG "Tried to run QC without valid BAM files that match the in Pattern, dying\n"; 
		die "no Input Files matching the pattern\n";
	}

	for my $file (@inputFiles){
		print "Processing $file\n";
		my $outfolder = $file;
		if($outfolder =~ m/(Tophat.+\/)/){
			print("outfolder\t$1\n");
			$outfolder = $1;
		}
		my $markedDupFile = markDup($file) unless $noMarkDup; 
	    $file = $markedDupFile unless !$markedDupFile;
		
	    print $RUNLOG "\n\nMARK DUP VALIDATION\n";
		tools::runval::checkIfExists( file=>$markedDupFile , RUNLOG=>$RUNLOG);
	    my $indexFile = $markedDupFile;
	    $indexFile =~ s/.bam$/.bai/;
		tools::runval::checkIfExists( file=>$indexFile , RUNLOG=>$RUNLOG);
		
		unless ( $noQC ){
			runAndLog("bam_stat.py -i $file 2> $file.bam_stat.tab");
			runRSEQC({
			file => $file,
			refgene => $refgene,
			});	
			print $RUNLOG "\n\nrunQC VALIDATION for $file\n";
			tools::runval::checkRunQC(
			outfolder => $file."_QC/",
			RUNLOG => $RUNLOG,
			);
		}
	}	
}


sub makeProcessedBam{
	my %args = @_;
	my $RUNLOG = $args{RUNLOG};
	print "Beginning to process marked files\n";
	my @markedFiles = glob "Tophat*/*marked_dup.bam";

	for my $file (@markedFiles){	
		print "processing\t$file";
		my $baiFile = $file;
		$baiFile =~ s/.bam$/.bai/;	
		print $RUNLOG "\n\nmarkedDup BAM index file validation\n";
		
		$file = filterBam($file);
		print "processed file:\t$file";
		buildBamIndex($file);
		$baiFile = $file;
		$baiFile =~ s/.bam$/.bai/;
		print $RUNLOG "\n\nprocessed BAM index file validation\n";
		tools::report::checkIfExists(file=>$baiFile, RUNLOG=>$RUNLOG);
	}
}



## this runs the convert_SAM_or_BAM_for_Genome_Browser.pl on a collection of bams
sub makeTracks{
	my %args = @_;
	my $analysisDir = $args{analysisDir};
	my $pattern1 = $args{pattern1};
	my $pattern2 = $args{pattern2};
	my $RUNLOG = $args{RUNLOG};
	if (!$pattern1){
		$pattern1 = "Tophat*/*marked_dup.bam";
	}
	if (!$pattern2){
		$pattern2 = "Tophat*/*processed.bam";
	}
	my @files = glob(${analysisDir}.$pattern1);
	
	if (@files == 0){
		die "Track Errror: Could not find ".$pattern1."\n"
	}
	print $RUNLOG "\n\nBROWSER TRACK VALIDATION\n";
	for my $file (@files){
		my $trackGenCommand = "convert_SAM_or_BAM_for_Genome_Browser.pl --nosort $file";
		print "$trackGenCommand\n";
		system($trackGenCommand) unless ($dry);
		tools::runval::checkBrowserTracks( file => $file , RUNLOG => $RUNLOG );
	}
	print "Finished making browser tracks for markedDup.bam's\n";

	@files = glob(${analysisDir}.$pattern2);
	if (@files == 0 ){
		"print no files with ".$pattern2."\n";
		return;
	}
	print $RUNLOG "\n\nBROWSER TRACK VALIDATION\n";
	for my $file (@files){
		my $trackGenCommand = "convert_SAM_or_BAM_for_Genome_Browser.pl --nosort $file";
		print "$trackGenCommand\n";
		system($trackGenCommand) unless ($dry);
		tools::runval::checkBrowserTracks( file => $file , RUNLOG => $RUNLOG );
	}
	print "Finished making browser tracks for processed.bam's\n";
}


## this runs USEQ's DRDS package
sub runDRDS{
	my %args 	= 	@_;
	my $bamID = $args{bamID};
	my $sampleHash = $args{sampleHash};
	my $analysisDir = $args{analysisDir};
	my $dry = $args{dry};
	my $species = $args{species};
	my $groupHash = {};
	my $analysisDirComplete = 0;

	my $USEQ_CONDdir = "${analysisDir}${bamID}_USEQ_CONDITIONS/";
	runAndLog("mkdir -p $USEQ_CONDdir");

	my $lsCommand = "ls ${analysisDir}Tophat*/accepted_hits_$bamID.bam";
	print $lsCommand."\n";
	my @bamFiles = `$lsCommand`;
	print "bamFiles = @bamFiles\n";
	for my $file (@bamFiles){
		chomp $file;
		print "file\t$file\n";
		
		if ( $file =~ qq/Tophat_(.+)\/(accepted_hits_${bamID}.bam)/){
			my $bamRoot = $1;
			print "bamRoot = $bamRoot\t\t";
			my $rootPattern = $bamRoot;
			$rootPattern =~ s/X\d+$//;
			print "rootPattern = $rootPattern\n";
			$groupHash->{$rootPattern}->{$bamRoot} = $bamRoot;
		} else{
			print "noMatch for $file\n";
			die "no Match for $file";
		}
	}

	for my $key ( keys ( %$groupHash ) ){
		my $mkdirCommand = "mkdir -p ${USEQ_CONDdir}$key";
		print "$mkdirCommand\n";
		runAndLog($mkdirCommand);
		for my $bamRoot ( keys %{ $groupHash->{ $key } } ){
			my $linkCommand = "ln -s -f ${analysisDir}Tophat_${bamRoot}/accepted_hits_$bamID.bam ${USEQ_CONDdir}$key/${bamRoot}${key}_accepted_hits_processed.bam";
			print "linkCommand = $linkCommand\n";
			runAndLog($linkCommand) unless ($dry);
		}
	}
	
	## now to setup the Useq command: 
	my $USEQ_GENES_MERGED_MM9="/work/Common/Data/Annotation/mouse/mm9/Mus_musculus.NCBIM37.67.clean.constitutive.table";
	my $USEQ_GENES_MERGED_HG19="/work/Common/Data/Annotation/human/Homo_sapiens.GRCh37.71.clean.constitutive.table";
	my $USEQ_GENES_MERGED_GG4="/work/Common/Data/Annotation/chicken/Gallus_gallus.Galgal4.72.table";
	my $MM9_VERS="M_musculus_Jul_2007"; # ← this only matters for hyperlinking
	my $HG19_VERS="H_sapiens_Feb_2009"; # ← this only matters for hyperlinking
	my $GG4_VERS="G_gallus_Nov_2011";
	my $UCSCGENES;
	my $GENVERSION;

	if (lc($species) eq "mouse" ){
		print "Species : mouse\n";
		$UCSCGENES=${USEQ_GENES_MERGED_MM9};
		$GENVERSION=${MM9_VERS};
	} elsif (lc($species) eq "human\n" ){
		print "Species : human";
		$UCSCGENES=${USEQ_GENES_MERGED_HG19};
		$GENVERSION=${HG19_VERS};
	} elsif (lc($species) eq "chicken\n" ){
		print "Species : human";
		$UCSCGENES=${USEQ_GENES_MERGED_GG4};
		$GENVERSION=${HG19_VERS};
	}else {
		die "Undefined species: script needs configuring : contact bioinformatician\n";
	}
	my $DRDSJAR = "/work/Apps/USeq_8.5.7/Apps/DefinedRegionDifferentialSeq";
	my $USEQOUT = "${analysisDir}${bamID}_USEQ_OUT/";
	system("mkdir -p $USEQOUT");
	my $GIGB=8;
	my $WHICHR = `which R`;
	chomp ($WHICHR);
	my $MINMAPPINGREADS=0;
	my $MIN_LOG10_SPACE_FDR=0;
	my $MIN_LOG2_SPACE_FOLD_CHANGE=0;
	my $CONDITIONDIR="${analysisDir}${bamID}_USEQ_CONDITIONS"; ## this has as SUB FOLDERS all the actual conditions!
	my $MAXBASEALIGNDEPTH=100000000;

	my $USEQCommand = "java -Xmx${GIGB}G -jar ${DRDSJAR}  -r ${WHICHR}  -s ${USEQOUT} ";
   $USEQCommand .= "-c ${CONDITIONDIR}  -g ${GENVERSION}  -x ${MAXBASEALIGNDEPTH} ";
   $USEQCommand .= "-f ${MIN_LOG10_SPACE_FDR}  -l ${MIN_LOG2_SPACE_FOLD_CHANGE} ";
   $USEQCommand .= "-e ${MINMAPPINGREADS}  -u ${UCSCGENES} -t";
   $USEQCommand .= "2>&1 | tee --append ${USEQOUT}useq.log.txt";
   print "$USEQCommand\n" ;
   runAndLog($USEQCommand) unless($dry);
}


## extracts all duplicate reads and puts in a new bam file
sub collectDups{
	my %args = @_;
	my $file = $args{file};
	my $dupfile = $file;
	$dupfile =~ s/marked_dup.bam$/dupes.bam/;
	my $gatherDupCommand = (
		qq{ samtools view -bhu -f 0x400 $file} ## only grab those marked as duplicates
		. qq{ | samtools view -hbu -F 0x004 - } ## remove unmapped reads
		. qq{ | samtools view -hbu -F 0x100 - } ## remove non-primary alignment
		. qq{ | samtools view -hbu -F 0x008 - } ## remove mate pair did not map
		. qq{ | samtools view -hb -F 0x200 -}   ## remove reads failed QC
		. qq{ > $dupfile });
	print $gatherDupCommand."\n";
	system($gatherDupCommand);
}


## this sub runs a host of the RSEQC programs.
sub runRSEQC{
	my %args = @_;
	print "args ". Dumper %args ."\n";
	my $outdir   = $args{outdir};

	
	my $file	= $args{file};
	my $refgene = $args{refgene};

	if( !$outdir ){
		$outdir = $file;
		$outdir .= "_QC/";
	}
	print ("outdir\t$outdir");
	my $outpath = $outdir;
	
	runAndLog("mkdir -p $outpath");

	geneCoverage({
		outpath => $outpath,
		file 	=> $file,
		refgene => $refgene,
		});

	innerDistance({
		outpath => $outpath,
		file 	=> $file,
		refgene => $refgene,
		});

	# readQuality({
	# 	outpath => $outpath,
	# 	file 	=> $file,
	# 	refgene => $refgene,
	# 	});

	readDuplication({
		outpath => $outpath,
		file 	=> $file,
		refgene => $refgene,
		});

	# readGC({
	# 	outpath => $outpath,
	# 	file 	=> $file,
	# 	refgene => $refgene,
	# 	});

	# readNVC({
	# 	outpath => $outpath,
	# 	file 	=> $file,
	# 	refgene => $refgene,
	# 	});

	RPKMsat({
		outpath => $outpath,
		file 	=> $file,
		refgene => $refgene,
		});

	# juncAnnot({
	# 	outpath => $outpath,
	# 	file 	=> $file,
	# 	refgene => $refgene,
	# 	});

	 juncSat({
	 	outpath => $outpath,
	 	file 	=> $file,
	 	refgene => $refgene,
	 	});

	readDist({
		outpath => $outpath,
		file 	=> $file,
		refgene => $refgene,
		});

}

sub geneCoverage{
	my ($args) = @_;
	my $outpath = $args->{outpath};
	my $file	= $args->{file};
	my $refgene = $args->{refgene};
	runAndLog("geneBody_coverage.py -i $file -r $refgene -o ${outpath}gene_body_coverage");
}

sub innerDistance{
	my ($args) = @_;
	my $outpath = $args->{outpath};
	my $file	= $args->{file};
	my $refgene = $args->{refgene};
	runAndLog("inner_distance.py -i $file -o ${outpath}inner_distance -r $refgene -u 700");
}

sub readQuality{
	my ($args) = @_;
	my $outpath = $args->{outpath};
	my $file	= $args->{file};
	my $refgene = $args->{refgene};
	runAndLog("read_quality.py -i $file -o ${outpath}read_quality");
}

sub readDuplication{
	my ($args) = @_;
	my $outpath = $args->{outpath};
	my $file	= $args->{file};
	my $refgene = $args->{refgene};
	# runAndLog("read_duplication.py -i $file -o ${outpath}read_duplication -u 20000");
	runAndLog("read_duplication.py -i $file -o ${outpath}read_duplication");
}

sub readGC{
	my ($args) = @_;
	my $outpath = $args->{outpath};
	my $file	= $args->{file};
	my $refgene = $args->{refgene};
	runAndLog("read_GC.py -i $file -o ${outpath}read_GC");
}

sub readNVC{
	my ($args) = @_;
	my $outpath = $args->{outpath};
	my $file	= $args->{file};
	my $refgene = $args->{refgene};
	runAndLog("read_NVC.py -i $file -o ${outpath}read_NVC -x");
}


sub RPKMsat{
	my ($args) = @_;
	my $outpath = $args->{outpath};
	my $file	= $args->{file};
	my $refgene = $args->{refgene};
	runAndLog("RPKM_saturation.py -i $file -o ${outpath}RPKM_saturation -r $refgene ");
}

sub juncAnnot{
	my ($args) = @_;
	my $outpath = $args->{outpath};
	my $file	= $args->{file};
	my $refgene = $args->{refgene};
	runAndLog("junction_annotation.py -i $file -r $refgene -o ${outpath}junction_annotation");
}

sub juncSat{
	my ($args) = @_;
	my $outpath = $args->{outpath};
	my $file	= $args->{file};
	my $refgene = $args->{refgene};
	runAndLog("junction_saturation.py -i $file -r $refgene -o ${outpath}junction_saturation");
}

sub readDist{
	my ($args) = @_;
	my $outpath = $args->{outpath};
	my $file	= $args->{file};
	my $refgene = $args->{refgene};
	runAndLog("read_distribution.py  -i $file -r $refgene > ${outpath}read_distribution.txt");
}	

sub buildBamIndex{
	# my $buildBamIndexPath = "/bioinformatics/installations/picard/picard-tools-1.89/";
	my $buildBamIndexPath = `which BuildBamIndex.jar`;
	chomp($buildBamIndexPath);
	my $file = shift;
	my $GIGABYTES_FOR_PICARD = 8;
	my $buildBamIndexCommand = (
		qq{java -Xmx${GIGABYTES_FOR_PICARD}g -jar ${buildBamIndexPath} }
		. qq{ INPUT=$file }
		);
	runAndLog($buildBamIndexCommand);
}

sub markDup{
	my $ULIMIT_RESULT = 1024; ## result of running the shell command ulimit -n. Since this is a shell built-in, it can, for some reason, not be run like a real command, so backticks don't work. `ulimit -n`; chomp($ULIMIT_RESULT);
	my $maxFileHandles = int(0.8 * $ULIMIT_RESULT); ## This is for the MAX_FILE_HANDLES_FOR_READ_ENDS_MAP parameter for MarkDuplicates: From the Picard docs: "Maximum number of file handles to keep open when spilling read ends to disk. Set this number a little lower than the per-process maximum number of file that may be open. This number can be found by executing the 'ulimit -n' command on a Unix system. Default value: 8000."
	# my $MARKDUPLICATES_PATH = `which MarkDuplicates.jar`; chomp($MARKDUPLICATES_PATH); ## Tries to find "MarkDuplicates.jar" from the Picard suite in the user's $PATH
	my $GIGABYTES_FOR_PICARD = 8;

	# my $MARKDUPLICATES_PATH = "/work/Apps/Bio/picard/picard-tools-1.90/MarkDuplicates.jar"; 
	my $MARKDUPLICATES_PATH = `which MarkDuplicates.jar`; 
	chomp($MARKDUPLICATES_PATH);
	# my $SORTSAM_PATH = "/bioinformatics/installations/picard/picard-tools-1.89/SortSam.jar";
	my $SORTSAM_PATH = `which SortSam.jar`;
	chomp($SORTSAM_PATH);
	my $inFile = shift;
	my $outFile = $inFile;
	$outFile =~ s/.bam$//;
	$outFile .= "_marked_dup.bam";
	my $dedupExtraMetricsFile = $outFile;
	$dedupExtraMetricsFile .= "_mark_dup_metrics.txt";
	
	my $dedupCmd = (qq{java -Xmx${GIGABYTES_FOR_PICARD}g -jar ${MARKDUPLICATES_PATH} }
		    . qq{ INPUT=$inFile } ## Picard SortSam.jar accepts both SAM and BAM files as input!
		    . qq{ REMOVE_DUPLICATES=FALSE }
		    . qq{ CREATE_INDEX=TRUE }
		    #. ((!$shouldSort) ? qq{ ASSUME_SORTED=TRUE } : qq { }) ## <-- if we use --nosort, then ASSUME they are sorted no matter what!
		    . qq{ MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=$maxFileHandles }
		    . qq{ OPTICAL_DUPLICATE_PIXEL_DISTANCE=10 } ## <-- 100 is default but our libraries have single integer coordinate values
		    . qq{ METRICS_FILE=$dedupExtraMetricsFile }
		    . qq{ OUTPUT=$outFile }
	);
    runAndLog($dedupCmd);
    return $outFile;
}

sub filterBam{
	my $file = shift;
	my $deDupFile = $file;
	$deDupFile =~ s/marked_dup.bam$//;
	$deDupFile .= "processed.bam";
	my $filterMappedOnlyCmd = (qq{     samtools view -hbu -F 0x004  $file } ## <-- only include the mapped (mapping flag is "4" apparently) reads, only include PRIMARY alignments, and skip any failing-QC reads
	   . qq{ | samtools view -hbu -F 0x100 - } ## remove "alignment is not primary"
	   . qq{ | samtools view -hbu -F 0x400 - } ## <-- remove "read is PCR or optical duplicate"
	   . qq{ | samtools view -hbu -F 0x008 - } ## <-- remove "mate pair did not map" reads
	   . qq{ | samtools view -hb  -F 0x200 - } ## <-- remove "read fails QC". Note this DOES NOT have the '-u' flag since it's the last filter!
	   . qq{ > $deDupFile });   ## <-- output file location		
	runAndLog($filterMappedOnlyCmd);
	return $deDupFile;
}

# Prints to standard out and to the log file
sub runAndLog{
	my $command = shift;
	my $time = localtime;
	print ("$time\t$command");
	system($command);
}
1;