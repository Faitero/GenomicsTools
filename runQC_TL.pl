#!/usr/bin/perl -w

use strict;
# use warnings;
use Getopt::Long;
use IO::Handle;
my $test = 0;


# testing block
if($test){
catQCtab();
die;
}


# This file is to run several of the QC metrics in the rseqc package

my $inPattern = "Tophat*/accepted_hits.bam";
my $startingDir = "./";
my $outDir		= "output/";
my $noMarkDup; 						## this flag will skip the dedup step.		
my $noClipping;						## this flag will skip the clipping profile
my $markedDupFile;
# This is the bd they recommend my $refgene = "/work/Common/Data/mouse/mm9/mm9_Ensembl_gene.bed";	## this is a bed annotation of the genes
# our Mus...clean.gtf -> Mus...clean.bed
my $refgene = "/work/Common/Data/Annotation/mouse/mm9/Mus_musculus.NCBIM37.67.fixed.bed";
my $now_string = localtime;
$now_string =~ s/\s/_/g;
my $logfile = $startingDir."runQC_".$now_string.".log";
open ERROR,  '>', $logfile  or die $!;
STDERR->fdopen( \*ERROR,  'w' ) or die $!;
# printL STDERR "printing in Main\n";

# testSTDERRprintInSub();


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
my $noCatTab;						## prevent cating together QC tables
my $species = "mouse";							## specify use to human bed file

GetOptions(
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
	);

my $outPath = $startingDir.$outDir;

if (lc($species) eq "human"){
	printL ("species: human\n");
	my $refgene = "/work/Common/Data/Annotation/human/2011_Archive/Homo_sapiens.GRCh37.56.chr.bed"
} elsif (lc($species) eq "mouse"){
	printL ("species: mouse\n");
	my $refgene = "/work/Common/Data/Annotation/mouse/mm9/Mus_musculus.NCBIM37.67.fixed.bed";
} else {
	die "species not recognized -- contact bioinformatician\n";
}

# my $SORTSAM_PATH = `which SortSam.jar`;
# my $MARKDUPLICATES_PATH = `which MarkDuplicates.jar`;

# runAndLog("mkdir -p $outDir");
printL("searching for inpattern = \t$inPattern");
my @inputFiles = glob $inPattern;

if (@inputFiles == 0 ){
	die "no Input Files matching the pattern\n";
}

for my $file (@inputFiles){
	printL("Processing $file");
	buildBamIndex($file);
	## make folder for each file
	my $outfolder = $file;
	if($outfolder =~ m/(Tophat.+\/)/){
		printL("outfolder\t$1\n");
		$outfolder = $1;
	}
	my $markedDupFile = markDup($file) unless $noMarkDup; 

    
    $file = $markedDupFile unless !$markedDupFile;
     

    runAndLog("bam_stat.py -i $file 2> $file.bam_stat.tab") unless $noBamStat;   
	runQC({
		file => $file,
		refgene => $refgene,
		}) unless $noQC;	

}	



# use Data::Dumper;
# printL Dumper @outTable;

#Filter the marked reads

# unless($noFilter && $noMarkDup){
	printL("Beginning to process marked files");
	my @markedFiles = glob "Tophat*/*marked_dup.bam";
	for my $file (@markedFiles){	
		printL("processing\t$file");
		## make folder for each file

		$file = filterBam($file);
		printL("processed file:\t$file");
		
		buildBamIndex($file);

		my $outfolder = $file;
		my $outFileRoot = $file;
		
		runQC({
			file => $file,
			refgene => $refgene,
		});
	}
# }


# catQCtab() unless ($noCatTab && $noBamStat);

sub runQC{
	my ($args) = @_;
	my $outdir   = $args->{outdir};
	my $outfile  = $args->{outfile};
	
	my $file	= $args->{file};
	my $refgene = $args->{refgene};

	if( !$outdir ){
		$outdir = $file;
		$outdir .= "_QC/";
	}
	printL("outdir\t$outdir");
	my $outpath = $outdir.$outfile;
	
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

	readQuality({
		outpath => $outpath,
		file 	=> $file,
		refgene => $refgene,
		});

	readDuplication({
		outpath => $outpath,
		file 	=> $file,
		refgene => $refgene,
		});

	readGC({
		outpath => $outpath,
		file 	=> $file,
		refgene => $refgene,
		});

	readNVC({
		outpath => $outpath,
		file 	=> $file,
		refgene => $refgene,
		});

	RPKMsat({
		outpath => $outpath,
		file 	=> $file,
		refgene => $refgene,
		});

	juncAnnot({
		outpath => $outpath,
		file 	=> $file,
		refgene => $refgene,
		});

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
	runAndLog("read_duplication.py -i $file -o ${outpath}read_duplication -u 20000");
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


sub runAndLog{
	my $command = shift;
	my $time = localtime;
	printL("$time\t$command");
	system($command);
}

# Prints to standard out and to the log file
sub printL{
	my $string = shift;
	print STDERR "$string\n\n";
	print "$string\n\n";
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
		    #. ((!$shouldSort) ? qq{ ASSUME_SORTED=TRUE } : qq { }) ## <-- if we use --nosort, then ASSUME they are sorted no matter what!
		    . qq{ MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=$maxFileHandles }
		    . qq{ OPTICAL_DUPLICATE_PIXEL_DISTANCE=10 } ## <-- 100 is default but our libraries have single integer coordinate values
		    . qq{ METRICS_FILE=$dedupExtraMetricsFile }
		    . qq{ OUTPUT=$outFile }
		    

	);
    runAndLog($dedupCmd);
    return $outFile;
}

# Tests log file printing in sub
sub testSTDERRprintInSub{
	printL("printing in Sub");
}

# pastes together QC tables
sub catQCtab{
	my ($args) = @_;
	my $baseFolder = $args->{baseFolder};
	
	my $bamFolder = $args->{bamFolder};

	my $qcTabPattern =$args->{qcTabPattern};
	
	my $markDupPattern = $args->{markDupPattern};


	$baseFolder = "./" unless $baseFolder;
	
	$bamFolder = "Tophat*/" unless $bamFolder;

	$qcTabPattern = "*.bam_stat.tab" unless $qcTabPattern;

	$markDupPattern = "*_mark_dup_metrics.txt" unless $markDupPattern;

	printL "\$qcTabPattern\t$qcTabPattern";
	# my @qcTables = glob $baseFolder.$qcTabPattern;
	my @folders = glob $baseFolder.$bamFolder;
	printL("folders\t@folders");
	my $arraySize = @folders;
	# printL "$arraySize\t entries\n";
	my @outTable; # = ( [ "sample" , "Total_Records" , "QC_failed" , "Optical_PCR_Duplicates" ,
	 # ] );
	## build file fields
	my @sampleCollection;

	my @qcTables = glob $baseFolder.$bamFolder.$qcTabPattern;
	my @dupStats = glob $baseFolder.$bamFolder.$markDupPattern;

	for ( my $i = 0 ; $i < $arraySize ; $i++){
		print "$i\t$qcTables[$i]\t$dupStats[$i]\n";
	}
	

	# my $dupHeader = `cat $dupStats[0] | grep LIBRARY`;
	# print "$dupHeader\n";
	my @dupHeadArray = split('\t', `grep LIBRARY $dupStats[0]`);
	# my @dupHeadArray = split('\t', `cat $dupStats[0] | grep Unknown`);

	for ( my $i = 0 ; $i < @dupHeadArray; $i++ ){
		print "$i\t$dupHeadArray[$i]\n";
	}


	
	open IN, $qcTables[0] or die $!;
	my $i = 1;
	$outTable[0][0] = "Sample:";
	## make header
	while(<IN>){
		if( $i == 4){
			$outTable[0][4] = $dupHeadArray[4];
			$outTable[0][5] = $dupHeadArray[5];
			$outTable[0][6] = $dupHeadArray[6];
			$outTable[0][7] = $dupHeadArray[7]; 
			$i = 8;
		}


		if( $_ =~ m/^(.+)\s+(\d+)$/ ){
				
				# printL "match 1: $1\\\n";
				# printL "match 2: $2\\\n";
				my $entry = $1;
				$entry =~ s/^\s*(.*?)\s*$/$1/;
				$outTable[0][$i] = $entry;
				$i++;
			}
	}

	# for( $i = 0 ; $i < @{$outTable[0]}; $i++){
	# 	print "$i\t$outTable[0][$i]\n";
	# }

	# die;
	close IN;
	open ALLTAB, "> $baseFolder/quality_summary.tab" or die $!;

	my $j = 1; 	## sample counter
	## get data for tables
	print @folders."\n";
	for my $folder (@folders){
		# my $pause = <STDIN>;
		# printL $table."\n";
		my $table = `ls ${folder}${qcTabPattern}`;
		my $dupTab = `ls ${folder}${markDupPattern}`;
		print "$table\n";
		print "$dupTab\n";
		# print $dupTab;
		# my $test = `grep Unknown $dupTab`;

		# print $test."\n";

		

		my @dupData = split('\t', `grep Unknown $dupTab`);

		for ($i = 0 ; $i < @dupData ; $i ++){

			print "$i\t$dupData[$i]\n";
		}

		print "processing\t$table\t$dupTab\n"; 

		open IN, "$table" or die $!;
		# my $sampleName = $table;
		my $sampleName;
	
		if( $table =~ m/(Tophat.*\/)/){
			$sampleName = $1;
		}
		printL("sampleName\t$sampleName");
	
		$outTable[$j][0] = $sampleName;
		$i = 1;		## element counter
		while(<IN>){
			if ($i == 4 ){
				$outTable[$j][4] = $dupData[4];
				$outTable[$j][5] = $dupData[5];
				$outTable[$j][6] = $dupData[6];
				$outTable[$j][7] = $dupData[7]; 
				$i = 8;
			}


			if( $_ =~ m/^(.+)\s+(\d+)$/ ){
			
			# printL "match 1: $1\\\n";
			# printL "match 2: $2\\\n";
				$outTable[$j][$i] = $2;
				$i++;
			}
		}
		for ($i = 0 ; $i < @{$outTable[$j]} ; $i++){
			# $pause = <STDIN>;
			print "$i\t$outTable[$j][$i]\n";
		}
		close IN;
		$j++;

	}	
	my $l = @outTable;
	printL "$l\n";
	my $k = @{$outTable[0]};
	printL "there are $k rows in the table\n";

	for ($i = 0 ; $i < $k ; $i++ ){			
			for ($j = 0 ; $j < $l; $j++ ){
				print ALLTAB "$outTable[$j][$i]\t";
				# print  "$outTable[$j][$i]\t";

				# printL "i:$i\tj:$j\t"
			}
		print "\n";

		print ALLTAB "\n";
	}
		close ALLTAB;
	printL "Finished making $baseFolder/bam_stat_all_files.tab\n";

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
sub printL{
	my $string = shift;
	print STDERR "$string\n\n";
	print "$string\n\n";
}
sub runAndLog{
	my $command = shift;
	my $time = localtime;
	printL("$time\t$command");
	system($command);
}

