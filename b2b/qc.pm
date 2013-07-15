package b2b::qc;
## test to save
use strict;
# use warnings;
use Getopt::Long;
use IO::Handle;
use diagnostics -verbose;


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

sub runQC{
	my ($args) = @_;
	my $outdir   = $args->{outdir};

	
	my $file	= $args->{file};
	my $refgene = $args->{refgene};

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

	print "\$qcTabPattern\t$qcTabPattern";
	# my @qcTables = glob $baseFolder.$qcTabPattern;
	my @folders = glob $baseFolder.$bamFolder;
	print "folders\t@folders";
	my $arraySize = @folders;
	# print "$arraySize\t entries\n";
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
				
				# print "match 1: $1\\\n";
				# print "match 2: $2\\\n";
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
		# print $table."\n";
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
		print "sampleName\t$sampleName";
	
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
			
			# print "match 1: $1\\\n";
			# print "match 2: $2\\\n";
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
	print "$l\n";
	my $k = @{$outTable[0]};
	print "there are $k rows in the table\n";

	for ($i = 0 ; $i < $k ; $i++ ){			
			for ($j = 0 ; $j < $l; $j++ ){
				print ALLTAB "$outTable[$j][$i]\t";
				# print  "$outTable[$j][$i]\t";

				# print "i:$i\tj:$j\t"
			}
		print "\n";

		print ALLTAB "\n";
	}
		close ALLTAB;
	print "Finished making $baseFolder/bam_stat_all_files.tab\n";

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