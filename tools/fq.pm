package b2b::tools;
use strict;
use Getopt::Long;
use Data::Dumper;
use File::Basename;
use File::Glob;
use IO::Handle;

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
	my $USEQ_GENES_MERGED_GG4="/work/Common/Data/Annotation/chickenGallus_gallus.Galgal4.72.tab"
	my $MM9_VERS="M_musculus_Jul_2007"; # ← this only matters for hyperlinking
	my $HG19_VERS="H_sapiens_Feb_2009"; # ← this only matters for hyperlinking
	my $GG4_VERS="galGal_Nov_2011";
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
	}elsif (lc($species) eq "chicken\n" ){
		print "Species : human";
		$UCSCGENES=${USEQ_GENES_MERGED_HG19};
		$GENVERSION=${GG4_VERS}; 
	else {
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

# this takes a experiment name and a reference to a sample hash and runs tophat on the relavent
# experiments
sub runTophat{
	my	 $HAT_ARGS="--min-anchor=5 --segment-length=25 --no-coverage-search --segment-mismatches=2 --splice-mismatches=2 --microexon-search --no-discordant --no-mixed ";
	my	 $NUM_THREADS=12;
	my	 $MM9_INDEX="/work/Common/Data/Bowtie_Indexes/mm9";
	my	 $HG19_INDEX="/work/Common/Data/Bowtie_Indexes/hg19";
	my   $GG4_INDEX="/work/Common/Data/Bowtie_Indexes/galGal4_ERCC92";
	my	 $HG_GTF="/work/Common/Data/Annotation/human/Homo_sapiens.GRCh37.71.fixed.gtf";
	my	 $MM_GTF="/work/Common/Data/Annotation/mouse/mm9/Mus_musculus.NCBIM37.67.fixed.gtf";
	my 	 $GG4_GTF="/work/Common/Data/Annotation/ercc/galGal4.72_ERCC92.gtf";
	my 	 $transcriptomeIndex = "/work/Common/Data/Annotation/bowtie_transcriptome_index";
	my	 $TOPGTF;
	my	 $TOPINDEX;
	my $PAIRED_MATE_INNER_DIST=250;
	my %args 	= 	@_;
	my $exp 	= 	$args{exp};
	my $sampleHash = $args{sampleHash};
	my $fastqDir  = $args{fastqDir};
	my $analysisDir = $args{analysisDir};
	my $dry = $args{dry};
	my $pairedEnd;
	my $read1;
	my $read2;
	my $sample1 = $exp;
	$sample1 =~ s/R//;
	$sample1 .= "X1";
	print "sample : $sample1\n";

	my $org = ${$sampleHash}{$sample1}{'Organism'};

	# print "Organism = $org\n";
	if ( !defined($species)) {die "No Species Defined\n"};
	if ( lc($species) eq "mouse" ){
		print "Species : Mouse\n";
		$TOPGTF = $MM_GTF;
		$TOPINDEX = $MM9_INDEX;
		$transcriptomeIndex .= "/mm9";
	} elsif ( lc($species) eq "human"){
		print "Sample : Human\n";
		$TOPGTF = $HG_GTF;
		$TOPINDEX = $HG19_INDEX;
		$transcriptomeIndex .= "/hg19";
	} elsif ( lc($species) eq "chicken"){
		print "Sample : Chicken\n";
		$TOPGTF = $GG4_GTF;
		$TOPINDEX = $GG4_INDEX;
		$transcriptomeIndex .= "/gg4";
	}
	else {
		die "No definition for ".${$sampleHash}{$sample1}{'Organism'}." contact comeone who can fix this";
	}

	$exp =~ s/R//;
	for my $key ( keys( %$sampleHash )  ){
		$read2 = "";	
		$pairedEnd = 0;
		unless ( defined($sampleHash->{$key}->{"used"}) ){
		if ( $key =~ m/^${exp}X\d+/  ){
			my $files = ${$sampleHash}{$key}{"Associated Files"};
			$files =~ s/"//g;
			$files =~ s/,//g;
			print "files : $files\n";
			my @files = split(" ", $files);
			
			print "Sequencing Read Type:\t".lc(${$sampleHash}{$key}{"Sequencing Read Type"})."\n";
			if ( @files == 2 ){
				print "paired end library\n";
				$files[0] =~ s/[,;]//;
				$files[1] =~ s/[,;]//;
				$pairedEnd=1;
				$read1 = `find $fastqDir -name $files[0]*`;
				$read2 = `find $fastqDir -name $files[1]*`;
				$sampleHash->{$key}->{"used"} = 1;
			} elsif ( lc(${$sampleHash}{$key}{"Sequencing Read Type"}) =~ m/paired/ ){
				die("Annotated as a paired end read but only has 1 file listed");
			}
			else{
				print "single read library\n";
				$pairedEnd=0;
				$read1 = `find $fastqDir -name $files[0]*`;
				$sampleHash->{$key}->{"used"} = 1;
			}
			chomp $read1;
			chomp $read2;
			print "read1 : $read1\n";
			if ( defined($read2)){
				print "read2 : $read2\n";
			}
			if( !defined (${$sampleHash}{$key}{"Bam Root"}) ) {die "Missing Bam Root field\n"};
			my $outpath = $analysisDir."Tophat_".${$sampleHash}{$key}{'Bam Root'};
			my $aboutfile = "${outpath}/RunLog-Tophat_".${$sampleHash}{$key}{'Bam Root'}.".log";
			print "about file = $aboutfile\n";
			print "outpath : $outpath\n";
			my $tophatCommand = "tophat -o ${outpath} ${HAT_ARGS} ";
			$tophatCommand .= "--GTF=${TOPGTF} --num-threads=${NUM_THREADS} ";
			$tophatCommand .= "--mate-inner-dist=250 " if ($pairedEnd);
			$tophatCommand .= "--transcriptome-index=${transcriptomeIndex} ";
			$tophatCommand .= "${TOPINDEX} $read1 ";
			$tophatCommand .= "$read2 " if ($pairedEnd);
			$tophatCommand .= "2>&1 | tee --append $aboutfile";
			print "tophatCommand: $tophatCommand\n";

			runAndLog($tophatCommand) unless($dry);
		}
	}
	}
	return 0;
}


# returns path to fastqs for a given experiment
sub findFastQPath{
	my %args 	= 	@_;
	print "Finding directory for experiment : $args{exp}\n";
	my $exp 	= 	$args{exp};
	my $dataPath  = $args{path};
	my $sampleHash = $args{sampleHash};
	$exp =~ s/R//;
	# print "findFastQPath experiment\t".$exp."\n";
	my $sample = $exp."X1";
	print "sample: $sample\n";
	# print "$%{$sampleHash}{$sample}{Associated Files}\n";
	if( !defined (${$sampleHash}{$sample}{"Associated Files"})) {die "Missing Associated Files field\n"};
	my $files = ${$sampleHash}{$sample}{"Associated Files"};
	$files =~ s/"//g;
	$files =~ s/,//g;
	print "files : $files\n";
	my @files = split(" ", $files);
	print "searching for file : $files[0] in $dataPath\n";
	my $dir = `find $dataPath -name *$files[0]*`;
	if (!defined($dir)){
		die "Could not find the correct directory for experiment : ${exp}\n";
	}
	chomp($dir);
	print "found $dir\n";
	$dir = dirname( $dir );
	# $dir =~ s/\/.+$//;
	$dir .= "/";
	print "returning $dir\n";
	return $dir;
}


sub runAndLog{
	my $command = shift;
	my $time = localtime;
	print "$time\t$command\n";
	system($command);
}



return 1;