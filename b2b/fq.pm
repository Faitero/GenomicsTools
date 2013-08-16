package b2b::tools;
use strict;
use Getopt::Long;
use Data::Dumper;
use File::Basename;
use File::Glob;
use IO::Handle;

# this runs fastQC on you sample fastq and puts the results in the 'fastqc' directory
# of the experiment
sub fastQC{
	my %args = @_;
	my $analysisDir = $args{analysisDir};
	my $fastqcDirCommand = "mkdir -p ${analysisDir}fastqc";
	system($fastqcDirCommand);
	if(-d "${analysisDir}fastqc"){
		print $RUNLOG "${analysisDir}fastqc created successfully";
	} else {
		die "unable to create ${analysisDir}fastqc -- Error\n$!";
	}
	my $fastQCcommand = "fastqc --outdir=${analysisDir}fastqc ${fastqDir}*" unless ($nofastqc);
	system("$fastQCcommand") unless($dry || $nofastqc);
}

# this takes a experiment name and a reference to a sample hash and runs tophat on the relavent
# experiments
sub runTophat{
	my	 $HAT_ARGS="--min-anchor=5 --segment-length=25 --no-coverage-search --segment-mismatches=2 --splice-mismatches=2 --microexon-search --no-discordant --no-mixed ";
	my	 $NUM_THREADS=12;
	my	 $MM9_INDEX="/work/Common/Data/Bowtie_Indexes/mm9";
	my	 $HG19_INDEX="/work/Common/Data/Bowtie_Indexes/hg19";
	my	 $HG_GTF="/work/Common/Data/Annotation/human/Homo_sapiens.GRCh37.71.fixed.gtf";
	my	 $MM_GTF="/work/Common/Data/Annotation/mouse/mm9/Mus_musculus.NCBIM37.67.fixed.gtf";
	my 	 $GG_GTF="/work/Common/Data/Annotation/ercc/galGal4.72_ERCC92.gtf";
	my 	 $GG_INDEX="/work/Common/Data/Bowtie_Indexes/galGal4_ERCC92";
	my 	 $transcriptomeIndex = "/work/Common/Data/Annotation/bowtie_transcriptome_index";
	my	 $TOPGTF;
	my	 $TOPINDEX;
	my $PAIRED_MATE_INNER_DIST=250;
	my %args 	= 	@_;
	my $exp 	= 	$args{exp};
	my $sampleHash = $args{sampleHash};
	my $fastqDir  = $args{fastqDir};
	my $analysisDir = $args{analysisDir};
	my $species = $args{species};
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
		$TOPGTF = $GG_GTF;
		$TOPINDEX = $GG_INDEX;
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

			system($tophatCommand) unless($dry);
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