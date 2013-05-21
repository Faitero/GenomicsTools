#!/usr/bin/perl -w


use strict;
# use warnings;
use Getopt::Long;

my $HAT_ARGS="--min-anchor=5 --segment-length=25 --no-coverage-search --segment-mismatches=2 --splice-mismatches=2 --microexon-search";
my $NUM_THREADS=6;

my $MM9_INDEX="/work/Common/Data/Bowtie_Indexes/mm9";
my $HG19_INDEX="/work/Common/Data/Bowtie_Indexes/hg19";
my $HG_GTF="/work/Common/Data/Annotation/human/Homo_sapiens.GRCh37.68.clean.gtf";
my $MM_GTF="/work/Common/Data/Annotation/mouse/mm9/Mus_musculus.NCBIM37.67.fixed.gtf";

my $INPUT_ARCHIVE_PATTERN="*.tar.gz";
my $INPUT_FASTQ_PATTERN="*.fq.bz2";
my $PAIRED_MATE_INNER_DIST=115;         # make sure to set this!!!

#for human
# TOPGTF=${HG_GTF}
# TOPINDEX=${HG19_INDEX}

#for mouse
my $TOPGTF=${MM_GTF};
my $TOPINDEX=${MM9_INDEX};


my @inputArchives = glob $INPUT_ARCHIVE_PATTERN;

print @inputArchives;

## makes the directories for the Tophat output
for my $archive (@inputArchives){
	if ( $archive =~ m/^(.+)_\d_sequence.tar.gz$/ ){
		my $folder = $1;
		print "$1\n";
		system ("mkdir -p Tophat-$1");
	}
}

my $root = ".";

opendir my $dh, $root
  or die "$0: opendir: $!";

my @dirs = grep {-d "$root/$_" && ! /^\.{1,2}$/} readdir($dh);

for my $dir (@dirs){
	# print "$dir\n";
	unzip_and_align($dir);

}




sub unzip_and_align{
	my $dir = shift;
	print $dir."\n";
	chdir($dir) or die "$!";
	my $filepattern;
	if ($dir =~ m/Tophat-(.+)$/){
		$filepattern = $1;
	}
	my @files = `ls ../$filepattern*tar.gz`;
	for my $file (@files){
		print "decompressing $file\n";
		system("tar -kxzpPf $file");
	}
	my $logfile = "${dir}.log";
	system("touch $logfile");
	my @fastq = `ls *sequence.txt`;
	my $read1 = $fastq[0];
	my $read2 = $fastq[1];
	chomp($read1);
	chomp($read2);
	print "aligning $read1 and $read2";
	system("unset BOWTIE_INDEXES");
	my $tophat_command = "tophat -o ./ --GTF=$TOPGTF $HAT_ARGS --num-threads=$NUM_THREADS --mate-inner-dist=$PAIRED_MATE_INNER_DIST --no-discordant --no-mixed $TOPINDEX $read1 $read2 2>&1 | tee --append $logfile";
	print($tophat_command."\n");
	system($tophat_command);
	# system('rm *sequence.txt');
	chdir('..') or die "$!";

}