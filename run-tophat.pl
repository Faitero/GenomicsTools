#!/usr/bin/perl -w

use strict;
use Getopt::Long;

my $species;
# my $compress;

my $usage_message = "run-tophat.pl -species=[m/h]";

GetOptions(
"species=s"		=>		\$species,
# "compress=s"	=>		\$compress,
	);



my $HAT_ARGS="--min-anchor=5 --segment-length=25 --no-coverage-search --segment-mismatches=2 --splice-mismatches=2 --microexon-search";
my $NUM_THREADS=6;
my $MM9_INDEX="/work/Common/Data/Bowtie_Indexes/mm9";
my $HG19_INDEX="/work/Common/Data/Bowtie_Indexes/hg19";
my $HG_GTF="/work/Common/Data/Annotation/human/Homo_sapiens.GRCh37.68.clean.gtf";
my $MM_GTF="/work/Common/Data/Annotation/mouse/mm9/Mus_musculus.NCBIM37.67.fixed.gtf";
my $TOPGTF;
my $TOPINDEX;

# if (!$compress){
# 	print $usage_message && die;
# }

if ($species == "m"){
	$TOPGTF = $MM_GTF;
	$TOPINDEX = $MM9_INDEX;
}
elsif ($species == "h"){
	$TOPGTF = $HG_GTF;
	$TOPINDEX = $HG19_INDEX;	
}
else {
	print $usage_message && die;
}

$INPUT_FASTQ_PATTERN="*.fq.bz2";
$PAIRED_MATE_INNER_DIST=250;

