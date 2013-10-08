use File::Copy;
use strict;
use warnings;

my $basedir = "/home/laurentt/exp/";

my @readCovs = glob($basedir."77R/Top*/*marked_dup*/*geneBodyCoverage.txt");

my $out = $basedir."geneBodyTest/geneBodyCovAll.txt";

open OUT , ">", $out;

for my $read (@readCovs){
	print $read."\n";
	my $outfile = $read;
	$outfile =~ s/\//_/g;
	$outfile =~ m/(\d+X\d+)/;
	my $id = $1."\n";
	open IN , "<" , $read;
	<IN>;
	<IN>;
	<IN>;
	while (my $line = <IN> ){
		chomp $line;
		print OUT $line."\t".$id;
	}
}

close OUT;