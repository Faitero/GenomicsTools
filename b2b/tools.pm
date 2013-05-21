#!/usr/bin/perl -w

use strict;
package b2b;

# parse sample sheet takes a path to a tsv sample sheet and returns a hash with the following structure:
# SampleID-><each property>-><each property's value>
sub parseSampleSheet{
	my $path = shift;
	open( IN, "<", $path ) or die "cannot open $path: $!\n";

	my $outHash = {};
	my $line = <IN>;
	chomp($line);
	my @header = split(/\t/, $line);

	my $i;
	my $sampleNumberIndex;
	for ( $i = 0 ; $i < @header ; $i++ ){
		if ($header[$i] == "Sample Number"){
			$sampleNumberIndex = $i;
		}
	}

	while ( $line = <IN> ){
		my @row = split(/\t/, $line);
		chomp($row[-1]);
		for ($i = 0 ; $i < @row; $i++) {
			if ($row[$i] != ""){
				print "$row[$i]\n";
				# $outHash->{$row[$sampleNumberIndex]}->{$header[$i]} = $row[$i];
			}
		}
	}
	return $outHash;
}

# this takes a experiment name and a reference to a sample hash and runs tophat on the relavent
# experiments
sub runTophat{
	my %args  	=	@_;
	my $sampleHash = $args{sampleHash};
	my $sample 	   = $args{sample};
	if ($sample =~ m/\d+R/ ){
		$sample =~ s/(\d+)R/$1/ ;
	}


}





return 1;