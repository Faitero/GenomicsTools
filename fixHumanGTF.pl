#!/usr/bin/perl -w

use strict;

my $gtfin = "Homo_sapiens.GRCh37.71.gtf";
my $gtfout = "Homo_sapiens.GRCh37.71.fixed.gtf";

open IN, '<', $gtfin;

open OUT, '>', $gtfout;
my $line;

while (<IN>) {
	$line = $_;
	# print $line."\n";
	if ( $line =~ s/^(\d+)(\t.+$)/chr$1$2/ ){
		print "Changed:\n$_\nto:\n$line\n";
		print OUT $line;
	} elsif ( $line =~ s/^([YX])(\t.+$)/chr$1$2/ ){
		print "Changed:\n$_\nto:\n$line\n";
		print OUT $line;
	} elsif ( $line =~ s/^(MT)(\t.+$)/chr$1$2/ ){
		print "Changed:\n$_\nto:\n$line\n";
		print OUT $line;
	} 	
}
close IN;
close OUT;