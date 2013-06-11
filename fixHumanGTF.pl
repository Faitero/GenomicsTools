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
	if ( $line =~ s/^.*(CHR\w+)\W.*\t(.+$)/$1\t$2/ ){
		print "Changed:\n$_\nto:\n$1\t$2\n";

	} else {
		# print "noMatch throwing out $line\n";
	}
}