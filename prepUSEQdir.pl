#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Cwd;

# This script takes a command line argument of a comma 
# separated list of patterns that can be used to partition
# sequence reads into conditions folders and then creates
# symlinks to the .bam files in these directories.

my $patterns;

GetOptions(

	"p=s"		=>		\$patterns,

	);


my @patterns = split( ',' , $patterns );

my $cwd = getcwd;
	

for my $pattern ( @patterns ){
	my $cond_dir = "USEQ_CONDITIONS/${pattern}/";
	system("mkdir -p $cond_dir");
	# print "$pattern\n";
	my @dirs = `ls -d Top*${pattern}*`;
	for my $dir (@dirs){
		chomp $dir;
		print "$pattern\t$dir\t$cond_dir\n";
		my $ln_command = "ln -s ${cwd}/${dir}/accepted_hits_processed.bam ${cond_dir}${dir}_accepted_hits_processed.bam";
		print "Executing\t$ln_command\n";
		system($ln_command);
	}
}




# for f in `ls -d Top*CTRL*`; do ln -s /work/b2b/exp_116R/116R/$f/accepted_hits.processed.bam USEQ_CONDITIONS/CTRL/$f.processed.bam ; done
# for f in `ls -d Top*HLHS*`; do ln -s /work/b2b/exp_116R/116R/$f/accepted_hits.processed.bam USEQ_CONDITIONS/HLHS/$f.prociessed.bam ; done
