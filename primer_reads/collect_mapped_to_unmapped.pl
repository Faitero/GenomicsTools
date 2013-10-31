#!/usr/bin/perl
use warnings;
use strict;
use File::Basename;
use Data::Dumper;

our $rootDir;
our $fh;


sub main{
    #print caller()."\n";
    my %args = @_;
    my $print = $args{print};
    if (!defined($args{rootDir})){
        $rootDir = "/Data01/gnomex/Analysis/experiment/";
    }
    my @ah = `find $rootDir -name accepted_hits.bam`;

    for my $ah ( @ah ){
        chomp $ah;
        #      print $ah."\n";
        my $dir = get_dir_from_path($ah);
        $fh->{$dir}->{accepted_hits} = -s $ah;
        $fh->{$dir}->{unmapped} = -s $dir."/unmapped.bam";
        my $total = $fh->{$dir}->{accepted_hits} + $fh->{$dir}->{unmapped};
        $fh->{$dir}->{percent_mapped} = ($fh->{$dir}->{accepted_hits}) / $total  ;
    }
    #print Dumper $fh;

    if ($print){
        print "#sample_dir\tpercent_mapped\taccepted_hits_size\tunmapped_size\n";
        my @outar;
        for my $samp (keys  %$fh ){
            my $line = $fh->{$samp}->{percent_mapped}."\t"; ## percentmapped
            $line .= $samp."\t" ;                           ## sample dir
            $line .= $fh->{$samp}->{accepted_hits}."\t";    ## accpeted hits size
            $line .= $fh->{$samp}->{unmapped}."\n";         ## unmapped
            push(@outar, $line);
        }
        for my $line (sort @outar){
            chomp $line;
           my  @flds = split ("\t", $line );
            print $flds[1]."\t".$flds[0]."\t".$flds[2]."\t".$flds[3]."\n";
        }
    } else {
        return $fh;
    }

}

sub get_dir_from_path{
    my $fdir = shift;
   return dirname($fdir);
}


unless (caller){
    main(print =>1);
}

