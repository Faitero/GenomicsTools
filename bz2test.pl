my $cmd = " fastq-mcf -S -o /Data01/gnomex/Analysis/experiment/171R/171X1/trimmed/LIB006971_CHS00011793_CTTGTA_L007_R1.fastq.bz2._MCF.gz /work/Common/Data/contaminants/contaminants.fa <(bunzip2 -c /Data01/gnomex/ExperimentData/2013/171R/upload_staging/LIB006971_CHS00011793_CTTGTA_L007_R1.fastq.bz2)";

system_bash ($cmd);


sub system_bash {
  my @args = ( "bash", "-c", shift );
  system(@args);
}
