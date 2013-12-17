use warnings;
use strict;
use GnomexAlign;
use tools::mongo;
use Test::More;
use Data::Dumper;

my $basedir = "/Data01/gnomex/Analysis/test/";
## 5R is a single end experiment, 
my @testExps = ("5R","91R");
## test that all the fastqc files are slurped correctly

# [laurentt@osono testdata]$ find /Data01/gnomex/Analysis/test/ -name fastqc_data.txt | sort
# /Data01/gnomex/Analysis/test/5R/5X1/fastqc/7940X13_81BNAABXX_8_sequence_fastqc/fastqc_data.txt
# /Data01/gnomex/Analysis/test/5R/5X2/fastqc/7940X14_81BNAABXX_8_sequence_fastqc/fastqc_data.txt
# /Data01/gnomex/Analysis/test/5R/5X3/fastqc/7940X15_81BNAABXX_8_sequence_fastqc/fastqc_data.txt
# /Data01/gnomex/Analysis/test/91R/91X1/fastqc/120607-FC437-L1-ATCACG--Bruneau--BK15-Tbx5Nkx2.5-DKO-d10--Mm--R1.fq_fastqc/fastqc_data.txt
# /Data01/gnomex/Analysis/test/91R/91X1/fastqc/120607-FC437-L1-ATCACG--Bruneau--BK15-Tbx5Nkx2.5-DKO-d10--Mm--R2.fq_fastqc/fastqc_data.txt
# /Data01/gnomex/Analysis/test/91R/91X2/fastqc/120607-FC437-L1-CGATGT--Bruneau--BK29-Tbx5Nkx2.5-DKO-d6--Mm--R1.fq_fastqc/fastqc_data.txt
# /Data01/gnomex/Analysis/test/91R/91X2/fastqc/120607-FC437-L1-CGATGT--Bruneau--BK29-Tbx5Nkx2.5-DKO-d6--Mm--R2.fq_fastqc/fastqc_data.txt
# /Data01/gnomex/Analysis/test/91R/91X3/fastqc/120607-FC437-L1-TTAGGC--Bruneau--BK30-Tbx5Nkx2.5-DKO-d6--Mm--R1.fq_fastqc/fastqc_data.txt
# /Data01/gnomex/Analysis/test/91R/91X3/fastqc/120607-FC437-L1-TTAGGC--Bruneau--BK30-Tbx5Nkx2.5-DKO-d6--Mm--R2.fq_fastqc/fastqc_data.txt

my $expected_paths = [
"/Data01/gnomex/Analysis/test/5R/5X1/fastqc/7940X13_81BNAABXX_8_sequence_fastqc/fastqc_data.txt",
"/Data01/gnomex/Analysis/test/5R/5X2/fastqc/7940X14_81BNAABXX_8_sequence_fastqc/fastqc_data.txt",
"/Data01/gnomex/Analysis/test/5R/5X3/fastqc/7940X15_81BNAABXX_8_sequence_fastqc/fastqc_data.txt",
"/Data01/gnomex/Analysis/test/91R/91X1/fastqc/120607-FC437-L1-ATCACG--Bruneau--BK15-Tbx5Nkx2.5-DKO-d10--Mm--R1.fq_fastqc/fastqc_data.txt",
"/Data01/gnomex/Analysis/test/91R/91X1/fastqc/120607-FC437-L1-ATCACG--Bruneau--BK15-Tbx5Nkx2.5-DKO-d10--Mm--R2.fq_fastqc/fastqc_data.txt",
"/Data01/gnomex/Analysis/test/91R/91X2/fastqc/120607-FC437-L1-CGATGT--Bruneau--BK29-Tbx5Nkx2.5-DKO-d6--Mm--R1.fq_fastqc/fastqc_data.txt",
"/Data01/gnomex/Analysis/test/91R/91X2/fastqc/120607-FC437-L1-CGATGT--Bruneau--BK29-Tbx5Nkx2.5-DKO-d6--Mm--R2.fq_fastqc/fastqc_data.txt",
"/Data01/gnomex/Analysis/test/91R/91X3/fastqc/120607-FC437-L1-TTAGGC--Bruneau--BK30-Tbx5Nkx2.5-DKO-d6--Mm--R1.fq_fastqc/fastqc_data.txt",
"/Data01/gnomex/Analysis/test/91R/91X3/fastqc/120607-FC437-L1-TTAGGC--Bruneau--BK30-Tbx5Nkx2.5-DKO-d6--Mm--R2.fq_fastqc/fastqc_data.txt"];


my $funname = "tools::mongo::get_fastqc_data";
my $fun = \&$funname;
#print $funname."\n";
#print &$fun."\n";
my $fqc_data = &$fun(basedir => $basedir);

my $i = 0;
#print Dumper $fqc_data;
for ( $i = 0 ; $i < scalar(@$fqc_data); $i++ ){

	ok(${$fqc_data}[$i] eq ${$expected_paths}[$i] , "$funname: $i test for collecting paths" );
	
}

my $path_to_samples = {
"/Data01/gnomex/Analysis/test/5R/5X1/fastqc/7940X13_81BNAABXX_8_sequence_fastqc/fastqc_data.txt" => "5X1",
"/Data01/gnomex/Analysis/test/5R/5X2/fastqc/7940X14_81BNAABXX_8_sequence_fastqc/fastqc_data.txt" => "5X2",
"/Data01/gnomex/Analysis/test/5R/5X3/fastqc/7940X15_81BNAABXX_8_sequence_fastqc/fastqc_data.txt" => "5X3",	
"/Data01/gnomex/Analysis/test/91R/91X1/fastqc/120607-FC437-L1-ATCACG--Bruneau--BK15-Tbx5Nkx2.5-DKO-d10--Mm--R1.fq_fastqc/fastqc_data.txt" => "91X1",
"/Data01/gnomex/Analysis/test/91R/91X1/fastqc/120607-FC437-L1-ATCACG--Bruneau--BK15-Tbx5Nkx2.5-DKO-d10--Mm--R2.fq_fastqc/fastqc_data.txt" => "91X1",
"/Data01/gnomex/Analysis/test/91R/91X2/fastqc/120607-FC437-L1-CGATGT--Bruneau--BK29-Tbx5Nkx2.5-DKO-d6--Mm--R1.fq_fastqc/fastqc_data.txt" => "91X2",
"/Data01/gnomex/Analysis/test/91R/91X2/fastqc/120607-FC437-L1-CGATGT--Bruneau--BK29-Tbx5Nkx2.5-DKO-d6--Mm--R2.fq_fastqc/fastqc_data.txt" => "91X2",
"/Data01/gnomex/Analysis/test/91R/91X3/fastqc/120607-FC437-L1-TTAGGC--Bruneau--BK30-Tbx5Nkx2.5-DKO-d6--Mm--R1.fq_fastqc/fastqc_data.txt" => "91X3",
"/Data01/gnomex/Analysis/test/91R/91X3/fastqc/120607-FC437-L1-TTAGGC--Bruneau--BK30-Tbx5Nkx2.5-DKO-d6--Mm--R2.fq_fastqc/fastqc_data.txt" => "91X3"};

$funname = "tools::mongo::get_sample_from_path";
$fun = \&$funname;

my $path_to_samp = &$fun(@$fqc_data);

#test for getting the sample from the fastq path
for my $path ( sort keys %$path_to_samp ) {
	ok($path_to_samp->{$path} eq $path_to_samples->{$path} , "$funname: $path maps to ".$path_to_samples->{$path} );
}

# test that the fascqc_path_to_file_regex works
$funname = "tools::mongo::fascqc_path_to_file_pattern";
$fun = \&$funname;
ok( &$fun( "" ) == 0 ,  "$funname: empty string returns 0 " );
my $testpath = "/Data01/gnomex/Analysis/test/5R/5X1/fastqc/7940X13_81BNAABXX_8_sequence_fastqc/fastqc_data.txt";
my $expected = "7940X13_81BNAABXX_8_sequence";
ok( &$fun( $testpath ) eq $expected ,  "$funname: $testpath string returns $expected" );
my $testpath = "/Data01/gnomex/Analysis/test/91R/91X2/fastqc/120607-FC437-L1-CGATGT--Bruneau--BK29-Tbx5Nkx2.5-DKO-d6--Mm--R2.fq_fastqc/fastqc_data.txt";
my $expected = "120607-FC437-L1-CGATGT--Bruneau--BK29-Tbx5Nkx2.5-DKO-d6--Mm--R2.fq";
ok( &$fun( $testpath ) eq $expected ,  "$funname: $testpath string returns $expected" );


## use the pipeline object to get the correct read from the samples
my $pipe = new Pipeline();
# print Dumper( $pipe->{samples}->{"5X1"});
# if (exists($pipe->{samples}->{"5X1"}->{files}->{fastqRead2}) ){
# 	print "exists\n";
# } else {
# 	print "doesn't exist\n";
# }
# <>;
$funname = "tools::mongo::get_read_from_path";
$fun = \&$funname;
my $srp = $pipe->$fun( $path_to_samp) ;

my $samp = "5X1";
my $read = "read1";
my $i = 0;
my $exp = "/Data01/gnomex/Analysis/test/5R/5X1/fastqc/7940X13_81BNAABXX_8_sequence_fastqc/fastqc_data.txt";
ok( ${$srp->{$samp}->{$read}}[$i] eq $exp , "$funname: $samp $read $i path $exp"  );

 $samp = "5X1";
 $read = "read2";
 $i = 0;
ok( !exists($srp->{$samp}->{$read}) , "$funname: $samp $read $i should be nonexistant"  );

 $samp = "91X1";
 $read = "read1";
 $i = 0;
 $exp = "/Data01/gnomex/Analysis/test/91R/91X1/fastqc/120607-FC437-L1-ATCACG--Bruneau--BK15-Tbx5Nkx2.5-DKO-d10--Mm--R1.fq_fastqc/fastqc_data.txt";
ok( ${$srp->{$samp}->{$read}}[$i] eq $exp , "$funname: $samp $read $i path $exp"  );

 $samp = "91X1";
 $read = "read2";
 $i = 0;
 $exp = "/Data01/gnomex/Analysis/test/91R/91X1/fastqc/120607-FC437-L1-ATCACG--Bruneau--BK15-Tbx5Nkx2.5-DKO-d10--Mm--R2.fq_fastqc/fastqc_data.txt";
ok( ${$srp->{$samp}->{$read}}[$i] eq $exp , "$funname: $samp $read $i path $exp"  );
# print Dumper($srp);
# print Dumper( $pipe->{samples}->{"5X1"});
# if (exists($pipe->{samples}->{"5X1"}->{files}->{fastqRead2} ) ){
# 	print "exists\n";
# } else {
# 	print "doesn't exist\n";
# }
# <>;
#### note need to test cases where there are more than 1 read1 and read2

$funname = "Pipeline::get_fastqc_data_paths";
$fun = \&$funname;
$pipe->$fun();
$read = "fastqRead1";
$samp = "91X1";

# print Dumper $pipe->{samples}->{$samp}->{files}->{fastqc};


ok(${$pipe->{samples}->{$samp}->{files}->{fastqc}->{$read}}[0] eq '/Data01/gnomex/Analysis/experiment/91R/91X1/fastqc/120607-FC437-L1-ATCACG--Bruneau--BK15-Tbx5Nkx2.5-DKO-d10--Mm--R1.fq_fastqc/fastqc_data.txt',
	"$funname: $samp $read fastqcPath element 0 should equal '/Data01/gnomex/Analysis/experiment/91R/91X1/fastqc/120607-FC437-L1-ATCACG--Bruneau--BK15-Tbx5Nkx2.5-DKO-d10--Mm--R1.fq_fastqc/fastqc_data.txt'");


$funname = "tools::mongo::process_fastqc_data";
$fun = \&$funname;
$pipe->$fun();


done_testing();