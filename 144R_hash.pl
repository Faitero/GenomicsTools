use strict;
use warnings;
use Data::Dumper;
use GnomexAlign;

my $sh = {

 "144X5" =>{
 	"type" => "MRNASEQ",
 	"species" => "chicken",
 	"files" => {
   'fastqRead1' => {
                                                                            'path' =>
 "/Data01/gnomex/ExperimentData/2013/144R/130308-FC589-L3-ACAGTG--Benoit--HH25-4--Gg--R1.fq.bz2"}, 
 'fastqRead2' => {
                                                                           'path' => 
 "/Data01/gnomex/ExperimentData/2013/144R/130308-FC589-L3-ACAGTG--Benoit--HH25-4--Gg--R2.fq.bz2"}}},
 "144X1" =>{
  	"type" => "MRNASEQ",
 	"species" => "chicken",
 "files" => {
   'fastqRead1' => {
                                                                            'path' =>
 "/Data01/gnomex/ExperimentData/2013/144R/130308-FC589-L3-ATCACG--Benoit--HH18-20-5--Gg--R1.fq.bz2"}, 
 'fastqRead2' => {
                                                                           'path' => 
 "/Data01/gnomex/ExperimentData/2013/144R/130308-FC589-L3-ATCACG--Benoit--HH18-20-5--Gg--R2.fq.bz2"}}},
 "144X2" =>{
  	"type" => "MRNASEQ",
 	"species" => "chicken",
 "files" => {
   'fastqRead1' => {
                                                                            'path' =>
 "/Data01/gnomex/ExperimentData/2013/144R/130308-FC589-L3-CGATGT--Benoit--HH19-1--Gg--R1.fq.bz2"}, 
 'fastqRead2' => {
                                                                           'path' => 
 "/Data01/gnomex/ExperimentData/2013/144R/130308-FC589-L3-CGATGT--Benoit--HH19-1--Gg--R2.fq.bz2"}}},
 "144X4" =>{
  	"type" => "MRNASEQ",
 	"species" => "chicken",
 "files" => {
   'fastqRead1' => {
                                                                            'path' =>
 "/Data01/gnomex/ExperimentData/2013/144R/130308-FC589-L3-TGACCA--Benoit--HH25-1--Gg--R1.fq.bz2"}, 
 'fastqRead2' => {
                                                                           'path' => 
 "/Data01/gnomex/ExperimentData/2013/144R/130308-FC589-L3-TGACCA--Benoit--HH25-1--Gg--R2.fq.bz2"}}},
 
 "144X3" =>{
  	"type" => "MRNASEQ",
 	"species" => "chicken",
 "files" => {
   'fastqRead1' => {
                                                                            'path' =>
 "/Data01/gnomex/ExperimentData/2013/144R/130308-FC589-L3-TTAGGC--Benoit--HH19-3--Gg--R1.fq.bz2"}, 
 'fastqRead2' => {
                                                                           'path' => 

 "/Data01/gnomex/ExperimentData/2013/144R/130308-FC589-L3-TTAGGC--Benoit--HH19-3--Gg--R2.fq.bz2"}}},
 "144X8" =>{
  	"type" => "MRNASEQ",
 	"species" => "chicken",
 "files" => {
   'fastqRead1' => {
                                                                            'path' =>
 "/Data01/gnomex/ExperimentData/2013/144R/130308-FC589-L6-ACTTGA--Benoit--HH30-5--Gg--R1.fq.bz2"}, 
 'fastqRead2' => {
                                                                           'path' => 
 "/Data01/gnomex/ExperimentData/2013/144R/130308-FC589-L6-ACTTGA--Benoit--HH30-5--Gg--R2.fq.bz2"}}},


 "144X7" =>{
  	"type" => "MRNASEQ",
 	"species" => "chicken",
 "files" => {
   'fastqRead1' => {
                                                                            'path' =>
 "/Data01/gnomex/ExperimentData/2013/144R/130308-FC589-L6-CAGATC--Benoit--HH30-4nn--Gg--R1.fq.bz2"}, 
 'fastqRead2' => {
                                                                           'path' => 
 "/Data01/gnomex/ExperimentData/2013/144R/130308-FC589-L6-CAGATC--Benoit--HH30-4nn--Gg--R2.fq.bz2"}}},

 "144X9" =>{
    "type" => "MRNASEQ",
    "species" => "chicken",
 "files" => {
   'fastqRead1' => {
                                                                            'path' =>
 "/Data01/gnomex/ExperimentData/2013/144R/130308-FC589-L6-GATCAG--Benoit--HH32-3--Gg--R1.fq.bz2"}, 
 'fastqRead2' => {
                                                                           'path' => 
 "/Data01/gnomex/ExperimentData/2013/144R/130308-FC589-L6-GATCAG--Benoit--HH32-3--Gg--R2.fq.bz2"}}},
 "144X6" =>{
  	"type" => "MRNASEQ",
 	"species" => "chicken",
 "files" => {
   'fastqRead1' => {
                                                                            'path' =>
 "/Data01/gnomex/ExperimentData/2013/144R/130308-FC589-L6-GCCAAT--Benoit--HH25-5--Gg--R1.fq.bz2" }, 
 'fastqRead2' => {
                                                                           'path' => 

 "/Data01/gnomex/ExperimentData/2013/144R/130308-FC589-L6-GCCAAT--Benoit--HH25-5--Gg--R2.fq.bz2"}}}

};



my $pl = new Pipeline( noload => 1 );

$pl->{samples}  = $sh ;
print Dumper $pl;
$pl->getLanes;
$pl->buildFastQCCommandHash;
$pl->multiprocessQ(%{$pl->{commands}->{fastqc}}) ;
$pl->make_fastq_MCF_commands;
$pl->multiprocessQ(%{$pl->{commands}->{trim}});
$pl->buildAlignCommandHash;
$pl->multiprocessQ(%{$pl->{commands}->{align}});

print Dumper $pl->{commands};








