use strict;
use warnings;
use GnomexAlign;
use Data::Dumper;


my $dry = 0;
my $pl = new Pipeline();

print Dumper $pl;

# $pl->findFileLocations;

$pl->getLanes;

$pl->buildFastQCCommandHash;



$pl->make_fastq_MCF_commands;


$pl->buildAlignCommandHash;

print Dumper $pl->{commands};
<>;



# $pl->{mdb}->addPipelineSamples( sh => $pl->{samples} );

$pl->multiprocessQ(%{$pl->{commands}->{fastqc}}) unless $dry;

$pl->multiprocessQ(%{$pl->{commands}->{trim}});

$pl->multiprocessQ(%{$pl->{commands}->{align}}) unless $dry;
$pl->logFieldReport;





# print Dumper($pl);