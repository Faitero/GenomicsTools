use strict;
use warnings;
use GnomexAlign;
use Data::Dumper;


my $dry = 0;
my $pl = new Pipeline();

# print Dumper $pl;

# $pl->findFileLocations;

$pl->getLanes;

$pl->buildFastQCCommandHash;

# print Dumper $pl->{commands};
 # <>;

# $pl->multiprocessQ(%{$pl->{commands}->{fastqc}}) unless $dry;

$pl->make_fastq_MCF_commands;

# print Dumper $pl->{commands};
# <>;

$pl->multiprocessQ(%{$pl->{commands}->{trim}});

$pl->buildAlignCommandHash;

# print Dumper $pl->{commands};
#  <>;

# $pl->multiprocessQ(%{$pl->{commands}->{align}}) unless $dry;

$pl->makeDeDupeCommandHash;

# print Dumper $pl->{commands};
#  <>;

# $pl->multiprocessQ( %{$pl->{commands}->{markDupes} } ) unless $dry;



$pl->makeRseqcCommandHash();

print Dumper $pl->{commands};
<>;

print Dumper $pl->{commands}->{rseqc};

$pl->multiprocessPerlQ(%{$pl->{commands}->{rseqc}});

# $pl->{mdb}->addPipelineSamples( sh => $pl->{samples} );









# $pl->logFieldReport;





# print Dumper($pl);