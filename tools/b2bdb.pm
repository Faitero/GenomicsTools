package tools::b2bdb;
use strict;
use warnings;
use IPC::System::Simple qw(system);
use IPC::System::Simple qw(capture);
use File::Basename;
use Data::Dumper;

our $dir = dirname(__FILE__);
my $b2bUN = $ENV{'b2bUN'};
our $b2bconn = $b2bUN.'@b2b.hci.utah.edu';
our $b2bssh = 'ssh '.$b2bconn;


## This sub connects to the b2b server via ssh and executes a query on the 
## gnomex database to get the Request #, the sample #, the associated files field
## the species and the type of experiment.
sub getSampMetaAssocFiles{
	my %args = @_;
	my $limit = $args{limit};
	# print $b2bssh."\n";
	my $qry = "select Request.number, Sample.number, PropertyEntry.valueString, organism , codeApplication from Request join Sample on Sample.idRequest = Request.idRequest join Property on Property.name='Associated Files' left join PropertyEntry on PropertyEntry.idSample=Sample.idSample and PropertyEntry.idProperty=Property.idProperty join Organism where Organism.idOrganism = Sample.idOrganism";
	if ($limit) {
		$qry .= " limit $limit";
	}
	system ("mkdir -p $dir/sql");
	open ( my $fh , '+>', "$dir/sql/AssocFiles.sql" ) || die $!;
	print $fh $qry;
	close $fh;
	my $scpComm = "rsync -e ssh $dir/sql/AssocFiles.sql $b2bconn:sql/AssocFiles.sql";
	print $scpComm."\n";
	system ($scpComm);
	my $sqlComm = " \"mysql -u ".$ENV{'b2bsqlUN'}." -p".$ENV{'b2bsqlPW'}." --database gnomex <sql/AssocFiles.sql\"";
	print $b2bssh.$sqlComm."\n";
	my $res = capture ($b2bssh.$sqlComm);
	my @rows = split("\n", $res);
	my $i = 0;
	my $sh = {};
	for my $row (@rows){
		# print $row."\n";
		my @flds = split ("\t", $row);
		# print Dumper(@flds);
		if ($i == 0){
			$i++;
			next;
		}
		# for my $fld (@flds){
		# 	print $_." field\n";
		# }
		$flds[2] =~ s/[,";]//g;
		my @files = split (' ', $flds[2]);
		if ($flds[2] ne '' && @files <= 2 && $flds[2] ne 'NULL' && !(lc($flds[2]) =~ m/.jpg$/) && !(lc($flds[2]) =~ m/.pptx$/)){
			## paired end sample
			if (@files == 2){
				$sh->{$flds[1]}->{files}->{read1}->{name} = $files[0];
				$sh->{$flds[1]}->{files}->{read2}->{name} = $files[1];
			} else {
				$sh->{$flds[1]}->{files}->{read1}->{name} = $files[0];
			}
			$sh->{$flds[1]}->{species} = $flds[3];
			$sh->{$flds[1]}->{type} = $flds[4];
		}
	}
	# print Dumper($sh);
	return $sh;
} 



1;