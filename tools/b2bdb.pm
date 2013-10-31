package tools::b2bdb;
use strict;
use warnings;
use IPC::System::Simple qw(system);
use IPC::System::Simple qw(capture);
use File::Basename;
use Data::Dumper;

our $debug = 1;
our $dir = dirname(__FILE__);
my $b2bUN = $ENV{'b2bUN'};
our $b2bconn = $b2bUN.'@b2b.hci.utah.edu';
our $b2bssh = 'ssh '.$b2bconn;
our $qry = <<QRY; 
select distinct Sample.number, organism, concat(propertyValue, year(Request.createDate), '/', filename), codeApplication, codeSampleFileType
    from Request join Sample on (Sample.idRequest = Request.idRequest)
    join Organism on Organism.idOrganism = Sample.idOrganism
    join SampleExperimentFile on SampleExperimentFile.idSample = Sample.idSample
    join ExperimentFile on ExperimentFile.idExperimentFile = SampleExperimentFile.idExperimentFile
    join PropertyDictionary on PropertyDictionary.propertyName='experiment_directory'
    where SampleExperimentFile.codeSampleFileType = 'fastqRead1' or SampleExperimentFile.codeSampleFileType = 'fastqRead2';
QRY





sub queryB2BSampleFiles{
	my %args = @_;
	print "Script version of B2BDB\n";
	my $limit = $args{limit};
	# print $b2bssh."\n";
	system ("mkdir -p $dir/sql");
	my $sqlPath = "sql/getReads.sql";
	open ( my $fh , '+>', "$dir/$sqlPath" ) || die $!;
    print $fh $qry;
    close $fh;
    my $scpComm = "rsync -e ssh $dir/$sqlPath $b2bconn:$sqlPath";
    print $scpComm."\n";
	system ($scpComm);
	my $sqlComm = " \"mysql -u ".$ENV{'b2bsqlUN'}." -p".$ENV{'b2bsqlPW'}." --database gnomex <$sqlPath\"";
	print $b2bssh.$sqlComm."\n";
	system ("mkdir -p $dir/sql/out");
	open SQLOUT, ">", "$dir/sql/out/res.out";
	my $res = capture ($b2bssh.$sqlComm);
	print SQLOUT $res;
	my @rows = split("\n", $res);
	my $i = 0;
	my $sh = {};
	for my $row (@rows){
		print $row."\n" if $debug;
		chomp $row;
		next if ($i++ == 0 ); ## prevent the header row from going into table
		my @flds = split("\t", $row);
		for my $fld (@flds){
			chomp $fld;
		}
		$sh->{$flds[0]}->{files}->{$flds[4]} .= $flds[2]."," ; ##if (defined $sh->{$flds[0]}->{files}->{$flds[4]} && !($sh->{$flds[0]}->{files}->{$flds[4]} =~ $flds[2]) );

		# push (@{$sh->{$flds[0]}->{files}->{$flds[4]}} , $flds[2] );
		$sh->{$flds[0]}->{species} = $flds[1];
		$sh->{$flds[0]}->{type} = $flds[3];
	}

	## there are duplicate rows
	print Dumper ($sh) if $debug;
	for my $samp (keys(%{$sh})){
		for my $file ( keys %{$sh->{$samp}->{files}} ){
			chop $sh->{$samp}->{files}->{$file};
		}
	}
	print Dumper ($sh) if $debug;
	print Dumper ($sh->{"37X5"}->{files}->{fastqRead2}) if $debug;
	for my $file ($sh->{"37X5"}->{files}->{fastqRead2}){
		print $file."\n";
	}
	return $sh;
}


unless ( caller ) {
	queryB2BSampleFiles;

}





## This sub connects to the b2b server via ssh and executes a query on the 
## gnomex database to get the Request #, the sample #, the associated files field
## the species and the type of experiment.
sub getSampMetaAssocFiles{
	my %args = @_;
	my $limit = $args{limit};
	# print $b2bssh."\n";
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
	system ("mkdir -p $dir/sql/out");
	open SQLOUT, ">", "$dir/sql/out/res.out";
	my $res = capture ($b2bssh.$sqlComm);
	print SQLOUT $res;
	my @rows = split("\n", $res);
	my $i = 0;
	my $sh = {};
	for my $row (@rows){
		print $row."\n" if $debug;
		chomp $row;
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
	print Dumper($sh) if $debug;
	return $sh;
} 



1;