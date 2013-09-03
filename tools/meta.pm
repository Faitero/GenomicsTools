package tools::meta;
use strict;
use Getopt::Long;
use Data::Dumper;
use File::Basename;
use File::Glob;
use IO::Handle;

## returns species of experiment
sub getSpecies{
	my %args = @_; 
	my $sampleHash = $args{sampleHash};
	my @samps = keys(%{ $sampleHash });
	for my $samp (@samps){
		if (exists($sampleHash->{$samp}->{'Organism'})){
			print "getSpecies\t" . $sampleHash->{$samp}->{'Organism'};
			return $sampleHash->{$samp}->{'Organism'};
		}
	}
	print "'Organism' field not defined in the sample sheet.\nEnter species (Mouse, Chicken, Human)\n>";
	
	my $return = <STDIN>;
	chomp $return;
	return $return;
}

# takes the sampleHash and the experiment and returns a
# has reference containing only sample from the chosen 
# experiment.
sub makeExpSampleHash{
	my %args = @_;
	my $samplesheet = $args{samplesheet};
	my $sampleHash = parseSampleSheet($samplesheet);
	my $expSampleHash = {};
	my $exp = $args{exp};
	$exp =~ s/R//;
	## check if experiment exists in this sample sheet
	if ( !defined ( $sampleHash )){
		die "exp: $exp not found in the sample sheet\n";
	}

	print "Building hash of experiment $exp samples\n";
	for my $sample ( keys( %$sampleHash ) ){
		if( $sample =~ m/^${exp}X\d+/ ){
			$expSampleHash->{$sample} = ${$sampleHash}{$sample};
		}
	}
	if (keys( %$expSampleHash ) == 0 ){
		die "Invalid experiment! try again\n";
	}
	print "Found samples:\n";
	for  my $key ( keys( %$expSampleHash ) ){
		if ( !defined( $expSampleHash->{$key}{"Bam Root"} ) ){
			die "Bam Root not defined for sample $key, please fix";
		}
		print "$key\n";
	}
	print "\n";
	return $expSampleHash;
}


# parse sample sheet takes a path to a tsv sample sheet and returns a hash with the following structure:
# SampleID-><each property>-><each property's value>
sub parseSampleSheet{
	my $path = shift;
	system("mac2unix $path");
	open( IN, "<", $path ) or die "cannot open $path: $!\n";
	my $outHash = {};
	my $line = <IN>;
	# print $line."\n";
	chomp($line);
	my @header = split(/\t/, $line);
	my $headerNum = @header;
	print "headerNum : $headerNum\n";
	# my $res = <stdin>;
	my $i;
	my $sampleNumberIndex;
	# print "test";
	# print __LINE__, "\n";
	for ( $i = 0 ; $i < @header ; $i++ ){
		# if ($debug) {print "$header[$i]\n";}
		if ($header[$i] eq "Sample Number" || $header[$i] eq "ID"){
			$sampleNumberIndex = $i;
		}
	}
	# $line = <IN>;
	# print $line."\n";
	# print "trace1\n";
	while ( $line = <IN> ){
		# my $res = <stdin>;
		print "line : $line\n";
		my @row = split(/\t/, $line);
		for ($i = 0 ; $i < @row; $i++) {
			if ($row[$i] ne "" ){
				$row[$i] =~ s/"//g;
				chomp($row[$i]);
				print "$row[$i]\t$i\n";
				print "adding $row[$sampleNumberIndex]->$header[$i] = $row[$i]\n";
				$outHash->{$row[$sampleNumberIndex]}{$header[$i]} = $row[$i];
			}
		}
		# readline;
	}
	close IN;
	return $outHash;
}
1;