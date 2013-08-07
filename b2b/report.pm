package b2b::report;
use strict;
use Data::Dumper;
# use db;

sub writeHeader{
	my %args = @_;
	my $fh = $args{fh};
	my $exp = $args{exp};
	my $curtime = $args{curtime};
	print $fh qq{
	<!DOCTYPE html>
	<html>
	<head>
		<title>Bench to Bassinet: QC for experiment $exp; $curtime</title>
	</head>
	<body>
		<h1>Bench to Bassinet: QC summary for experiment $exp</h1>
		<p>$curtime</p>
	};
}

sub makeLegend {
	my %args= @_;
	my $exp = $args{exp};
	my $fh = $args{fh};
	my $sh = b2b::db::getSampleMap(exp => $exp);
	print $fh qq{
		<h3>Experiment $exp: group to sample table</h3>
		<table>
			<tr>
				<th>Group</th> <th>SampleID</th>
			</tr>
	};
	for my $group (sort(keys(%$sh))){
		# print Dumper($sh);
		print $group."\n";
		for my $sample (sort(keys( %{$sh->{$group}} ) ) ){
			print $sample."\n";
			print $fh qq{
				<tr>
					<td>$group</td><td>$sample</td>
				</tr>
			};
		}
	}
	print $fh qq{
		</table>
	};
}

sub writeBamStat{
	my %args = @_;
	my $exp = $args{exp};
	my $fh = $args{fh};
	my $sh = b2b::db::getBamStat(exp => $exp);
	print Dumper($sh);
	my @flds = @{$sh->{fields}};
	my @samps = sort( keys ( %{ $sh->{samples} } ) );
	my $html = "<h2>Bam Stat Results </h2><br>";
	
	$html .= "<table>
				<tr><th></th>";
	for my $samp (@samps){
		print $samp."\n";
		$html .= "<th>$samp</th>";
	}
	$html .= "</tr>";
	print "fields\t".@flds."\n";
	for my $fld (@flds){
		$html .= "<tr><td>$fld</td>";
		for my $samp ( @samps ){
			$html .= "<td>".$sh->{samples}->{$samp}->{$fld}."</td>";
		}
		$html .= "</tr>";
	}
	$html .= "</table>";
	$html .= "<p>Created with RSeQC Sequence data analysis package</p>"
	print $fh $html;
}

sub writeFooter{
	my %args = @_;
	my $fh = $args{fh};
	print $fh qq{
		</body>
		</html>
	};
}

1;