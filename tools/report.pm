package b2b::report;
use strict;
use Data::Dumper;
# use db;

sub commify {
  local $_ = shift;
  s{(?:(?<=^)|(?<=^-))(\d{4,})}
   {my $n = $1;
    $n=~s/(?<=.)(?=(?:.{3})+$)/,/g;
    $n;
   }e;
  return $_;
}

sub sortSamps{
	my ($samps) = @_;
	# print Dumper($samps);
	my @samps = @$samps;
	# print "sortSamps\t@samps\n";
	my @fixSamps;
	my @fixedSamps;
	# print @samps."\n";
	for my $samp ( @samps ){
		# print "beforeSub\t".$samp."\n";
		if ($samp =~ m/^(\d+X)(\d)$/){
			$samp = $1.'0'.$2;
		} 
		# print "afterSub\t".$samp."\n";
		push(@fixSamps, $samp);
	}
	@fixSamps = sort(@fixSamps);
	for my $samp ( @fixSamps ){
		# print "afterSort".$samp."\n";
		$samp =~ s/X0/X/;
		# print "afterSubSort".$samp."\n";
		push (@fixedSamps, $samp);
	}
	print @fixedSamps."\n";
	return @fixedSamps;
}

sub writeHeader{
	my %args = @_;
	my $fh = $args{fh};
	my $exp = $args{exp};
	my $curtime = $args{curtime};
	print $fh qq{
	<!DOCTYPE html>
	<html>
	<head>
		<style>
		table, th, td { border: 1px solid black; padding : 5px 10px 5px 10px;}
		table { border-collapse:collapse;}
		td, th {text-align : center; vertical-align:middle}
		th { background-color: #AAAAAA}
		.ref { color : #444444 ; font-size : .7em}
		.num { text-align: right}
		</style>
		<title>Bench to Bassinet: QC for experiment $exp; $curtime</title>
	</head>
	<body>
		<h1>Bench to Bassinet: QC summary for experiment $exp</h1>
		<p>$curtime</p>
	};
}

sub writeLegend {
	my %args= @_;
	my $exp = $args{exp};
	my $fh = $args{fh};
	my $sh = b2b::db::getSampleMap(exp => $exp);
	print $fh qq{
		<h2>Legend</h2>
		<h3>Experiment $exp: group to sample table</h3>
		<table>
			<tr>
				<th>Group</th> <th>SampleID</th>
			</tr>
	};
	for my $group (sort((keys(%$sh)))){
		# print Dumper($sh);
		print $group."\n";
		my @samps = keys( %{$sh->{$group}});
		for my $sample (sortSamps( \@samps )){
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

sub makeTableString{
	my %args = @_;
	my $sh = $args{sh};
	my $title = $args{title};
	my $ref = $args{ref};
	my $totfld = $args{totfld};
	my @flds = @{$sh->{fields}};
	# print "make Table String ". keys ( %{ $sh->{samples} } )."\n";
	my @samps = keys ( %{ $sh->{samples} } ) ;
	@samps = sortSamps(\@samps);
	my $html = "<h2>$title</h2><br>
		<table>
				<tr><th></th>";
	for my $samp (@samps){
		print $samp."\n";
		$html .= "<th>$samp</th>";
	}
	$html .= "</tr>";
	for my $fld (@flds){
		my $printfld = $fld;
		$printfld =~ s/_/ /g;
		$printfld =~ s/LT/ < /g;
		$printfld =~ s/GTE/ >= /g;
		$printfld =~ s/cal-PCR /cal\/PCR /g;
		$html .= "<tr><td>$printfld</td>";
		for my $samp ( @samps ){
			if ( $sh->{samples}->{$samp}->{$fld} > 1 ){
				$sh->{samples}->{$samp}->{$fld} = commify($sh->{samples}->{$samp}->{$fld});
			} elsif ($sh->{samples}->{$samp}->{$fld} >0 ) {
				$sh->{samples}->{$samp}->{$fld} = sprintf '%.6f',$sh->{samples}->{$samp}->{$fld};
			}
			$html .= "<td class=\"num\">".$sh->{samples}->{$samp}->{$fld}."</td>";
		}
		$html .= "</tr>";
	}	
	$html .= "</table>";
	$html .= "<p class=\"ref\">Data generated using $ref</p>";
	return $html;
}

sub writeBamStat{
	my %args = @_;
	my $exp = $args{exp};
	my $fh = $args{fh};
	my $sh = b2b::db::getBamStat(exp => $exp);
	my $html = makeTableString( title => "Bam Statistics", ref =>"RSeQC Bamstat", sh =>$sh );	
	print $fh $html;
}

sub writeMarkDup{
	my %args = @_;
	my $exp = $args{exp};
	my $fh = $args{fh};
	my $sh = b2b::db::getMarkDup(exp => $exp);
	my $html = makeTableString( title => "Duplication Statistics" ,  ref => "PicardTools MarkDuplicates", sh=> $sh);
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