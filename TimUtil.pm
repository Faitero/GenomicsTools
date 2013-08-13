use strict;
## miscellanious utilities

package TimUtil;

##	generates a current time of the form  "YYYY-MM-DD-HH-MM-SS-"
sub getCurrentTime{
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	$mon = $mon + 1;
	$mon = prependZero($mon);
	$mday = prependZero($mday);
	$hour = prependZero($hour);
	$min = prependZero($min);
	$sec = prependZero($sec);

	my $curTime = $year + 1900 ."-".( $mon)  ."-". $mday . "-" . $hour . "-" . $min . "-" . $sec ."-";
	return $curTime;
}

sub prependZero{
	my $num = shift;
	if (length ($num ) == 1 ){
		$num = "0" . $num;
	}
	return $num;
}


1;