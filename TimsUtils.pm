use strict;

package TimsUtils;


sub getCurrentTime{
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	my $curTime = $year + 1900 ."-".( $mon + 1)  ."-". $mday . "-" . $hour . "-" . $min . "-" . $sec ."-";
	return $curTime;
}




1;