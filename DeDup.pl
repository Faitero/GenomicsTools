package DeDup;

use strict;
use warnings;
use IPC::System::Simple qw(system);
use IPC::System::Simple qw(capture);
use File::Basename;
use Proc::Queue size => 4, debug => 1;
use POSIX ":sys_wait_h";

our $debug = 1;

my $deDup = DeDup->new;

$deDup->getNonDedupedBams;

$deDup->makeDedupCmdHash;

$deDup->multiprocessQ(%{$deDup->{dedupCmds}});


sub makeDedupCmdHash{
	my $self = shift;
	my $ULIMIT_RESULT = 1024; ## result of running the shell command ulimit -n. Since this is a shell built-in, it can, for some reason, not be run like a real command, so backticks don't work. `ulimit -n`; chomp($ULIMIT_RESULT);
	my $maxFileHandles = int(0.8 * $ULIMIT_RESULT/4); ## This is for the MAX_FILE_HANDLES_FOR_READ_ENDS_MAP parameter for MarkDuplicates: From the Picard docs: "Maximum number of file handles to keep open when spilling read ends to disk. Set this number a little lower than the per-process maximum number of file that may be open. This number can be found by executing the 'ulimit -n' command on a Unix system. Default value: 8000."
	# my $MARKDUPLICATES_PATH = `which MarkDuplicates.jar`; chomp($MARKDUPLICATES_PATH); ## Tries to find "MarkDuplicates.jar" from the Picard suite in the user's $PATH
	my $GIGABYTES_FOR_PICARD = 6;

	# my $MARKDUPLICATES_PATH = "/work/Apps/Bio/picard/picard-tools-1.90/MarkDuplicates.jar"; 
	my $MARKDUPLICATES_PATH = `which MarkDuplicates.jar`; 
	chomp($MARKDUPLICATES_PATH);
	# my $SORTSAM_PATH = "/bioinformatics/installations/picard/picard-tools-1.89/SortSam.jar";
	my $SORTSAM_PATH = `which SortSam.jar`;
	chomp($SORTSAM_PATH);
	
	for my $inFile ( sort keys %{$self->{toDeDup}}){
		my $outFile = $inFile;
		$outfile =~ s/.bam//;
		$outfile .= "_marked_dup.bam";
		my $dedupExtraMetricsFile = $outFile;
		$dedupExtraMetricsFile .= "_mark_dup_metrics.txt";
		my $dedupCmd = (qq{java -Xmx${GIGABYTES_FOR_PICARD}g -jar ${MARKDUPLICATES_PATH} }
			    . qq{ INPUT=$inFile } ## Picard SortSam.jar accepts both SAM and BAM files as input!
			    . qq{ REMOVE_DUPLICATES=FALSE }
			    . qq{ CREATE_INDEX=TRUE }
			    #. ((!$shouldSort) ? qq{ ASSUME_SORTED=TRUE } : qq { }) ## <-- if we use --nosort, then ASSUME they are sorted no matter what!
			    . qq{ MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=$maxFileHandles }
			    . qq{ OPTICAL_DUPLICATE_PIXEL_DISTANCE=10 } ## <-- 100 is default but our libraries have single integer coordinate values
			    . qq{ METRICS_FILE=$dedupExtraMetricsFile }
			    . qq{ OUTPUT=$outFile }
		);	
	
		$self->{dedupCmds}->{$dedupCmd} = 1;	
	}
}


## this creates a multiprocessed queue taking an array of system calls
sub multiprocessQ{
	my $self = shift;
	my %args = @_;
	for my $call (sort keys %args){
		my $f = fork;
		if (defined ($f) and $f==0 ){
			
			system("nice $call");
			exit(0);
		}
		1 while waitpid(-1, WNOHANG) > 0;
	}
}

sub new{
	my $class = shift;
	my $dir = shift;
	my $self = {};
	$dir = "/Data01/gnomex/Analysis/experiment" unless defined $dir;
	$self->{dir} = $dir;
	bless $self, $class;
	return $self;
}

sub getNonDedupedBams{
	my $self = shift;
	my $cmd = "find ".$self->{dir}." -type f -name \"accepted_hits.bam\"";
	my @accHits = capture($cmd);
	
	for my $hit ( @accHits ){
		if (-s $hit < 104857600 ){
			next;
		}
		chomp $hit;
		my $dedupName = dirname($hit)."/accepted_hits_marked_dup.bam";
		chomp $dedupName;
		unless ( -e $dedupName ){
			print "adding $hit to toDeDup\n";
			$self->{toDeDup}->{$hit} = 1;
		} else {
			print "adding $dedupName to Deduped\n";
			$self->{DeDuped}->{$dedupName} = 1;
		}
	}
}