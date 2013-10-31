package Pipeline;

use File::Basename;
use strict;
use warnings;
use tools::b2bdb;
use Data::Dumper;
use tools::fq;
use IPC::System::Simple qw(system);
use IPC::System::Simple qw(capture);
use Getopt::Long;
use tools::mongo;
use Proc::Queue size => 20, debug => 1;
use POSIX ":sys_wait_h";
use Switch;
use Scalar::Util qw(blessed);

our $dry = 0;
our $debug = 1;

my $pl = new Pipeline();

$pl->findFileLocations;

$pl->getLanes;
# $pl->buildFastQCCommandHash;
$pl->buildAlignCommandHash;
$pl->{mdb}->addPipelineSamples( sh => $pl->{samples} );

# $pl->multiprocessQ(%{$pl->{fastqc}->{commands}}) unless $dry;
$pl->multiprocessQ(%{$pl->{align}->{commands}}) unless $dry;
$pl->logFieldReport;


print Dumper($pl);
# print Dumper($pl->{fastqc});

# print Dumper($pl->{align});

## this creates a new blessed pipeline object and populates it with data from Gnomex database
sub new{
	my $class = shift;
	my $self = {
		time => getCurrentTime(),
		samples => tools::b2bdb::getSampMetaAssocFiles,
		analysisDir=> "/Data01/gnomex/Analysis/experiment",
		logDir => "/Data01/gnomex/Analysis/log"
	};
	system("mkdir -p ".$self->{logDir});
	$self->{errLog} =  $self->{logDir}."/error.log";
	$self->{runLog} = $self->{logDir}."/run.log";
	bless $self, $class;
	$self->runLog($self->{time}."\tNEW RUN #######################");
	$self->errLog($self->{time}."\tNEW RUN #######################");
	$self->{db}->{chicken}->{version} = "Galgal4.72";
	$self->{db}->{mouse}->{version} = "mm9_NCBI37.67";
	$self->{db}->{human}->{version} = "hg19_GRCh37.71";
	$self->{db}->{zebrafish}->{version} = "Zv9.73";
	$self->{db}->{chicken}->{gtf} = "/work/Common/Data/Annotation/chicken/Gallus_gallus.Galgal4.72.gtf";
	$self->{db}->{mouse}->{gtf} = "/work/Common/Data/Annotation/mouse/mm9/Mus_musculus.NCBIM37.67.fixed.gtf";
	$self->{db}->{human}->{gtf} = "/work/Common/Data/Annotation/human/Homo_sapiens.GRCh37.71.fixed.gtf";
	$self->{db}->{zebrafish}->{gtf} = "/work/Common/Data/Annotation/zebrafish/Danio_rerio.Zv9.73.gtf";
	$self->{db}->{chicken}->{index} = "/work/Common/Data/Bowtie_Indexes/galGal4.72_ERCC92";
	$self->{db}->{mouse}->{index} = "/work/Common/Data/Bowtie_Indexes/mm9";
	$self->{db}->{human}->{index} = "/work/Common/Data/Bowtie_Indexes/hg19";
	$self->{db}->{zebrafish}->{index} = "/work/Common/Data/Bowtie_Indexes/danRer_7V9.73";
	$self->{db}->{chicken}->{refgene} = "/work/Common/Data/Annotation/chicken/Gallus_gallus.Galgal4.72.bed";
	$self->{db}->{mouse}->{refgene} = "/work/Common/Data/Annotation/mouse/mm9/Mus_musculus.NCBIM37.67.fixed.bed";
	$self->{db}->{human}->{refgene} = "/work/Common/Data/Annotation/human/Homo_sapiens.GRCh37.71.fixed.bed";
	$self->{db}->{zebrafish}->{refgene} = "/work/Common/Data/Annotation/zebrafish/Danio_rerio.Zv9.73.bed";
	$self->{mdb} = tools::mongo->new ;
	print "self->{bdb} = ".blessed($self->{mdb}->{db})."\n";
	return $self;
}

## makes a report of the log field
sub logFieldReport{
	my $self = shift;
	open LOG, ">", $self->{logDir}."/stepLog.txt";
	my @steps = keys ( %{  $self->{log} });
	my @samps = keys ( %{ $self->{log}->{$steps[0]} });
	my $header = "sample";
	for my $step (@steps){
		$header .= "\t$step";
	}
	$header .= "\n";
	print LOG $header;
	for my $samp (@samps){
		my $row = $samp;
		for my $step (@steps){
			$row .= "\t";
			$row .= $self->{$step}->{$samp}->{success}.";";
			$row .= $self->{$step}->{$samp}->{err};
		}
		print LOG $row;
	}
}


## This write the argument to the runlog
sub runLog{
	my $self = shift;
	my $msg = shift;
	print $msg."\n";
	system ("echo \"$msg\" >> ".$self->{runLog});
}

## This write the argument to the errlog
sub errLog{
	my $self = shift;
	my $msg = shift;
	print $msg."\n";
	system ("echo \"$msg\" >> ".$self->{errLog});
}

## This adds msg to err log field
sub logErr{
	my $self = shift;
	my $msg = shift;
	$self->logWrite($msg, "err");
}

## This adds msg to succ log field
sub logSuccess{
	my $self = shift;
	my $msg = shift;
	$self->logWrite($msg, "success");
}

## This adds write to succ
sub logWrite{
	my $self = shift;
	my $msg = shift;
	my $fld = shift;
	my $samp = $self->{curSamp};
	my $step = $self->{curStep};
	$self->{log}->{$step}->{$samp}->{$fld} = $msg;
}

sub logMsg{
	my $msg = shift;
	my $log = shift;
	system("echo \"$msg\" >>$log");
}


## this adds FS locations to each file field
sub findFileLocations {
	my $self = shift;
	$self->{curStep} = "findFileLocations";

	my $findFilesLog = $self->{logDir}."/findFiles.log";
	logMsg($self->{time}, $findFilesLog);
	for my $samp ( sort keys %{$self->{samples}} ) {
		$self->{curSamp} = $samp;
		my $succMsg;
		my $errMsg;
		for my $read ( sort keys %{$self->{samples}->{$samp}->{files}} ){
			my $logmsg;
			my $fname = $self->{samples}->{$samp}->{files}->{$read}->{name};
			my $file = $self->{mdb}->checkSampFileDB(
					samp => $samp,
					read => $read
				);
			## check if in the database
			if (defined($file) && $file ne "0"){
				my $msg = "$file found in DB; ";
				print "$msg\n" if ($debug);
				my $succMsg .= $msg;
				$self->{samples}->{$samp}->{files}->{$read}->{name} = $file;
				$logmsg = "$samp\t$read\t$msg";
				logMsg($logmsg, $findFilesLog);
			} else{
				## check if path can be found on fs
				my $exp = $samp;
				$exp =~ s/X.*$//;
				$exp .= "R";
				if (my $path = findFastQPathSample(file => $fname, exp => $exp) ){
					my $msg = "$fname found at $path";
					print "$msg\n" if ($debug);
					$succMsg .= "$msg; ";
					# $path =~ s/ /\\ /g;
					if ($self->isFastQ($path) ){
						$msg = "$path is a fastq file";
						print "$msg\n";
						$succMsg .= "$msg; ";
						$self->{samples}->{$samp}->{files}->{$read}->{name} = $path;
						$logmsg = "$samp\t$read\t$msg";
						logMsg($logmsg, $findFilesLog);
						
					} else {
						$msg = "$samp : $path is not a fastq file";
						print "$msg\n";
						# clear out read name information
						$self->{samples}->{$samp}->{files}->{$read}->{name} = 0;
						$self->{samples}->{$samp}->{files}->{$read}->{error} = $msg;
						$self->errLog ($msg);
						$errMsg .= "$msg; ";
						$logmsg = "$samp\t$read\t$msg";
						logMsg($logmsg, $findFilesLog);
					}
				} else {
					## could not be found on fs
					my $msg = "$samp : $fname could not be found on the filesystem\n";
					$self->{samples}->{$samp}->{files}->{$read}->{name} = 0;
					$self->{samples}->{$samp}->{files}->{$read}->{error} = $msg;
					$self->errLog("Error: findFastQPathSample: $msg");
					print $msg."\n" if ($debug); 
					$errMsg .= "$msg. ";
					$logmsg = "$samp\t$read\t$msg";
					logMsg($logmsg, $findFilesLog);
				}
			}
		}
		if ( defined ( $succMsg ) ){
			$self->logSuccess($succMsg);
		}
		if ( defined ( $errMsg ) ){
			$self->logErr($errMsg);
		}
	}
}

## this checks if a purported fastq file is a fastq file and prints to the error file if it detects an
## error or if the file doesn't appear to be a fastq
sub isFastQ{
	my $self = shift;
	my $file = shift;
	chomp $file;
	if ( -e $file){
		my $IN = $self->openFastQ($file);
		if (!($IN)){
			$self->errLog("Error: isFastQ: $file could not be opened");
			return;
		}
		#  check for error messages while 
		my $line1;
		my $line2;
		my $line5;
			$line1 = <$IN>;
			$line2 = <$IN>;
			<$IN>;
			<$IN>;
			$line5 = <$IN>;
		#  check for errors reading file
		if( $@ ){
			$self->errLog("Error: isFastQ: Could not read from $file - $@");
			return 0;
		}
		if ( $line1 =~ m/^@/ && $line5 =~ m/^@/ && $line2 =~ m/[ACGTN]/){
			return 1;
		} else {
			$self->errLog("Error: isFastQ: $file doesn't appear to be a valid FastQ file");
		}
	} else {
		$self->errLog("Error: isFastQ: $file doesn't exist on the filesystem");
	}
}

##this gets a lane for each fastq in the pipeline object and adds it to the lane field for the onject
sub getLanes{
	my $self = shift;
	for my $samp (keys (%{$self->{samples}} )){
		my $name = $self->{samples}->{$samp}->{files}->{read1}->{name};
		if ($name && -e $name){
			my $IN = $self->openFastQ($name);
			my $line = <$IN>;
			my @flds = split (":", $line);
			my $lane;
			if ( $line =~ m/\D\d$/){
				$lane = $flds[1];
			} else {
				$lane = $flds[3];
			}
			$self->{samples}->{$samp}->{lane} = $lane;
		}
	}
}


##this opens a fastq file
sub openFastQ{
	my $self = shift;
	my $file = shift;
	my $IN;
	if ($file =~ m/.bz2$/){
		open $IN , "bzcat \"$file\" |" or $self->errLog("Error: isFastQ: Could not open $file - $!"); 
	} elsif ( $file =~ m/.gz$/ ){
		open $IN , "zcat \"$file\" |" or $self->errLog("Error: isFastQ: Could not open $file - $!");
	} else {
		open $IN , '<' , $file or $self->errLog("Error: isFastQ: Could not open $file - $!");
	}
	if($IN){
		return $IN;
	} else {
		return 0;
	}
}


## this creates a multiprocessed queue taking an array of system calls
sub multiprocessQ{
	my $self = shift;
	my %args = @_;
	for my $call (sort keys %args){
		my $f = fork;
		$self->runLog($call);
		if (defined ($f) and $f==0 ){
			$self->runLog($call);
			system("nice $call");
			exit(0);
		}
		1 while waitpid(-1, WNOHANG) > 0;
	}
}

## this builds a hash reference of system commands to run fastqc on the sample fastqs
##  will be referred to as 'fastqc' in log field
sub buildFastQCCommandHash{
	my $self = shift;
	$self->runLog("Building FastQC Command Hash");
	$self->{curStep} = "fastqc";
	my $analysisDir = $self->{analysisDir};
	my $log = $self->{logDir}."/buildFastQCCommands.log";
	logMsg($self->{time}, $log);
	for my $samp ( sort keys %{ $self->{samples} } ){
		$self->{curSamp} = $samp;
		my $exp = getExp($samp);
		if (!exists($self->{samples}->{$samp}->{files}->{read1}->{name}) ){
			my $msg = "Error: buildFastQCCommandHash: Read 1 doesn't exist for sample $samp";
			$self->errLog( $msg );
			$self->logErr($msg);
			logMsg("$samp\tRead1\t$msg", $log);
			next;
		}
		unless ($self->isPairedEnd && exists($self->{samples}->{$samp}->{files}->{read2}->{name}) ){
			my $msg = "Error: buildFastQCommandHash: Read 2 is paired end but it doesn't exits for $samp";
			$self->errLog($msg);
			$self->logErr($msg);
			logMsg("$samp\tRead2\t$msg", $log);
			next;
		}
		my $succMsg;
		my $outdir = $analysisDir."/".$exp."/fastqc";
		system ("mkdir -p $outdir");
		my $read1 = $self->{samples}->{$samp}->{files}->{read1}->{name};
		my $read1base = basename($read1);
		if ($read1base =~ m/([^.]+)\./){
			$read1base = $1;
		}
		if (-e $outdir."/".$read1base."_fastqc.zip"){
			$succMsg = "Success: $outdir/${read1base}_fastqc.zip exists;";
			
		} else {
			my $command = "fastqc --outdir=$outdir ".$read1;
			$self->{fastqc}->{commands}->{$command} = $samp;
			$succMsg = $command;
			
		}
		if ($self->isPairedEnd){
			my $read2 = $self->{samples}->{$samp}->{files}->{read2}->{name};
			my $read2base = basename($read2);
			if ($read2base =~ m/([^.]+)\./){
				$read2base = $1;
			}	
			if ( -e  $outdir."/".$read2base."_fastqc.zip") {
				$succMsg .= " $outdir/${read1base}_fastqc.zip exists;";
			} else {
				my $command = "fastqc --outdir=$outdir ".$read2;
				$self->{fastqc}->{commands}->{$command} = $samp;
				$succMsg .= " ; $command";
			}
		}
		logMsg("$samp\t$succMsg", $log);
		$self->runLog($succMsg);
		$self->logSuccess($succMsg);
	}
}



## this builds a hash reference of system commands to align the sample fastqs
## will be referred to as 'align' in log field
sub buildAlignCommandHash{
	my $self = shift;
	print Dumper ($self);
	$self->runLog("Building alignment command hash");
	$self->{curStep} = 'align';
	my $log= $self->{logDir}."/align.log";
	logMsg($self->{time}, $log);
	my $analysisDir = $self->{analysisDir};
	for my $samp ( sort keys %{ $self->{samples} }){
		$self->{surSamp} = $samp;
		my $exp = getExp($samp);
		my $type = $self->{samples}->{$samp}->{type};
		my $tool;
		my $pattern1;
		my $pattern2; 
		## check if read 1 exists
		if (!exists($self->{samples}->{$samp}->{files}->{read1}->{name}) ){
			my $msg = "Error: buildAlignCommandHash: Read 1 doesn't exist for sample $samp";
			$self->errLog( $msg );
			$self->logErr($msg);
			logMsg("$samp\t$msg", $log);
			next;
		}
		
		if (!(-e $self->{samples}->{$samp}->{files}->{read1}->{name})){
			my $msg = "$samp doesn't exist on the filesystem";
			$self->errLog($msg);
			$self->logErr($msg);
			logMsg("$samp\t$msg", $log);
			next;
		}
		if ($type eq "CHIPSEQ" || $type eq "DNASEQ" ){
			## usebowtie
			$tool = "bowtie";
			$pattern1 = "Bowtie";
			my $msg = "Bowtie alignment currently not supported";
			$self->logErr($msg);
			next;
		} elsif ($type eq "MRNASEQ" || $type eq "DMRNASEQ"){
			## use tophat
			$tool = "tophat";
			$pattern1 = "Tophat";
			$pattern2 = "accepted_hits"
		} elsif( $type eq "TDNASEQ"){
			## use BWA
			$tool = "bwa";
			my $msg = "BWA alignment currently not supported";
			$self->logErr($msg);
			next;
		}else {
			my $msg = "Error: buildAlignCommandHash: $samp has invalid type : $type";
			$self->errLog($msg);
			$self->logErr($msg);
			next;
		}
		
		# Check if the force flag is ineffect
		unless ( exists ($self->{align}->{force}) && $self->{align}->{force} ) {
			$self->runLog("checking for alignment file");
			if (my @files = glob("$analysisDir/$exp/$pattern1*--$samp/$pattern2*.bam")){
				## alignment has been done already
				if (-s $files[0] < 104857600 ){
					logMsg("$samp\t".$files[0]."is less than 100BM -- rerunning", $log);
				} else{
					my $msg = "$samp already has a bam file";
					$self->runLog($msg);
					$self->logSuccess($msg);
					logMsg("$samp\t$msg ".$files[0], $log);
					next;
				}
			} 
		}
		if ($tool eq "bowtie"){
			next;
		} elsif ($tool eq "bwa"){
			next;
		} elsif ( $tool eq "tophat"){
			my $command = $self->generateTophatCommand($samp);
			if (defined($command) && $command ne ""){
				$self->{align}->{commands}->{$command} = $samp;
				$self->logSuccess("created command: $self->generateTophatCommand($samp)");
				logMsg("$samp\t$command", $log);
			} else {
				logMsg($samp."\tthere was an error creating the command", $log);
			}
		}		
	}
}

## this takes a sample and generates a new Tophat command
sub generateTophatCommand{
	my $self = shift;
	my $samp = shift;
	my $species = $self->getSpecies($samp);
	my $paired = $self->isPairedEnd;
	my $HAT_ARGS="--min-anchor=5 --segment-length=25 --no-coverage-search --segment-mismatches=2 --splice-mismatches=2 --microexon-search --no-mixed ";
	my $PAIRED_MATE_INNER_DIST=250;
	my $NUM_THREADS = 12;
	if ( !( $species ) ){
		return;
	}
	if (!($paired) && $self->{samples}->{$samp}->{files}->{read1}->{name} =~ m/R[12]/){
		my $msg = $self->{samples}->{$samp}->{files}->{read1}->{name}." looks like a paired end read but $samp marked as paired end check database";
		print $msg."\n";

		return;
	}
	$self->runLog("making Tophat command for $samp");
	my $TOPGTF = $self->{db}->{$species}->{gtf};
	my $TOPINDEX = $self->{db}->{$species}->{index};
	my $exp = getExp($samp);
	my $outpath = $self->{analysisDir}."/$exp/Tophat_B2B$exp--$samp";
	system ("mkdir -p $outpath");
	my $aboutfile = "$outpath/RunLog-Tophat_B2B$exp--$samp.log";
	my $read1 = $self->{samples}->{$samp}->{files}->{read1}->{name};
	my $tophatCommand = "tophat -o ${outpath} ${HAT_ARGS} ";
	$tophatCommand .= "--GTF=${TOPGTF} --num-threads=${NUM_THREADS} ";
	$tophatCommand .= "--mate-inner-dist=250 " if ($paired);
	$tophatCommand .= "--rg-sample $samp ";
	$tophatCommand .= "--rg-id $samp:lane_".$self->{samples}->{$samp}->{lane}." ";
	$tophatCommand .= "${TOPINDEX} $read1 ";
	$tophatCommand .= $self->{samples}->{$samp}->{files}->{read2}->{name}." " if ($paired);
	$tophatCommand .= "2>&1 | tee --append $aboutfile";
	return $tophatCommand;
}

## returns true if a paired end sample
sub isPairedEnd{
	my $self = shift;
	my $samp = $self->{curSamp};
	if (exists( $self->{samples}->{$samp}->{files}->{read2}->{name}) && $self->{samples}->{$samp}->{files}->{read2}->{name} ){
		return 1;
	}
	return 0;
}

### takes a sample and returns an exeriment;
sub getExp{
	my $samp = shift;
	my $exp = $samp;
	$exp =~ s/X.*//;
	$exp .= "R";
	return $exp;
}

## this returns a species from the sample
sub getSpecies{
	my $self = shift;
	my $samp = shift;
	my $rawSpecies = $self->{samples}->{$samp}->{species};
	switch ( lc($rawSpecies )){
		case /mouse/	{ return "mouse" }
		case /human/	{ return "human" }
		case /chicken/  { return "chicken" }
		case /zebrafish/ {return "zebrafish"}
	}
	$self->errLog("Error: getSpecies: $rawSpecies not supported");
	return 0;
}

# returns hash ref with path to fastqs or error for a given sample;
sub findFastQPathSample{
	my %args 	= 	@_;
	my $file = $args{file};
	my $exp = $args{exp};
	print "Finding path for file \t $file\n";
	my $dataPath  = $args{path};
	$dataPath = "/Data01/gnomex/ExperimentData/" unless defined $dataPath;
	my $cmd = "find \"$dataPath\" -type f -path \"*$exp*/*$file*\"";
	print "Executing\t$cmd\n";
	my $path = capture($cmd);
	print "path=\t$path" if ($debug);
	if ($path eq '') {
		return 0;
	} else {
		chomp $path;
		return $path;
	}
}

## what it says
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