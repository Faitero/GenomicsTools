package tools::mongo;
use strict;
use warnings;
use Data::Dumper;
use MongoDB::MongoClient;
use Scalar::Util qw(blessed);
use String::Diff;
our @ISA = qw(MongoDB::MongoClient);
our $debug = 1;


## gets the list of fastqc data files and returns a hashref to them
sub get_fastqc_data {
	my %args = @_;
	my $basedir = $args{basedir};
	if (!defined ($basedir) ){
		$basedir = "/Data01/gnomex/Analysis/experiment/";
	}
	my @oh = `find $basedir -name fastqc_data.txt | sort`;
	for my $entry (@oh){
		chomp $entry;
	}
	return \@oh;
}


## gets the sample from the a collection of fastqc paths
## returns hashRef->$path->$sample 
sub get_sample_from_path {
	my @paths = @_;
	my $oh = {};
	for my $path ( @paths ){
		$path =~ m/\/(\d+X\d+)\//;
		$oh->{$path} = $1;
	}
	return $oh;
}

## takes a string of a fastqc path 
## returns string that can be used to search for the fastq
sub fascqc_path_to_file_pattern{
	my $path  = shift;
	if (  $path =~ m/\/([^\/]+)_fastqc\// ){
		print caller()." $1\n" if $debug;
		return $1;
	}	
	else {return 0};
}


## get the read1 and read 2 from the sample
## takes a path-samp hash
## returns a samp->read1[path array]
## 			samp->read2[path array]
sub get_read_from_path {
	my $pipe = shift;
	## first invert hash
	my $p2s = shift;
	# print Dumper $pipe;
	## make out hash
	my $oh = {};
	print Dumper $pipe->{samples} if $debug;
	#make sample files Ref hash reference
	for my $path ( sort keys %$p2s){
		my $samp = $p2s->{$path};
		#make sample files Ref hash reference
		my $fr = $pipe->{samples}->{$samp}->{files}; 
		my $pattern = fascqc_path_to_file_pattern($path);
		my $i = 0 ; 
		for my $read1 ( $fr->{fastqRead1}->{path} ){
			print $read1."\n" if $debug;
			if ( $read1 =~ $pattern ){
				if (!exists($oh->{$samp}->{read1} )){
					$oh->{$samp}->{read1} = ();
				}
				${$oh->{$samp}->{read1}}[$i] = $path;
				# print Dumper $oh;
				# <>;
				next; 
			}
			$i++;
		}
		$i = 0;
		if (exists ($fr->{fastqRead2} )  ){
			if ($samp eq "5X1"){
				print "$samp in Read2 loop\n";
				<>;
			}
			for my $read2 ( $fr->{fastqRead2}->{path} ){
				print $read2."\n" if $debug;
				if ( $read2 =~ $pattern ){
					if (!exists($oh->{$samp}->{read2} )){
						$oh->{$samp}->{read2} = ();
					}
					${$oh->{$samp}->{read2}}[$i] = $path;
					next; 
				}
				$i++;
			}
		}	
	}
	return $oh;
}

## makes objects from fastqc data output
sub process_fastqc_data{
	my $pipe = shift;
	my $mh = ();
	my $measColl = $pipe->{mdb}->{db}->get_collection('measurements');
	for my $samp (sort keys ( %{ $pipe->{samples} } ) ){
		for my $read ( sort keys ( %{ $pipe->{samples}->{$samp}->{files}->{fastqc} } ) ){
			my $i ;
			my $rh = $pipe->{samples}->{$samp}->{files}->{fastqc}->{$read};
			for ( $i = 0 ; $i <  scalar ( @{$rh} ) ; $i++ ){
				
				my $file = $$rh[$i];
				if ( -e $file ){
					my $key = "${samp}  ${read}  No. ".eval($i+1);
					## get per sequence quality
					my $data = get_per_sequence_qual($key, $file, $samp) ;
					my $id = $data->{measurement}."__".$key;
					$id =~ s/ /_/g;
					$measColl->update({'_id' => $id }, {'$set' => { %{$data} }}, {'upsert' => 1, 'safe' => 1});
					## get per sequence GC content
					$data = get_per_sequence_gc($key, $file, $samp);
					my $id = $data->{measurement}."__".$key;
					$id =~ s/ /_/g;
					$measColl->update({'_id' => $id }, {'$set' => { %{$data} }}, {'upsert' => 1, 'safe' => 1});

				}
			}
		}
	}
}

## makes a measurement object from fastqc file for per sequence quality scores
sub get_per_sequence_gc{
	my $key = shift;
	my $file = shift;
	my $samp = shift;
	my $out = {};
	my $meas = "Per sequence GC content";
	open my $IN,  "<" , $file ;
	until(  (my $line = <$IN>) =~ qq/$meas/  ){};
	my $line = <$IN>;
	$line =~ s/#//;
	chomp $line;
	my @flds = split (" ", $line);
	# print $flds[0]."\t".$flds[1]."\n";
	# $out->{_id} = $_id;
	$out->{key} = $key;
	$out->{file} = $file;
	$out->{x_label} = $flds[0];
	$out->{y_label} = "Proportion of Reads";
	$out->{values} = ();
	$out->{samp} = $samp;
	$out->{measurement} = $meas; 
	until ( ($line = <$IN>) =~ m/>>END_MODULE/ ){
		my $val = {};
		chomp $line ;
		my @vals = split("\t", $line);
		$vals[1] =~ s/E/e/;
		$val->{x} = $vals[0];
		$val->{y} = $vals[1];
		$out->{y_total} += $vals[1];
		push ( @{ $out->{values} } , $val );
	}
	for my $val ( @{ $out->{values} } ){
		$val->{y} /= $out->{y_total};
	}
	return $out; 
}


## makes a measurement object from fastqc file for per sequence quality scores
sub get_per_sequence_qual{
	my $key = shift;
	my $file = shift;
	my $samp = shift;
	my $out = {};
	my $meas = "Per sequence quality scores";
	open my $IN,  "<" , $file ;
	until(  (my $line = <$IN>) =~ qq/$meas/  ){};
	my $line = <$IN>;
	$line =~ s/#//;
	chomp $line;
	my @flds = split (" ", $line);
	# print $flds[0]."\t".$flds[1]."\n";
	# $out->{_id} = $_id;
	$out->{key} = $key;
	$out->{file} = $file;
	$out->{x_label} = $flds[0];
	$out->{y_label} = "Proportion of Reads";
	$out->{values} = ();
	$out->{samp} = $samp;
	$out->{measurement} = $meas; 
	until ( ($line = <$IN>) =~ m/>>END_MODULE/ ){
		my $val = {};
		chomp $line ;
		my @vals = split(" ", $line);
		$vals[1] =~ s/E/e/;
		$val->{x} = $vals[0];
		$val->{y} = $vals[1];
		$out->{y_total} += $vals[1];
		push ( @{ $out->{values} } , $val );
	}
	for my $val ( @{ $out->{values} } ){
		$val->{y} /= $out->{y_total};
	}
	return $out; 
}






# ## this takes an array ref like created from get_read_from_path
# ## and outputs a data object suitable for insertion into the mondodb
# ## collection
# sub parse_fastqc_data{
# 	my $sh = shift;
# 	## get list 
# 	my $oh = {};
# 	for my $samp ( sort keys %$sh ){
# 		for my $read ( sort keys %{$sh->{$samp} }){
# 			my $i = 0;
# 			for ( $i = 0 ; $i < scalar @{$sh->{$samp}->{$read:while () {
# 				# body...
# 			}}} ; $i++ ){
# 				my $key = "$samp $read $i";
# 				open my $IN , "<" , ${$sh->{$samp}->{$read}}[$i];
# 				my $version = <$IN>;
# 				chomp $version;
# 				<>;
# 				<>;
# 				my $line = <$IN>;
# 				chomp $line; 
# 				my $filename = @{split( "\s" , $line )}[1];	
# 				print $filename;
# 			}
# 		}
# 	}
# }






# returns a handle to the samples collection.
sub getDB{
	
	my $client = MongoDB::MongoClient->new;
	my $db = $client->get_database( 'restTest' );
	
	return $db;
}

sub new{
	my $class = shift;
	my $self = {};
	$self->{db} = getDB;
	print "db->{db} = ". blessed($self->{db})."\n";
	bless ($self, $class);
	return $self;
}

# Adds entries from sample hashes to the sample collection.
# Reqd. parameters:
# sh=>$sampleHash
sub addSamps{
	my $db = shift;
	my %args = @_;
	my $sh = $args{sh};
	my $sampColl = $db->{db}->get_collection( 'samples' );
	my $expColl = $db->{db}->get_collection( 'experiments' );
	# for testing DB
	# $sampColl->drop;
	# $expColl->drop;
	for my $key (keys %$sh){
	    my $expAtts = $sh->{$key};
		my $exp = $key;
		$exp =~ s/(\d+)X\d+/$1/;
		$exp .= "R";
	
	    print "Adding sample $key to the databse\n";
	    $sampColl->update({"_id"=>$key}, {'$set' => { "meta"=>{%$expAtts} }},{"upsert"=> 1, "safe"=>1});
	    $expColl->update({"_id"=>$exp}, {'$set' => {} },{"upsert"=> 1, "safe"=>1});
	}
	## add organism field and directory to experiments collection
	my $expCurs = $expColl->find;
	while (my $doc = $expCurs->next){
    	my $exp = $doc->{'_id'};
    	my $sampCur = $sampColl->query({"experiment_id"=>$exp}, {limit => 1});
    	my $samp = $sampCur->next;
    	my $org = $samp->{"Organism"};
    	if( $org ){
    		$expColl->update({"_id"=>$exp}, { '$set'=>{"Organism"=>$org}}, {"safe" => 1});
    	} else {
    		$expColl->update({"_id"=>$exp}, { '$set'=>{"Organism"=>"NA"}}, {"safe" => 1});
    	}
    	## add experiment directory
    	if (!defined( $doc->{"Fastq Dir"} )){
    		my $files = $samp->{'Associated Files'};
    		my $dir;
    		if ($files){
	    		$files =~ s/[",]//g;
	    		my @files = split (" ", $files);
    			my $dataPath = "/Data01/gnomex/ExperimentData/";
    			$dir = `find $dataPath -name *$files[0]*`;

    			if ($dir){
    				chomp $dir;
    			} else {
    				$dir = "Not found";
    			}
    		} else {
    			$dir = "No Files on Sample Sheet";
    		}
    		$expColl->update({"_id"=>$doc->{'_id'}}, {'$set' => { "Fastq Dir" => $dir }}, {"safe" => 1});
    		print "added Dir\t$dir\n";
    	}
   	}
}



## this adds the pipeline->Samples object to the samples databse
## required Params : sh => hashRef of Samples 
sub lineSamples{
	my $db = shift;
	my %args = @_;
	my $sh = $args{sh};
	my $sampColl = $db->{db}->get_collection( 'samples' );
	for my $key ( keys %$sh ){
		print "Adding sample $key to the database\n";
		# print Dumper($sh->{$key});
		$sampColl->update({"_id"=>$key},
		 { '$set'=>{%{$sh->{$key}}} }, 
		 {'upsert'=>1, "safe"=>1} );
	}
}

##checks the database for files returning the path if present
## and 0 if not present;
sub checkSampFileDB{
	my $db = shift;
	my %args= @_;
	my $samp = $args{samp};
	my $read = $args{read};
	print "mongo::checkSampFileDB\t$samp\t$read\n";
	my $coll = $db->{db}->get_collection('samples');
	my $file = $coll->find( { "_id" => $samp} );
	my $fileObj = $file->next;
	my $filename = $fileObj->{files}->{$read}->{name};
	# print "Filename=\t$filename\n";
	if ($filename && -e $filename){
		print "FileName in database\t".$filename."\n";
		if (-e $filename){
			print "File exists\n";
			return $filename;
		}
	} else {
		print "$filename not in database\n";
		return 0;
	}
}




# adds field to the given collection
# req params:
# coll => $collection,
# qfld => $queryField
# qval => $queryVal
# fld => $field,
# val => $value,
# db => $database,
sub addField{
	my $db = shift;
	my %args = @_;
	my $coll = $args{coll};
	my $qfld = $args{qfld};
	my $qval = $args{qval};
	my $fld = $args{fld};
	my $val = $args{val};
	$coll = $db->{db}->get_collection($coll);
	$coll->update(
		{$qfld=>$qval},
		{'$set'=> {$fld => $val} });
}

#adds the fastq dir to the experiment collection
# requred parameters :
# exp => $exp,
# dir => $dir,
sub addFastQDir {
	my $db = shift;
	my %args = @_;
	my $exp = $args{exp};
	my $dir = $args{dir};
	$db->addField(
		coll => "experiments",
		qfld => "_id",
		qval => $exp,
		fld => "fastqdir",
		val => $dir,
		);
}

## add rows to meas collection
sub addBamStat{
	my $db = shift;
	my $sampColl = $db->{db}->get_collection('samples');
	my $measColl = $db->{db}->get_collection('meas');
	$measColl->drop;
	my %args = @_;
	my $exp = $args{exp};
	my $file = $args{file};
	my $expDir = $args{expdir};
	my $expPat = $exp;
	$expPat =~ s/R//;
	$expPat.="X";
	# print $exp."\n";
	$file = "accepted_hits_marked_dup.bam" unless ( defined ( $file ) );
	$expDir = "/Data01/gnomex/Analysis/experiment/$exp" unless ( defined ($expDir ) );
	my $samps = $sampColl->find( { "_id" => qr/^$expPat/} );
	

	while ( my $sampObj = $samps->next){
		my @orderedFlds = ();
		my $samp = $sampObj->{'_id'};
		my $bsTab = $expDir."/Tophat_".$sampObj->{"meta"}->{'Bam Root'}."/$file.bam_stat.tab";
		my $rh = {};	 # row hash	
		$rh->{meas} = "bamstat";
		open IN, "< $bsTab" or die "cannot open $bsTab $!\n";		
		while (my $line = <IN> ) {
			if( $line =~ /\d+$/){
				## data row
				# print $line."\n";
				if ( $line =~ m/Non Primary Hits\s+(\d+)/) {
					$rh->{'data'}->{Non_Primary_Hits} = $1;
					push(@orderedFlds,"Non_Primary_Hits");
				} else {
					my @elems = split ( ":", $line );
					for my $elem ( @elems ){
						chomp($elem);
						$elem =~ s/^\s+//;
					}
					$elems[1] =~ s/^\s+//;
					$rh->{'data'}->{$elems[0]} = $elems[1];
					push(@orderedFlds, $elems[0]);
				}
			}
		}
		$rh->{file}= $file;
		$rh->{samples_id} = $samp; 
		# print "adding to sample $samp\n";
		$rh->{fields} =\@orderedFlds;
		## check if the bamstat is set yet
		# $measColl->update({'_id' => "bamstat_".$samp."_".$file}, { '$set' => { "file" => $file, "samples_id" => $samp, 'type'=>'bamstat', 'data' =>{%$rh}} }, {'upsert' => 1, 'safe' => 1});
		# $measColl->update({'_id' => "bamstat_".$samp."_".$file}, { '$addToSet' => { "fields" => {'$each' =>  [@orderedFlds] } } }, {'upsert' => 1, 'safe' => 1});
		$db->insertMeasColl(
			samp => $samp,
			sh => $rh,
			file => $file,
		);
		# print Dumper($rh);
	}
}


## this takes an experiment and a file agrument and returns a hash reference:
## returnHash-> {id} -> {fldName} ->{value}
##		
sub getBamStat{
	my $db = shift;
	my %args = @_;
	return $db->queryMeasColl(
		exp => $args{exp},
		file => $args{file},
		expDir => $args{expDir},
		meas => "bamstat",
	);
}

sub addMarkDuplicates{
	my $db = shift;
	my $sampColl = $db->{db}->get_collection('samples');
	my $measColl = $db->{db}->get_collection('meas');
	my %args = @_;
	my $exp = $args{exp};
	my $file = $args{file};
	my $expDir = $args{expdir};
	my $expPat = $exp;
	$expPat =~ s/R//;
	$expPat.="X";
	# print $exp."\n";
	$file = "accepted_hits_marked_dup.bam" unless ( defined ( $file ) );
	$expDir = "/Data01/gnomex/Analysis/experiment/$exp" unless ( defined ($expDir ) );
	my $samps = $sampColl->find( { "_id" => qr/^$expPat/} );
	
	while ( my $sampObj = $samps->next){
		my @orderedFlds = ();
		my $sh = {}; ##samp hash
		$sh->{meas} = "markDuplicates";
		my $samp = $sampObj->{'_id'};
		my $bsTab = $expDir."/Tophat_".$sampObj->{"meta"}->{'Bam Root'}."/${file}_mark_dup_metrics.txt";
		open IN, "< $bsTab" or die "cannot open $bsTab $!\n";		
		<IN>;
		my $command = <IN>;
		chomp $command;
		$command =~ s/^#\s+//;
		# print "command ".$command."\n";
		<IN>;
		my $date = <IN>;
		chomp $date;
		$date =~ s/^#\s+//;
		# print "date ".$date."\n";
		$sh->{command} = $command;
		$sh->{date} = $date;
		$sh->{file} = $file;
		$sh->{samples_id} = $samp;
		
		<IN>;
		<IN>;
		my @flds = split("\t",<IN>);
		my @vals = split("\t",<IN>);
		for (my $i = 0 ; $i < @flds ; $i ++ ){
			my $fld = $flds[$i];
			my $val = $vals[$i];
			chomp $fld;
			chomp $val;
			push(@orderedFlds, $fld);
			# print $fld."\t".$val."\n";
			$sh->{data}->{$fld} = $val;
		}
		## this will ensure the order of the fields
		$sh->{fields} = \@orderedFlds;
		# $measColl->update({'_id' => "markdup_".$samp."_".$file}, { '$set' => {"file" => $file, "samples_id" => $samp, 'type'=>'markdup', "command"=>$command, 'data' =>{%$rh}} }, {'upsert' => 1, 'safe' => 1});
		# $measColl->update({'_id' => "markDuplicates_".$samp."_".$file}, {'$set' => { %{$sh} }}, {'upsert' => 1, 'safe' => 1});
		$db->insertMeasColl(
			samp => $samp,
			sh => $sh,
			file => $file,
		);
	}
}

sub getMarkDuplicates {
	my $db = shift;
	my %args = @_;
	my $exp = $args{exp};
	my $file = $args{file};
	my $expDir = $args{expdir};

	return $db->queryMeasColl (
		exp => $exp,
		file => $file,
		expDir =>$expDir,
		meas => "markDuplicates",
	);
}

sub addInnerDistance {
	my $db = shift;
	my $sampColl = $db->{db}->get_collection('samples');
	my %args = @_;
	my $exp = $args{exp};
	my $file = $args{file};
	my $expDir = $args{expdir};
	my $expPat = $exp;
	$expPat =~ s/R//;
	$expPat.="X";
	# print $exp."\n";
	$file = "accepted_hits_marked_dup.bam" unless ( defined ( $file ) );
	$expDir = "/Data01/gnomex/Analysis/experiment/$exp" unless ( defined ($expDir ) );
	my $samps = $sampColl->find( { "_id" => qr/^$expPat/} );
	
	while ( my $sampObj = $samps->next){
		my @orderedFlds = ();
		my $sh = {}; ##samp hash
		my $samp = $sampObj->{'_id'};
		my $bsTab = $expDir."/Tophat_".$sampObj->{"meta"}->{'Bam Root'}."/${file}_QC/inner_distance.inner_distance_freq.txt";
		open IN, "< $bsTab" or die "cannot open $bsTab $!\n";		
		my @bin1 = ();
		my @bin2 = ();
		my @freq = ();
		while ( my $line = <IN>){
			my @flds = split(" ", $line);
			for ( my $i = 0 ; $i < @flds ; $i++ ){
				chomp $flds[$i];
			}
			push (@bin1, $flds[0]);
			push (@bin2, $flds[1]);
			push (@freq, $flds[2]);
		}
		$sh->{samples_id} = $samp;
		$sh->{file} = $file;
		$sh->{meas} = "inner_distance";
		$sh->{data}->{bin1} = \@bin1;
		$sh->{data}->{bin2} = \@bin2;
		$sh->{data}->{freq} = \@freq;
		$db->insertMeasColl(
			samp => $samp,
			sh => $sh,
			file => $file,
		);
	}
}

sub getInnerDistance {
	my $db = shift;
	my %args = @_;
	return $db->queryMeasColl(
		exp => $args{exp},
		file => $args{file},
		expDir => $args{expDir},
		meas => "inner_distance",
	);
}

sub addGeneBodyCov {
	my $db = shift;
	my $sampColl = $db->{db}->get_collection('samples');
	my %args = @_;
	my $exp = $args{exp};
	my $file = $args{file};
	my $expDir = $args{expdir};
	my $expPat = $exp;
	# my $meas = "markDuplicates";
	$expPat =~ s/R//;
	$expPat.="X";
	# print $exp."\n";
	$file = "accepted_hits_marked_dup.bam" unless ( defined ( $file ) );
	$expDir = "/Data01/gnomex/Analysis/experiment/$exp" unless ( defined ($expDir ) );
	my $samps = $sampColl->find( { "_id" => qr/^$expPat/} );
	
	while ( my $sampObj = $samps->next){
		my @orderedFlds = ();
		my $sh = {}; ##samp hash
		my $samp = $sampObj->{'_id'};
		my $bsTab = $expDir."/Tophat_".$sampObj->{"meta"}->{'Bam Root'}."/${file}_QC/gene_body_coverage.geneBodyCoverage.txt";
		open IN, "< $bsTab" or die "cannot open $bsTab $!\n";		
		my $line  = <IN>;
		chomp $line;
		my @totReads = split(":", $line);
		$totReads[1] =~ s/^\s+//g;
		$sh->{data}->{$totReads[0]} = $totReads[1];
		$line = <IN>;
		chomp $line;
		my @fragNo = split(":", $line);
		$fragNo[1] =~ s/^\s+//g;
		$sh->{data}->{$fragNo[0]} = $fragNo[1];
		<IN>;
		my @percentile = ();
		my @count = ();
		while (my $line = <IN>){
			chomp $line;
			my @flds = split(" ", $line);
			push(@percentile, $flds[0]);
			push(@count, $flds[1]);
		}
		$sh->{data}->{percentile} = \@percentile;
		$sh->{data}->{count} =  \@count;
		$sh->{samples_id} = $samp;
		$sh->{file} = $file;
		$sh->{meas} = "gene_body_coverage";
		# print Dumper($sh);
		$db->insertMeasColl(
			samp => $samp,
			sh => $sh,
			file => $file,
		);
	}
}

sub getGeneBodyCov {
	my $db = shift;
	my %args = @_;
	return queryMeasColl(
		exp => $args{exp},
		file => $args{file},
		expDir => $args{expDir},
		meas => "gene_body_coverage",
	);
}






sub insertMeasColl{
	my $db = shift;
	my $measColl = $db->{db}->get_collection('meas');
	my %args = @_;
	my $file = $args{file};
	my $sh = $args{sh};
	my $samp = $args{samp};
	my $meas = $sh->{meas};
	# print Dumper($sh);
	my $_id = $meas."__".$samp."__".$file;
	print "inserting\t$_id\n";	

	$measColl->update({'_id' => $meas."__".$samp."__".$file}, {'$set' => { %{$sh} }}, {'upsert' => 1, 'safe' => 1});
}

sub queryMeasColl{
	my $db = shift;
	my $measColl = $db->{db}->get_collection('meas');
	my %args = @_;
	my $exp = $args{exp};
	my $file = $args{file};
	my $expDir = $args{expdir};
	my $meas = $args{meas};
	my $expPat = $exp;
	$expPat =~ s/R//;
	$expPat.="X";
	# print $exp."\n";
	$file= "accepted_hits_marked_dup.bam" unless ( defined ( $file ) );
	$expDir = "/Data01/gnomex/Analysis/experiment/$exp" unless ( defined ($expDir ) );
	my $samps = $measColl->find( { "samples_id" => qr/^$expPat/, 'file' => $file , 'type' => $meas} );
	my $rethash = {};
	while (my $sampObj = $samps->next){
		$rethash->{$sampObj->{samples_id} } = $sampObj;	
	}
	# print "subDumper\t".Dumper($rethash)."\n";
	return $rethash;
}


1;