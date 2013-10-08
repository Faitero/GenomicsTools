package tools::mongo;
use strict;
use warnings;
use Data::Dumper;
use MongoDB::MongoClient;


# returns a handle to the samples collection.
sub getDB{
	my $client = MongoDB::MongoClient->new;
	my $db = $client->get_database( 'b2bPipeline' );
	
	return $db;
}

# Adds entries from sample hashes to the sample collection.
# Reqd. parameters:
# sh=>$sampleHash
sub addSamps{
	my %args = @_;
	my $sh = $args{sh};
	my $db = getDB();
	my $sampColl = $db->get_collection( 'samples' );
	my $expColl = $db->get_collection( 'experiments' );
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
sub addPipelineSamples{
	my %args = @_;
	my $sh = $args{sh};
	my $db = getDB;
	my $sampColl = $db->get_collection( 'samples' );
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
	my %args= @_;
	my $samp = $args{samp};
	my $read = $args{read};
	print "mongo::checkSampFileDB\t$samp\t$read\n";
	my $db = getDB();
	my $coll = $db->get_collection('samples');
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
	my %args = @_;
	my $coll = $args{coll};
	my $qfld = $args{qfld};
	my $qval = $args{qval};
	my $fld = $args{fld};
	my $val = $args{val};
	my $db = getDB();
	$coll = $db->get_collection($coll);
	$coll->update(
		{$qfld=>$qval},
		{'$set'=> {$fld => $val} });
}

#adds the fastq dir to the experiment collection
# requred parameters :
# exp => $exp,
# dir => $dir,
sub addFastQDir {
	my %args = @_;
	my $exp = $args{exp};
	my $dir = $args{dir};
	addField(
		coll => "experiments",
		qfld => "_id",
		qval => $exp,
		fld => "fastqdir",
		val => $dir,
		);
}

## add rows to meas collection
sub addBamStat{
	my $db = getDB();
	my $sampColl = $db->get_collection('samples');
	my $measColl = $db->get_collection('meas');
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
		insertMeasColl(
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
	my %args = @_;
	return queryMeasColl(
		exp => $args{exp},
		file => $args{file},
		expDir => $args{expDir},
		meas => "bamstat",
	);
}

sub addMarkDuplicates{
	my $db = getDB();
	my $sampColl = $db->get_collection('samples');
	my $measColl = $db->get_collection('meas');
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
		insertMeasColl(
			samp => $samp,
			sh => $sh,
			file => $file,
		);
	}
}

sub getMarkDuplicates {
	my %args = @_;
	my $exp = $args{exp};
	my $file = $args{file};
	my $expDir = $args{expdir};

	return queryMeasColl (
		exp => $exp,
		file => $file,
		expDir =>$expDir,
		meas => "markDuplicates",
	);
}

sub addInnerDistance {
	my $db = getDB();
	my $sampColl = $db->get_collection('samples');
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
		insertMeasColl(
			samp => $samp,
			sh => $sh,
			file => $file,
		);
	}
}

sub getInnerDistance {
	my %args = @_;
	return queryMeasColl(
		exp => $args{exp},
		file => $args{file},
		expDir => $args{expDir},
		meas => "inner_distance",
	);
}

sub addGeneBodyCov {
	my $db = getDB();
	my $sampColl = $db->get_collection('samples');
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
		insertMeasColl(
			samp => $samp,
			sh => $sh,
			file => $file,
		);
	}
}

sub getGeneBodyCov {
	my %args = @_;
	return queryMeasColl(
		exp => $args{exp},
		file => $args{file},
		expDir => $args{expDir},
		meas => "gene_body_coverage",
	);
}






sub insertMeasColl{
	my $db = getDB();
	my $measColl = $db->get_collection('meas');
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
	my $db = getDB();
	my $measColl = $db->get_collection('meas');
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