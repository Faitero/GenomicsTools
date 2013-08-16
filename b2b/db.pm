package b2b::db;
use strict;
use DBI;
## connects to the b2b database
sub dbConnect{
	my $dbp = "DBI:mysql:b2b";
	my $user = "b2b";
	my $pw = "cvdc";
	my $dbh = DBI->connect($dbp, $user, $pw) or die "Couldn't connect to database : ".DBI->errstr;
	# $dbh->{RaiseError} = 1;
	return $dbh;
}

## gets grouping information for each sample and returns a hash
## of group->sampleID
sub getSampleMap{
	my %args = @_;
	my $exp = $args{exp};
	my $sh = {};
	my $dbh = dbConnect();
	my $ps = qq{
		select `group` , sampleID from samples where `experiment` = "$exp"
	};
	print $ps."\n";
	my $sth = $dbh->prepare($ps);
	$sth->execute();
	while ( my @row = $sth->fetchrow_array ){
		print "@row\n";
		print "row0 ".$row[0]."\n";
		print "row1 ".$row[1]."\n";		
		$sh->{$row[0]}->{$row[1]} = 1;
	}
	return $sh;
}

## returns a hash of the information in the bamstat analysis for an
## experiment
sub getBamStat{
	my %args = @_;
	my $exp = $args{exp};
	my $dbh = dbConnect();
	my $rethash = {};
	my @flds = ("Total_records", "Optical-PCR_duplicate", "Non_Primary_Hits", "Unmapped_Reads",
		 "mapqLTmapq_cut_non_unique", "mapqGTEmapq_cut_unique", "Read-1", "Read-2", "Reads_mapped_to_+",
		 "Reads_mapped_to_-", "Non-splice_reads", "Splice_reads", "Reads_mapped_in_proper_pairs"
		);
	print @flds."\n";
	$rethash->{fields} = ["Total_records", "Optical-PCR_duplicate", "Non_Primary_Hits", "Unmapped_Reads",
		 "mapqLTmapq_cut_non_unique", "mapqGTEmapq_cut_unique", "Read-1", "Read-2", "Reads_mapped_to_+",
		 "Reads_mapped_to_-", "Non-splice_reads", "Splice_reads", "Reads_mapped_in_proper_pairs"
		];
	my $ps = "select `sampleID`, ";
	for my $fld (@flds){
		$ps .= "`$fld`";
		if ($fld  ne $flds[-1]){
			$ps .= ", ";
		}
	}
	$ps .= " from samples natural join bamstat where `experiment` = '$exp' ";

	print $ps."\n";

	my $sth = $dbh->prepare($ps);

	$sth->execute();
	my $totfld = "Total_records";
	while (my $row = $sth->fetchrow_hashref ){
		# print "rowID ". $row->{sampleID}."\n";
		for my $fld (@flds){
			# print "$fld\t".$row->{$fld}."\n";
			if ($fld eq $totfld){
				$rethash->{samples}->{ $row->{sampleID} }->{ $fld } = $row->{$fld} ;
			}	else {
				my $val = $row->{$fld} / $row->{$totfld};
				$rethash->{samples}->{ $row->{sampleID} }->{ $fld } = $val ;
			}
		}
	}
	return $rethash;
} 

## returns hash_ref of the information the the markDup table
sub getMarkDup{
	my %args = @_;
	my $exp = $args{exp};
	my $dbh = dbConnect();
	my $rethash = {};
	my @flds = ("UNPAIRED_READS_EXAMINED", "READ_PAIRS_EXAMINED", "UNPAIRED_READ_DUPLICATES", 
			"READ_PAIR_DUPLICATES", "READ_PAIR_OPTICAL_DUPLICATES", "PERCENT_DUPLICATION", "ESTIMATED_LIBRARY_SIZE");
	$rethash->{fields} = ["UNPAIRED_READS_EXAMINED", "READ_PAIRS_EXAMINED", "UNPAIRED_READ_DUPLICATES", 
			"READ_PAIR_DUPLICATES", "READ_PAIR_OPTICAL_DUPLICATES", "PERCENT_DUPLICATION", "ESTIMATED_LIBRARY_SIZE"];
	my $ps = "select `sampleID`, ";
	for my $fld (@flds){
		$ps .= "`$fld`";
		if ($fld  ne $flds[-1]){
			$ps .= ", ";
		}
	}
	$ps .= " from samples natural join markDup where `experiment` = '$exp' ";			
	my $sth = $dbh->prepare($ps);
	my $totfld1 = "READ_PAIRS_EXAMINED";
	my $totfld2 = "UNPAIRED_READS_EXAMINED";
	$sth->execute();
	while (my $row = $sth->fetchrow_hashref ){
		# print "rowID ". $row->{sampleID}."\n";
		for my $fld (@flds){
			# print "$fld\t".$row->{$fld}."\n";
			if ($fld eq $totfld1 || $fld eq 'ESTIMATED_LIBRARY_SIZE' || $fld eq $totfld2 || $fld eq 'PERCENT_DUPLICATION'){
				$rethash->{samples}->{ $row->{sampleID} }->{ $fld } = $row->{$fld} ;
			} else{
				my $val = $row->{$fld} / ($row->{$totfld1} + $row->{$totfld2});
				$rethash->{samples}->{ $row->{sampleID} }->{ $fld } = $val ;
			}
		}
	}
	return $rethash;
}

## add rows to the readdisttibution table 
sub addReadDistribution {
	my %args = @_;
	my $dbh = dbConnect();
	my $sh = $args{sh};
	my $expDir = $args{expdir};
	my $type = $args{type};
	my @samps = keys( %$sh );
	$type = "accepted_hits_marked_dup" unless ( defined ( $type ) );
	my $exp = $samps[0];
	$exp =~ s/(\d+).+/$1/;
	$exp .= "R";
	unless ( defined ($expDir ) ){
		$expDir = "/Data01/gnomex/Analysis/experiment/$exp";
	}
	for my $samp ( keys ( %$sh ) ){
		my $RDTab = $expDir."/Tophat_".$sh->{$samp}->{'Bam Root'}."/$type.bam_QC/read_distribution.txt";
		my $rh = {};
		open IN, "< $RDTab" or die "cannot open $RDTab, $!\n";
		while ( my $line = <IN> ){
			chomp $line;
			if ($line =~ m/^Total/ || $line=~ m/^=/ || $line =~ /^Group/){
				next;
			}
			my @elems = split(/\s+/, $line);
			$elems[0] =~ s/'//;
			print $elems[0]."\t".$elems[3]."\n";
			$rh->{$elems[0]} = $elems[3];
		}
		my $ps = qq{
			insert into readDistribution (bam, sampleID, CDS_Exons, `Exons_5'UTR`, `Exons_3'UTR`,
				Introns, TSS_up_1kb, TSS_up_5kb, TSS_up_10kb, TES_down_1kb, TES_down_5kb,
				TES_down_10kb) values (?,?,?,?,?,?,?,?,?,?,?,?) on duplicate key update
				`CDS_Exons` = ?, `Exons_5'UTR` = ? , `Exons_3'UTR` = ? , `Introns` = ?, `TSS_up_1kb` = ?,
				`TSS_up_5kb` = ? , `TSS_up_10kb` = ? , `TES_down_1kb` = ? , `TES_down_5kb` = ? ,
				`TES_down_10kb` = ? 
		};
		print $ps."\n";
		my $sth = $dbh->prepare($ps);
		$sth->execute( $type, $samp, $rh->{CDS_Exons}, $rh->{'5UTR_Exons'}, $rh->{'3UTR_Exons'},
			$rh->{Introns}, $rh->{TSS_up_1kb}, $rh->{TSS_up_5kb}, $rh->{TSS_up_10kb},
			$rh->{TES_down_1kb}, $rh->{TES_down_5kb}, $rh->{TES_down_10kb}, $rh->{CDS_Exons}, $rh->{'5UTR_Exons'}, $rh->{'3UTR_Exons'},
			$rh->{Introns}, $rh->{TSS_up_1kb}, $rh->{TSS_up_5kb}, $rh->{TSS_up_10kb},
			$rh->{TES_down_1kb}, $rh->{TES_down_5kb}, $rh->{TES_down_10kb}
			);
	}
}

## add rows to bamstat table
sub addBamStat{
	my $dbh = dbConnect();
	my %args = @_;
	my $sh = $args{sh};
	my $expDir = $args{expdir};
	my $type = $args{type};
	my @samps = keys( %$sh );
	$type = "accepted_hits_marked_dup" unless ( defined ( $type ) );
	my $exp = $samps[0];
	$exp =~ s/(\d+).+/$1/;
	$exp .= "R";
	unless ( defined ($expDir ) ){
		$expDir = "/Data01/gnomex/Analysis/experiment/$exp";
	}
	# my @bamStats = glob $expdir."/Tophat*/$type.bam.bam_stat.tab";
	for my $samp ( keys (%$sh) ){
		my $bsTab = $expDir."/Tophat_".$sh->{$samp}->{'Bam Root'}."/$type.bam.bam_stat.tab";
		print "bamStatTab\t".$bsTab."\n";
		my $rh = {};	 # row hash	
		open IN, "< $bsTab" or die "cannot open $bsTab $!\n";
		while (my $line = <IN> ) {
			if( $line =~ /\d+$/){
				## data row
				print $line."\n";
				if ( $line =~ m/Non Primary Hits\s+(\d+)/) {
					$rh->{Non_Primary_Hits} = $1;
				} else {
					my @elems = split ( ":", $line );
					for my $elem ( @elems ){
						chomp($elem);
					}
					chomp($elems[1]);
					$elems[0] =~ s/ /_/g;
					$elems[0] =~ s/\//-/g;
					$elems[0] =~ s/</LT/g;
					$elems[0] =~ s/>=/GTE/g;
					$elems[0] =~ s/'//g;
					$elems[0] =~ s/\(//g;
					$elems[0] =~ s/\)//g;
					$elems[1] =~ s/\s+//g;
					print $elems[0]."\n";
					print $elems[1]."\n";
					# print length($elems[1])."\n";
					$rh->{$elems[0]} = $elems[1];
				}
			}
		}

		my $ps = qq{
			insert into bamstat (bam, sampleID, Total_records, QC_failed, `Optical-PCR_duplicate`, Non_Primary_Hits,
				Unmapped_Reads, mapqLTmapq_cut_non_unique, mapqGTEmapq_cut_unique, `Read-1`, `Read-2`, `Reads_mapped_to_+`,
				`Reads_mapped_to_-`, `Non-splice_reads`, Splice_reads, Reads_mapped_in_proper_pairs )
				VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?) on duplicate key update
				Total_records = ?, QC_failed = ? , `Optical-PCR_duplicate` = ? , Non_Primary_Hits= ? , Unmapped_Reads = ? , mapqLTmapq_cut_non_unique = ? ,
				mapqGTEmapq_cut_unique = ? , `Read-1` = ? , `Read-2` = ?, `Reads_mapped_to_+` = ? , `Reads_mapped_to_-` = ?, `Non-splice_reads`=?, Splice_reads = ?,
				Reads_mapped_in_proper_pairs=? 
		};
		print $ps."\n";
		my $sth = $dbh->prepare($ps);
		$sth->execute( $type, $samp, $rh->{Total_Records}, $rh->{QC_failed}, $rh->{'Optical-PCR_duplicate'}, $rh->{Non_Primary_Hits}, $rh->{Unmapped_reads},
			$rh->{'mapq_LT_mapq_cut_non-unique'}, $rh->{mapq_GTE_mapq_cut_unique}, $rh->{'Read-1'}, $rh->{'Read-2'}, $rh->{'Reads_map_to_+'},
			$rh->{'Reads_map_to_-'}, $rh->{'Non-splice_reads'}, $rh->{'Splice_reads'}, $rh->{'Reads_mapped_in_proper_pairs'}, 
			$rh->{Total_records}, $rh->{QC_failed}, $rh->{'Optical-PCR_duplicate'}, $rh->{Non_Primary_Hits}, $rh->{Unmapped_reads},
			$rh->{'mapq_LT_mapq_cut_non-unique'}, $rh->{mapq_GTE_mapq_cut_unique}, $rh->{'Read-1'}, $rh->{'Read-2'}, $rh->{'Reads_map_to_+'},
			$rh->{'Reads_map_to_-'}, $rh->{'Non-splice_reads'}, $rh->{'Splice_reads'}, $rh->{'Reads_mapped_in_proper_pairs'}
		);
	}
}



## add rows to mark dup table
sub addMarkDup{
	my $dbh = dbConnect();
	my %args = @_;
	my $sh = $args{sh};
	my $expDir = $args{expDir};
	my @samps = keys( %$sh );
	my $exp = $samps[0];
	my $type = $args{type};
	$type = "accepted_hits_marked_dup" unless ( defined ( $type ) );
	$exp =~ s/(\d+).+/$1/;
	$exp .= "R";
	unless ( defined ($expDir ) ){
		$expDir = "/Data01/gnomex/Analysis/experiment/$exp";
	}
	my @markDupFiles = glob $expDir."/Tophat*/*dup_metrics.txt";
	for my $samp (@samps){
		my $dupTab = $expDir."/Tophat_".$sh->{$samp}->{'Bam Root'}."/accepted_hits_marked_dup.bam_mark_dup_metrics.txt";
		print "dupTab = $dupTab\n";
		if (! $dupTab ){
			die "Cannot find dupe metrics file for sample $samp\n";		
		} else {
			my @header = split('\t', `grep LIBRARY $dupTab`);
			my @data = split('\t', `grep Unknown $dupTab`);
			my $ps = qq{
				INSERT INTO markDup (bam, sampleID, $header[1], $header[2], $header[3], $header[4], $header[5], $header[6], $header[7], $header[8]) 
				VALUES (?,?,?,?,?,?,?,?,?,?) on duplicate key update
				$header[1]=?, 
				$header[2]=?,
				$header[3]=?,
				$header[4]=?, 
				$header[5]=?,
				$header[6]=?, 
				$header[7]=?, 
				$header[8]=?
			};
			print $ps;
			my $sth = $dbh->prepare($ps);
			$sth->execute( $type, $samp, $data[1], $data[2], $data[3], $data[4], $data[5], $data[6], $data[7], $data[8],$data[1], $data[2], $data[3], $data[4], $data[5], $data[6], $data[7], $data[8]);
		}
	}
}


## add rows to b2b bam table
sub addBams{
	my $dbh = dbConnect();
	my %args = @_;
	my $sh = $args{sh}; ##
	my $type = $args{type};

	for my $samp ( keys %$sh ){
		my $ps = qq {
			insert into bam (bam, sampleID)
			VALUES
			(?,?) on duplicate key update bam =?, sampleID=?
		};
		my $sth = $dbh->prepare($ps);
		$sth->execute($type , $samp,$type , $samp);
	}
}


## add rows to db sample table
sub addSamples{
	my $dbh = dbConnect();
	my %args = @_;
	my $sh = $args{sh};  ## sample hash object
	my $fastqDir= $args{fastqdir};
	for my $samp ( keys ( %$sh ) ){
		my $sampID = $samp;
		my $species = $sh->{$samp}->{Organism};
		my $files = $sh->{$samp}->{"Associated Files"};
		print "Files\t$files\n";
		my @files = split(" ", $files);
		my ($read1, $read2);
		if ( @files == 2 ){
			$files[0] =~ s/[,;]//;
			print "file1\t".$files[0]."\n";
			$files[1] =~ s/[,;]//;
			print "file2\t".$files[1]."\n";
			$read1 = `find $fastqDir -name $files[0]*`;
			$read2 = `find $fastqDir -name $files[1]*`;
		} else {
			$read1 = `find $fastqDir -name $files[0]*`;
			$read2 = "NULL";
		}
		print "read1\t".$read1."\n";
		print "read2\t$read2\n";
		my $exp = $samp;
		$exp =~ s/(\d+).+/$1/;
		$exp .= "R";
		my $bamroot = $sh->{$samp}->{'Bam Root'};
		chomp($bamroot);
		my $group = $bamroot;
		$group =~ s/(.+)--\d+X\d+$/$1/;
		print "group\t$group\n";
		print "experiment : $exp\n";
		my $ps = qq{
		INSERT INTO samples 
		(`sampleID`, `read1`, `read2`, `experiment`, `species`, `bam_root`, `group`)
		VALUES
		(?,?,?,?,?,?,?)
		on duplicate key update 
		read1 = ?,
		read2 = ?,
		experiment = ?,
		species = ?,
		bam_root = ? ,
		`group` = ?
		};
		print "prepared statement: ".$ps."\n";
		my $sth = $dbh->prepare($ps);
		$sth->execute($sampID, $read1, $read2, $exp, $species, $bamroot, $group, $read1, $read2, $exp, $species, $bamroot, $group);
		# $sth->execute($sampID, $files[0], defined $files[1] ? $files[1] : 0 , $exp, $species, $bamroot, $group, $read1, $read2, $exp, $species, $bamroot, $group);
	}
}

## database build command
sub buildDatabase{
	my %args = @_;
	my $dbh = $args{dbh};
	# -- Create Table: bam
	# --------------------------------------------------------------------------------
	my $sth = $dbh->prepare(qq{
	CREATE TABLE IF NOT EXISTS bam
	(
	`ID` INT NOT NULL AUTO_INCREMENT
	,`bam` VARCHAR(250) NOT NULL 
	,`sampleID` VARCHAR(250) NOT NULL 
	, PRIMARY KEY (bam, sampleID)
	, KEY (ID)
	)
	ENGINE=INNODB});
	$sth->execute();

	# -- Create Table: samples
	# --------------------------------------------------------------------------------
	$sth = $dbh->prepare(qq{
	CREATE TABLE IF NOT EXISTS samples
	(
		`ID` INT NOT NULL AUTO_INCREMENT
		,`sampleID` VARCHAR(250) NOT NULL 
		,PRIMARY KEY (sampleID)
		,`read1` VARCHAR(250) NOT NULL 
		,`read2` VARCHAR(250)  NULL 
		,`experiment` VARCHAR(250) NOT NULL 
		,`species` VARCHAR(250) NOT NULL 
		, `bam_root` VARCHAR(250) NOT NULL
		, `group` VARCHAR(250) NOT NULL
		, KEY (ID)
	)
	ENGINE=INNODB
	});
	$sth->execute();
	# -- Create Table: bamstat
	# --------------------------------------------------------------------------------
	$sth = $dbh->prepare(qq{
	CREATE TABLE IF NOT EXISTS bamstat
	(
	`ID` INT NOT NULL AUTO_INCREMENT
	, KEY(ID)
	,`bam` VARCHAR(250) NOT NULL 
	,`sampleID` VARCHAR(250) NOT NULL
	, PRIMARY KEY (bam, sampleID)
	,`Total_records` INT NOT NULL 
	,`QC_failed` INT NOT NULL 
	,`Optical-PCR_duplicate` INT NOT NULL 
	,`Non_Primary_Hits` INT NOT NULL 
	,`Unmapped_Reads` INT NOT NULL 
	,`mapqLTmapq_cut_non_unique` INT NOT NULL 
	,`mapqGTEmapq_cut_unique` INT NOT NULL 
	,`Read-1` INT NOT NULL 
	,`Read-2` INT NULL 
	,`Reads_mapped_to_+` INT NOT NULL 
	,`Reads_mapped_to_-` INT NOT NULL 
	,`Non-splice_reads` INT NOT NULL 
	,`Splice_reads` INT NOT NULL 
	,`Reads_mapped_in_proper_pairs` INT NOT NULL 
	,`Proper-paired_reads_map_to_different_chromosome` INT  NULL 

	)
	ENGINE=INNODB
	});
	$sth->execute();
	# -- Create Table: markDup
	# --------------------------------------------------------------------------------
	$sth = $dbh->prepare(qq{
	CREATE TABLE IF NOT EXISTS markDup
	(
		`ID` INT  NOT NULL AUTO_INCREMENT
		, KEY(ID)
		,`bam` VARCHAR(250) NOT NULL 
		, `sampleID` VARCHAR(250) NOT NULL
		, PRIMARY KEY (bam, sampleID)
		,`UNPAIRED_READS_EXAMINED` INT NOT NULL 
		,`READ_PAIRS_EXAMINED` INT NOT NULL 
		,`UNMAPPED_READS` INT NOT NULL 
		,`UNPAIRED_READ_DUPLICATES` INT NOT NULL 
		,`READ_PAIR_DUPLICATES` INT NOT NULL 
		,`READ_PAIR_OPTICAL_DUPLICATES` INT NOT NULL 
		,`PERCENT_DUPLICATION` FLOAT NOT NULL 
		,`ESTIMATED_LIBRARY_SIZE` INT NOT NULL 
	)
	ENGINE=INNODB
	});
	$sth->execute();
	# create Table :readDistrib
	$sth = $dbh->prepare(qq{
		CREATE TABLE IF NOT EXISTS readDistribution
		(
			`ID` INT NOT NULL AUTO_INCREMENT
			, KEY(ID)
			, `bam` VARCHAR(250) NOT NULL
			, `sampleID` VARCHAR(250) NOT NULL
			, PRIMARY KEY (bam, sampleID)
			, `CDS_Exons` float NOT NULL
			, `Exons_5'UTR` float NOT NULL
			, `Exons_3'UTR` float NOT NULL
			, `Introns` float NOT NULL
			, `TSS_up_1kb` float NOT NULL
			, `TSS_up_5kb` float NOT NULL
			, `TSS_up_10kb` float NOT NULL
			, `TES_down_1kb` float NOT NULL
			, `TES_down_5kb` float NOT NULL
			, `TES_down_10kb` float NOT NULL
			)
			ENGINE=INNODB
	});

	# -- Create Foreign Key: bam.sampleID -> samples.sampleID
	$sth = $dbh->prepare(qq{ALTER TABLE bam ADD FOREIGN KEY (sampleID) REFERENCES samples(sampleID)});
	$sth->execute();

	# -- Create Foreign Key: bamstat.bam -> bam.bam
	$sth = $dbh->prepare(qq{ALTER TABLE bamstat ADD FOREIGN KEY (bam, sampleID ) REFERENCES bam(bam, sampleID)});
	$sth->execute();
	# $sth = $dbh->prepare(qq{ALTER TABLE bamstat ADD FOREIGN KEY (sampleID ) REFERENCES bam(bam)});
	# $sth->execute();
	# -- Create Foreign Key: markDup.bam -> bam.bam
	$sth = $dbh->prepare(qq{ALTER TABLE markDup ADD FOREIGN KEY (bam, sampleID) REFERENCES bam(bam, sampleID)});
	$sth->execute();
	# $sth = $dbh->prepare(qq{ALTER TABLE markDup ADD FOREIGN KEY (sampleID) REFERENCES bam(bam)});
	# $sth->execute();


	print "tables built!\n";
}



1;


