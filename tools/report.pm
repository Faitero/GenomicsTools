package tools::report;
use strict;
use Data::Dumper;
use IO::Handle;
use File::Basename;
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
	# my ($samps) = @_;
	# print Dumper($samps);
	my @samps = @_;#$samps;
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
	my $sh = tools::db::getSampleMap(exp => $exp);
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
	print "sh dump\t".Dumper($sh);
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
	print "fields\t".Dumper(@flds);
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
	my $sh = tools::db::getBamStat(exp => $exp);
	my $html = makeTableString( title => "Bam Statistics", ref =>"RSeQC Bamstat", sh =>$sh );	
	print $fh $html;
}

sub writeMarkDup{
	my %args = @_;
	my $exp = $args{exp};
	my $fh = $args{fh};
	my $sh = tools::db::getMarkDup(exp => $exp);
	my $html = makeTableString( title => "Duplication Statistics" ,  ref => "PicardTools MarkDuplicates", sh=> $sh);
	print $fh $html;
}

sub writeFigTables{
	my %args = @_;
	my $analysisDir = $args{analysisDir};
	my $fh = $args{fh};
	my $sampleHash = $args{sampleHash}; 
	my $repdir = $args{repdir};
	my $qcfile;
	my $plotdir = $analysisDir.$repdir."plots/";
	
	system ("mkdir $plotdir");
	my @files = ( "gene_body_coverage.geneBodyCoverage.pdf", "inner_distance.inner_distance.txt", "junction_saturation.junctionSaturation_plot.r",  "read_GC.GC_plot.pdf",    "read_quality.qual.boxplot.pdf",   "RPKM_saturation.saturation.r",
"gene_body_coverage.geneBodyCoverage_plot.r",  "junction_annotation.junction_plot.r", "read_distribution.txt", "read_GC.GC_plot.r", 
"gene_body_coverage.geneBodyCoverage.txt", "junction_annotation.junction.xls", "read_duplication.DupRate_plot.pdf", "read_GC.GC.xls", "read_quality.qual.r",
"inner_distance.inner_distance_freq.txt",  "read_duplication.DupRate_plot.r", "read_NVC.NVC_plot.pdf", "RPKM_saturation.eRPKM.xls",
"inner_distance.inner_distance_plot.pdf",  "read_duplication.pos.DupRate.xls", "read_NVC.NVC_plot.r", "RPKM_saturation.rawCount.xls",
"inner_distance.inner_distance_plot.r", "junction_saturation.junctionSaturation_plot.pdf", "read_duplication.seq.DupRate.xls", "read_NVC.NVC.xls", "RPKM_saturation.saturation.pdf");
	for my $id (sortSamps(keys (%$sampleHash))){
		# print "plotdir\t$plotdir\n";
		my $outfolder = getPathFromID($analysisDir,  $sampleHash,$id);
		$outfolder .= "/accepted_hits_marked_dup.bam_QC";
		# print "outfolder\t".$outfolder."\n";
		for my $file (@files) {
			my $checkfile = $outfolder."/".$file; 
			if ($checkfile =~ m/(Tophat.*)$/){
				# print "\$1 = $checkfile\n";
				$qcfile = $1;
				$qcfile =~ s/\//\./g;
				print $plotdir.$qcfile."\n";
				my $cpStr = "cp $checkfile $plotdir$qcfile";
				print $cpStr."\n";
				system ( "cp $checkfile $plotdir$qcfile");
			}
		}
	}
	for my $file (@files) {
		print "$file\n";
		if ($file =~ m/\.pdf$/){
			print "match\n";
			writeFigureTable(
				file => $file,
				sh 	=> $sampleHash,
				fh =>	$fh,
				repdir => $repdir,
				plotdir => $plotdir,
				cols => 3,
				);
		} else {

			print "no match\n";
		}
	}
	
}


sub writeFigureTable{
	print "in writeFigure Table\n";
	my %args = @_;
	my $ft = $args{file};
	my $sh = $args{sh};
	my $fh = $args{fh};
	my $cols = $args{cols};
	my $repdir = $args{repdir};
	my $plotdir = $args{plotdir};
	# print "sh\t".Dumper(%$sh);
	my @samps = sortSamps( keys ( %$sh ) );
	my $title;
	my $tabStr;
	if ( $ft =~ m/^(\w+)\.\w+/ ){
		my $title = $1;
		$title =~ s/_/ /g;
		# print $title."\n";
		$tabStr = "<h3>$title</h3>";
	}
	$tabStr .= "<table id='$ft'>";
	my $i = 1;
	for my $samp (@samps){
		print "samp\t$samp\n";
		# opens row
		if ($i == 1) {
			$tabStr .= "<tr>";
		} 
		my $globStr = "$plotdir*--$samp.*$ft";
		print "globstr\t$globStr\n";
		my @pdf = glob( "$plotdir*--$samp.*$ft" );
		print "pdf\t".$pdf[0]."\n";
		my $pdf = $pdf[0];
		$pdf =~ m/^(.*)(.{3})$/;
		my $base = $1;
		my $jpg = $base."jpg";
		system ("convert $pdf $jpg");
		$tabStr .= "<td><p>Sample $samp</p>";
		my $bn = basename($jpg);
		$tabStr .= "<img src='plots/$bn'></td>";

		# closes row
		if ($i == $cols){
			$i = 1;
			$tabStr .= "</tr>";
		} else {
			$i++;
		}
	}
	$tabStr =~ m/(.{5})$/;
	if ($1 ne "</tr>") {
		$tabStr .= "</tr>";
	}

	$tabStr .= "</table>";
	print $fh $tabStr;
}

sub writeFooter{
	my %args = @_;
	my $fh = $args{fh};
	print $fh qq{
		</body>
		</html>
	};
}

# ## takes the fastqDir, the analysis dir and the filehandle to the runlog file as 
# ## inputs and checks and validates that the fastqc dir has been created and that 
# ## directories and zip files were created for each fastq file and writes the results of theis to 
# ## the runlog file.
# sub checkFastQC{
# 	my %args = @_;
# 	my $RUNLOG = $args{RUNLOG};
# 	my $analysisDir = $args{analysisDir};
# 	my $fastqcDir = $analysisDir."fastqc/";
# 	if (-d $fastqcDir){
# 		print $RUNLOG "$fastqcDir created successfully\n";
# 	} else {
# 		print $RUNLOG "$fastqcDir was not created -- Error\n";
# 		# die "$fastqcDir was not created";
# 	}
# 	print $fastqcDir."*.zip\n";
# 	my @files = glob($fastqcDir."*.zip");
# 	print Dumper(@files);
# 	for my $file (@files){
# 		chomp($file);
# 		print "file = $file\n";
# 		my $fileDir = $fastqcDir.$file;
# 		my $zipfile = $fastqcDir.$file."_fastqc.zip";
# 		if( -e $zipfile){
# 			my $size = -s $zipfile;
# 			print $RUNLOG $zipfile." created successfully\n";
# 			print $RUNLOG "size\t$size\n\n";
# 		} else {
# 			print $RUNLOG $fastqcDir.$file."_fastqc.zip not created -- Error\n";
# 		}
# 		if ( -d $fileDir ){
# 			print $RUNLOG "$fileDir created successfully\n";
# 		} else {
# 			print $RUNLOG "$fileDir could not be created -- Error\n";
# 		}
# 	## this could be expanded to look for each file but the .zip and folder should be sufficient
# 	## indicia of a successful run.
# 	}
# }

## takes the sample hash, the analysis dir and the filehandle to the runlog as inputs
## and writes validations for successful tophat runs to the the runlog file
sub checkTophat{
	my %args = @_;
	my $sampleHash = $args{sampleHash};
	my $analysisDir = $args{analysisDir};
	my $RUNLOG = $args{RUNLOG};
	for my $key ( sortSamps(keys (%$sampleHash)) ){
		my $outpath = getPathFromID($analysisDir,$sampleHash,$key);
		checkIfDirExists( file=>$outpath , RUNLOG=>$RUNLOG);
		checkIfExists( file=>"$outpath/accepted_hits.bam" , RUNLOG=>$RUNLOG);
	}
}

## get the tophat dir from analysis dir and sample id 
sub getPathFromID{
	# my %args = @_;
	my $anDir = shift;
	
	my $sh = shift;
	# print Dumper($sh);
	my $samp = shift;
	print $anDir."Tophat_".$sh->{$samp}->{'Bam Root'}."\n";;
	return $anDir."Tophat_".$sh->{$samp}->{'Bam Root'};
}

## this checks if the file exists and writes if it does to  the logfile
## along with the file size.
sub checkIfExists{
	my %args = @_;
	my $file = $args{file};
	my $RUNLOG = $args{RUNLOG};
	if( -e $file ){
    	my $size = -s $file;
    	print $RUNLOG "$file\tsize\t$size\tcreated successfully\n";
    } else {
    	print $RUNLOG "$file not created -- Error\n";
    }
}

sub checkIfDirExists{
	my %args = @_;
	my $file = $args{file};
	my $RUNLOG = $args{RUNLOG};
	if( -e $file ){
    	my $size = -s $file;
    	print $RUNLOG "$file\tsize\t$size\tcreated successfully\n";
    } else {
    	print $RUNLOG "$file not created -- Error\n";
    }
}

sub checkRunQC {
	my %args = @_;
	my $sampleHash = $args{sampleHash};
	my $analysisDir = $args{analysisDir};
	my $RUNLOG = $args{RUNLOG};
	my @files = ( "gene_body_coverage.geneBodyCoverage.pdf", "inner_distance.inner_distance.txt", "junction_saturation.junctionSaturation_plot.r",  "read_GC.GC_plot.pdf",    "read_quality.qual.boxplot.pdf",   "RPKM_saturation.saturation.r",
"gene_body_coverage.geneBodyCoverage_plot.r",  "junction_annotation.junction_plot.r", "read_distribution.txt", "read_GC.GC_plot.r", "read_quality.qual.heatmap.pdf",
"gene_body_coverage.geneBodyCoverage.txt", "junction_annotation.junction.xls", "read_duplication.DupRate_plot.pdf", "read_GC.GC.xls", "read_quality.qual.r",
"inner_distance.inner_distance_freq.txt", "junction_annotation.splice_events.pdf", "read_duplication.DupRate_plot.r", "read_NVC.NVC_plot.pdf", "RPKM_saturation.eRPKM.xls",
"inner_distance.inner_distance_plot.pdf", "junction_annotation.splice_junction.pdf", "read_duplication.pos.DupRate.xls", "read_NVC.NVC_plot.r", "RPKM_saturation.rawCount.xls",
"inner_distance.inner_distance_plot.r", "junction_saturation.junctionSaturation_plot.pdf", "read_duplication.seq.DupRate.xls", "read_NVC.NVC.xls", "RPKM_saturation.saturation.pdf");
	for my $id (sortSamps(keys (%$sampleHash))){
		my $outfolder = getPathFromID($analysisDir,  $sampleHash,$id);
		$outfolder .= "/accepted_hits_marked_dup.bam_QC";
		for my $file (@files) {
			my $checkFile = $outfolder."/".$file; 
			checkIfExists( RUNLOG=>$RUNLOG , file=>$checkFile);
		}
	}
}

sub checkBrowserTracks{
	my %args = @_;
	my $file = $args{file};
	my $RUNLOG = $args{RUNLOG};
	$file =~ s/\//_/g;
	$file = "Browser.Track.Descriptions.".$file;
	my $bwFile = $file.".scaled_by_1.000.bw";
	checkIfExists( RUNLOG=>$RUNLOG , file=>$bwFile );
}

sub checkDRDS{
	my %args = @_;
	my $sampleHash = $args{sampleHash};
	my $analysisDir = $args{analysisDir};
	my $USEQ_CONDdir = "${analysisDir}USEQ_CONDITIONS/";
}


1;