package b2b::runval;
use strict;
use Getopt::Long;
use IO::Handle;

## takes the fastqDir, the analysis dir and the filehandle to the runlog file as 
## inputs and checks and validates that the fastqc dir has been created and that 
## directories and zip files were created for each fastq file and writes the results of theis to 
## the runlog file.
sub checkFastQC{
	my %args = @_;
	my $RUNLOG = $args{RUNLOG};
	my $fastqDir = $args{fastqDir};
	my $analysisDir = $args{analysisDir};
	my $fastqcDir = $analysisDir."fastqc/";
	if (-d $fastqcDir){
		print $RUNLOG "$fastqcDir created successfully\n";
	} else {
		print $RUNLOG "$fastqcDir was not created -- Error\n";
		# die "$fastqcDir was not created";
	}
	my @files = <$fastqDir>;
	for my $file (@files){
		chomp($file);
		my $fileDir = $fastqcDir.$file;
		my $zipfile = $fastqcDir.$file."_fastqc.zip";
		if( -e $zipfile){
			my $size = -s $zipfile;
			print $RUNLOG $zipfile." created successfully\n";
			print $RUNLOG "size\t$size\n\n";
		} else {
			print $RUNLOG $fastqcDir.$file."_fastqc.zip not created -- Error\n";
		}
		if ( -d $fileDir ){
			print $RUNLOG "$fileDir created successfully\n";
		} else {
			print $RUNLOG "$fileDir could not be created -- Error\n";
		}
	## this could be expanded to look for each file but the .zip and folder should be sufficient
	## indicia of a successful run.
	}
}

## takes the sample hash, the analysis dir and the filehandle to the runlog as inputs
## and writes validations for successful tophat runs to the the runlog file
sub checkTophat{
	my %args = @_;
	my $sampleHash = $args{sampleHash};
	my $analysisDir = $args{analysisDir};
	my $RUNLOG = $args{RUNLOG};
	for my $key ( keys (%$sampleHash) ){
		my $outpath = $analysisDir."Tophat_".${$sampleHash}{$key}{'Bam Root'};
		checkIfDirExists( file=>$outpath , RUNLOG=>$RUNLOG);
		checkIfExists( file=>"$outpath/accepted_hits.bam" , RUNLOG=>$RUNLOG);
	}
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
	my $outfolder = $args{outfolder};
	my $RUNLOG = $args{RUNLOG};
	my @files = ( "gene_body_coverage.geneBodyCoverage.pdf", "inner_distance.inner_distance.txt", "junction_saturation.junctionSaturation_plot.r",  "read_GC.GC_plot.pdf",    "read_quality.qual.boxplot.pdf",   "RPKM_saturation.saturation.r",
"gene_body_coverage.geneBodyCoverage_plot.r",  "junction_annotation.junction_plot.r", "read_distribution.txt", "read_GC.GC_plot.r", "read_quality.qual.heatmap.pdf",
"gene_body_coverage.geneBodyCoverage.txt", "junction_annotation.junction.xls", "read_duplication.DupRate_plot.pdf", "read_GC.GC.xls", "read_quality.qual.r",
"inner_distance.inner_distance_freq.txt", "junction_annotation.splice_events.pdf", "read_duplication.DupRate_plot.r", "read_NVC.NVC_plot.pdf", "RPKM_saturation.eRPKM.xls",
"inner_distance.inner_distance_plot.pdf", "junction_annotation.splice_junction.pdf", "read_duplication.pos.DupRate.xls", "read_NVC.NVC_plot.r", "RPKM_saturation.rawCount.xls",
"inner_distance.inner_distance_plot.r", "junction_saturation.junctionSaturation_plot.pdf", "read_duplication.seq.DupRate.xls", "read_NVC.NVC.xls", "RPKM_saturation.saturation.pdf");
	for my $file (@files) {
		my $checkFile = $outfolder.$file; 
		checkIfExists( RUNLOG=>$RUNLOG , file=>$checkFile);
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