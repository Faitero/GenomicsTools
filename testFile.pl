
for my $file (@markedFiles){	
	print "processing\t$file";
	## make folder for each file
	$file = filterBam($file);
	print "processed file:\t$file";
	
	buildBamIndex($file);

	my $baiFile = $file;
	$baiFile =~ s/.bam$/.bai/;

	print $RUNLOG "\n\nBAM index file validation\n";
	b2b:report::checkIfExists(file=>$baiFile, RUNLOG=>$RUNLOG);

	runQC({
		file => $file,
		refgene => $refgene,
	});

	print $RUNLOG "\n\nrunQC VALIDATION for $file\n";
	b2b::report::checkRunQC(
		outfolder => $file."_QC/",
		RUNLOG => $RUNLOG,
		);
}

## not that useful
# b2b::qc::catQCtab() unless (($noCatTab && $noBamStat) || $dry || $repOnly);	





print $RUNLOG "\n\nQUALITY SUMMARY TABLE VALIDATION\n";
b2b::report::checkIfExists( file=>"quality_summary.tab" , RUNLOG=>$RUNLOG );


TRACK:
## make Genome tracks

my @files = glob("${analysisDir}Tophat*/*processed.bam");

print $RUNLOG "\n\nBROWSER TRACK VALIDATION\n";
for my $file (@files){
	my $trackGenCommand = "convert_SAM_or_BAM_for_Genome_Browser.pl --nosort $file";
	print "$trackGenCommand\n";
	b2b::tools::runAndLog($trackGenCommand) unless ($dry);
	b2b::report::checkBrowserTracks( file => $file , RUNLOG => $RUNLOG );
}

print "Finished making browser tracks\n";

## check for track only
if ($TRACKonly) {
	
	goto FINISH;
}
