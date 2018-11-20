#!/usr/bin/perl

$fastqFile = "myData.fastq";
$fastqOutFile = "out/myData";
$barcodeStartPos = 3; # barcode starts in 3rd position
$barcodeLength = 4; # barcode length

open(FASTQ, '<', $fastqFile) or die "Could not open file $fastqFile";
my %outFiles;


while(my @lines=readRecord(\*FASTQ))
{
	if(scalar @lines == 0) { print "Finished. @lines"; last; }
	if(scalar @lines != 4) { print "Warning: fastq file doesn't end in quartets"; last; }
	
	my $id = $lines[0];
	my $seq = $lines[1];
	my $sep = $lines[2];
	my $qual = $lines[3];
	
	# Barcode starts in $barcodeStartPos position, and has length of $barcodeLength
	my $barcode = substr($seq, $barcodeStartPos, $barcodeLength);
	
	if(!defined($outFiles{$barcode}))
	{
		open($outFiles{$barcode}, ">".$fastqOutFile."_".$barcode.".fastq") || die;
	}
	print {$outFiles{$barcode}} $id;
	print {$outFiles{$barcode}} $seq;
	print {$outFiles{$barcode}} $sep;
	print {$outFiles{$barcode}} $qual;
	
}
close FASTQ;

foreach my $barcode (keys %outFiles)
{
	close $outFiles{$barcode};
}



sub readRecord
{
	my $fh = shift;
	my @lines;
	$i = 0;
	while($i < 4 && (my $line=<$fh>))
	{
		push(@lines, $line);
		$i++;
	}
	return(@lines);
}