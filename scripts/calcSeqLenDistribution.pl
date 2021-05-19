#!/usr/bin/perl

use strict;
use warnings;
use v5.10;

my $usage = "usage: perl calcSeqLenDistribution.pl inputFastaFile outlenghtfile outsummaryfile\n";

my $infile = shift or die $usage;
my $outlengthfile = shift or die $usage;
my $outsummaryfile = shift or die $usage;
my $count = my $nosgacount = my $bigdeletecount = 0;
my $name = "";
my (@names, %nameseq, @lens);
print "\n=== processing $infile ===\n";
open IN, $infile or die "couldn't open $infile: $!\n";
while (my $line = <IN>) {
	$line =~ s/\R$//;
	next if $line =~ /^\s*$/;
	if ($line =~ /^>(\S+)/) {
		$name = $1;
		unless ($name =~ /HXB2/) {
			push @names, $name;
			++$nosgacount;
		}
		++$count;
	}else {
		$line =~ s/\-//g;
		$line =~ s/\*//;
		$nameseq{$name} .= $line;
	}
}
close IN;

foreach my $name (@names) {
	my $seq = $nameseq{$name};
	die "$name has * in sequence\n" if ($seq =~ /\*/);
	my $len = length($nameseq{$name});
	push @lens, $len;
}

my @ids = split /_/, $infile;
my $sampleid = $ids[0]."_".$ids[1]."_".$ids[2]."_".$ids[3]."_".$ids[4];
@lens = sort{$a <=> $b} @lens;
open OUT, ">>", $outlengthfile or die "could't open $outlengthfile: $!\n";
foreach my $len (@lens) {
	print OUT "$sampleid,$len\n";
}
close OUT;

my @stats = DoStatistics(@lens);
my $median = $stats[4];

foreach my $len (@lens) {
	if ($len < $median * 0.8) {
		++$bigdeletecount;
	}
}

open OUT, ">>", $outsummaryfile or die "could't open $outsummaryfile: $!\n";
my $fraction = int($bigdeletecount / $nosgacount * 10000 + 0.5) / 10000;
print OUT "$sampleid,$nosgacount,$bigdeletecount,$fraction,",join(",", @stats),"\n";
close OUT;

print "total $count sequences, $nosgacount no SGA.\n";


sub DoStatistics {
	my @diverArray = @_;	
	my ($avg, $std, $sem, $median, $q1, $q3, $min, $max, $arraySize);
	my @resultArray = ();
	$arraySize = scalar @diverArray;
	
	# calculate average
	my $diver = 0;
	foreach (@diverArray) {
		$diver += $_;
	}
	$avg = int($diver / $arraySize + 0.5);
	
	# calculate standard deviation
	my $squareSum = 0;
	foreach (@diverArray) {
		$squareSum += ($_ - $avg)*($_ - $avg);
	}
	if($arraySize == 1) {
		$std = 0;
	}else {
		$std = sqrt($squareSum / ($arraySize - 1));
	}
	$sem = $std / sqrt($arraySize);
	$std = int($std * 10000 + 0.5) / 10000;
	
	# calculate median, first quartile, third quartile
	if($arraySize == 1) {
		$median = $min = $q1 = $max = $q3 = $avg;
	}elsif($arraySize == 2) {
		$median = $avg;
		$min = $q1 = $diverArray[0];
		$max = $q3 = $diverArray[1];
	}elsif($arraySize == 3) {
		$median = $diverArray[1];
		$min = $q1 = $diverArray[0];
		$max = $q3 = $diverArray[2];
	}else {
		my $middle = int($arraySize / 2);
		my @firstHalf;
		my @secondHalf;
		my $secondHalfStart;
		
		if($arraySize % 2 == 0) {
			$secondHalfStart = $middle;					
		}else {
			$secondHalfStart = $middle + 1;
		}
		
		for(my $i = 0; $i < $middle; $i++) {
			push(@firstHalf, $diverArray[$i]);
		}
		for(my $i = $secondHalfStart; $i < $arraySize; $i++) {
			push(@secondHalf, $diverArray[$i]);
		}
		
		$median = GetMedian(@diverArray);
		$q1 = GetMedian(@firstHalf);
		$q3 = GetMedian(@secondHalf);
		$min = $diverArray[0];
		$max = $diverArray[$arraySize-1];
	}
	return ($avg, $std, $min, $q1, $median, $q3, $max);
}

sub GetMedian {
	my @inArray = @_;
	my $arraySize = scalar @inArray;
	my $median;
	my $middle = int($arraySize / 2);
	
	if($arraySize % 2 == 0) {
		$median = ($inArray[$middle-1] + $inArray[$middle]) / 2;
	}else {
		$median = $inArray[$middle];
	}
	$median = int($median + 0.5);
	return $median;
}
