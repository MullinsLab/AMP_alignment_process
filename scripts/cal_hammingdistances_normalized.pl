#!/usr/bin/perl

################################################################################
# Program: cal_hammingdistances_normalized.pl
# Purpose: calculate normalized hamming distance of an alignment
# Normalize hamming distance: 1. pairwise hamming distance divide by the shorter
# sequence length of the pair; 2. sum all normalized pariwise hamming distances
# and then divided by number of pairwise comparison
# Author: Wenjie Deng
# Date: 2020-11-18
################################################################################

use strict;
use warnings;
use v5.10;

my $usage = "usage: perl cal_hammingdistances_normalized.pl inAlignmentFastaFile outputPairwiseHammingDistanceCSVfile\n";
my $inFile = shift or die $usage;
my $outfile = shift or die $usage;
my $seqCount = my $len = 0;
my $seqName = '';
my @seqNames = my %nameSeq = ();
open IN, $inFile or die "couldn't open $inFile: $!\n";
while (my $line = <IN>) {
	$line =~ s/\R$//;
	next if $line =~ /^\s*$/;
	if ($line =~ /^>(\S+)/) {
		if ($seqCount) {
			my $seqLen = length $nameSeq{$seqName};
			if ($seqCount == 1) {
				$len = $seqLen;
			}else {
				if ($seqLen != $len) {
					die "error: length of sequences are not same, your sequences are probably not aligned.\n";
				}
			}
		}
		$seqName = $1;
		push @seqNames, $seqName;
		$seqCount++;
	}else {
		$nameSeq{$seqName} .= uc $line;
	}
}
# last sequence
my $seqLen = length $nameSeq{$seqName};
if ($seqLen != $len) {
	die "error: length of sequences are not same, your sequences are probably not aligned.\n";
}
close IN;

my $distCount = my $nhdistSum = my $nhdistAvg = 0;
open(OUT, ">", $outfile) or die "couldn't Open $outfile: $!\n";
print OUT "Id1,Id2,distance\n";
while (@seqNames) {	
	my $firstname = shift @seqNames;
	my $seq = $nameSeq{$firstname};
	my @seqNas = split //, $seq;
	$seq =~ s/\-//g;
	my $seqlen = length $seq;			
	foreach my $restname (@seqNames) {
		++$distCount;
		my $hdist = 0;
		my $restseq = $nameSeq{$restname};
		my @restNas = split //, $restseq;
		for (my $i = 0; $i < $len; $i++) {
			if ($seqNas[$i] =~ /[A-Z]/ and $restNas[$i] =~ /[A-Z]/) {
				if ($seqNas[$i] ne $restNas[$i]) {
					++$hdist;
				}
			}
		}
		$restseq =~ s/\-//g;
		my $minlen = length $restseq;
		if ($seqlen < $minlen) {
			$minlen = $seqlen;
		}
		my $nhdist = $hdist / $minlen;
		$nhdistSum += $nhdist;
		print OUT "$firstname,$restname,$nhdist\n";
	}
}
if ($distCount) {
	$nhdistAvg = int ($nhdistSum / $distCount * 1000000 + 0.5) / 1000000;
}
print OUT "\nAverage hamming distance,,$nhdistAvg\n";
print "Total $seqCount sequences, alignment length $len bp, average normalized hamming distance for the alignment is $nhdistAvg\n";

