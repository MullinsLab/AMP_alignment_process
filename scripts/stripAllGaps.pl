#!/usr/bin/perl -w

#########################################################################################################
# Program: stripAllGaps.pl
# Purpose: parse the fasta alignment file, remove the columns which all gaps
# Input: sequence alignment fasta file
# Output: cleaned sequence alignment file
# Author: Wenjie Deng
# Date: 2020-07-15
###########################################################################################################

use strict;
use warnings;
use v5.10;
use Getopt::Long;

my $usage = "\nusage: stripAllGaps.pl infile outfile\n";
my $infile = shift or die $usage;
my $outfile = shift or die $usage;
my $flag = my $totalSeq = my $rmCount = my $len = my $seqcount = 0;
my $seqName = my $seq = "";
my (%nameSeq, @seqNames, %namestatus, %removeSite);

open IN, $infile or die "couldn't open $infile: $!\n";
while (my $line = <IN>) {
	$line =~ s/\R$//;
	next if ($line =~ /^\s*$/);
	if ($line =~ /^>(\S+)/) {
		if ($totalSeq) {
			my $seqLen = length $seq;
			if ($totalSeq == 1) {
				$len = $seqLen;
			}else {
				if ($seqLen != $len) {
					die "error: length of sequences are not same ($seqLen vs $len), your sequences are probably not aligned.\n";
				}
			}
			@{$nameSeq{$seqName}} = split //, $seq;
		}
		$seqName = $1;
		if (!$namestatus{$seqName}) {
			$namestatus{$seqName} = 1;
			unless ($seqName =~ /HXB2/) {
				push @seqNames, $seqName;
				++$seqcount;
			}			
		}else {
			die "duplicate names: $seqName\n";
		}		
		$seq = "";	
		++$totalSeq;	
	}else {
		$seq .= uc $line;
	}
}
my $seqLen = length $seq;
if ($seqLen != $len) {
	die "error: length of sequences are not equal ($seqLen vs $len), your sequences are probably not aligned.\n";
}
@{$nameSeq{$seqName}} = split //, $seq;
close IN;

for (my $i = 0; $i < $len; $i++) {
	my $gapCount = 0;
	foreach my $seqName (@seqNames) {
		if ($nameSeq{$seqName}[$i] eq '-') {
			++$gapCount;
		}
	}
	if ($gapCount == $seqcount) {
		$removeSite{$i} = 1;
		++$rmCount;
	}
}

open OUT, ">", $outfile or die "Couldn't open $outfile: $!\n";
foreach my $seqName (@seqNames) {
	print OUT ">$seqName\n";
	for (my $i = 0; $i < $len; $i++) {
		unless ($removeSite{$i}) {
			print OUT $nameSeq{$seqName}[$i];
		}
	}
	print OUT "\n";
}
close OUT;

print "total $seqcount sequences, alignment length $len. $rmCount all gap columns removed.\n";


