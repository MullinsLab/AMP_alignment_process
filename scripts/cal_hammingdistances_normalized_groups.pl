#!/usr/bin/perl

################################################################################
# Program: cal_hammingdistances_normalized_groups.pl
# Purpose: calculate normalized hamming distance of an alignment of group info in sequence name
# Normalize hamming distance: 1. pairwise hamming distance divide by the shorter
# sequence length of the pair; 2. sum all normalized pariwise hamming distances
# and then divided by number of pairwise comparison
# Outputs: pairwise distance .csv files for all sequences, within and between groups
# Author: Wenjie Deng
# Date: 2020-12-14
################################################################################

use strict;
use warnings;
use v5.10;

my $usage = "usage: perl cal_hammingdistances_normalized_groups.pl inAlignmentFastaFile\n";
my $infile = shift or die $usage;
my $outdir = "hdist_outputs";
my $outfile = $outdir."/".$infile;
$outfile =~ s/\.fasta/_hdist.csv/;
unless (-e $outdir) {
	mkdir $outdir;
}
my $seqCount = my $len = 0;
my $seqName = '';
my (@seqNames, %nameSeq);
print "\n=== Processing $infile ===\n";
open IN, $infile or die "couldn't open $infile: $!\n";
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

my $distCount = my $nhdistSum = my $nhdistAvg = my $groupdistcount = 0;
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
#print OUT "\nAverage hamming distance,,$nhdistAvg\n";
close OUT;

my (%tpsline, @tpss, %tpsstatus);
open IN, "<", $outfile or die "couldn't open $outfile: $!\n";
while (my $line = <IN>) {
	$line =~ s/\R$//;
	next if ($line =~ /^\s*$/ or $line =~ /^Id1,/);
	my ($id1, $id2, $hdist) = split /,/, $line;
	my $tp1 = my $tp2 = "";
	if ($id1 =~ /^(.*?)_\d+_(.*?)_/) {
		$tp1 = $2;
	}else {
		die "id1 not formatted: $id1\n";
	}
	if ($id2 =~ /^(.*?)_\d+_(.*?)_/) {
		$tp2 = $2;
	}else {
		die "id2 not formatted: $id2\n";
	}
	if ($tp1 gt $tp2) {
		my $tmp = $tp1;
		$tp1 = $tp2;
		$tp2 = $tmp;
	}
	my $tps = $tp1."-".$tp2;
	if (!$tpsstatus{$tps}) {
		$tpsstatus{$tps} = 1;
		push @tpss, $tps;
	}
	push @{$tpsline{$tps}}, $line;
}
close IN;

foreach my $tps (@tpss) {
	my $groupoutfile = $outfile;
	$groupoutfile =~ s/\.csv/_$tps.csv/;
	open OUT, ">", $groupoutfile or die "couldn't open $groupoutfile: $!\n";
	print OUT "Id1,Id2,hdist\n";
	foreach my $line (@{$tpsline{$tps}}) {
		print OUT $line,"\n";
		++$groupdistcount;
	}
	close OUT;
}

if ($groupdistcount != $distCount) {
	die "dist count is different: groupdistcount $groupdistcount vs. distCount $distCount\n";
}

print "Total $seqCount sequences, alignment length $len bp, $distCount pairwise comparisons. done.\n";

