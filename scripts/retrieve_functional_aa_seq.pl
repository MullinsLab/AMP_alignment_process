#!/usr/bin/perl

# from a protein sequence alignment, retrieve functional amino acid sequences (stop codon 
# within 3' end of 5% of median alignment length and < 20% of deletions of median length
# and output aa alignment without defectives.

use strict;
use warnings;
use v5.10;

my $usage = "perl retrieve_functional_aa_seq.pl inAaAlignment outsummaryfile\n";
my $infile = shift or die $usage;
my $summaryfile = shift or die $usage;
my $outfile = my $defectivefile = $infile;
$outfile =~ s/\.fasta/_functional.fasta/;
$defectivefile =~ s/\.fasta/_defective.fasta/;
my $count = my $seqcount = my $hxb2alignlen = my $functionalcount = my $deletioncount = my $prematurecount = my $missing5count = my $maxalignlen = my $notMcount = 0;
my $name = "";
my (@names, %nameseq, @lens, @alignlens, %namestatus);
open IN, $infile or die "couldn't open $infile: $!\n";
while (my $line = <IN>) {
	$line =~ s/\R$//;
	next if $line =~ /^\s*$/;
	if ($line =~ /^>(.*)/) {
		++$count;
		$name = $1;
		push @names, $name;
	}else {
		$line = uc $line;
		$nameseq{$name} .= $line;
	}
}
close IN;

foreach my $name (@names) {
	my $seq = $nameseq{$name};
	if ($name =~ /HXB2/) {
		$hxb2alignlen = length $seq;
	}else {
		unless ($seq =~ /^\-+$/) {
			++$seqcount;
			my $alignlen = length $seq;
			if ($alignlen > $maxalignlen) {
				$maxalignlen = $alignlen;
			}
			push @alignlens, $alignlen;
			$seq =~ s/\-//g;
			my $len = length $seq;
			push @lens, $len;
		}
	}
}
if ($hxb2alignlen > $maxalignlen) {
	$maxalignlen = $hxb2alignlen;
}
my @sortedlens = sort {$a <=> $b} @lens; 
my @sortedalignlens = sort {$a <=> $b} @alignlens; 
my $mediancount = int ($seqcount / 2);
my $medianlen = $sortedlens[$mediancount];
my $medianalignlen = $sortedalignlens[$mediancount];

foreach my $name (@names) {
	my $seq = my $gapstripseq = $nameseq{$name};
	my $seqlen = length $seq;
	$gapstripseq =~ s/\-//g;
	if ($name =~ /HXB2/) {
		$namestatus{$name} = 1;
		for (my $i = $seqlen; $i < $maxalignlen; $i++) {
			$nameseq{$name} .= "-";
		}
	}else {
		if ($seq =~ /^\-/) {
			++$missing5count;
		}elsif ($seq !~ /^M/) {
			++$notMcount;
		}elsif ($seqlen >= $medianalignlen * 0.95) {
			if (length $gapstripseq >= 0.8 * $medianlen) {
				$namestatus{$name} = 1;
				++$functionalcount;
				for (my $i = $seqlen; $i < $maxalignlen; $i++) {
					$nameseq{$name} .= "-";
				}
			}else {
				++$deletioncount;
			}
		}else {
			++$prematurecount;
		}
	}
}

open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
open DFCT, ">", $defectivefile or die "couldn't open $defectivefile: $!\n";
foreach my $name (@names) {
	my $seq = $nameseq{$name};
	if ($namestatus{$name}) {
		print OUT ">$name\n$seq\n";
	}else {
		print DFCT ">$name\n$seq\n";
	}
}
close OUT;
close DFCT;

open SUMMARY, ">>", $summaryfile or die "couldn't open $summaryfile: $!\n";
print SUMMARY "$infile,$seqcount,$functionalcount,$prematurecount,$deletioncount,$missing5count,$notMcount,$medianlen\n";
close SUMMARY;

print "Processing total $seqcount sequences, $functionalcount functional, $missing5count missing 5' end, $notMcount not start with 'M', $deletioncount with big deletions (> 20%), $prematurecount premature stop (< 95%), maxlen: $maxalignlen\n";

