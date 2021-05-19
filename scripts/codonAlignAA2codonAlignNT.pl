#!/usr/bin/perl -w

#########################################################################################################
# Program: codonAlignAA2codonAlignNT.pl
# Purpose: from the original codon aligned NT alignment and its corresponding AA alignment, if AA alignment
# changes, output the adjusted codon NT alignment
# Input: codon NT and AA alignment fasta files, and updated AA alignment fasta file
# Output: updated NT codon alignment fasta file
# Author: Wenjie Deng
# Date: 2020-12-17
###########################################################################################################

use strict;
use warnings;
use v5.10;
use Getopt::Long;

my $usage = "\nusage: codonAlignAA2codonAlignNT.pl inntfile inaafile inupdateaafile outupdatentfile\n";
my $inntfile = shift or die $usage;
my $inaafile = shift or die $usage;
my $inupdateaafile = shift or die $usage;
my $outupdatentfile = shift or die $usage;
my $flag = my $aaseqcount = my $ntseqcount = my $ntlen = my $aalen = my $updateaaseqcount = my $updateaalen = 0;
my $seqName = my $aaseq = my $ntseq = "";
my (%nameaas, %naments, @seqNames, %namestatus, %namecodonaas, %namecodonnts, %nameupdateaas);

open AA, $inaafile or die "couldn't open $inaafile: $!\n";
while (my $line = <AA>) {
	$line =~ s/\R$//;
	next if ($line =~ /^\s*$/);
	if ($line =~ /^>(\S+)/) {
		if ($aaseqcount) {
			my $seqLen = length $aaseq;
			if ($aaseqcount == 1) {
				$aalen = $seqLen;
			}else {
				if ($seqLen != $aalen) {
					die "error: length of sequences are not same ($seqLen vs $aalen), your sequences are probably not aligned.\n";
				}
			}
			@{$nameaas{$seqName}} = split //, $aaseq;
		}
		$seqName = $1;
		if (!$namestatus{$seqName}) {
			$namestatus{$seqName} = 1;
			push @seqNames, $seqName;
		}else {
			die "duplicate name: $seqName\n";
		}		
		$aaseq = "";	
		++$aaseqcount;	
	}else {
		$line =~ s/\#/X/g;
		$aaseq .= uc $line;
	}
}
my $seqlen = length $aaseq;
if ($seqlen != $aalen) {
	die "error: length of sequences are not equal ($seqlen vs $aalen), your sequences are probably not aligned.\n";
}
@{$nameaas{$seqName}} = split //, $aaseq;
close AA;

$seqName = "";
open NT, $inntfile or die "couldn't open $inntfile: $!\n";
while (my $line = <NT>) {
	$line =~ s/\R$//;
	next if ($line =~ /^\s*$/);
	if ($line =~ /^>(\S+)/) {
		if ($ntseqcount) {
			my $seqLen = length $ntseq;
			if ($ntseqcount == 1) {
				$ntlen = $seqLen;
			}else {
				if ($seqLen != $ntlen) {
					die "error: length of sequences are not same ($seqLen vs $ntlen), your sequences are probably not aligned.\n";
				}
			}
			@{$naments{$seqName}} = split //, $ntseq;
		}
		$seqName = $1;
		if (!$namestatus{$seqName}) {
			die "No $seqName in input AA file\n";
		}
		$ntseq = "";	
		++$ntseqcount;	
	}else {
		$ntseq .= uc $line;
	}
}
my $seqLen = length $ntseq;
if ($seqLen != $ntlen) {
	die "error: length of sequences are not equal ($seqLen vs $ntlen), your sequences are probably not aligned.\n";
}
@{$naments{$seqName}} = split //, $ntseq;
close NT;

if ($ntseqcount != $aaseqcount) {
	die "sequence count in NT and AA alignments are not same, probably not corresponded\n";
}

if ($ntlen != $aalen * 3) {
	die "alignment length in NT and AA alignemnts are not corresponed, $ntlen vs $aalen\n";
}

foreach my $name (@seqNames) {
	for (my $i = 0; $i < $aalen; $i++) {
		if ($nameaas{$name}[$i] ne "-") {
			push @{$namecodonaas{$name}}, $nameaas{$name}[$i];
		}
	}
	for (my $i = 0; $i < $ntlen; $i += 3) {
		my $codon = $naments{$name}[$i].$naments{$name}[$i+1].$naments{$name}[$i+2];
		if ($codon ne "---") {
			push @{$namecodonnts{$name}}, $codon;
		}
	}
	if (scalar @{$namecodonaas{$name}} != scalar @{$namecodonnts{$name}}) {
		die "Codon AA and Nt are not corresponded\n";
	}
}

open AA, $inupdateaafile or die "couldn't open $inupdateaafile: $!\n";
while (my $line = <AA>) {
	$line =~ s/\R$//;
	if ($line =~ /^>(\S+)/) {
		if ($updateaaseqcount) {
			my $seqLen = length $aaseq;
			if ($updateaaseqcount == 1) {
				$updateaalen = $seqLen;
			}else {
				if ($seqLen != $updateaalen) {
					die "error: length of sequences are not same in updated AA sequences ($seqLen vs $updateaalen), your sequences are probably not aligned in $inupdateaafile.\n";
				}
			}
			@{$nameupdateaas{$seqName}} = split //, $aaseq;
		}
		$seqName = $1;
		if (!$namestatus{$seqName}) {
			die "No $seqName in input AA file $inaafile\n";
		}		
		$aaseq = "";	
		++$updateaaseqcount;	
	}else {
		$line =~ s/\#/X/g;
		$aaseq .= uc $line;
	}
}
if ($updateaalen != length $aaseq) {
	die "error: length of sequences are not equal (length $aaseq vs $aalen), your sequences are probably not aligned in $inupdateaafile.\n";
}
@{$nameupdateaas{$seqName}} = split //, $aaseq;
close AA;

open OUT, ">", $outupdatentfile or die "couldn't open $outupdatentfile: $!\n";
foreach my $name (@seqNames) {
	my $idx = 0;
	my $ntseq = "";
	if ($nameupdateaas{$name}) {
		foreach my $aa (@{$nameupdateaas{$name}}) {
			if ($aa eq "-") {
				$ntseq .= "---";
			}else {
				if ($aa eq $namecodonaas{$name}[$idx]) {
					$ntseq .= $namecodonnts{$name}[$idx];
					++$idx;
				}else {
					die "AA are not same. idx: $idx, aa: $aa, namecodonaas{$name}[$idx]: $namecodonaas{$name}[$idx]\n";
				}				
			}
		}
		if (length $ntseq != $updateaalen * 3) {
			die "updated NT sequence length not correponding to updated AA sequence lenght of $updateaalen\n";
		}
		my $nogapsntseq = join('', @{$naments{$name}});
		$nogapsntseq =~ s/\-//g;
		my $nogapsupdatentseq = $ntseq;
		$nogapsupdatentseq =~ s/\-//g;
		if ($nogapsntseq eq $nogapsupdatentseq) {
			print OUT ">$name\n$ntseq\n";
		}else {
			die "updated NT sequence is not same as input NT sequence for $name: $nogapsntseq\n$nogapsupdatentseq\n";
		}		
	}else {
		print "* Missing $name in updated AA alignment file *\n"
	}
}
close OUT;

print "total $ntseqcount sequences, AA alignment length $aalen, NT alignment length $ntlen. $updateaaseqcount updated AA sequences, AA alignment length: $updateaalen\n";


