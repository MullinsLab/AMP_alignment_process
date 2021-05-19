#!/usr/bin/perl

# translates nucleotide codon alignment into amino acid codon alignment

use strict;
use warnings;
use v5.10;

my $usage = "perl codonAlignNT2codonAlignAA.pl inNtAlignment outAAAlignment\n";
my $infile = shift or die $usage;
my $outfile = shift or die $usage;
my $count = my $gaps = 0;
my $name = "";
my (@names, %nameseq, @codonnts, @codongaps);
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

open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
foreach my $name (@names) {	
	my $seq = $nameseq{$name};
	my $len = length $seq;
	my @nts = split //, $seq;
	my $aaseq = my $partialaaseq = "";
	my $flag = 0;
	my (@codonnts, @codongaps);
	for (my $i = 0; $i < $len; $i += 3) {
		my $codon = $nts[$i].$nts[$i+1].$nts[$i+2];
		my $aa = "";
		if ($codon eq "---" or $codon =~ /[A-Z]{3}/) {
			$aa = translation($codon);
		}else {
			$aa = "X";
		}
		$aaseq .= $aa;
	}
	print OUT ">$name\n$aaseq\n";
}
close OUT;

print "Total $count sequences.\n";


sub translation {
	my $codon = shift;
	my %codon2aa = (
		"---" => "-",
		"ATT" => "I",
		"ATC" => "I",
		"ATA" => "I",
		"CTT" => "L",
		"CTC" => "L",
		"CTA" => "L",
		"CTG" => "L",
		"TTA" => "L",
		"TTG" => "L",
		"GTT" => "V",
		"GTC" => "V",
		"GTA" => "V",
		"GTG" => "V",
		"TTT" => "F",
		"TTC" => "F",
		"ATG" => "M",
		"TGT" => "C",
		"TGC" => "C",
		"GCT" => "A",
		"GCC" => "A",
		"GCA" => "A",
		"GCG" => "A",
		"GGT" => "G",
		"GGC" => "G",
		"GGA" => "G",
		"GGG" => "G",
		"CCT" => "P",
		"CCC" => "P",
		"CCA" => "P",
		"CCG" => "P",
		"ACT" => "T",
		"ACC" => "T",
		"ACA" => "T",
		"ACG" => "T",
		"TCT" => "S",
		"TCC" => "S",
		"TCA" => "S",
		"TCG" => "S",
		"AGT" => "S",
		"AGC" => "S",
		"TAT" => "Y",
		"TAC" => "Y",
		"TGG" => "W",
		"CAA" => "Q",
		"CAG" => "Q",
		"AAT" => "N",
		"AAC" => "N",
		"CAT" => "H",
		"CAC" => "H",
		"GAA" => "E",
		"GAG" => "E",
		"GAT" => "D",
		"GAC" => "D",
		"AAA" => "K",
		"AAG" => "K",
		"CGT" => "R",
		"CGC" => "R",
		"CGA" => "R",
		"CGG" => "R",
		"AGA" => "R",
		"AGG" => "R",
		"TAA" => "*",
		"TAG" => "*",
		"TGA" => "*",
	);
	return $codon2aa{$codon};
}
