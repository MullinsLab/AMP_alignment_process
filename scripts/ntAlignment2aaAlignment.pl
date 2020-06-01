#!/usr/bin/perl

# translates any nucleotide alignment into amino acid alignment
# xxxx------xxxxx --> X--XX
# xxxxx------xxxx --> XX--X

use strict;
use warnings;
use v5.10;

my $usage = "perl ntAlignment2aaAlignment.pl inNtAlignment\n";
my $infile = shift or die $usage;
my $outfile = $infile;
$outfile =~ s/_NT/_AA/;
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
	my (@codonnts, @codongaps);
	for (my $i = 0; $i < $len; $i += 3) {
		my $ntcount = scalar @codonnts;
		if ($nts[$i] =~ /[ACTG]/) {
			push @codonnts, $nts[$i];
		}elsif ($nts[$i] eq "-") {
			push @codongaps, $nts[$i];
		}
		if ($nts[$i+1] =~ /[ACTG]/) {
			push @codonnts, $nts[$i+1];
		}elsif ($nts[$i+1] eq "-") {
			push @codongaps, $nts[$i+1];
		}
		if ($nts[$i+2] =~ /[ACTG]/) {
			push @codonnts, $nts[$i+2];
		}elsif ($nts[$i+2] eq "-") {
			push @codongaps, $nts[$i+2];
		}
		
		if (scalar @codongaps >= 3) {
			shift @codongaps;
			shift @codongaps;
			shift @codongaps;
			$partialaaseq .= "-";
			if ($ntcount == 0) {
				$aaseq .= "-";
				$partialaaseq = "";
			}
		}
		
		if (scalar @codonnts >= 3) {
			my $codon = shift @codonnts;
			$codon .= shift @codonnts;
			$codon .= shift @codonnts;
			my $aa = translation($codon);
			if ($ntcount == 0) {
				$aaseq .= $aa;
				if ($aa eq "*") {
					last;
				} 
			}elsif ($ntcount == 1) {
				$partialaaseq .= $aa;
				if ($aa eq "*") {
					$aaseq .= $aa;
					last;
				}else {
					$aaseq .= $partialaaseq;
				}
			}elsif ($ntcount == 2) {
				$partialaaseq = $aa.$partialaaseq;
				if ($aa eq "*") {
					$aaseq .= $aa;
					last;
				}else {
					$aaseq .= $partialaaseq;
				}
			}else {
				die "impossible: ntcount = $ntcount\n";
			}
			$partialaaseq = "";
		}
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
