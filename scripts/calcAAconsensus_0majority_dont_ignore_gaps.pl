#!/usr/bin/perl

# calculate consensus sequence from an AA sequence alignment
# majority character (including "-") gets consensus. If two or more AA share the majority, "X" will be consensus.
# If an AA and "-" share the majority, the "AA" will be consensus

use strict;
use warnings;
use v5.10;

my $usage = "perl calcAAconsensus_0majority_dont_ignore_gaps.pl infile outconsensusfile\n";
my $infile = shift or die $usage;
my $outconsfile = my $outaafile = ""; 
my $subject = "";
if ($infile =~ /^(V70\d)_\d\d\d\d_(.*?)_(.*?)_AA_/) {
	my @fields = split /\_/, $infile;
	$outconsfile = $fields[0]."_".$fields[3]."_AA_functional_consensus.fasta";
	$outaafile = $fields[0]."_".$fields[3]."_AA_functional.fasta";
	$subject = $fields[0]."_".$fields[1]."_".$fields[2]."_".$fields[3];
}else {
	die "file name is not formatted: $infile\n";
}
my $name = '';
my $count = my $processed = 0;
my (@names, %nameSeq, %nameAAs, @consAAs);
open IN, "<", $infile or die "couldn't open $infile: $!\n";
while (my $line = <IN>) {
	$line =~ s/\R$//;
	next if ($line =~ /^\s*$/);
	if ($line =~ /^>(\S+)/) {
		$name = $1;
		unless ($name =~ /HXB2_/) {
			++$count;
			push @names, $name;
		}		
	}else {
		unless ($name =~ /HXB2_/) {
			$nameSeq{$name} .= $line;
		}		
	}
}
close IN;

my $alignlen = 0;
foreach my $name (@names) {
	++$processed;
	my $seq = $nameSeq{$name};
	my $len = length $seq;
	if ($processed == 1) {
		$alignlen = $len;
	}
	if ($len != $alignlen) {
		die "sequence not aligned: $name\n";
	}
	@{$nameAAs{$name}} = split //, $seq;
}

for (my $i = 0; $i < $alignlen; $i++) {
	my %posAAcount = ();
	foreach my $name (@names) {
		if ($name =~ /_(\d+)$/) {
			$posAAcount{$nameAAs{$name}[$i]} += $1;
		}else {
			die "name not formatted: $name\n";
		}
	}
	
	my @sortedaas = sort {$posAAcount{$b} <=> $posAAcount{$a}} keys %posAAcount;
	if (scalar @sortedaas == 1) {
		push @consAAs, $sortedaas[0];
	}elsif ($posAAcount{$sortedaas[0]} > $posAAcount{$sortedaas[1]}) {
		push @consAAs, $sortedaas[0];
	}elsif ($posAAcount{$sortedaas[0]} == $posAAcount{$sortedaas[1]}) {
		if ($posAAcount{$sortedaas[0]} ne "-" and $posAAcount{$sortedaas[1]} ne "-") {
			push @consAAs, "X";
		}else {
			if (scalar @sortedaas == 2 or $posAAcount{$sortedaas[1]} > $posAAcount{$sortedaas[2]}) {
				if ($posAAcount{$sortedaas[0]} ne "-") {
					push @consAAs, $sortedaas[0];
				}else {
					push @consAAs, $sortedaas[1];
				}
			}elsif ($posAAcount{$sortedaas[1]} == $posAAcount{$sortedaas[2]}) {
				push @consAAs, "X";
			}else {
				die "impossible! posAAcount{$sortedaas[1]}: $posAAcount{$sortedaas[1]} < posAAcount{$sortedaas[2]}: $posAAcount{$sortedaas[2]}\n";
			}			
		}
	}else {
		die "impossible! posAAcount{$sortedaas[0]}: $posAAcount{$sortedaas[0]} < posAAcount{$sortedaas[1]}: $posAAcount{$sortedaas[1]}\n";
	}
}

my $outdir = "calcAAconsensus_0majority_dont_ignore_gaps_outputs";
$outconsfile = "$outdir/$outconsfile";
$outaafile = "$outdir/$outaafile";
unless (-e $outdir) {
	mkdir $outdir;
}
open CONS, ">>", $outconsfile or die "couldn't open $outconsfile: $!\n";
print CONS ">$subject","_AA_functional_consensus\n";
print CONS join('', @consAAs), "\n";
close CONS;
my $aaalignfile = $infile;
$aaalignfile =~ s/\.fasta/_withCons.fasta/;
$aaalignfile = "$outdir/$aaalignfile";
open ALIGN, ">", $aaalignfile or die "couldn't open $aaalignfile: $!\n";
print ALIGN ">$subject","_AA_functional_consensus\n";
print ALIGN join('', @consAAs), "\n";
open AA, ">>", $outaafile or die "couldn't open $outaafile: $!\n";
foreach my $name (@names) {
	print ALIGN ">$name\n";
	print ALIGN join('', @{$nameAAs{$name}}), "\n";
	print AA ">$name\n";
	print AA join('', @{$nameAAs{$name}}), "\n";
}
close AA;
close ALIGN;

print "In $subject: total $count AA sequences, processed $processed AA sequences\n";