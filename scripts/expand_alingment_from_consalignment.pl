#!/usr/bin/perl -w

#########################################################################################################
# Program: expand_alignment_from_consalignment.pl
# Purpose: from an alignment of consensus sequences, expand the alignment to sequences of those consensus
# Author: Wenjie Deng
# Date: 2020-10-12
###########################################################################################################

use strict;
use v5.10;
use File::Copy;

my $usage = "perl expand_alignment_from_consalignment.pl consensus_alignment expand_name all_seq outfile\n";
my $consalignfile = shift or die $usage;
my $expandname = shift or die $usage;
my $seqfile = shift or die $usage;
my $outfile = shift or die $usage;
my $subjectcount = my $flag = 0;
my $name = my $AAalignseq = "";
my %CONnameAAs = ();
open NT, "<", $seqfile or die "couldn't open $seqfile: $!\n";
while (my $line = <NT>) {
	$line =~ s/\R$//;
	next if $line =~ /^\s*$/;
	if ($line =~ /^>(\S+)/) {
		$name = $1;
	}else {
		my @nts = split //, $line;
		my $len = length $line;
		for (my $i = 0; $i < $len; $i++) {
			$CONnameAAs{$name}{$i} = $nts[$i];
		}
	}
}
close NT;

print "\n*** expanding alignment based on consensus alignment ***\n";
my $ntalignfile = $seqfile;
$ntalignfile =~ s/\.fasta/_aln.fasta/;
open IN, "<", $consalignfile or die "couldn't open $consalignfile: $!\n";
while (my $line = <IN>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	if ($line =~ /^>/) {
		if ($line =~ /^>$expandname/) {
			$flag = 1;
		}else {
			$flag = 0;
		}
	}elsif ($flag) {
		$AAalignseq .= $line;
	}
}
close IN;

copy($consalignfile, $outfile);

my $seqcount = 0;
open OUT, ">>", $outfile or die "couldn't open $outfile: $!\n";
my @aas = split //, $AAalignseq;
my $len = length $AAalignseq;
foreach my $name (sort {$a cmp $b} keys %CONnameAAs) {
	++$seqcount;
	my $idx = -1;
	my @nts = ();
	for (my $i = 0; $i < $len; $i++) {
		if ($aas[$i] eq "-") {
			push @nts, "-";
		}else {
			++$idx;
			push @nts, $CONnameAAs{$name}{$idx};
		}
	}
	print OUT ">$name\n";
	print OUT join('', @nts),"\n";
}	
close OUT;

print "Expanded total $seqcount sequences for $expandname\n";

