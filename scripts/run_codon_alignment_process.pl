#!/usr/bin/perl -w

##########################################################################################
# Program: run_codon_alignment_process.pl
# Purpose: collapse sequences into unique sequences, add HXB2 sequence, muscle align 
# collapsed sequences and HXB2, codon align the collapsed alignment
# Author: Wenjie Deng
# Date: 2021-04-06
##########################################################################################

use strict;
use warnings;
use v5.10;
use Getopt::Long;
use File::Basename;

my %option = (
	'i' => '',
	'r' => '',
);

my $usage = "\nusage: run_codon_alignment_process.pl [-option value]

options:  
-i		input nucleotide sequence fasta file
-r		reference sequence fasta file
";

GetOptions (\%option, 'i=s', 'r=s');

my $infile = $option{'i'} or die $usage;
my $reffile = $option{'r'} or die $usage;
my $scriptspath = dirname(__FILE__);
# collapse sequences
print "\n=== Collapse sequences in $infile ===\n";
system("perl $scriptspath/collapse_seqs.pl -if $infile");
my $collapsedfile = my $collapsedWithHXB2File = my $alignedfile = my $HXB2topalignedfile = my $codonalignfile = $infile;
$collapsedfile =~ s/\.fasta/_collapsed.fasta/;
$collapsedWithHXB2File =~ s/\.fasta/_collapsed_withHXB2.fasta/;
$alignedfile =~ s/\.fasta/_collapsed_withHXB2_aligned.fasta/;
$HXB2topalignedfile =~ s/\.fasta/_collapsed_withHXB2_top_aligned.fasta/;
$codonalignfile =~ s/\.fasta/_collapsed_codon_aligned.fasta/;
# add HXB2 sequence
open OUT, ">", $collapsedWithHXB2File or die "couldn't open $collapsedWithHXB2File: $!\n";
open HXB2, "<", $reffile or die "couldn't open $reffile: $!\n";
while (my $line = <HXB2>) {
	$line =~ s/\R//;
	next if $line =~ /^\s*$/;
	print OUT "$line\n";
} 
close HXB2;
open IN, "<", $collapsedfile or die "couldn't open $collapsedfile: $!\n";
while (my $line = <IN>) {
	$line =~ s/\R//;
	next if $line =~ /^\s*$/;
	print OUT "$line\n";
} 
close IN;
close OUT;
# align collapsed sequences and HXB2
print "=== Align collapsed file $collapsedWithHXB2File ===\n";
system("muscle -quiet -in $collapsedWithHXB2File -out $alignedfile");
# guaranty that HXB2 sequence at the top of alignment
my $name = "";
my (@names, %nameseq);
open IN, "<", $alignedfile or die "couldn't open $alignedfile: $!\n";
while (my $line = <IN>) {
	$line =~ s/\R//;
	next if $line =~ /^\s*$/;
	if ($line =~ /^>(\S+)/) {
		$name = $1;
		if ($name =~ /HXB2/) {
			unshift @names, $name;
		}else {
			push @names, $name;
		}
	}else {
		$nameseq{$name} .= $line;
	}
} 
close IN;
open OUT, ">", $HXB2topalignedfile or die "couldn't open $HXB2topalignedfile: $!\n";
foreach my $name (@names) {
	print OUT ">$name\n$nameseq{$name}\n";
}
close OUT;
# run codon alignment via CodonAlign.pl
print "== Codon alignment for $HXB2topalignedfile ==\n";
system("$scriptspath/codonalign_uw.pl $HXB2topalignedfile 0 1 15 .");
