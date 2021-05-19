#!/usr/bin/perl -w

##########################################################################################
# Program: left_aling_muscle.pl
# Purpose: make a left aligned sequence alignment via muscle by reversing sequences, 
# running muscle, reversing back sequences
# Author: Wenjie Deng
# Date: 2020-09-30
##########################################################################################

use strict;
use warnings;
use v5.10;
use Getopt::Long;
use File::Basename;

my $usage = "\nusage: left_align_muscle.pl inputfastafile\n";

my $infile = shift or die $usage;
my $outfile = my $reversefile = my $revesealignfile = $infile;
$outfile =~ s/\.fasta/_leftAligned.fasta/;
$reversefile =~ s/\.fasta/_reversed.fasta/;
$revesealignfile =~ s/\.fasta/_reverseAligned.fasta/;
my @names = my %nameSeq = ();
my $name = "";
my $count = 0;
print "\n=== implementing $infile ===\n";
open IN, $infile or die "couldn't open $infile: $!\n";
while (my $line = <IN>) {
	$line =~ s/\R$//;
	next if $line =~ /^\s*$/;
	if ($line =~ /^>/) {
		$name = $line;
		push @names, $name;
		++$count;		
	}else {
		$nameSeq{$name} .= $line;
	}
}
close IN;

open REV, ">", $reversefile or die "couldn't open $reversefile: $!\n";
foreach my $name (@names) {
	my $seq = reverse $nameSeq{$name};
	print REV "$name\n$seq\n";
}
close REV;

system("muscle -quiet -in $reversefile -out $revesealignfile");

$name = "";
%nameSeq = ();
open IN, $revesealignfile or die "couldn't open $revesealignfile: $!\n";
open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
while (my $line = <IN>) {
	$line =~ s/\R$//;
	next if $line =~ /^\s*$/;
	if ($line =~ /^>/) {
		$name = $line;
	}else {
		$nameSeq{$name} .= $line;
	}
}
close IN;

open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
foreach my $name (@names) {
	my $seq = reverse $nameSeq{$name};
	print OUT "$name\n$seq\n";
}
close OUT;

print "total $count sequences\n";

