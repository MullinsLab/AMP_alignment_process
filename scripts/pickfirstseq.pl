#!/usr/bin/perl -w


use strict;
use warnings;
use v5.10;

my $usage = "perl pickfirstseq.pl infasta outfasta\n";
my $inFasta = shift or die $usage;
my $outfasta = shift or die $usage;
my $count = 0;
open IN, $inFasta or die "couldn't open $inFasta: $!\n";
open OUT, ">>", $outfasta or die "couldn't open $outfasta: $!\n";
while (my $line = <IN>) {
	$line =~ s/\R$//;
	next if $line =~ /^\s*$/;	
	if ($line =~ />(.*)/) {
		++$count;
		last if $count > 1;	
	}
	print OUT "$line\n";
}
close IN;
close OUT;
