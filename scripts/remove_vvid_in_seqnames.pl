#!/usr/bin/perl

# remove Viroverse ID in sequence names 
# original name: \d+\.\d|(RV|V)\d+_
# new name: (RV|V)\d+_

use strict;
use warnings;
use v5.10;

my $usage = "perl remove_vvid_in_seqnames.pl inFastaFile\n";
my $infasta = shift or die $usage;
my $outdir = "remove_vvid_outputs";
my $outfile = "$outdir/$infasta";
unless (-e $outdir) {
	mkdir $outdir;
}
my $name = "";
my $count = my $vvflag = my $dflag = my $vvidcount = 0;
my (@names, %nameseq);
open IN, $infasta or die "couldn't open $infasta: $!\n";
while (my $line = <IN>) {
	$line =~ s/\R$//;
	next if $line =~ /^\s*$/;	
	if ($line =~ />(\S+)/) {
		$name = $1;
		if ($name =~ /\|/) {
			my @fields = split /\|/, $name;
			$name = $fields[1];
			$vvflag = 1;
			++$vvidcount;
		}
		if ($name =~ /_\d$/) {
			print "$name -> ";
			$name =~ s/_(\d)$/\-$1/;
			$dflag = 1;
			print "$name\n";
		}
		push @names, $name;
		++$count;
	}else {
		$nameseq{$name} .= $line;
	}
}
close IN;

open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
foreach my $name (@names) {
	print OUT ">$name\n$nameseq{$name}\n";
}
close OUT;
print "\nOrigianl file: $infasta, renamed file: $outfile, total $count sequences, $vvidcount with vvid\n";



