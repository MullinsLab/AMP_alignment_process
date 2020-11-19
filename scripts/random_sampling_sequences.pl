#!/usr/bin/perl

##########################################################################################
# Program: random_sampling_sequences.pl
# Purpose: In a sequence fasta file, randomly sampling certain number of sequences 
# Author: Wenjie Deng
# Date: 2020-11-17
##########################################################################################

use strict;
use warnings;
use v5.10;

my $usage = "Usage: random_sampling_sequences.pl infile numberOfSampling numberOfSequences";

my $infile = shift or die $usage;
my $sampling = shift || 10;
my $seqs = shift || 20;
my $name = "";
my $count = 0;
my (@names, %nameSeq);
open IN, $infile or die "couldn't open $infile: $!\n";
while (my $line = <IN>) {
	$line =~ s/\R$//;
	next if ($line =~ /^\s*$/);
	if ($line =~ /^>(.*)/) {
		++$count;
		$name = $1;
		push @names, $name;
	}else {
		$nameSeq{$name} .= $line;
	}
}
close IN;

for (my $i = 1; $i <= $sampling; $i++) {
	my $outfile = $infile;
	$outfile =~ s/\.fasta/_sampling_$i.fasta/;
	my (%seqStatus, @seqnames);
	my $seqcount = 0;
	while (1) {
        if ($seqcount == $seqs) {
			open(OUT, ">", $outfile) or die "couldn't Open $outfile: $!\n";
			foreach my $name (@seqnames) {
				print OUT ">$name\n$nameSeq{$name}\n";
			}
			close OUT;
            last;
        }else {
			my $idx = int(rand($count));
			if (!$seqStatus{$idx}) {
                $seqStatus{$idx} = 1;
				++$seqcount;
				push @seqnames, $names[$idx];
            }            
		}        
    }   
}

say "total $count sequences, generated $sampling samplings with $seqs random sequences\n";