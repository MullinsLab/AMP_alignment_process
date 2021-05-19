#!/usr/bin/perl -w

################################################################################
# Program: collapse_seqs_by_groups.pl
# Purpose: collapse group of sequences to unique sequences, output unique sequence files 
# and corresponding naming file
# Author: Wenjie Deng
# Date: 2021-02-08
################################################################################

use strict;
use warnings;
use v5.10;
use Getopt::Long;
use File::Basename;

my $usage = "\nusage: collapse_seqs.pl infastafile\n";
my $infile = shift or die $usage;
my $count = 0;
my $seqName = my $outfile = my $namefile = '';
my (%nameSeq, %uniqDup, %seqCount, %idseqnames, %idseqcount);
if ($infile =~ /(.*?)_uncollapse/) {
	$outfile = $1."_collapsed_by_timepoint.fasta";
	$namefile = $1."_name_by_timepoint.txt";
}else {
	die "file name is not formatted: $infile\n";
}
open INFASTA, $infile or die "couldn't open $infile: $!\n";
while (my $line = <INFASTA>) {
	$line =~ s/\R$//;
	next if ($line =~ /^\s*$/);
	if ($line =~ /^>(\S+)/) {
		$seqName = $1;
		if ($seqName =~ /\|/) {
			my @fields = split /\|/, $seqName; 
			$seqName = $fields[1];
		}
		if ($seqName =~ /^(RV\d+)_(\d+)_(.*?)_(GP|REN)_/ or $seqName =~ /^(V\d+)_(\d+)_(.*?)_(GP|REN)_/) { # allow multiple time points file: (.*?) could be \d+ or \d+-\d+ or \d+-\d+-\d+ and so on
			my $id = $1."_".$2."_".$3."_".$4;
            push @{$idseqnames{$id}}, $seqName;
			++$count;
        }else {
        	die "sequence name is not formatted: $seqName\n";
        }
	}else {
		$nameSeq{$seqName} .= uc $line;
	}  	
}
close INFASTA;

foreach my $id (keys %idseqnames) {
	foreach my $name (@{$idseqnames{$id}}) {
		my $seq = $nameSeq{$name};
		if (!$seqCount{$id}{$seq}) {
			$seqCount{$id}{$seq} = 0;
		}
		push @{$uniqDup{$id}{$seq}}, $name;
		++$seqCount{$id}{$seq};
		++$idseqcount{$id};
	}
	
}

open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
open NAME, ">", $namefile or die "couldn't open $namefile: $!\n";
print "Total $count sequences in $infile.\n";
foreach my $id (sort {$a cmp $b} keys %seqCount) {
	my $uniqueCount = 0;
	foreach my $seq (sort {$seqCount{$id}{$b} <=> $seqCount{$id}{$a}} keys %{$seqCount{$id}}) {
		$uniqueCount++;
		my $name = $id."_NT_".$uniqueCount."_".$seqCount{$id}{$seq};
		print OUT ">$name\n$seq\n";
		print NAME $name, "\t", join(',', @{$uniqDup{$id}{$seq}}), "\n";
	}
	print "$idseqcount{$id} sequences, $uniqueCount unique sequences in $id.\n";
}
print "\n";
close OUT;
close NAME;

