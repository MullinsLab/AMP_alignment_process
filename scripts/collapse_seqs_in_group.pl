#!/usr/bin/perl -w

################################################################################
# Program: collapse_seqs_in_group.pl
# Purpose: collapse one group of sequences to unique sequences, output unique sequence files 
# and corresponding naming file
# Author: Wenjie Deng
# Date: 2020-11-04
################################################################################

use strict;
use warnings;
use v5.10;
use Getopt::Long;
use File::Basename;

my %option = (
	'if'    => '',
	'grpid' => '',
);

my $usage = "\nusage: collapse_seqs.pl [-option value]

options:  
-if		input sequence fasta file
-grpid  group id for the group of sequences collapsed (e.g. V704_0123_456_GP)

";

GetOptions (\%option, 'if=s', 'grpid=s');

my $infile = $option{'if'} or die $usage;
my $grpid = $option{'grpid'} or die $usage;
my $nametag = $grpid."_NT";

my $count = my $flag = 0;
my $seq = my $seqName = '';
my (%seqCount);
my %seqLen;
my %nameTag;
my (@seqnames, %nameSeq, @uniqSeqs, $uniqDup);
open INFASTA, $infile or die "couldn't open $infile: $!\n";
while (my $line = <INFASTA>) {
	$line =~ s/\R$//;
	next if ($line =~ /^\s*$/);
	if ($line =~ /^>(\S+)/) {
		$seqName = $1;
		$flag = 0;
		if ($seqName =~ /$grpid/) {
			$flag = 1;
            push @seqnames, $seqName;
			++$count;
        }
	}elsif ($flag) {
		$nameSeq{$seqName} .= uc $line;
	}  	
}
close INFASTA;

foreach my $name (@seqnames) {
	my $seq = $nameSeq{$name};
	if (!$seqCount{$seq}) {
		$seqCount{$seq} = 0;
	}
	push @{$uniqDup->{$seq}}, $seqName;
	++$seqCount{$seq};
}

my $uniqueCount = 0;
my $outfile = $infile;
$outfile =~ s/\.fasta$/_collapsed.fasta/;
my $namefile = $infile.".name";
open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
open NAME, ">", $namefile or die "couldn't open $namefile: $!\n";
foreach my $seq (sort {$seqCount{$b} <=> $seqCount{$a}} keys %seqCount) {
	$uniqueCount++;
	my $name = $nametag."_".$uniqueCount."_".$seqCount{$seq};
	print OUT ">",$name,"\n",$seq,"\n";
	print NAME $name, "\t", join(',', @{$uniqDup->{$seq}}), "\n";
}
close OUT;
close NAME;
print "Total $count sequences in $infile. $uniqueCount unique sequences.\n";
