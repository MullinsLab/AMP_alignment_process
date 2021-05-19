#!/usr/bin/perl -w

################################################################################
# Program: collapse_seqs.pl
# Purpose: collapse sequences to unique sequences, output unique sequence files 
# and corresponding naming file
# Author: Wenjie Deng
# Date: 2020-01-28
# Modified: 2020-04-02
################################################################################

use strict;
use warnings;
use v5.10;
use Getopt::Long;
use File::Basename;

my %option = (
	'if' => '',
);

my $usage = "\nusage: collapse_seqs.pl [-option value]

options:  
-if		input sequence fasta file

";

GetOptions (\%option, 'if=s');

my $infile = $option{'if'} or die $usage;
my $nametag = $infile;
$nametag = basename($nametag);
if ($nametag =~ /(.*)\.fasta$/) {
	$nametag = $1;
	$nametag =~ s/\s+/_/g;
}else {
	die "Not correct fasta file extension, must be '.fasta'\n";
}

my $count = 0;
my $seq = my $seqName = my $hxb2name = my $hxb2seq = '';
my (%seqCount, @names, %nameSeq);
my %seqLen;
my %nameTag;
my (@uniqSeqs, $uniqDup);
open INFASTA, $infile or die "couldn't open $infile: $!\n";
while (my $line = <INFASTA>) {
	$line =~ s/\R$//;
	next if ($line =~ /^\s*$/);
	if ($line =~ /^>(.*)/) {
		$seqName = $1;
		if ($seqName =~ /HXB2/) {
			$hxb2name = $seqName;
		}else {
			push @names, $seqName;
		}	
	}else {
		$nameSeq{$seqName} .= uc $line;
	}		
}
close INFASTA;

foreach my $name (@names) {
	$seq = $nameSeq{$name};
	if (!$seqCount{$seq}) {
		$seqCount{$seq} = 0;
	}
	push @{$uniqDup->{$seq}}, $name;
	++$seqCount{$seq};			
	++$count;	
}

my $uniqueCount = 0;
my $outfile = $infile;
$outfile =~ s/\.fasta$/_collapsed.fasta/;
my $namefile = $infile.".name";
open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
open NAME, ">", $namefile or die "couldn't open $namefile: $!\n";
#print OUT ">$hxb2name\n$nameSeq{$hxb2name}\n";
foreach my $seq (sort {$seqCount{$b} <=> $seqCount{$a}} keys %seqCount) {
	$uniqueCount++;
	my $name = $nametag."_".$uniqueCount."_".$seqCount{$seq};
	print OUT ">",$name,"\n",$seq,"\n";
	print NAME $name, "\t", join(',', @{$uniqDup->{$seq}}), "\n";
}
close OUT;
close NAME;
print "Total $count sequences in $infile. $uniqueCount unique sequences.\n";
