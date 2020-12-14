#!/usr/bin/perl -w


use strict;
use warnings;
use v5.10;

my $usage = "perl verify_seq_origin.pl inSeqFastaFile directoryForOriginalFiles originalFileNameSuffix(including file extension)\n";
my $inFasta = shift or die $usage;
my $dir = shift or die $usage;
my $suffix = shift or die $usage;
my $originfile = "";

if ($inFasta =~ /^V(\d+)_(\d+)_(.*?)_/) {
    $originfile = $dir."/V".$1."_".$2."_".$3."_".$suffix;
}elsif ($inFasta =~ /^RV(\d+)_(\d+)_(.*?)_/) {
    $originfile = $dir."/RV".$1."_".$2."_".$3."_".$suffix;
}else {
	die "file name not formatted: $inFasta\n";
}
print "\n*** $inFasta vs $originfile ***\n";
my $name = my $originname = "";
my $count = my $origincount = 0;
my (%nameSeq, %originnameSeq, @names);
open ORIGIN, $originfile or die "couldn't open $originfile: $!\n";
while (my $line = <ORIGIN>) {
	$line =~ s/\R$//;
	next if $line =~ /^\s*$/;
	if ($line =~ /^>(.*)/) {
		++$origincount;
		$originname = $1;
	}else {
		$line = uc $line;
		$line =~ s/\-//g;
		$originnameSeq{$originname} .= $line;
	}
}
close ORIGIN;

open IN, $inFasta or die "couldn't open $inFasta: $!\n";
while (my $line = <IN>) {
	$line =~ s/\R$//;
	next if $line =~ /^\s*$/;	
	if ($line =~ />(.*)/) {
		++$count;
		$name = $1;
	}else {
		$line = uc $line;
		$line =~ s/\-//g;
		$nameSeq{$name} .= $line;
	}
}
close IN;

if ($count != $origincount) {
	print "number of sequences are not same. count $count, original count $origincount\n";
}

my $exactcount = my $containcount = 0;
foreach my $originalname (keys %originnameSeq) {
	if ($nameSeq{$originalname}) {
		my $seq = $nameSeq{$originalname};
		my $originalseq = $originnameSeq{$originalname};
		$originalseq =~ s/\*/Z/;
		if ($seq eq $originalseq) {
			++$exactcount;
		}else {
			my $idx = index ($originalseq, $seq);
			if ($idx >= 0) {
				++$containcount;
			}else {
				die "sequence not in original sequence! nameSeq{$originalname}: $nameSeq{$originalname}\noriginnameSeq{$originalname}: $originnameSeq{$originalname}\n";
			}
		}
	}else {
		print "$originalname missing in $inFasta file\n";
	}
}

print "total $count sequences, $origincount original sequences, $exactcount exact match, $containcount contained in original sequences, all verified.\n";
