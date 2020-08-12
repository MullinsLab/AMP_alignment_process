#!/usr/bin/perl -w


use strict;
use warnings;
use v5.10;

my $usage = "perl verify_seq_origin.pl inSeqFastaFile\n";
my $inFasta = shift or die $usage;
my $originfile = "";
if ($inFasta =~ /^V(\d+)_(\d+)_(.*?)_(.*?)_NT_/) {
	$originfile = "../../NameNoriginalFiles/V".$1."_".$2."_".$3."_".$4."_NT.fasta";
#	$originfile = "../../batch2_submission/REN/V".$1."_".$2."_".$3."_".$4."_NT_uncollapsed.fasta";
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
		if ($name =~ /\|/) {
			my @fields = split /\|/, $name;
			$name = $fields[1];
		}
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

foreach my $originalname (keys %originnameSeq) {
	if ($nameSeq{$originalname}) {
		my $seq = $nameSeq{$originalname};
		my $originalseq = $originnameSeq{$originalname};
#		my $idx = index ($originalseq, $seq);
#		if ($idx < 0) {
#			die "sequence not in original sequence! nameSeq{$originalname}: $nameSeq{$originalname}\noriginnameSeq{$originalname}: $originnameSeq{$originalname}\n";
#		}
		unless ($nameSeq{$originalname} eq $originnameSeq{$originalname}) {
			print "sequences not same! nameSeq{$originalname}: $nameSeq{$originalname}\noriginnameSeq{$originalname}: $originnameSeq{$originalname}\n";
		}
	}else {
		print "$originalname missing in $inFasta file\n";
	}
}

print "total $count sequences, all verified.\n";