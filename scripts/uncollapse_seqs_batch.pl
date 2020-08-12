#!/usr/bin/perl -w

# un-collapse unique sequences into original sequences
# input: unique sequence fasta file, naming file of the relationship between unique name and it original names
# output: un-collapse sequence fasta file with original sequence names

use strict;
use warnings;
use v5.10;

my $usage = "perl uncollapse.pl inCollapsedSeqFastaFile\n";
my $inFasta = shift or die $usage;
my $namefile = "";
my $outfile = $inFasta;
$outfile =~ s/\.fasta|\.fas/_uncollapsed.fasta/;
if ($inFasta =~ /^V(\d+)_(\d+)_(.*?)_(.*?)_NT_/) {
	$namefile = "../../NameNoriginalFiles/V".$1."_".$2."_".$3."_".$4."_NT.fasta.name";
}else {
	die "file name not formatted: $inFasta\n";
}
print "\n*** uncollapsing $inFasta ***\n";
my $name = "";
my $count = my $uniqcount = 0;
my (%nameNames, %nameSeq, @names);
open NAME, $namefile or die "couldn't open $namefile: $!\n";
while (my $line = <NAME>) {
	$line =~ s/\R$//;
	next if $line =~ /^\s*$/;
	my ($uniqname, $originalnames) = split /\t/, $line;
	@{$nameNames{$uniqname}} = split /,/, $originalnames;
}
close NAME;

open IN, $inFasta or die "couldn't open $inFasta: $!\n";
while (my $line = <IN>) {
	$line =~ s/\R$//;
	next if $line =~ /^\s*$/;	
	if ($line =~ />(\S+)/) {
		$name = $1;
		push @names, $name;
		++$uniqcount;
	}else {
		$nameSeq{$name} .= uc $line;
	}
}
close IN;
open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
foreach my $name (@names) {
	if (exists $nameNames{$name}[0]) {
		if ($name =~ /(.*)_(\d+)$/) {
			my $duplicates = $2;
			my $members = scalar @{$nameNames{$name}};
			if ($duplicates == $members) {
				foreach my $originalname (@{$nameNames{$name}}) {
					++$count;
					print OUT ">$originalname\n$nameSeq{$name}\n";
				}
			}else {
				die "$name, nameNames{$name}: @{$nameNames{$name}}, $members, $duplicates, duplicates not same\n";
			}
		}		
	}else {
		++$count;
		print OUT ">$name\n$nameSeq{$name}\n";
	}
}
close OUT;

print "total $uniqcount unique sequences, expanded to $count sequences\n";