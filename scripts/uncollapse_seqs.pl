#!/usr/bin/perl -w

# un-collapse unique sequences into original sequences
# input: unique sequence fasta file, naming file of the relationship between unique name and it original names
# output: un-collapse sequence fasta file with original sequence names

use strict;

my $usage = "perl uncollapse.pl inCollapsedSeqFastaFile inNameFile inNameFile1 (inNameFile2 inNameFile3 ...)\n";
my $inFasta = shift @ARGV or die $usage;
if (scalar @ARGV == 0) {
	die $usage;
}
my $outfile = $inFasta;
$outfile =~ s/\.fasta/_uncollapse.fasta/;
my $name = "";
my $count = my $uniqcount = 0;
my (%nameNames, %nameSeq, @names);
foreach my $inName (@ARGV) {
	open NAME, $inName or die "couldn't open $inName: $!\n";
	while (my $line = <NAME>) {
		chomp $line;
		next if $line =~ /^\s*$/;
		my ($uniqname, $originalnames) = split /\t/, $line;
		@{$nameNames{$uniqname}} = split /,/, $originalnames;
	}
	close NAME;
}

open IN, $inFasta or die "couldn't open $inFasta: $!\n";
while (my $line = <IN>) {
	chomp $line;
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