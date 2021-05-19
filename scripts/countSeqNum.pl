#!/usr/bin/perl -w

use strict;

my $usage = "usage: perl countSeqNum.pl inputFastaFile outfile\n";

my $infile = shift or die $usage;
my $outfile = shift or die $usage;
my $count = 0;
open IN, $infile or die "couldn't open $infile: $!\n";
while (my $line = <IN>) {
	chomp $line;
	if ($line =~ /^>(\S+)/) {
		my $name = $1;
		unless ($name =~ /HXB2/) {
			++$count;
		}
	}
}
close IN;
open OUT, ">>", $outfile or die "could't open $outfile: $!\n";
print OUT "$infile,$count\n";

print "In $infile: total $count sequences.\n";

