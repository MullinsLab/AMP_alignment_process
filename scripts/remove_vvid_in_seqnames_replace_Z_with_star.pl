#!/usr/bin/perl

# remove Viroverse ID in sequence names 
# original name: \d+\.\d|(RV|V)\d+_
# new name: (RV|V)\d+_

use strict;
use warnings;
use v5.10;

my $usage = "perl remove_vvid_in_seqnames.pl inFastaFile\n";
my $infasta = shift or die $usage;
my $outdir = "outputs";
my $outfile = "$outdir/$infasta";
unless (-e $outdir) {
	mkdir $outdir;
}
my $name = "";
my $count = my $vvflag = my $dflag = my $vvidcount = 0;
my (@names, %nameseq);
open IN, $infasta or die "couldn't open $infasta: $!\n";
open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
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
		print OUT ">$name\n";
		++$count;
	}else {
		if ($name =~ /HXB2/) {
			if ($line =~ /\-+$/) {
				$line =~ s/(\-+)$/\*$1/;
				$line =~ s/\-$//;
			}else {
				die "HXB2 $infasta\n";
			}
		}else {
			$line =~ s/Z/\*/g;
		}
		if ($line =~ /[A-Z]$/) {
			print "file: $infasta\n";
		}
		print OUT "$line\n";
	}
}
close IN;
close OUT;

if ($vvflag or $dflag) {
	print "\nOrigianl file: $infasta, renamed file: $outfile, total $count sequences, $vvidcount with vvid\n";
}


