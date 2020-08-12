#!/usr/bin/perl -w


use strict;
use warnings;
use v5.10;

my $usage = "perl extract_timepoint_sequences.pl inSeqFastaFile\n";
my $infile = shift or die $usage;
print "\n*** Processing $infile ***\n";
if ($infile =~ /^V(\d+)_(\d+)_(.*?)_(.*)/) {
	my $vid = $1;
	my $sid = $2;
	my $tps = $3;
	my $rest = $4;
	my @tpss = split /\-/, $tps;
	if (scalar @tpss == 1) {
		print "Only one time point\n";
		exit;
	}else {
		my $firsttp = $tpss[0];
		my $outfile = "V".$1."_".$2."_".$firsttp."_".$rest;
		my $tpreg = "_".$firsttp."_";
		my $name = "";
		my $count = my $tpcount = my $flag = 0;
		open IN, $infile or die "couldn't open $infile: $!\n";
		open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
		while (my $line = <IN>) {
			$line =~ s/\R$//;
			next if $line =~ /^\s*$/;	
			if ($line =~ />(.*)/) {
				++$count;
				$name = $1;
				if ($name =~ /$tpreg/) {
					++$tpcount;
					$flag = 1;
					print OUT ">$name\n";
				}else {
					$flag = 0;
				}
			}elsif ($flag) {
				print OUT "$line\n";
			}
		}
		close IN;
		close OUT;
		print "total $count sequences, $tpcount with first time point of $firsttp.\n";
	}	
}else {
	die "file name not formatted: $infile\n";
}

