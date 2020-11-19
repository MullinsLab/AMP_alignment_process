#!/usr/bin/perl


use strict;
use warnings;
use v5.10;

my $usage = "perl separate_timepoint_sequences.pl inSeqFastaFile\n";
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
		foreach my $tp (@tpss) {
			my $outfile = "V".$vid."_".$sid."_".$tp."_".$rest;
			my $tpreg = "_".$tp."_";
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
			print "total $count sequences, $tpcount with time point of $tp.\n";
		}
	}	
}else {
	die "file name not formatted: $infile\n";
}

