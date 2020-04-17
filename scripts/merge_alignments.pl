#!/usr/bin/perl -w

##########################################################################################
# Program: merge_alignments.pl
# Purpose: In a directory with alignment fasta files, merges alignments based on file name
# Author: Wenjie Deng
# Date: 2020-04-07
##########################################################################################

use strict;
use Getopt::Long;

my %option = (
	'id' => '.',
);

my $usage = "\nusage: unique_reads.pl [-option value]

options:  
-id		input directory with fasta files (default: . )

";

GetOptions (\%option, 'id=s');

my $indir = $option{'id'} or die $usage;
my $name = "";
my (%sampleRegionNameSeq, %sampleRegionNameFlag, %sampleRegionTP, %sampleRegionTpStatus, %regionStatus, @regions);
opendir DIR, $indir or die "couldn't open $indir: $!\n";
while (my $file = readdir DIR) {
	if ($file =~ /(.*)\.fasta$/) {
		$file = $indir."/".$file;
		if ($file =~ /V(\d+)_(\d+)_(\d+)_([A-Z]+)/) {
			my $id = "V".$1."_".$2;
			my $tp = $3;
			my $region = $4;
			if (!$regionStatus{$region}) {
				$regionStatus{$region} = 1;
				push @regions, $region;
			}
			if (!$sampleRegionTpStatus{$id}{$region}{$tp}) {
				$sampleRegionTpStatus{$id}{$region}{$tp} = 1;
				push @{$sampleRegionTP{$id}{$region}}, $tp;
			}			
			open IN, $file or die "couldn't open $file: $!\n";
			while (my $line = <IN>) {
				chomp $line;
				next if ($line =~ /^\s*$/);
				if ($line =~ /^>(\S+)/) {
					$name = $1;
					if ($name =~ /_rpt_/) {
						$name =~ s/_rpt_/_/;
					}
					if (!$sampleRegionNameFlag{$id}{$region}{$name}) {
						$sampleRegionNameFlag{$id}{$region}{$name} = 1;
					}else {
						++$sampleRegionNameFlag{$id}{$region}{$name};
						$name = $name."_".$sampleRegionNameFlag{$id}{$region}{$name};
					}
				}else {
					$line =~ s/\-//g;
					$sampleRegionNameSeq{$id}{$region}{$name} .= $line;
				}
			}
			close IN;
		}
	}
}
closedir DIR;

foreach my $rg (@regions) {
	my $outdir = $indir."/".$rg;
	unless (-e $outdir) {
		mkdir $outdir;
	}
}

foreach my $id (sort {$a cmp $b} keys %sampleRegionTP) {
	foreach my $rg (sort {$a cmp $b} keys %{$sampleRegionTP{$id}}) {
		my @tps = sort @{$sampleRegionTP{$id}{$rg}};
		my $tp = join("-", @tps);
		my $outfile = $indir."/".$rg."/".$id."_".$tp."_".$rg."_NT.fasta";
		my $count = 0;
		open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
		foreach my $name (sort {$a cmp $b} keys %{$sampleRegionNameSeq{$id}{$rg}}) {
			print OUT ">$name\n$sampleRegionNameSeq{$id}{$rg}{$name}\n";
			++$count;
		}
		close OUT;
		print "Write file $outfile, total $count sequences\n";
	}
}



