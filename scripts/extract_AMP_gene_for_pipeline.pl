#!/usr/bin/perl -w

######################################################################################
# program: extract_align_region_w_ref_pos.pl
# purpose: extract alignment portion between user defined biggning and ending postion 
# of reference sequence which is at the top of alignment
# Author: Wenjie Deng
# Date: 2010-09-01
# Modified: 2011-04-22
# Changes: added one option to only extract the sequences that cover the portion
# Modified: 2014-06-04
# Changes: save memory usage by implementing reads one by one rather than putting into 
# hash table, it's good for very large dataset
# Modified : 2020-04-02
# Changes: For AMP alignments, add one option of region, change output file name and 
# sequence names according to region 
# Modified (2020-04-13): For AMP alignment pipeline, added option '-oa', removed '-rg'
# Modified (2020-04-23): HXB2 doesn't have to be at the top of alignment
######################################################################################

use strict;
use warnings;
use v5.10;
use Getopt::Long;

my $usage = "\nusage: extract_align_portion_w_ref_pos.pl [-option value]

options:  
-ia		input sequence alignment fasta file
-oa		output extracted sequence alignment fasta file
-rs		start position in reference
-re		ending postion in reference
-cp		cover portion flag (default: false). if set, will output only the sequences that cover the entire portion.
		if not, will output all the sequences in the portion
-rf		flag for removing reference sequence in output file (default: false)
-gf		flag for removing gaps in output file, i.e. output sequences are no longer aligned. (default: false)

";

my $inFile = '';
my $outFile = '';
my $refStart = 0;
my $refEnd = 0;
my $cpFlag = '';
my $rFlag = '';
my $gFlag = '';

GetOptions ('ia=s' => \$inFile, 'oa=s' => \$outFile, 'rs=i' => \$refStart, 're=i' => \$refEnd, 'cp' => \$cpFlag, 'rf' => \$rFlag, 'gf' => \$gFlag);
die $usage unless ($inFile && $outFile && $refStart && $refEnd);

if ($refEnd <= $refStart) {
	die "Ending position must be greater than start position\n";
}

my $refFlag = my $count = my $extractCount = 0;
my $refName = my $refSeq = my $seq = my $seqName = '';

open IN, $inFile or die "couldn't open $inFile: $!\n";
while (my $line = <IN>) {
	$line =~ s/\R$//;
	next if $line =~ /^\s*$/;
	if (!$refFlag and $line =~ /^>(\S*HXB2\S*)/) {	# reference sequences must have HXB2 sequence
		$refName = $1;	
		$refFlag = 1;
	}elsif ($refFlag) {
		if ($line =~ /^>/) {
			$refFlag = 0;
			last;
		}else {
			$refSeq .= $line;
		}		
	}
}
close IN;

if (!$refName or !$refSeq) {
	die "No HXB2 reference sequence\n";
}

my $refseqnogaps = $refSeq;
$refseqnogaps =~ s/\-//g;
my $reflen = length $refseqnogaps;
if ($refEnd > $reflen) {
	$refEnd = $reflen;
}

my @refNas = split //, $refSeq;
my $idx = my $start = my $end = 0;
for (my $i = 0; $i < @refNas; $i++) {
	if ($refNas[$i] =~ /[A-Za-z]/) {
		$idx++;
		if ($idx == $refStart) {
			$start = $i;
		}elsif ($idx == $refEnd) {
			$end = $i;
			last;
		}
	}
}

open OUT, ">", $outFile or die "couldn't open $outFile: $!\n";
open IN, $inFile or die "couldn't open $inFile: $!\n";
while (my $line = <IN>) {
	$line =~ s/\R$//;
	next if $line =~ /^\s*$/;
	if ($line =~ />(\S+)/) {
		if ($count >= 1) {	# extract region
			my $naStart = my $naEnd = 0;
			if ($cpFlag) {					
				my @nas = split //, $seq;
				my $len = scalar @nas;
				for (my $i = 0; $i < $len; $i++) {
					if ($nas[$i] =~ /[A-Za-z]/) {
						$naStart = $i;
						last;
					}
				}
				for (my $i = $len - 1; $i >= 0; $i--) {
					if ($nas[$i] =~ /[A-Za-z]/) {
						$naEnd = $i;
						last;
					}
				}						
			}
			my $portion = substr($seq, $start, $end-$start+1);
			if ($gFlag) {
				$portion =~ s/\-//g;
			}
			if ($cpFlag) {
				if ($naStart <= $start && $naEnd >= $end) {
					print OUT ">$seqName\n$portion\n";
					++$extractCount;
				}
			}else {
				print OUT ">$seqName\n$portion\n";
				++$extractCount;
			}	
		}
		$seqName = $1;
		$seq = '';
		++$count;			
	}else {
		$seq .= $line;
	}
}
close IN;

# last sequence
my $naStart = my $naEnd = 0;
if ($cpFlag) {					
	my @nas = split //, $seq;
	my $len = scalar @nas;
	for (my $i = 0; $i < $len; $i++) {
		if ($nas[$i] =~ /[A-Za-z]/) {
			$naStart = $i;
			last;
		}
	}
	for (my $i = $len - 1; $i >= 0; $i--) {
		if ($nas[$i] =~ /[A-Za-z]/) {
			$naEnd = $i;
			last;
		}
	}						
}
my $portion = substr($seq, $start, $end-$start+1);
if ($gFlag) {
	$portion =~ s/\-//g;
}
if ($cpFlag) {
	if ($naStart <= $start && $naEnd >= $end) {
		print OUT ">$seqName\n$portion\n";
		++$extractCount;
	}
}else {
	print OUT ">$seqName\n$portion\n";
	++$extractCount;
}	
close OUT;

print "total $count sequences, extracted $extractCount sequences in defined region.\n";