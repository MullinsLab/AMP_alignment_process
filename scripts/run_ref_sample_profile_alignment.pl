#!/usr/bin/perl -w

################################################################################
# Program: run_ref_sample_profile_alignment.pl
# Purpose: In a directory with reference sequence alignment and sample alignment fasta 
# files, do profile alignment between reference and sample alignment
# Author: Wenjie Deng
# Date: 2020-04-10
################################################################################

use strict;
use Getopt::Long;
use File::Basename;

my %option = (
	'id' => '.',
	'ref' => '',
);

my $usage = "\nusage: run_ref_sample_profile_alignment.pl [-option value]

options:  
-id     input directory with fasta files (default: . )
-ref    name of reference alignment fasta file in the input directory

";

GetOptions (\%option, 'id=s', 'ref=s');

my $indir = $option{'id'} or die $usage;
my $reffile = $option{'ref'} or die $usage;
$reffile = $indir."/".$reffile;

my $scriptspath = dirname(__FILE__);

opendir DIR, $indir or die "couldn't open $indir: $!\n";
while (my $file = readdir DIR) {
	if ($file =~ /_NT_collapse\.fasta$/) {
		$file = $indir."/".$file;
		my $outfile = $file;
		$outfile =~s/\.fasta/_withRef.fasta/;
		print "=== Profile alignment on $reffile and $file ===\n";
		system("muscle -profile -quiet -in1 $reffile -in2 $file -out $outfile");
	}
}
closedir DIR;


