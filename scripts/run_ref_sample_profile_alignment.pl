#!/usr/bin/perl -w

################################################################################
# Program: run_ref_sample_profile_alignment.pl
# Purpose: In a directory with reference sequence alignment and sample alignment fasta 
# files, clean alignment by stripping columns with all gaps, do profile alignment between
# reference and sample alignment
# Author: Wenjie Deng
# Date: 2020-04-10
# Modified: 2020-10-08
################################################################################

use strict;
use warnings;
use v5.10;
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
my $outdir = $indir."/Alignments_withRef";
if (-e $outdir) {
	unlink $outdir;
}
mkdir $outdir;

my $scriptspath = dirname(__FILE__);

opendir DIR, $indir or die "couldn't open $indir: $!\n";
while (my $file = readdir DIR) {
	if ($file =~ /\.fasta$/) {
		my $gapstripfile = my $outfile = $outdir."/".$file;
		$gapstripfile =~ s/\.fasta/_stripgaps.fasta/;
		$outfile =~s/\.fasta/_withRef.fasta/;
		my $infile = $indir."/".$file;		
		print "\n=== Strip all gap columns in $infile ===\n";
		system ("$scriptspath/stripAllGaps.pl $infile $gapstripfile");		
		print "=== Profile alignment on $reffile and $gapstripfile ===\n";
		system("muscle -profile -quiet -in1 $gapstripfile -in2 $reffile -out $outfile");
	}
}
closedir DIR;


