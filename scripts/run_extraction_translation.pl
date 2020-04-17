#!/usr/bin/perl -w

##########################################################################################
# Program: run_extraction_translation.pl
# Purpose: In a directory of sequence alignments with reference (HXB2) at the top of 
# alignment, extract the region of alignment according to the positions of reference, 
# translate the nucleotide sequences of the region into amino acid sequences
# Author: Wenjie Deng
# Date: 2020-04-13
##########################################################################################

use strict;
use Getopt::Long;
use File::Basename;

my %option = (
	'id' => '.',
	'rg' => '',
	'rs' => 0,
	're' => 0,
);

my $usage = "\nusage: run_extraction_translation.pl [-option value]

options:  
-id     input directory with fasta files (default: . )
-rg     name of extracted region (Gag, Pol, Rev, Vpu, Env, ...)
-rs     reference start position for extraction
-re     reference end position for extraction

";

GetOptions (\%option, 'id=s', 'rg=s', 'rs=i', 're=i');

my $indir = $option{'id'} or die $usage;
my $region = $option{'rg'} or die $usage;
my $rs = $option{'rs'} or die $usage;
my $re = $option{'re'} or die $usage;
my $scriptspath = dirname(__FILE__);
my $outdir = $indir."/".$region;
unless (-e $outdir) {
	mkdir $outdir;
}

opendir DIR, $indir or die "couldn't open $indir: $!\n";
while (my $file = readdir DIR) {
	if ($file =~ /_withRef\.fasta$/) {
		my $outfile = $outdir."/".$file;
		$file = $indir."/".$file;		
		if ($outfile =~ /\d+_(GP|REN)_NT/) {
			$outfile =~ s/GP|REN/$region/;
		}else {
			die "file name does not include 'GP or REN'\n";
		}
		print "\n=== Extract $region at reference position of $rs to $re in $file ===\n";
		system("$scriptspath/extract_AMP_gene_for_pipeline.pl -ia $file -oa $outfile -rs $rs -re $re");
		print "\n=== Translate nucleotide to amino acid sequences ===\n";
		system("$scriptspath/ntAlignment2aaAlignment.pl $outfile");
	}
}
closedir DIR;


