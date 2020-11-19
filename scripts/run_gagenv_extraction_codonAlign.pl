#!/usr/bin/perl -w

##########################################################################################
# Program: run_extraction_codonAlign_fuseBack.pl
# Purpose: In a directory of GP or REN sequence alignments with reference (HXB2) at the  
# top of alignment, extract gag/env of alignment according to the Gag/Rev start position  
# in HXB2, align extract sequences using muscle to make right aligned, codonAlign to 
# do codon alignment and insert codon aligned sequences into original alignment
# Author: Wenjie Deng
# Date: 2020-05-18
##########################################################################################

use strict;
use warnings;
use v5.10;
use Getopt::Long;
use File::Basename;

my %option = (
	'id' => '.',
	'sp' => 1,
);

my $usage = "\nusage: run_extraction_translation.pl [-option value]

options:  
-id     input directory with fasta files (default: . )
-sp     start position of Gag or Rev in HXB2 within reference alignment (default: 1, means alignment right starts HXB2's Gag or Rev)

";

GetOptions (\%option, 'id=s', 'sp=i');

my $indir = $option{'id'} or die $usage;
my $start = $option{'sp'} or die $usage;
my $scriptspath = dirname(__FILE__);
my $gagstart = 4;
my $gagend = 1506;
my $envstart = 256;
my $envend = 2826;

opendir DIR, $indir or die "couldn't open $indir: $!\n";
while (my $file = readdir DIR) {
	if ($file =~ /_withRef\.fasta$/) {
		if ($file =~ /\d+_(GP|REN)_NT/) {
			my $region = $1;
			my $gene = "";
			my $genestart = my $geneend = 0;
			my $infile = $indir."/".$file;
			
			if ($region eq "GP") {
				$gene = "Gag";
				$genestart = $start + $gagstart - 1;
				$geneend = $start + $gagend - 1;
			}elsif ($region eq "REN") {
				$gene = "Env";
				$genestart = $start + $envstart - 1;
				$geneend = $start + $envend - 1;
			}
			my $outdir = "$indir/$gene";
			unless (-e $outdir) {
				mkdir $outdir;
			}
			my $outfile = $outdir."/".$file;
			if ($outfile =~ /\d+_GP|REN_NT/) {
				$outfile =~ s/GP|REN/$gene/;
			}else {
				die "file name does not include GP or REN\n";
			}
			print "\n=== Extract $gene at HXB2 position of $genestart to $geneend in $infile ===\n";
			system("$scriptspath/extract_AMP_gene_for_pipeline.pl -ia $infile -oa $outfile -rs $genestart -re $geneend");
			my $muscleoutfile = my $codoninfile = $outfile;
			$muscleoutfile =~ s/\.fasta/_muscle.fasta/;
			print "== muscle aligning ==\n";
			system("muscle -quiet -in $outfile -out $muscleoutfile");
			# collapse sequences and move HXB2 at the top of alignment
			my $name = my $refname = my $refseq = "";
			my $noreffile = $muscleoutfile;
			$noreffile =~ s/_withRef//;
			my $collapsedfile = $noreffile;
			$collapsedfile =~ s/\.fasta/_collapsed.fasta/;
			open(IN, "<", $muscleoutfile) or die "couldn't Open $muscleoutfile: $!\n";
			open(OUT, ">", $noreffile) or die "couldn't Open $noreffile: $!\n";
			while (my $line = <IN>) {
                $line =~ s/\R$//;
				next if $line =~ /^\s*$/;
				if ($line =~ /^>(.*)/) {
					$name = $1;
					if ($name =~ /HXB2/) {
						$refname = $name;
					}else {
						print OUT ">$name\n";
					}
				}else {
					if ($name =~ /HXB2/) {
						$refseq .= $line;
					}else {
						print OUT "$line\n";
					}					
				}
            }
            close IN;
			close OUT;
			print "== collapse sequences in $noreffile to $collapsedfile ==\n";
			system("$scriptspath/collapse_seqs.pl -if $noreffile");
			
			$codoninfile =~ s/\.fasta/_collapsed.fasta/;
			my (@names, %nameseq);
			open(OUT, ">", $codoninfile) or die "could't open $codoninfile: $!\n";
			print OUT ">$refname\n$refseq\n";
			open IN, $collapsedfile or die "could't open $collapsedfile: $!\n";
			while (my $line = <IN>) {
				$line =~ s/\R$//;
				next if $line =~ /^\s*$/;
				print OUT "$line\n";
			}
			close IN;
			close OUT;
			print "== Codon alignment for $codoninfile ==\n";
			system("$scriptspath/codonalign_uw.pl $codoninfile 0 1 15 .");
		}		
	}
}
closedir DIR;


