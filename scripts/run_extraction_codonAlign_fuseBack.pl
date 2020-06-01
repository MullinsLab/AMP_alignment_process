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
	'sp' => 0,
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
my $gagstart = 1;
my $gagend = 1503;
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
			print "\n=== muscle aligning ==\n";
			system("muscle -quiet -in $outfile -out $muscleoutfile");
			# move HXB2 at the top of alignment
			$codoninfile =~ s/\.fasta/_muscle_input.fasta/;
			my $name = "";
			my (@names, %nameseq);
			open IN, $muscleoutfile or die "could't open $muscleoutfile: $!\n";
			while (my $line = <IN>) {
				$line =~ s/\R$//;
				next if $line =~ /^\s*$/;
				if ($line =~ /^>(.*)/) {
					$name = $1;
					if ($name =~ /HXB2/) {
						unshift @names, $name;
					}else {
						push @names, $name;
					}
				}else {
					$nameseq{$name} .= $line;
				}
			}
			close IN;
			open OUT, ">", $codoninfile or die "couldn't open $codoninfile: $!\n";
			foreach my $name (@names) {
				print OUT ">$name\n$nameseq{$name}\n";
			}
			close OUT;
			print "=== Codon alignment ===\n";
			system("$scriptspath/codonalign_uw.pl $codoninfile 0 1 15 .");
			# fuse back codon alignment to original alignment
			my %codonnameseq = ();
			$name = "";
			my $genecodonfile = $codoninfile;
			$genecodonfile =~ s/\.fasta/_codon_aligned_nt.fasta/;
			open CODON, $genecodonfile or die "couldn't open $genecodonfile: $!\n";
			while (my $line = <CODON>) {
				$line =~ s/\R$//;
				if ($line =~ /^>(.*)/) {
					$name = $1;
				}else {
					$codonnameseq{$name} .= $line;
				}
			}
			close CODON;
			my (%originalnameseq, @originalnames); 
			$name = "";
			open IN, $infile or die "couldn't open $infile: $!\n";
			while (my $line = <IN>) {
				$line =~ s/\R$//;
				if ($line =~ /^>(.*)/) {
					$name = $1;
					push @originalnames, $name;
				}else {
					$originalnameseq{$name} .= $line;
				}
			}
			close IN;
			my $start = my $end = 0;
			foreach my $name (@originalnames) {
				if ($name =~ /HXB2/) {
					my @refNas = split //, $originalnameseq{$name};
					my $idx = 0;
					for (my $i = 0; $i < @refNas; $i++) {
						if ($refNas[$i] =~ /[A-Za-z]/) {
							$idx++;
							if ($idx == $genestart) {
								$start = $i;
							}elsif ($idx == $geneend) {
								$end = $i;
								last;
							}
						}
					}
					last;
				}
			}			
			my $replacelength = $end - $start + 1;
			my $codonoutfile = $infile;
			$codonoutfile =~ s/\.fasta/_codon.fasta/;
			open OUT, ">", $codonoutfile or die "couldn't open $codonoutfile: $!\n";
			foreach my $name (@originalnames) {
				my $originalseq = $originalnameseq{$name};
				my $replaceseq = $codonnameseq{$name};
				my $finalseq = substr($originalseq, $start, $replacelength, $replaceseq);
				print OUT ">$name\n$originalseq\n";
			}
			close OUT;
		}		
	}
}
closedir DIR;


