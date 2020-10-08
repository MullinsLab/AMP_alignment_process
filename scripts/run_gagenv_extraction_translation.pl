#!/usr/bin/perl -w

##########################################################################################
# Program: run_gagenv_extraction_translation.pl
# Purpose: In a directory of GP or REN sequence alignments with reference (HXB2),   
# extract the gag/env region of alignment according to the positions of 
# reference, translate the nucleotide sequences of the region into amino acid sequences
# Author: Wenjie Deng
# Date: 2020-04-13
# Modified: 2020-04-17
# only provides the parameter of Gag/Rev start position of HXB2 in reference alignment,
# program will calculate the start and end position in alignment for each gene
# Modified: 2020-09-20
# extract region from start position only to the end of sequences
# retrieve functional protein sequences (sequence aligned length up to first stop >= 95% 
# of median alignment length and sequence length >= 80% of median sequence length in the
# alignment)
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

my $usage = "\nusage: run_gagenv_extraction_translation.pl [-option value]

options:  
-id     input directory with fasta files (default: . )
-sp     start position of Gag or Rev in HXB2 within reference alignment (default: 1, means alignment right starts HXB2's Gag or Rev)

";

GetOptions (\%option, 'id=s', 'sp=i');

my $indir = $option{'id'} or die $usage;
my $start = $option{'sp'} or die $usage;
my $scriptspath = dirname(__FILE__);
#my @GPgenes = qw(Gag Pol Prot);
#my @RENgenes = qw(Rev1 Rev2 Vpu Env Nef);
my @GPgenes = qw(Gag);
my @RENgenes = qw(Env);
my %GPgeneStart = (
	"Gag"  => 1,
	"Pol"  => 1296,
	"Prot" => 1464,
);
my %GPgeneEnd = (
	"Gag"  => 1503,
	"Pol"  => 4307,
	"Prot" => 1761,
);
my %RENgeneStart = (
	"Rev1" => 1,
	"Rev2" => 2410,
	"Vpu"  => 93,
	"Env"  => 256,
	"Nef"  => 2828,
);
my %RENgeneEnd = (
	"Rev1" => 76,
	"Rev2" => 2684,
	"Vpu"  => 341,
	"Env"  => 2826,
	"Nef"  => 3448,
);

my $filecount = 0;
opendir DIR, $indir or die "couldn't open $indir: $!\n";
while (my $file = readdir DIR) {
	if ($file =~ /_withRef\.fasta$/) {
		++$filecount;
		if ($file =~ /\d+_(GP|REN)_NT/) {
			my $region = $1;
			if ($region eq "GP") {
				foreach my $gene (@GPgenes) {
					my $outdir = "$indir/$gene";
					unless (-e $outdir) {
						mkdir $outdir;
					}
					my $genesummaryfile = "$outdir/$gene"."_summary.csv";
					my $genestart = $start + $GPgeneStart{$gene} - 1;
					my $geneend = $start +$GPgeneEnd{$gene} - 1;
					my $outfile = $outdir."/".$file;		
					if ($outfile =~ /\d+_GP_NT/) {
						$outfile =~ s/GP/$gene/;
					}else {
						die "file name does not include GP\n";
					}
					print "\n=== Extract $gene at HXB2 position of $genestart to the end of sequences in $file ===\n";
					#system("$scriptspath/extract_AMP_gene_for_pipeline.pl -ia $file -oa $outfile -rs $genestart -re $geneend");
					system("$scriptspath/extract_AMP_gene_for_pipeline_from_start_position.pl -ia $file -oa $outfile -rs $genestart");
					# translation
					my $aaalignfile = $outfile;
					$aaalignfile =~ s/_NT/_AA/;
					print "== Translate nucleotide to amino acid sequences for $outfile ==\n";
					system("$scriptspath/ntAlignment2aaAlignment.pl $outfile");
					# retrieve functional protein sequences and write summary
					if ($filecount == 1) {
						open GSF, ">", $genesummaryfile or die "couldn't open $genesummaryfile: $!\n";
						print GSF "file,sequences,functional,premature_stop,big_deletion,missing_5',not_start_M,median_sequence_length\n";
						close GSF;
					}
					print "== Retrieve functional amino acid sequences in $aaalignfile ==\n";
					system("$scriptspath/retrieve_functional_aa_seq.pl $aaalignfile $genesummaryfile");
					# remove HXB2 and strip all gap columns
					my $gapstripfile = my $functionalfile = $aaalignfile;
					$functionalfile =~ s/\.fasta/_functional.fasta/;
					$gapstripfile =~ s/_uncollapsed(.*?)withRef\.fasta/_functional.fasta/;
					print "== Remove HXB2 and strip all gap columns in $functionalfile ==\n";
					system("$scriptspath/stripAllGaps.pl $functionalfile $gapstripfile");
					# collapse functional protein sequence alignment
					print "== Collapse sequences in $gapstripfile ==\n";
					system("$scriptspath/collapse_seqs.pl -if $gapstripfile");
				}
			}elsif ($region eq "REN") {
				foreach my $gene (@RENgenes) {
					my $outdir .= "$indir/$gene";
					unless (-e $outdir) {
						mkdir $outdir;
					}
					my $genesummaryfile = "$outdir/$gene"."_summary.csv";
					my $genestart = $start + $RENgeneStart{$gene} - 1;
					my $geneend = $start +$RENgeneEnd{$gene} - 1;
					my $outfile = $outdir."/".$file;		
					if ($outfile =~ /\d+_REN_NT/) {
						$outfile =~ s/REN/$gene/;
					}else {
						die "file name does not include REN\n";
					}
					print "\n=== Extract $gene at HXB2 position of $genestart to the end of sequences in $file ===\n";
					#system("$scriptspath/extract_AMP_gene_for_pipeline.pl -ia $file -oa $outfile -rs $genestart -re $geneend");
					system("$scriptspath/extract_AMP_gene_for_pipeline_from_start_position.pl -ia $file -oa $outfile -rs $genestart");
					# translation
					my $aaalignfile = $outfile;
					$aaalignfile =~ s/_NT/_AA/;
					print "== Translate nucleotide to amino acid sequences for $outfile ==\n";
					system("$scriptspath/ntAlignment2aaAlignment.pl $outfile");
					# retrieve functional protein sequences and write summary
					if ($filecount == 1) {
						open GSF, ">", $genesummaryfile or die "couldn't open $genesummaryfile: $!\n";
						print GSF "file,sequences,functional,premature_stop,big_deletion,missing_5',not_start_M,median_sequence_length\n";
						close GSF;
					}
					print "== Retrieve functional amino acid sequences in $aaalignfile ==\n";
					system("$scriptspath/retrieve_functional_aa_seq.pl $aaalignfile $genesummaryfile");
					# remove HXB2 and strip all gap columns
					my $gapstripfile = my $functionalfile = $aaalignfile;
					$functionalfile =~ s/\.fasta/_functional.fasta/;
					$gapstripfile =~ s/_uncollapsed(.*?)withRef\.fasta/_functional.fasta/;
					print "== Remove HXB2 and strip all gap columns in $functionalfile ==\n";
					system("$scriptspath/stripAllGaps.pl $functionalfile $gapstripfile");
					# collapse functional protein sequence alignment
					print "== Collapse sequences in $gapstripfile ==\n";
					system("$scriptspath/collapse_seqs.pl -if $gapstripfile");
				}
			}
		}		
	}
}
closedir DIR;


