#!/usr/bin/perl -w

##########################################################################################
# Program: run_extraction_translation.pl
# Purpose: In a directory of GP or REN sequence alignments with reference (HXB2) at the  
# top of alignment, extract the regions of alignment according to the positions of 
# reference, translate the nucleotide sequences of the region into amino acid sequences
# Author: Wenjie Deng
# Date: 2020-04-13
# Modified: 2014-04-17
# only provides the parameter of Gag/Rev start position of HXB2 in reference alignment,
# program will calculate the start and end position in alignment for each gene
##########################################################################################

use strict;
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
my @GPgenes = qw(Gag Pol Prot);
my @RENgenes = qw(Rev1 Rev2 Vpu Env Nef);
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

opendir DIR, $indir or die "couldn't open $indir: $!\n";
while (my $file = readdir DIR) {
	if ($file =~ /_withRef\.fasta$/) {
		if ($file =~ /\d+_(GP|REN)_NT/) {
			my $region = $1;
			if ($region eq "GP") {
				foreach my $gene (@GPgenes) {
					my $outdir = "$indir/$gene";
					unless (-e $outdir) {
						mkdir $outdir;
					}
					my $genestart = $start + $GPgeneStart{$gene} - 1;
					my $geneend = $start +$GPgeneEnd{$gene} - 1;
					my $outfile = $outdir."/".$file;		
					if ($outfile =~ /\d+_GP_NT/) {
						$outfile =~ s/GP/$gene/;
					}else {
						die "file name does not include GP\n";
					}
					print "\n=== Extract $gene at HXB2 position of $genestart to $geneend in $file ===\n";
					system("$scriptspath/extract_AMP_gene_for_pipeline.pl -ia $file -oa $outfile -rs $genestart -re $geneend");
					print "=== Translate nucleotide to amino acid sequences ===\n";
					system("$scriptspath/ntAlignment2aaAlignment.pl $outfile");
				}
			}elsif ($region eq "REN") {
				foreach my $gene (@RENgenes) {
					my $outdir .= "$indir/$gene";
					unless (-e $outdir) {
						mkdir $outdir;
					}
					my $genestart = $start + $RENgeneStart{$gene} - 1;
					my $geneend = $start +$RENgeneEnd{$gene} - 1;
					my $outfile = $outdir."/".$file;		
					if ($outfile =~ /\d+_REN_NT/) {
						$outfile =~ s/REN/$gene/;
					}else {
						die "file name does not include REN\n";
					}
					print "\n=== Extract $gene at HXB2 position of $genestart to $geneend in $file ===\n";
					system("$scriptspath/extract_AMP_gene_for_pipeline.pl -ia $file -oa $outfile -rs $genestart -re $geneend");
					unless ($gene eq "Rev1") {
						if ($gene eq "Rev2") {
							# merge Rev1 and Rev2 nucleotide sequences and translate to amino acid sequences
							my $revdir = "$indir/Rev";
							unless (-e $revdir) {
								mkdir $revdir;
							}
							my $rev1file = "$indir/Rev1/".$file;
							my $rev2file = "$indir/Rev2/".$file;
							$outfile = "$revdir/$file";
							$rev1file =~ s/_REN_/_Rev1_/;
							$rev2file =~ s/_REN_/_Rev2_/;
							$outfile =~ s/_REN_/_Rev_/;
							my $count1 = my $count2 = 0;
							my $name = "";
							my (@names, %rev1nameseq, %rev2nameseq);
							open REV1, $rev1file or die "couldn't open $rev1file: $!\n";
							while (my $line = <REV1>) {
								chomp $line;
								next if $line =~ /^\s*$/;
								if ($line =~ /^>(\S+)/) {
									$name = $1;
									push @names, $name;
									++$count1;
								}else {
									$rev1nameseq{$name} .= $line;
								}
							}
							close REV1;
							open REV2, $rev2file or die "couldn't open $rev2file: $!\n";
							while (my $line = <REV2>) {
								chomp $line;
								next if $line =~ /^\s*$/;
								if ($line =~ /^>(\S+)/) {
									$name = $1;
									++$count2;
								}else {
									$rev2nameseq{$name} .= $line;
								}
							}
							close REV2;
							if ($count1 != $count2) {
								die "number of sequences in $rev1file $count1 and $rev2file $count2 are not equal\n";
							}
							open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
							foreach my $name (@names) {
								if ($rev1nameseq{$name} and $rev2nameseq{$name}) {
									my $revseq = $rev1nameseq{$name}.$rev2nameseq{$name};
									print OUT ">$name\n$revseq\n";
								}else {
									die "missing one of Rev sequences for $name: $file\n";
								}
							}
						}
						print "=== Translate nucleotide to amino acid sequences ===\n";
						system("$scriptspath/ntAlignment2aaAlignment.pl $outfile");
					}					
				}
			}
		}		
	}
}
closedir DIR;


