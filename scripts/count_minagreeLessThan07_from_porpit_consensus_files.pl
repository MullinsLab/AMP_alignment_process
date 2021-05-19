#!/usr/bin/perl -w

################################################################################
# Program: count_minagreeLessThan07_from_porpid_consensus_files.pl
# In the propid directory, run the script to go through each samples directory
# to read the consensus fasta file, from the sequence names in fasta file to
# count the number of consensus that min_agreement < 0.7
# Author: Wenjie Deng
# Date: 2020-11-30
################################################################################

use strict;
use warnings;
use v5.10;
use Getopt::Long;
use File::Basename;

my %option = (
	'id' => '.',   	
);

my $usage = "\nusage: count_minagreeLessThan07_from_porpid_consensus_files.pl [-option value]

options:  
-id     input directory of porpid's postproc output consensus files (default: ./)

";

GetOptions (\%option, 'id=s');

my $indir = $option{'id'} or die $usage;
my $outfile = "min_agreement_summary.csv";
my $dircount = my $conscount = 0;
open(OUT, ">", $outfile) or die "couldn't Open $outfile: $!\n";
print OUT "File,min_agree>=0.7,min_agree<0.7,total\n";
opendir INDIR, $indir or die "couldn't open $indir: $!\n";
while (my $dir = readdir INDIR) {
	unless ($dir =~ /^\./) {
		if (-d $dir and $dir =~ /^20/ and $dir !~ /_old/) {
			++$dircount;
			my $sampledir = $indir."/".$dir;
			print "dir: $dir, sampledir: $sampledir\n";
			opendir SDIR, $sampledir or die "couldn't open $sampledir: $!\n";
			while (my $sdir = readdir SDIR) {
				unless ($sdir =~ /^\./) {
#					$sdir = $sampledir."/$sdir";
					print "sdir: $sdir\n";
					if ($sdir eq "consensus") {
						my $consdir = $sampledir."/$sdir";
						print "sdir: $sdir, consdir: $consdir\n";
						opendir CDIR, $consdir or die "couldn't open $consdir: $!\n";
						while (my $file = readdir CDIR) {
							unless ($file =~ /^\./) {
								my $infile = $consdir."/".$file;
	#							die "infile: $infile\n";
								if ($infile =~ /\.gz$/) {
									system("gunzip $infile");
									$infile =~ s/\.gz$//;
								}
								if ($infile =~ /\.fasta$/) {
									my $templatecount = my $less07count = my $great07count = 0;								
									open(IN, $infile) or die "couldn't open $infile: $!\n";
									while (my $line = <IN>) {
										$line =~ s/\R$//;
										if ($line =~ /^>(.*?)min_agreement=(\S+)/) {										
											my $minagree = $2;
											++$templatecount;
											if ($minagree < 0.7) {
												++$less07count;
											}elsif ($minagree >= 0.7) {
												++$great07count;
											}
										}   
									}
									++$conscount;
									close IN;
									print OUT "$file,$great07count,$less07count,$templatecount\n";
								}
							}							    
						}
						closedir CDIR;
					}
				}
			}		
			closedir SDIR;
		}
    }
}
closedir INDIR;
close OUT;

print "\nTotal $dircount sequencing runs, $conscount samples sequenced.\n";
