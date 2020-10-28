#!/usr/bin/perl -w

################################################################################
# Program: count_CCS_reads_from_template_consensus_file.pl
# In the postproc directory, run the script to go through each samples directory
# to read the consensus fasta file, from the sequence names in fasta file to
# count the number of CCS reads and templates
# Author: Wenjie Deng
# Date: 2020-10-28
################################################################################

use strict;
use warnings;
use v5.10;
use Getopt::Long;
use File::Basename;

my %option = (
	'id' => '.',   	
);

my $usage = "\nusage: merge_alignments_by_muscle_profile.pl [-option value]

options:  
-id     input directory of porpid's postproc output consensus files (default: ./)

";

GetOptions (\%option, 'id=s');

my $indir = $option{'id'} or die $usage;
my $outfile = "template_ccs_summary.csv";
my $count = 0;
open(OUT, ">", $outfile) or die "couldn't Open $outfile: $!\n";
print OUT "File,Templates,CCS\n";
opendir INDIR, $indir or die "couldn't open $indir: $!\n";
while (my $dir = readdir INDIR) {
	unless ($dir =~ /^\./) {
		if (-d $dir) {
			++$count;
			$dir = $indir."/".$dir;
			opendir DIR, $dir or die "couldn't open $dir: $!\n";
			while (my $file = readdir DIR) {
				if ($file =~ /\.fasta$/) {
					my $templatecount = my $ccscount = 0;
					my $infile = $dir."/".$file;
                    open(IN, $infile) or die "couldn't open $infile: $!\n";
					while (my $line = <IN>) {
                        if ($line =~ /^>\S+\snumCCS=(\d+)/) {
                            $ccscount += $1;
							++$templatecount;
                        }   
                    }
					print OUT "$file,$templatecount,$ccscount\n";
				}	    
			}
			closedir DIR;
		}
    }
}
closedir INDIR;
close OUT;

print "\nTotal $count consensus files\n";
