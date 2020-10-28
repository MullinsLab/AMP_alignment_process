#!/usr/bin/perl -w

################################################################################
# Program: remove_sequences_from_list.pl
# Purpose: In a directory of sequence alignment, remove sequences in fasta files
# from a list
# Author: Wenjie Deng
# Date: 2020-10-26
################################################################################

use strict;
use warnings;
use v5.10;
use Getopt::Long;
use File::Basename;

my %option = (
	'id' => '.',
	'ls' => '',
);

my $usage = "\nusage: remove_sequences_from_list.pl [-option value]

options:  
-id     input directory with fasta files (default: . )
-ls     name of list file listing the sequence to be deleted

";

GetOptions (\%option, 'id=s', 'ls=s');

my $indir = $option{'id'} or die $usage;
my $lsfile = $option{'ls'} or die $usage;

my (%fileSeqnameStatus);
open(IN, $lsfile) or die "couldn't open $lsfile: $!\n";
while (my $line = <IN>) {
    $line =~ s/\R$//;
	next if ($line =~ /^\s*$/ or $line =~ /^SeqID/);
	my ($name, $seq) = split /,/, $line;
	if ($name =~ /^>(\S+)/) {
        $name = $1;
		my @fields = split /_/, $name;
		my $id = $fields[0]."_".$fields[1]."_".$fields[2]."_".$fields[3];
		$fileSeqnameStatus{$id}{$name} = 1;
    }else {
		die "sequence name $name not formatted\n";		
	}    
}
close IN;

opendir DIR, $indir or die "couldn't open $indir: $!\n";
while (my $file = readdir DIR) {
	if ($file =~ /\.fasta$/) {
		my @fields = split /_/, $file;
		my $id = $fields[0]."_".$fields[1]."_".$fields[2]."_".$fields[3];
		if ($fileSeqnameStatus{$id}) {
			print "\n=== processing $file ===\n";
            my $outfile = $file;
			$outfile =~ s/\.fasta/_updated.fasta/;
			my $flag = 1;
			my $count = my $removecount = my $keepcount = 0;
			open(IN, $file) or die "couldn't open $file: $!\n";
			open(OUT, ">", $outfile) or die "couldn't open $outfile: $!\n";
			while (my $line = <IN>) {
                $line =~ s/\R$//;
				if ($line =~ /^>(\S+)/) {
					++$count;
					$flag = 1;
					my $name = $1;
					if ($name =~ /\|/) {
                        my @fields = split /\|/, $name;
						$name = $fields[1];
                    }                   
                    if ($fileSeqnameStatus{$id}{$name}) {
                        $flag = 0;
						++$removecount;
                    }else {
						++$keepcount;
						print OUT ">$name\n";	
					}                  
                }elsif ($flag) {
					print OUT "$line\n";	
				}                
            }
			close IN;
			close OUT;
			print "Total $count sequences, removed $removecount, kept $keepcount sequences\n";
        }
	}
}
closedir DIR;


