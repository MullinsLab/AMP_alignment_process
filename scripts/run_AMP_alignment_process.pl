#!/usr/bin/perl -w

##########################################################################################
# Program: run_AMP_alignment_process.pl
# Purpose: In a directory with alignment fasta files, renames file and sequence names,
# merges alignments based on file name (different time points/repeat runs of sequencing),
# collapses sequences into unique sequences, muscle align collapsed sequences, refines
# collapsed alignment
# Author: Wenjie Deng
# Date: 2020-04-08
# Modified (2020-04-28): reversed sequences before muscle alignment and then reversed back
# to make sequences left aligned
##########################################################################################

use strict;
use Getopt::Long;
use File::Basename;

my %option = (
	'id' => '.',
);

my $usage = "\nusage: run_AMP_alignment_process.pl [-option value]

options:  
-id		input directory with fasta files (default: . )
";

GetOptions (\%option, 'id=s');

my $indir = $option{'id'} or die $usage;
my $scriptspath = dirname(__FILE__);

print "\n=== Rename file and sequence names ===\n";
opendir DIR, $indir or die "couldn't open $indir: $!\n";
while (my $file = readdir DIR) {
	if ($file =~ /^H\d+_(.*?)\.fasta$/) {
		$file = $indir."/".$file;
		system("perl $scriptspath/rename_file_seq.pl $file");
	}
}
closedir DIR;

print "\n=== Merge alignments ===\n";
system("perl $scriptspath/merge_alignments.pl -id $indir");

opendir DIR, $indir or die "couldn't open $indir: $!\n";
while (my $name = readdir DIR) {
	if (-d $name and ($name eq "GP" or $name eq "REN")) {
		my $subdir = $indir."/".$name;
		opendir SUBDIR, $subdir or die "couldn't open $subdir: $!\n";
		while (my $file = readdir SUBDIR) {
			if ($file =~ /_NT\.fasta$/) {
				$file = $subdir."/".$file;
				print "\n=== Collapse sequences in $file ===\n";
				system("perl $scriptspath/collapse_seqs.pl -if $file");
				my $collapsedfile = my $reversedfile = my $alignedfile = my $orderedalignfile = $file;
				$collapsedfile =~ s/\.fasta/_collapsed.fasta/;
				$reversedfile =~ s/\.fasta/_collapsed_reversed.fasta/;
				$alignedfile =~ s/\.fasta/_collapsed_reversed_aligned.fasta/;
				$orderedalignfile =~ s/\.fasta/_collapse.fasta/;
				# reverse collapsed sequences to make left aligned after MSA
				open IN, $collapsedfile or die "couldn't open $collapsedfile: $!\n";
				open OUT, ">", $reversedfile or die "couldn't open $reversedfile: $!\n";
				while (my $line = <IN>) {
					chomp $line;
					next if $line =~ /^\s*$/;
					if ($line =~ /^>/) {
						print OUT "$line\n";
					}else {
						my $rvseq = reverse $line;
						print OUT "$rvseq\n";
					}
				}
				close IN;
				close OUT;
				print "=== Align collapsed file $collapsedfile ===\n";
				system("muscle -quiet -in $reversedfile -out $alignedfile");
				my $name = "";
				my %idxName = my %nameSeq = ();
				open IN, $alignedfile or die "couldn't open $alignedfile: $!\n";
				while (my $line = <IN>) {
					chomp $line;
					next if $line =~ /^\s*$/;
					if ($line =~ /^>(\S+)/) {
						$name = $1;
						if ($name =~ /_(\d+)_\d+$/) {
							$idxName{$1} = $name;
						}else {
							die "sequence name not formatted, maybe not collapsed\n";
						}
					}else {
						$nameSeq{$name} .= $line;
					}
				}
				close IN;
				open OUT, ">", $orderedalignfile or die "couldn't open $orderedalignfile: $!\n";
				foreach my $idx (sort {$a <=> $b} keys %idxName) {
					my $name = $idxName{$idx};
					my $seq = reverse $nameSeq{$name};
					print OUT ">$name\n$seq\n";
				}
				close OUT;
				print "=== Verify sequences ===\n";
				system("perl $scriptspath/verify_seq_origin.pl $orderedalignfile $collapsedfile");
			}
		}
		closedir SUBDIR;
	}	
}
closedir DIR;
exit(0);



