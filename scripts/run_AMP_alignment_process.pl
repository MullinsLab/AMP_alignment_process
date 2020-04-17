#!/usr/bin/perl -w

##########################################################################################
# Program: run_AMP_alignment_process.pl
# Purpose: In a directory with alignment fasta files, renames file and sequence names,
# merges alignments based on file name (different time points/repeat runs of sequencing),
# collapses sequences into unique sequences, muscle align collapsed sequences, refines
# collapsed alignment
# Author: Wenjie Deng
# Date: 2020-04-08
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
				my $collapsedfile = my $alignedfile = my $orderedalignfile = my $refinedfile = $file;
				$collapsedfile =~ s/\.fasta/_collapsed.fasta/;
				$alignedfile =~ s/\.fasta/_collapsed_aligned.fasta/;
				$orderedalignfile =~ s/\.fasta/_collapsed_aligned_ordered.fasta/;
				$refinedfile =~ s/\.fasta/_collapse.fasta/;
				print "=== Align collapsed file $collapsedfile ===\n";
				system("muscle -quiet -in $collapsedfile -out $alignedfile");
				print "=== Order sequences in file $alignedfile ===\n";
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
					print OUT ">$name\n$nameSeq{$name}\n";
				}
				close OUT;				
				print "=== Refine alignment $orderedalignfile ===\n";
				system("perl $scriptspath/refineAlignment_left.pl -ia $orderedalignfile -uf");
				print "=== Verify sequences ===\n";
				system("perl $scriptspath/verify_seq_origin.pl $refinedfile $collapsedfile");
			}
		}
		closedir SUBDIR;
	}	
}
closedir DIR;
exit(0);



