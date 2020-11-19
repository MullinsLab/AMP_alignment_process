#!/usr/bin/perl -w

#########################################################################################################
# Program: retrieve_only_HXB2_and_individual_seqs_from_individual_AA_alignment.pl
# Purpose: from an alignment of references and individual sequence alignment, retrieve only the HXB2 and 
# each individual sequences 
# Author: Wenjie Deng
# Date: 2020-11-02
###########################################################################################################

use strict;
use v5.10;
use File::Copy;

my $usage = "perl retrieve_only_HXB2_and_individual_seqs_from_individual_AA_alignment.pl reference_alignment individualAlignment\n";
my $refalignfile = shift or die $usage;
my $inalignfile = shift or die $usage;
print "$inalignfile\n";
my $outfile = $inalignfile;
$outfile =~ s/_withRef_jm/_withHXB2/;
my $outdir = "WithHXB2";
my $subjectcount = my $flag = 0;
my $id = my $name = my $consname = my $AAalignseq = my $hxb2name = "";
my (%refstatus);
open REF, "<", $refalignfile or die "couldn't open $refalignfile: $!\n";
while (my $line = <REF>) {
	$line =~ s/\R$//;
	next if $line =~ /^\s*$/;
	if ($line =~ /^>(\S+)/) {
		my $refname = $1;
		unless ($refname =~ /HXB2/) {
			$refstatus{$refname} = 1;
		}	
	}
}
close REF;

unless (-e $outdir) {
	mkdir $outdir;
}
$outfile = $outdir."/".$outfile;
open IN, "<", $inalignfile or die "couldn't open $inalignfile in 41: $!\n";
open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
while (my $line = <IN>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	if ($line =~ /^>(\S+)/) {
		$flag = 0;
		$name = $1;
		unless ($refstatus{$name}) {
			print OUT ">$name\n";
		}
	}elsif (!$refstatus{$name}) {
		print OUT "$line\n";
	}
}
close IN;
close OUT;
