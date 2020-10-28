#!/usr/bin/perl -w

#########################################################################################################
# Program: expand_consalignment_to_individualalignments.pl
# Purpose: from an alignment of references and individual consensus sequences, expand the alignment to 
# each individual sequence alignment with references
# Author: Wenjie Deng
# Date: 2020-10-19
###########################################################################################################

use strict;
use v5.10;
use File::Copy;

my $usage = "perl expand_consalignment_to_individualalignments.pl reference_alignment consensus_alignment all_seq_file subdirectoryForIndividualSequenceAlignments\n";
my $refalignfile = shift or die $usage;
my $consalignfile = shift or die $usage;
my $seqfile = shift or die $usage;
my $outdir = shift || "Outputs";
my $subjectcount = my $flag = 0;
my $id = my $name = my $consname = my $AAalignseq = my $hxb2name = "";
my (%refstatus, @refnames, @ids, %idnames, %idstatus, @names, %nameAAs, %consnameSeq);
open REF, "<", $refalignfile or die "couldn't open $refalignfile: $!\n";
while (my $line = <REF>) {
	$line =~ s/\R$//;
	next if $line =~ /^\s*$/;
	if ($line =~ /^>(\S+)/) {
		my $refname = $1;
		$refstatus{$refname} = 1;
		if ($refname =~ /HXB2/) {
			$hxb2name = $refname;
		}else {
			push @refnames, $refname;
		}		
	}
}
close REF;
push @refnames, $hxb2name;

open NT, "<", $seqfile or die "couldn't open $seqfile: $!\n";
while (my $line = <NT>) {
	$line =~ s/\R$//;
	next if $line =~ /^\s*$/;
	if ($line =~ /^>(\S+)/) {
		$name = $1;
		if ($name =~ /^(.*?)_Env_AA/) {
			 my $id = $1;
			 unless ($idstatus{$id}) {
			 	$idstatus{$id} = 1;
			 	push @ids, $id;
			 }
			 push @{$idnames{$id}}, $name;
		}else {
			die "unformated name: $name\n";
		}
	}else {
		my @nts = split //, $line;
		my $len = length $line;
		for (my $i = 0; $i < $len; $i++) {
			$nameAAs{$name}{$i} = $nts[$i];
		}
	}
}
close NT;

open IN, "<", $consalignfile or die "couldn't open $consalignfile: $!\n";
while (my $line = <IN>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	if ($line =~ /^>(\S+)/) {
		$consname = $1;
	}else {
		$consnameSeq{$consname} .= $line;
	}
}
close IN;

unless (-e $outdir) {
	mkdir $outdir;
}

foreach my $id (@ids) {
	print "\n*** generating $id individual sequence alignment ***\n";
	my $outfile = $outdir."/".$id."_Env_AA_functional_collapsed_withRef.fasta";
	open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
	foreach my $ref (@refnames) {
		if ($consnameSeq{$ref}) {
			print OUT ">$ref\n$consnameSeq{$ref}\n";
		}else {
			die "No $ref in consensus alignment $consalignfile\n";
		}
	}
	my $AAalignseq = "";
	my $flag = 0;
	foreach my $consname (keys %consnameSeq) {
		if ($consname =~ /^$id/) {
			$AAalignseq = $consnameSeq{$consname};
			$flag = 1;
			last;
		}
	}
	if ($flag) {
		my $seqcount = 0;
		my @aas = split //, $AAalignseq;
		my $len = length $AAalignseq;
		foreach my $name (@{$idnames{$id}}) {
			++$seqcount;
			my $idx = -1;
			my @nts = ();
			for (my $i = 0; $i < $len; $i++) {
				if ($aas[$i] eq "-") {
					push @nts, "-";
				}else {
					++$idx;
					push @nts, $nameAAs{$name}{$idx};
				}
			}
			print OUT ">$name\n";
			print OUT join('', @nts),"\n";
		}
		print "expanded $seqcount sequences\n";	
	}else {
		print "!!!No consensus name with $id!!!\n";
	}	
	close OUT;	
}
