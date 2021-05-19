#!/usr/bin/perl -w


use strict;
use warnings;
use v5.10;

my $usage = "perl compare_alignments_for_batchinput.pl inSeqAlignmentFastaFile outputCSVFile directoryForOriginalFiles originalFileNameSuffix(including file extension)\n";
my $inFasta = shift or die $usage;
my $outfile = shift or die $usage;
my $dir = shift or die $usage;
my $suffix = shift or die $usage;
my $originfile = "";
my $studyid = my $subjectid = my $tpid = "";
if ($inFasta =~ /^V(\d+)_(\d+)_(.*?)_/) {
	$studyid = "V".$1;
	$subjectid = $2;
	$tpid = $3;
    $originfile = $dir."/V".$1."_".$2."_".$3."_".$suffix;
}elsif ($inFasta =~ /^RV(\d+)_(\d+)_(.*?)_/) {
	$studyid = "RV".$1;
	$subjectid = $2;
	$tpid = $3;
    $originfile = $dir."/RV".$1."_".$2."_".$3."_".$suffix;
}else {
	die "file name not formatted: $inFasta\n";
}

my $name = my $originname = "";
my $count = my $origincount = 0;
my (%nameSeq, %originnameSeq);
if (-e $originfile) {
	print "\n*** current $inFasta vs original $originfile ***\n";
	open ORIGIN, $originfile or die "couldn't open $originfile: $!\n";
	while (my $line = <ORIGIN>) {
		$line =~ s/\R$//;
		next if $line =~ /^\s*$/;
		if ($line =~ /^>(.*)/) {
			++$origincount;
			$originname = $1;
		}else {
			$line = uc $line;
			$originnameSeq{$originname} .= $line;
		}
	}
	close ORIGIN;
}else {
	my $id = $studyid."_".$subjectid."_".$tpid;
	opendir DIR, $dir or die "couldn't open dir $dir: $!\n";
	while (my $file = readdir DIR) {
		if ($file =~ /$suffix/) {
			if ($file =~ /$id/) {
				$originfile = $dir."/".$file;
				print "\n*** current $inFasta vs original $originfile ***\n";
				open ORIGIN, $originfile or die "couldn't open $originfile: $!\n";
				while (my $line = <ORIGIN>) {
					$line =~ s/\R$//;
					next if $line =~ /^\s*$/;
					if ($line =~ /^>(.*)/) {
						++$origincount;
						$originname = $1;
					}else {
						$line = uc $line;
						$originnameSeq{$originname} .= $line;
					}
				}
				close ORIGIN;
			}
		}
	}
}

open IN, $inFasta or die "couldn't open $inFasta: $!\n";
while (my $line = <IN>) {
	$line =~ s/\R$//;
	next if $line =~ /^\s*$/;	
	if ($line =~ />(.*)/) {
		++$count;
		$name = $1;
	}else {
		$line = uc $line;
		$nameSeq{$name} .= $line;
	}
}
close IN;

if ($count != $origincount) {
	print "number of sequences are not same. count $count, original count $origincount\n";
}

open OUT, ">>", $outfile or die "couldn't open $outfile: $!\n";
if (-z $outfile) {
	print OUT "original,count,current,count,unchanged,modified,missing,adding,flag\n";
}

my $missing = my $adding = my $modified = my $unchanged = my $flag = 0;
foreach my $originalname (keys %originnameSeq) {
	if ($nameSeq{$originalname}) {
		my $seq = $nameSeq{$originalname};
		my $originalseq = $originnameSeq{$originalname};
		if ($seq eq $originalseq) {
			++$unchanged;
		}else {
			++$modified;
		}
	}else {
		print "$originalname missing in $inFasta file\n";
		++$missing;
	}
}

foreach my $name (keys %nameSeq) {
	unless ($originnameSeq{$name}) {
		++$adding;
	}
}

if ($count == $origincount and $unchanged == $count) {
	$flag = 1;
}
print OUT "$originfile,$origincount,$inFasta,$count,$unchanged,$modified,$missing,$adding,$flag\n";
close OUT;

print "total $count sequences, $origincount original sequences, $unchanged unchanged, $modified modified, $missing missing, $adding added.\n";
