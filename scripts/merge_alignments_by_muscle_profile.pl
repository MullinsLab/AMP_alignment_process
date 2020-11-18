#!/usr/bin/perl -w

################################################################################
# Program: merge_alignments_by_muscle_profile.pl
# Purpose: merge two alignments by muscle's profile-profile alignment
# Author: Wenjie Deng
# Date: 2020-10-27
################################################################################

use strict;
use warnings;
use v5.10;
use Getopt::Long;
use File::Basename;

my %option = (
	'id1' => '',
	'id2' => '',
	'sf'  => '',
);

my $usage = "\nusage: merge_alignments_by_muscle_profile.pl [-option value]

options:  
-id1     input directory of first alignments
-ld2     input directory of second alignments
-sf      outfile suffix (including file extension)

";

GetOptions (\%option, 'id1=s', 'id2=s', 'sf=s');

my $indir1 = $option{'id1'} or die $usage;
my $indir2 = $option{'id2'} or die $usage;
my $suffix = $option{'sf'} or die $usage;
my $count = 0;
my (@filenames, %idfilename);
opendir DIR, $indir1 or die "couldn't open $indir1: $!\n";
while (my $file = readdir DIR) {
	if ($file =~ /\.fasta$/) {
        push @filenames, $file;
    }
}
closedir DIR;

opendir DIR, $indir2 or die "couldn't open $indir2: $!\n";
while (my $file = readdir DIR) {
	if ($file =~ /\.fasta$/) {
        my @fields = split /_/, $file;
        my $id = $fields[0]."_".$fields[1];
        $idfilename{$id} = $file;
    }
}

foreach my $filename (@filenames) {
    ++$count;
    my @fields = split /_/, $filename;
    my $id = $fields[0]."_".$fields[1];
    my $firstTP = $fields[2];
    my $region1 = $fields[3];
    if ($idfilename{$id}) {
        my $firstfile = $indir1."/".$filename;
        my $secondfile = $indir2."/".$idfilename{$id};
        my @parts = split /_/, $idfilename{$id};
        my $secondTP = $parts[2];
        my $region2 = $parts[3];
        if ($region1 eq $region2) {
            print "\n=== Profile-profile alignment for $firstfile and $secondfile ===\n";
            my $outfile = $id."_".$firstTP."-".$secondTP."_".$region1."_".$suffix;
            system("muscle -profile -quiet -in1 $firstfile -in2 $secondfile -out $outfile");
        }else {
            die "region is different: $region1 vs. $region2\n";
        }      
    }else {
        print "No corresponding second alignment for $filename\n";
    }
    
}
closedir DIR;
print "\nTotal profile aligned $count files\n";

