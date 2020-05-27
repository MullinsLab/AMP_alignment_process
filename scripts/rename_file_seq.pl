#!/usr/bin/perl -w

# rename file and sequence names 
# file: Hxxx_xxxx_xxx_REN_p.fasta --> Vxxx_xxxx_xxx_REN_NT.fasta
# sequence: Hxxx_xxxx_xxx_REN_pxxxxxxxx --> Vxxx_xxxx_xxx_REN_pblibxxxxxxxx
#           Hxxx_xxxx_xxx_REN_pe\d+ --> Vxxx_xxxx_xxx_REN_pbsga\d+
#           Hxxx_xxxx_xxx_REN_e\d+ --> Vxxx_xxxx_xxx_REN_sga\d+

use strict;

my $usage = "perl rename_file_seq.pl inFastaFile\n";
my $infasta = shift or die $usage;
my $outfile = $infasta;
$outfile =~ s/H(\d+)/V$1/;
my $wname = "";
my $count = 0;
my (%nameNames, %nameSeq, @names);

open IN, $infasta or die "couldn't open $infasta: $!\n";
open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
while (my $line = <IN>) {
	chomp $line;
	next if $line =~ /^\s*$/;	
	if ($line =~ />(.*)/) {
		$wname = $1;
		my @fields = split /\s+/, $wname;
		my $fname = shift @fields;
		my $name = "";
		if ($fname =~ /H(\d+)_(\d+)_(\d+)_(.*?)_pe(\d+)/) {
			$name = "V".$1."_".$2."_".$3."_".$4."_pbsga".$5;
		}elsif ($fname =~ /H(\d+)_(\d+)_(\d+)_(.*?)_p([A-Z]{8}_\d+)/) {
			$name = "V".$1."_".$2."_".$3."_".$4."_pblib".$5;
		}elsif ($fname =~ /H(\d+)_(\d+)_(\d+)_(.*?)_p([A-Z]{8})/) {
			$name = "V".$1."_".$2."_".$3."_".$4."_pblib".$5;
		}elsif ($fname =~ /H(\d+)_(\d+)_(\d+)_(.*?)_e(\d+)/) {
			$name = "V".$1."_".$2."_".$3."_".$4."_sga".$5;
		}else {
			$name = $fname;
		}
		if (@fields) {
			$name = $name." ".join(' ', @fields);
		}
		print OUT ">$name\n";
		++$count;
	}else {
		print OUT "$line\n";
	}
}
close IN;
close OUT;

print "Origianl file: $infasta, renamed file: $outfile, total $count sequences.\n";


