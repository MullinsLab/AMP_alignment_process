#!/usr/bin/perl
#   Responsible Author: bkg@lanl.gov
#   Author: Brian Gaschen

=pod
*************************************************************************************
This script has been approved by LANL/LANS for open source release.
Include the following copyright and disclaimer. [JM 2012]

Copyright 2012.  Los Alamos National Security, LLC. This material was produced under
U.S. Government contract DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL),
which is operated by Los Alamos National Security, LLC for the U.S. Department of
Energy. The U.S. Government has rights to use, reproduce, and distribute this software.
NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS
OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is
modified to produce derivative works, such modified software should be clearly marked,
so as not to confuse it with the version available from LANL.

Additionally, this program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free Software
Foundation; version 2.0 of the License. Accordingly, this program is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details.
*************************************************************************************
=cut

# Modified by Wenjie Deng: takes .fasta as input sequence alignment file, writes NT codon
# alignment and AA alignment to files.

use strict;
use warnings;
use v5.10;

# File - alignment file (table format)
my $file  = $ARGV[0];
my $MasterSeqnum = $ARGV[1];
my $frame = $ARGV[2];
my $compensatingNum = $ARGV[3];
my $dir = $ARGV[4]; #tmp dir already created in cgi script - passed along to write log file to it
my ($id, $seq, @ids, @seqs, $ref_seq, $seqnum,$seqs,$codon_aligned,$aaseq);

open(FILE,"$file")|| die "Cannot open $file: $!\n";
my $log_file = $dir . "/log.txt";
open(LOG,">$log_file")|| die "Cannot open $log_file: $!\n";

# initialize codon hash

my ($codon_r,$aux) = init_codon_table(1);
my $seqcount = 0;
while (my $line = <FILE>) {
   $line =~ s/\R$//;
   if ($line =~ /^>(.*)/) {
   		push @ids, $1;
   		if ($seqcount) {
   			push @seqs, $seq;
   		}
   		$seq = "";   		
   		++$seqcount;
   }else {
   		$seq .= $line;
   }
#   ($id,$seq) = (/(\S+)\t(\S+)/);
#   push @ids, $id;
#   push @seqs, $seq;
}
push @seqs, $seq;
close FILE;

for( my $y = 0; $y < @ids; $y++){
    print LOG "\tIN    " . $seqs[$y] . "\n";
}

$seqnum=$MasterSeqnum;

#$frame = findFrame($alignment_ref,$ref_seq);
$seqs = \@seqs;
#$seqs = adjustFrame($seqs,$frame) if ($frame-1);
$seqs = adjustFrame($seqs, $frame);
$codon_aligned = codon_align($seqs, $seqnum);
$codon_aligned = pull_unaligned_left($seqnum, $codon_aligned, \@ids);

my $codon_aligned_aa_file = my $codon_aligned_nt_file = $file;
$codon_aligned_aa_file =~ s/\.fasta/_codon_aligned_aa.fasta/;
$codon_aligned_nt_file =~ s/\.fasta/_codon_aligned_nt.fasta/;
open AA, ">", $codon_aligned_aa_file or die "couldn't open $codon_aligned_aa_file: $!\n";
open NT, ">", $codon_aligned_nt_file or die "couldn't open $codon_aligned_nt_file: $!\n";
for (my $i = 0; $i <= $#ids; $i++) {
   $id = $ids[$i];	
   $seq = $codon_aligned -> [$i];
   $seq =~ s/Z/-/g;
   my $seq_len = length $seq;
 	if ($seq_len % 3 == 2){
		$seq = $seq."-";
	}elsif ($seq_len % 3 ==1){
		$seq = $seq."--";
	}
	
   $aaseq = na2aa(uc $seq,$codon_r);
   print NT ">$id\n$seq\n";
   print AA ">$id\n$aaseq\n";
#   print "$id\t$seq\t$aaseq\n";
}
close AA;
close NT;
close LOG;

print "Codon aligned $seqcount sequences\n";


############# Sub ############

#sub findFrame {
#	use strict;
#	my ($alignment_ref, $ref_seq) = @_;
#	my ($offset, $frame);
#	$alignment_ref = uc $alignment_ref;
#	$ref_seq = uc $ref_seq;
#	$ref_seq =~ s/-//g;
#	$alignment_ref =~ s/-//g;
#	$offset = index($ref_seq, $alignment_ref);
	
#    $frame = $offset % 3;
#    return($frame);
#}	 

#Takes an array of sequences and for each sequence:
#  1. Replace all leading gaps with 'Z'
#  2. Depending on frame # prepend 1 or 2 Z's to each seq
#Returns array of modified sequences.
sub adjustFrame {
	use strict;
	my ($seqs, $frame) = @_;
	my ($seq, $i);
	
	for $i (0.. $#{$seqs}) {
		$seq = $seqs -> [$i];
		#Replace all leading gaps with Z
   		$seq =~ s/^(-+)/"Z" x length($1)/eg;

		if ($frame == 2){
			#$seq =~ s/^-/Z/;
			#$seq = substr ($seq,0,1)."ZZ".substr ($seq,1);
			$seq = "ZZ".substr ($seq,0,1).substr ($seq,1);
		}
		if ($frame == 3){
			#$seq =~ s/^--/ZZ/;
			#$seq =~ s/^-/Z/;	
			#$seq = substr ($seq,0,2)."Z".substr ($seq,2);
			$seq = "Z".substr ($seq,0,2).substr ($seq,2);
		}
		#$seq =~ s/^(-*)/$1Z/ if $frame == 3;
		#$seq =~ s/^(-*)/$1ZZ/ if $frame == 2;
		#$seq = "--".$seq if $frame == 2;
		$seqs -> [$i] = $seq;
	}
        
	return $seqs;
}


sub init_codon_table {
   use strict;

   my ($translate_ambiguity) = @_;
   my (%codon,%aux);
   
   
   $codon{TTT} = "F";
   $codon{TTC} = "F";
   $codon{TTA} = "L";
   $codon{TTG} = "L"; 
   $codon{CTT} = "L";
   $codon{CTC} = "L";
   $codon{CTA} = "L";
   $codon{CTG} = "L"; 
   $codon{ATT} = "I";
   $codon{ATC} = "I";
   $codon{ATA} = "I";
   $codon{ATG} = "M"; 
   $codon{GTT} = "V";
   $codon{GTC} = "V";
   $codon{GTA} = "V";
   $codon{GTG} = "V"; 
   $codon{TCT} = "S";
   $codon{TCC} = "S";
   $codon{TCA} = "S";
   $codon{TCG} = "S"; 
   $codon{CCT} = "P";
   $codon{CCC} = "P";
   $codon{CCA} = "P";
   $codon{CCG} = "P"; 
   $codon{ACT} = "T";
   $codon{ACC} = "T";
   $codon{ACA} = "T";
   $codon{ACG} = "T"; 
   $codon{GCT} = "A";
   $codon{GCC} = "A";
   $codon{GCA} = "A";
   $codon{GCG} = "A"; 
   $codon{TAT} = "Y";
   $codon{TAC} = "Y";
   $codon{TAA} = "\*";
   $codon{TAG} = "\*"; 
   $codon{CAT} = "H";
   $codon{CAC} = "H";
   $codon{CAA} = "Q";
   $codon{CAG} = "Q"; 
   $codon{AAT} = "N";
   $codon{AAC} = "N";
   $codon{AAA} = "K";
   $codon{AAG} = "K"; 
   $codon{GAT} = "D";
   $codon{GAC} = "D";
   $codon{GAA} = "E";
   $codon{GAG} = "E"; 
   $codon{TGT} = "C";
   $codon{TGC} = "C";
   $codon{TGA} = "\*";
   $codon{TGG} = "W"; 
   $codon{CGT} = "R";
   $codon{CGC} = "R";
   $codon{CGA} = "R";
   $codon{CGG} = "R"; 
   $codon{AGT} = "S";
   $codon{AGC} = "S";
   $codon{AGA} = "R";
   $codon{AGG} = "R"; 
   $codon{GGT} = "G";
   $codon{GGC} = "G";
   $codon{GGA} = "G";
   $codon{GGG} = "G"; 

   $aux{A} = "A";
   $aux{C} = "C";
   $aux{T} = "T";
   $aux{G} = "G";

   $aux{R} = "AG";
   $aux{Y} = "TC";
   $aux{W} = "AT";
   $aux{S} = "GC";
   $aux{M} = "AC";
   $aux{K} = "GT";
   $aux{H} = "ATC";
   $aux{B} = "GCT";
   $aux{V} = "GAC";
   $aux{D} = "GAT";
   $aux{N} = "GATC";

   $codon{'---'} = "-";
   $codon{'???'} = "?";
   $codon{'555'} = "#";
   $codon{'666'} = "X"; 
   $codon{'AAR'} = "K" if $translate_ambiguity == 1;
   $codon{'AAY'} = "N" if $translate_ambiguity == 1;
   $codon{'TCR'} = "T" if $translate_ambiguity == 1;
   $codon{'ACY'} = "T" if $translate_ambiguity == 1;
   $codon{'ACK'} = "T" if $translate_ambiguity == 1;
   $codon{'ACM'} = "T" if $translate_ambiguity == 1;
   $codon{'ACB'} = "T" if $translate_ambiguity == 1;
   $codon{'ACR'} = "T" if $translate_ambiguity == 1;
   $codon{'ACD'} = "T" if $translate_ambiguity == 1;
   $codon{'ACV'} = "T" if $translate_ambiguity == 1;
   $codon{'ACW'} = "T" if $translate_ambiguity == 1;
   $codon{'ACH'} = "T" if $translate_ambiguity == 1;
   $codon{'ACN'} = "T" if $translate_ambiguity == 1;
   $codon{'ACS'} = "T" if $translate_ambiguity == 1;
   $codon{'ATY'} = "I" if $translate_ambiguity == 1;
   $codon{'ATM'} = "I" if $translate_ambiguity == 1;
   $codon{'ATW'} = "I" if $translate_ambiguity == 1;
   $codon{'ATH'} = "I" if $translate_ambiguity == 1;
   $codon{'AGR'} = "R" if $translate_ambiguity == 1;
   $codon{'AGY'} = "S" if $translate_ambiguity == 1;
   $codon{'GGR'} = "G" if $translate_ambiguity == 1;
   $codon{'GGY'} = "G" if $translate_ambiguity == 1;
   $codon{'GGK'} = "G" if $translate_ambiguity == 1;
   $codon{'GGM'} = "G" if $translate_ambiguity == 1;
   $codon{'GGB'} = "G" if $translate_ambiguity == 1;
   $codon{'GGD'} = "G" if $translate_ambiguity == 1;
   $codon{'GGV'} = "G" if $translate_ambiguity == 1;
   $codon{'GGW'} = "G" if $translate_ambiguity == 1;
   $codon{'GGH'} = "G" if $translate_ambiguity == 1;
   $codon{'GGN'} = "G" if $translate_ambiguity == 1;
   $codon{'GGS'} = "G" if $translate_ambiguity == 1;
   $codon{'GCR'} = "A" if $translate_ambiguity == 1;
   $codon{'GCY'} = "A" if $translate_ambiguity == 1;
   $codon{'GCK'} = "A" if $translate_ambiguity == 1;
   $codon{'GCM'} = "A" if $translate_ambiguity == 1;
   $codon{'GCB'} = "A" if $translate_ambiguity == 1;
   $codon{'GCD'} = "A" if $translate_ambiguity == 1;
   $codon{'GCV'} = "A" if $translate_ambiguity == 1;
   $codon{'GCW'} = "A" if $translate_ambiguity == 1;
   $codon{'GCH'} = "A" if $translate_ambiguity == 1;
   $codon{'GCN'} = "A" if $translate_ambiguity == 1;
   $codon{'GCS'} = "A" if $translate_ambiguity == 1;
   $codon{'GAR'} = "E" if $translate_ambiguity == 1;
   $codon{'GAY'} = "D" if $translate_ambiguity == 1;
   $codon{'GTR'} = "V" if $translate_ambiguity == 1;
   $codon{'GTY'} = "V" if $translate_ambiguity == 1;
   $codon{'GTK'} = "V" if $translate_ambiguity == 1;
   $codon{'GTM'} = "V" if $translate_ambiguity == 1;
   $codon{'GTB'} = "V" if $translate_ambiguity == 1;
   $codon{'GTD'} = "V" if $translate_ambiguity == 1;
   $codon{'GTV'} = "V" if $translate_ambiguity == 1;
   $codon{'GTW'} = "V" if $translate_ambiguity == 1;
   $codon{'GTH'} = "V" if $translate_ambiguity == 1;
   $codon{'GTN'} = "V" if $translate_ambiguity == 1;
   $codon{'GTS'} = "V" if $translate_ambiguity == 1;
   $codon{'CGR'} = "R" if $translate_ambiguity == 1;
   $codon{'CGY'} = "R" if $translate_ambiguity == 1;
   $codon{'CGK'} = "R" if $translate_ambiguity == 1;
   $codon{'CGM'} = "R" if $translate_ambiguity == 1;
   $codon{'CGB'} = "R" if $translate_ambiguity == 1;
   $codon{'CGD'} = "R" if $translate_ambiguity == 1;
   $codon{'CGV'} = "R" if $translate_ambiguity == 1;
   $codon{'CGW'} = "R" if $translate_ambiguity == 1;
   $codon{'CGH'} = "R" if $translate_ambiguity == 1;
   $codon{'CGN'} = "R" if $translate_ambiguity == 1;
   $codon{'CGS'} = "R" if $translate_ambiguity == 1;
   $codon{'MGR'} = "R" if $translate_ambiguity == 1;
   $codon{'MGA'} = "R" if $translate_ambiguity == 1;
   $codon{'MGG'} = "R" if $translate_ambiguity == 1;
   $codon{'CCR'} = "P" if $translate_ambiguity == 1;
   $codon{'CCY'} = "P" if $translate_ambiguity == 1;
   $codon{'CCK'} = "P" if $translate_ambiguity == 1;
   $codon{'CCM'} = "P" if $translate_ambiguity == 1;
   $codon{'CCB'} = "P" if $translate_ambiguity == 1;
   $codon{'CCD'} = "P" if $translate_ambiguity == 1;
   $codon{'CCV'} = "P" if $translate_ambiguity == 1;
   $codon{'CCW'} = "P" if $translate_ambiguity == 1;
   $codon{'CCH'} = "P" if $translate_ambiguity == 1;
   $codon{'CCN'} = "P" if $translate_ambiguity == 1;
   $codon{'CCS'} = "P" if $translate_ambiguity == 1;
   $codon{'CAR'} = "Q" if $translate_ambiguity == 1;
   $codon{'CAY'} = "H" if $translate_ambiguity == 1;
   $codon{'CTR'} = "L" if $translate_ambiguity == 1;
   $codon{'CTY'} = "L" if $translate_ambiguity == 1;
   $codon{'CTK'} = "L" if $translate_ambiguity == 1;
   $codon{'CTM'} = "L" if $translate_ambiguity == 1;
   $codon{'CTB'} = "L" if $translate_ambiguity == 1;
   $codon{'CTD'} = "L" if $translate_ambiguity == 1;
   $codon{'CTV'} = "L" if $translate_ambiguity == 1;
   $codon{'CTW'} = "L" if $translate_ambiguity == 1;
   $codon{'CTH'} = "L" if $translate_ambiguity == 1;
   $codon{'CTN'} = "L" if $translate_ambiguity == 1;
   $codon{'CTS'} = "L" if $translate_ambiguity == 1;
   $codon{'YTR'} = "L" if $translate_ambiguity == 1;
   $codon{'YTA'} = "L" if $translate_ambiguity == 1;
   $codon{'YTG'} = "L" if $translate_ambiguity == 1;
   $codon{'TGY'} = "C" if $translate_ambiguity == 1;
   $codon{'TCR'} = "S" if $translate_ambiguity == 1;
   $codon{'TCY'} = "S" if $translate_ambiguity == 1;
   $codon{'TCK'} = "S" if $translate_ambiguity == 1;
   $codon{'TCM'} = "S" if $translate_ambiguity == 1;
   $codon{'TCB'} = "S" if $translate_ambiguity == 1;
   $codon{'TCD'} = "S" if $translate_ambiguity == 1;
   $codon{'TCV'} = "S" if $translate_ambiguity == 1;
   $codon{'TCW'} = "S" if $translate_ambiguity == 1;
   $codon{'TCH'} = "S" if $translate_ambiguity == 1;
   $codon{'TCN'} = "S" if $translate_ambiguity == 1;
   $codon{'TCS'} = "S" if $translate_ambiguity == 1;
   $codon{'TAR'} = "*" if $translate_ambiguity == 1;
   $codon{'TAY'} = "Y" if $translate_ambiguity == 1;
   $codon{'TTR'} = "L" if $translate_ambiguity == 1;
   $codon{'TTY'} = "F" if $translate_ambiguity == 1;

   return (\%codon,\%aux);
 }


 
sub replace_char
{
   use strict;
	my ($seq,$pos,$val) = @_;
	
	my $old_char;
	my $str1;
	my $str2;
	
	
	$old_char = substr($seq,$pos,1);
	$str1 = substr($seq,0,$pos);
	$pos++;
	$str2 = substr($seq,$pos);
	$seq = $str1.$val.$str2;
	return ($seq, $old_char);
}
	
   
   
   
sub scan_seq
{
   use strict;
	my ($pos,$seq) = @_;
	
	my $i;
	my @seq = split (//,$seq);
	
	for ($i = $pos; $i <= $#seq; $i++)
	{
		last if $seq[$i] ne "-";
	}
	$i--;
	return $i;
}	  
	

sub get_codon
{
   use strict;
   my ($start,$seq) = @_;
   my $len;
   my $codon;
 
   $len = length $seq;
  
   if ($len > $start)
   {
      $codon = substr($seq,$start,3);
   }
   else
   {
     $codon = "";
   }
   
   return $codon;
}

			
sub replace_codon
{
   use strict;
	my ($codon,$start,$seq) = @_;
	my ($str1, $str2, $codon_end);
		
	$str1 = substr($seq,0,$start);
	$codon_end = $start + 3;
	$str2 = substr($seq,$codon_end);
	$seq = $str1.$codon.$str2;
	
	return $seq;
}			
	
				
	
  
sub make_proteins 
{
   use strict;
	my ($seq_r,$seq_start,$num_codons) = @_;
	my ($i,$seq,$frame,$len,$start,$codon_count,$str1,$str2);
	my ($codon, $next_start, $next_codon, $next_codon_len, $gaps,  $orig_codon_len );
	
        for ($i = 0; $i <= $#{$seq_r}; $i++) 	
	{
		$seq = $seq_r -> [$i];
		$len = length $seq;
		
		for ($start = $seq_start; $start < $len; $start += 3)
		{
			$codon = get_codon($start,$seq);
                        $codon_count = 0;

			if ($codon =~ /\-/ and $codon =~ /\w/)
			{
                                $next_start =$start + 3;
                                $next_codon = get_codon($next_start,$seq);
                                $next_codon = uc $next_codon;
                                $codon_count++ if $next_codon =~ /\w/; 
                                while ($next_codon =~ /\w{3}|-{3}/ and $codon_count <= $num_codons and $next_start < $len) 
                                {
                                     $next_start+= 3; 
                                     $next_codon = get_codon($next_start,$seq);
                                     $next_codon = uc $next_codon;
                                     $codon_count++ if $next_codon =~ /\w/;

                                 }
                                 if ($next_codon =~ /\w/ and $next_codon =~ /\-/) 
				 {
                                    $frame = squeeze_gaps($start,$next_start,$seq);
                                    $str1 = substr($seq,0,$start);
                                    $next_start += 3;
                                    $str2 = substr($seq,$next_start);
                                    $seq = $str1.$frame.$str2;
                                 }
                         }
                  }
                                   
                                           
                  $seq_r -> [$i] = $seq;
               }
	
	return $seq_r;
}
sub clean_ends
{
   use strict;
	my ($seq_r,$seq_start) = @_;
	my ($i,$start,$seq, $name, $codon,  $next_codon, $next_codon_len, $next_start, $gaps, $len, $orig_codon_len );
	
	for ($i = 0; $i <= $#{$seq_r}; $i++) 
	{
           $seq = $seq_r -> [$i];
	   $len = length $seq;
		
           for ($start = $seq_start; $start < $len; $start += 3)
	   {
			$codon = substr($seq,$start,3);
			if ($codon =~ /\-/ and $codon =~ /\w/)
			{
                                $next_start =$start + 3;
                                $next_codon = get_codon($next_start,$seq);
                                while ($next_codon !~ /\w/ and $next_start < $len) 
                                {
                                     $next_start+= 3; 
                                     $next_codon = get_codon($next_start,$seq);
            
                                 }
                                 $next_codon .= "-" while ((length $next_codon) < 3);
                                 if ($next_codon =~ /\w/ and $next_codon =~ /\-/) 
				 				{
                                        if ($next_codon =~ /\-\-\w/)
                                        {
                                           ($codon,$next_codon) = move_nucs($codon,$next_codon,1);
                                        }
                                        elsif ($next_codon =~ /\-\w\-/)
                                        {
                                           ($codon,$next_codon) = move_nucs($codon,$next_codon,1);
                                        }
                                        elsif ($next_codon =~ /\w\-\-/)
                                        {
                                           ($codon,$next_codon) = move_nucs($codon,$next_codon,1);
                                        }
                                        elsif ($next_codon =~ /\-\-\w/)
                                        {
                                           ($codon,$next_codon) = move_nucs($codon,$next_codon,1);
                                        }
                                        elsif ($next_codon =~ /\-\w\w/)
                                        {
                                           ($codon,$next_codon) = move_nucs($codon,$next_codon,2);
                                        }
                                        elsif ($next_codon =~ /\w\-\w/)
                                        {
                                           ($codon,$next_codon) = move_nucs($codon,$next_codon,1) if ($codon =~ /\w\w\-/ or $codon =~ /\w\-\w/ or $codon =~ /\-\w\w/);
                                           ($codon,$next_codon) = move_nucs($codon,$next_codon,2) if ($codon =~ /\-\-\w/ or $codon =~ /\-\w\-/ or $codon =~ /\w\-\-/);
                                        }
                                        elsif ($next_codon =~ /\w\w\-/)
                                        {
                                           ($codon,$next_codon) = move_nucs($codon,$next_codon,1) if ($codon =~ /\w\w\-/ or $codon =~ /\w\-\w/ or $codon =~ /\-\w\w/);
                                           ($codon,$next_codon) = move_nucs($codon,$next_codon,2) if ($codon =~ /\-\-\w/ or $codon =~ /\-\w\-/ or $codon =~ /\w\-\-/);
                                        }


                                        $seq = replace_codon($codon,$start,$seq);
                                        $seq = replace_codon($next_codon,$next_start,$seq);
                                        $start -= 3 if $codon =~ /\-/ and $next_start < ($len - 3);
                                 }
              }
          }
                                 
                                           
          $seq_r -> [$i] = $seq;
    }
	
	return $seq_r;
}
 	
sub move_nucs
{
   use strict;
   my ($codon,$next_codon,$action) = @_;
   my ($needed_nuc,$rem,$donor,$donor_len,$recipient,$recipient_len,$gaps,$needed_nucs,$nucs);

      $codon =~ s/-//g;
      $next_codon =~ s/-//g;
      $donor = $next_codon if $action == 1;
      $donor = $codon if $action == 2;
      $recipient = $next_codon if $action == 2;
      $recipient = $codon if $action == 1;
      $donor_len = length $donor;
      $recipient_len = length $recipient;
      $needed_nuc = 3 - $recipient_len;
      if ($donor_len > $needed_nuc)
      {
         $rem = $donor_len - $needed_nuc if $action == 2;
         $nucs = substr($donor,$rem,$needed_nuc) if $action == 2;
         $donor = substr($donor,0,$rem) if $action == 2;
         $nucs = substr($donor,0,$needed_nuc) if $action == 1;
         $donor = substr($donor,$needed_nuc) if $action == 1;
         $recipient = $recipient.$nucs if $action == 1;
         $recipient = $nucs.$recipient if $action == 2;
   
      }
      else
      {
         $recipient = $recipient.$donor if $action ==1;
         $recipient = $donor.$recipient if $action ==2;
         $donor = "---";
      }
      $recipient_len = length $recipient;
      $donor_len = length $donor;
      $gaps = 3-$recipient_len;
      $recipient = $recipient . "-" x $gaps if $gaps > 0;
      $gaps = 3 - $donor_len;
      $donor  = $donor . "-" x $gaps if $gaps > 0;
      if ($action == 1)
      {
         $next_codon = $donor;
         $codon = $recipient;
      }
      else
      {
         $codon = $donor;
         $next_codon = $recipient;
      }

      return ($codon,$next_codon);

}
sub insert_gap
{
   use strict;
   my ($seq, $pos) = @_;
   my ($str1, $str2, $end);

   $end = $pos + 1;
   $str1 = substr($seq, 0, $end);
   $str2 = substr($seq, $end);
   $seq = $str1 . "-" . $str2;

   return $seq;
}

sub na2aa {
 use strict;

my ($fullseq,$codon_r) = @_;

my ($i,$len_fullseq,$codon,$aa,$seq);

   $len_fullseq = length $fullseq;
   $seq = "";

   for ($i = 0; $i < $len_fullseq; $i += 3) {
      $codon = substr($fullseq,$i,3);
      $codon = '???' if $codon =~ /\?/;
      $codon.= "-" while ((length $codon) < 3);
      $codon = "666" if (!defined $codon_r -> {$codon}) and $codon !~ /-/;
      $codon = '555' if (!defined $codon_r -> {$codon}) and $codon =~ /-/;
      $aa = $codon_r -> {$codon};
      $seq = $seq.$aa;
   }

   return $seq;
}

sub pull_unaligned_left{
    use strict;
    my ($seqnum, $seq_r, $ids_r) = @_;
    print LOG "\n\n--->SUB PULL_UNALIGNED_LEFT with seqnum " . $seqnum . "\n"; 
    my ($seq,$stdseq,$i,$chr,$end_dash,$stdseq_name);
    my (@id,@seq,@stdseq);
    $stdseq = $seq_r -> [$seqnum];
    @stdseq = split("",$stdseq);
    $stdseq_name = $ids_r -> [$seqnum];
    print LOG "stdseq_name: $stdseq_name, stdseq: $stdseq\n";

    #for every character in the reference sequence
    for ($i = 0; $i <= $#stdseq; $i++) {
        $chr = $stdseq[$i];
        if ($chr eq "-") {
            $end_dash = find_end(\@stdseq,$i);
            print LOG "\t\tpos $i ($stdseq_name); end_dash $end_dash\t";
#            $seq_r = pull_left($seq_r,$i,$end_dash,$ids_r);
            $i = $end_dash;
        }
    }
    for( my $x = 0; $x < @$seq_r; $x++){
        print LOG "\tOUT:  " . $ {$seq_r}[$x] . "\n"; 
    }
    print LOG "\n";
	$seq_r = gap_squeeze($seq_r);
    for( my $x = 0; $x < @$seq_r; $x++){
        print LOG "\tOUT2: " . $ {$seq_r}[$x] . "\n"; 
    }
    return($seq_r);

}

sub pull_left {
   use strict;
   my ($seq_r,$start_dash,$end_dash,$ids_r) = @_;

   my ($len,$seq,$subseq,$sublen,$indx, $next_codon);

   $len = $end_dash - $start_dash + 1;
   print LOG "--->PULL_LEFT with length $len\t"; 
   #process every sequence
   for ($indx = 0; $indx <= $#{$seq_r}; $indx++) {
       $seq = $seq_r -> [$indx];
       $next_codon = substr($seq, $end_dash+1, 3);
       my $seq_name = $ids_r -> [$indx];
       print LOG "\t\tseq $indx ($seq_name) with next_codon: '$next_codon'\t";

       $len=$len+3 if ($next_codon=~/-/);
       $subseq = substr($seq, $start_dash, $len);
       print LOG "subseq: '$subseq'\t";
       $subseq =~ s/-//g;
       print LOG "subseq no '-': '$subseq'\t";
       $sublen = length $subseq;

       while ($sublen < $len) {
          $subseq = $subseq . "-";
          $sublen++;
       }
       print LOG "new subseq: '$subseq' and replacing with start=" . $start_dash . " and len=" . $len . "\t";       
       #replace the substring inside the currently processed sequence
       substr($seq, $start_dash, $len, $subseq);
       $seq_r -> [$indx] = $seq;
   }

   return $seq_r;
}

sub find_end {
   use strict;

   my ($stdseq_r,$indx) = @_;
   my ($char);

   $char = $stdseq_r -> [$indx];
   while ($char eq "-" and $indx <= $#{$stdseq_r}) {
       $indx++;
       $char = $stdseq_r -> [$indx];
   }
   $indx--;
  
   return $indx;
}
       

sub gap_squeeze {
use strict;

my ($seq_r) = @_;

my ($bases,$seq,$dash,$len,$longest,$shortest,$i,$j,$codon,$seq2);
my (@repseq);


$longest = 0;
$shortest = 10000000000000000000000000000000000000000000000;

for $i (0 .. $#{$seq_r}) {
    $seq = $seq_r -> [$i];
    $len = length $seq;
    if ($len > $longest) {
       $longest = $len;
    }
    if ($shortest > $len) {
        $shortest = $len;
    }
}

if ($longest != $shortest) {
     for $i (0 .. $#{$seq_r}) {
        $seq = $seq_r -> [$i];
        $len = length $seq;
        while ($longest > $len) {
             $seq = $seq."-";
             $len++;
        }
        $seq_r -> [$i] = $seq;
    }
}


for ($i = 0; $i <= $longest; $i++) {
    $repseq[$i] = 0;
}

for ($j = 0; $j <= $#{$seq_r}; $j++) {
   $seq = $seq_r -> [$j];
   for ($i = 0; $i < $longest; $i = $i+3) {
       $codon = substr($seq,$i,3);
       $repseq[$i] = 1 if $codon ne "---";
   }
}

$bases = 0;
for ($j = 0; $j <= $#{$seq_r}; $j++) {
   $seq = $seq_r -> [$j];
   for ($i = 0; $i < $longest; $i = $i+3) {
     $bases = $repseq[$i];
     if ($bases == 1) {
         $codon = substr($seq,$i,3);
         $seq2 = $seq2.$codon;
     }
   }
   $seq_r -> [$j] = $seq2;
   $seq2 = "";
}
    return $seq_r;
}
sub codon_align
{
   print LOG "--->SUB CODON_ALIGN\n"; 
   use strict;
   my ($seq_r,$seqnum) = @_; 
   my $name;
   my $seq;
   my $seq_length;
   my $start;
   my $offset;
   my $old_char;
   my $next_pos;
   my %gaps = ();
   my ($pos1, $pos2, $pos3);
   my ($char1,$char2,$char,$char3);
   my $char_val;
   my @cut_stdseq;
   my @seq_num = ();
   my %index = ();
   my $i;
   my $codon;
   my $insert;
   my $num = 0;
   my $stdseq_length;
   my $cut_stdseq;

   $offset = 0;
   $cut_stdseq = $seq_r -> [$seqnum];
   $stdseq_length = length $cut_stdseq;

   for ($start = 0; $start < $stdseq_length; $start= $start + 3){
       @cut_stdseq = ();
       @cut_stdseq = split(//,$cut_stdseq);
       $pos1 = $start;
       $pos2 = $start + 1;
       $pos3 = $start + 2;
       my $c1 = $cut_stdseq[$pos1];
       my $c2 = $cut_stdseq[$pos2];
       my $c3 = $cut_stdseq[$pos3];
       
       #Situation: -nn
       if ($cut_stdseq[$pos1] eq "-" and $cut_stdseq[$pos2] ne "-" and $cut_stdseq[$pos3] ne "-"){
           print LOG "\tSituation: -nn: " . $cut_stdseq[$pos1] . $cut_stdseq[$pos2] . $cut_stdseq[$pos3] . "\n";;
           ($cut_stdseq) = insert_gap($cut_stdseq, $pos1);
           ($cut_stdseq) = insert_gap($cut_stdseq, $pos1);
           @cut_stdseq = ();
           @cut_stdseq = split(//, $cut_stdseq);
           
           for ($i= 0; $i <= $#{$seq_r}; $i++){  
               $seq = $seq_r -> [$i];
               ($seq) = insert_gap($seq,$pos1);
               ($seq) = insert_gap($seq,$pos1);
               $seq_r -> [$i] = $seq;
           }
           $stdseq_length = length $cut_stdseq;
       }
       #Situation: --n
       if ($cut_stdseq[$pos1] eq "-" and $cut_stdseq[$pos2] eq "-" and $cut_stdseq[$pos3] ne "-"){
           print LOG "\tSituation: --n: " . $cut_stdseq[$pos1] . $cut_stdseq[$pos2] . $cut_stdseq[$pos3] . "\n";;
           ($cut_stdseq) = insert_gap($cut_stdseq,$pos2);
           @cut_stdseq = ();
           @cut_stdseq = split(//,$cut_stdseq);
           
           for ($i= 0; $i <= $#{$seq_r}; $i++){  
               $seq = $seq_r -> [$i];
               ($seq) = insert_gap($seq,$pos2);
               $seq_r -> [$i] = $seq;
           }
           $stdseq_length = length $cut_stdseq;
       }
       #Situation: n--       
       if ($cut_stdseq[$pos1] ne "-" and $cut_stdseq[$pos2] eq "-" and $cut_stdseq[$pos3] eq "-"){
           print LOG "\tSituation: n--: " . $cut_stdseq[$pos1] . $cut_stdseq[$pos2] . $cut_stdseq[$pos3] . "\n";;
           $codon = get_codon($pos1,$cut_stdseq);
           $insert = 0;
           ($cut_stdseq,$old_char) = replace_char($cut_stdseq,$pos1,"-");
           $next_pos = scan_seq($pos3,$cut_stdseq);
           $insert = 1 if ($next_pos == $pos3);
           ($cut_stdseq,$old_char) = replace_char($cut_stdseq,$next_pos,$old_char);
           
           if ($insert == 1){
               ($cut_stdseq,$old_char) = replace_char($cut_stdseq,$pos3,"-");
               ($cut_stdseq,$old_char) = replace_char($cut_stdseq,$pos1,$old_char);
               $next_pos = scan_seq($pos3,$cut_stdseq);
               $next_pos++;
               ($cut_stdseq,$old_char) = replace_char($cut_stdseq,$next_pos,"-");
               ($cut_stdseq,$old_char) = replace_char($cut_stdseq,$pos2,$old_char);
               $next_pos = scan_seq($pos3,$cut_stdseq);
               $next_pos++;
               ($cut_stdseq,$old_char) = replace_char($cut_stdseq,$next_pos,"-");
               ($cut_stdseq,$old_char) = replace_char($cut_stdseq,$pos3,$old_char);
           }
 
           @cut_stdseq = ();
           @cut_stdseq = split(//,$cut_stdseq);
           
           for ($i= 0; $i <= $#{$seq_r}; $i++){
               $seq = $seq_r -> [$i];
               $char3 = substr($seq, $pos3, 1);
               $char2 = substr($seq, $pos2, 1);
               $char1 = substr($seq, $pos1, 1);
               
               if ($char2 eq "-" and $char3 eq "-" and $char1 ne "-"){
                   if ($insert == 1){
                       $next_pos = scan_seq($pos3,$seq);
                       $next_pos++;
                       ($seq,$old_char) = replace_char($seq,$next_pos,"-");
                       ($seq,$old_char) = replace_char($seq,$pos2,$old_char);
                       $next_pos = scan_seq($pos3,$seq);
                       $next_pos++;
                       ($seq,$old_char) = replace_char($seq,$next_pos,"-");
                       ($seq,$old_char) = replace_char($seq,$pos3,$old_char);
                   }
                   else{ 
                       ($seq,$old_char) = replace_char($seq,$pos1,"-");
                       $next_pos = scan_seq($pos3,$seq);
                       ($seq,$old_char) = replace_char($seq,$next_pos, $old_char);
                   }
               }
               $seq_r -> [$i] = $seq;
           }
           $stdseq_length = length $cut_stdseq;
       }
       #Situation: nn-       
       if ($cut_stdseq[$pos1] ne "-" and $cut_stdseq[$pos2] ne "-" and $cut_stdseq[$pos3] eq "-"){
           print LOG "\tSituation: nn-: " . $cut_stdseq[$pos1] . $cut_stdseq[$pos2] . $cut_stdseq[$pos3] . "\n";;
           $codon = get_codon($pos1,$cut_stdseq);
           $next_pos = scan_seq($pos3,$cut_stdseq);
           $next_pos++;
           ($cut_stdseq,$old_char) = replace_char($cut_stdseq,$next_pos,"-");
           ($cut_stdseq,$old_char) = replace_char($cut_stdseq,$pos3,$old_char);
           @cut_stdseq = ();
           @cut_stdseq = split(//,$cut_stdseq);
           
           for ($i= 0; $i <= $#{$seq_r}; $i++){
               $seq = $seq_r -> [$i];
               $char3 = substr ($seq, $pos3, 1);
               $char2 = substr ($seq, $pos2, 1);
               $char1 = substr ($seq, $pos1, 1);
               
               if ($char2 ne "-" and $char3 eq "-" and $char1 ne "-"){
                   $next_pos = scan_seq($pos3,$seq);
                   $next_pos++;
                   ($seq,$old_char) = replace_char($seq,$next_pos,"-");
                   ($seq,$old_char) = replace_char($seq,$pos3,$old_char);
               }
               
               $seq_r -> [$i] = $seq;
           }
           $stdseq_length = length $cut_stdseq;
       }
       #Situation: -n-
       if ($cut_stdseq[$pos1] eq "-" and $cut_stdseq[$pos2] ne "-" and $cut_stdseq[$pos3] eq "-"){
           print LOG "\tSituation: -n-: " . $cut_stdseq[$pos1] . $cut_stdseq[$pos2] . $cut_stdseq[$pos3] . "\n";;
           $insert = 0;
           $codon = get_codon($pos1,$cut_stdseq);
           
           ($cut_stdseq,$old_char) = replace_char($cut_stdseq,$pos2,"-");
           $next_pos = scan_seq($pos3,$cut_stdseq);
           $insert = 1 if ($next_pos == $pos3);
           ($cut_stdseq,$old_char) = replace_char($cut_stdseq,$next_pos,$old_char);
           
           if ($insert == 1){
               ($cut_stdseq,$old_char) = replace_char($cut_stdseq,$pos3,"-");
               ($cut_stdseq,$old_char) = replace_char($cut_stdseq,$pos1,$old_char);
               $next_pos = scan_seq($pos3,$cut_stdseq);
               $next_pos++;
               ($cut_stdseq,$old_char) = replace_char($cut_stdseq,$next_pos,"-");
               ($cut_stdseq,$old_char) = replace_char($cut_stdseq,$pos2,$old_char);
               $next_pos = scan_seq($pos3,$cut_stdseq);
               $next_pos++;
               ($cut_stdseq,$old_char) = replace_char($cut_stdseq,$next_pos,"-");
               ($cut_stdseq,$old_char) = replace_char($cut_stdseq,$pos3,$old_char);
           }
           
           @cut_stdseq = ();
           @cut_stdseq = split(//,$cut_stdseq);
           
           for ($i= 0; $i <= $#{$seq_r}; $i++){
               $seq = $seq_r -> [$i];
               $char3 = substr ($seq, $pos3, 1);
               $char2 = substr ($seq, $pos2, 1);
               $char1 = substr ($seq, $pos1, 1);
               
               if ($char1 eq "-" and $char2 ne "-" and $char3 eq "-"){
                   if ($insert == 1){
                       ($seq,$old_char) = replace_char($seq,$pos2,"-");
                       ($seq,$old_char) = replace_char($seq,$pos1,$old_char);
                       $next_pos = scan_seq($pos3,$seq);
                       $next_pos++;
                       ($seq,$old_char) = replace_char($seq,$next_pos,"-");
                       ($seq,$old_char) = replace_char($seq,$pos2,$old_char);
                       $next_pos = scan_seq($pos3,$seq);
                       $next_pos++;
                       ($seq,$old_char) = replace_char($seq,$next_pos,"-");
                       ($seq,$old_char) = replace_char($seq,$pos3,$old_char);
                   }else{
                       ($seq,$old_char) = replace_char($seq,$pos2,"-");
                       $next_pos = scan_seq($pos3,$seq);
                       ($seq,$old_char) = replace_char($seq,$next_pos,$old_char);
                   }
               }
               $seq_r -> [$i] = $seq;
           }
           $stdseq_length = length $cut_stdseq;
       }
       #Situation: n-n       
       if ($cut_stdseq[$pos1] ne "-" and $cut_stdseq[$pos2] eq "-" and $cut_stdseq[$pos3] ne "-"){
           print LOG "\tSituation: n-n: " . $cut_stdseq[$pos1] . $cut_stdseq[$pos2] . $cut_stdseq[$pos3] . "\n";;
           $codon = get_codon($pos1,$cut_stdseq);
           ($cut_stdseq,$old_char) = replace_char($cut_stdseq,$pos3,"-");
           ($cut_stdseq,$old_char) = replace_char($cut_stdseq,$pos2,$old_char);
           $next_pos = scan_seq($pos3,$cut_stdseq);
           $next_pos++;
           ($cut_stdseq,$old_char) = replace_char($cut_stdseq,$next_pos,"-");
           ($cut_stdseq,$old_char) = replace_char($cut_stdseq,$pos3,$old_char);
           @cut_stdseq = ();
           @cut_stdseq = split(//,$cut_stdseq);
           
           for ($i= 0; $i <= $#{$seq_r}; $i++){
               $seq = $seq_r -> [$i];
               $char3 = substr ($seq, $pos3, 1);
               $char2 = substr ($seq, $pos2, 1);
               $char1 = substr ($seq, $pos1, 1);
               
               if ($char1 ne "-" and $char2 eq "-" and $char3 ne "-"){
                   ($seq,$old_char) = replace_char($seq,$pos3,"-");
                   ($seq,$old_char) = replace_char($seq,$pos2,$old_char);
                   $next_pos = scan_seq($pos3,$seq);
                   $next_pos++;
                   ($seq,$old_char) = replace_char($seq,$next_pos,"-");
                   ($seq,$old_char) = replace_char($seq,$pos3,$old_char);
               }
               $seq_r -> [$i] = $seq;       
           }
           $stdseq_length = length $cut_stdseq;
       }    
   }

   for( my $x = 0; $x < @$seq_r; $x++){
       print LOG "\tOUT:  " . $ {$seq_r}[$x] . "\n"; 
   }
   $seq_r = clean_ends($seq_r, 0);
   $seq_r = make_proteins($seq_r,0,$compensatingNum);
   
   return $seq_r;
}

 
sub squeeze_gaps
{
   use strict;
   my ($start,$next_start,$seq) = @_;
   my ($diff,$i,$from,$frame_len,$frame,$gaps);

   $next_start = $next_start + 3;
   $diff = $next_start - $start;
   $frame = substr($seq,$start,$diff);
   $frame =~ s/-//g;
   $frame_len = length $frame;
   $gaps = $diff - $frame_len;
   $frame = $frame."-"x$gaps;
   return $frame;
}
