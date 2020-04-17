#!/usr/bin/perl -w

################################################################################
# Author: Wenjie Deng
# modifiy date: 12-18-2009
# changes: handle three input file formats (Mac, Unix and Dos)
# modify date: 2010-08-18
# changes: add option of refining alignment of unique sequences with duplicte
# infomation at the end of sequence name.
# modify date: 2010-10-07
# changes: align homopolymers left to right
# modify date: 2011-02-28
# changes: modified the function of calculating frequency, now the frequency 
# calculation is based on previous one. save memory and much fast
# modified: create output refined alignment file based on input file name to 
# run batch inputs
# modified (2020-04-09): For AMP scripts_for_alignments_pipeline
# modified (2020-04-16): Fix a bug by not reverse sequence at beginning
################################################################################

use strict;
use Getopt::Long;

my %option = (
	'ia' => '',
	'dt' => 'nt',
	'ws' => 10,
	'ss' => 5,
	'ms' => 10,
	'mm' => -9,
	'gp' => -15,
	'bp' => 1,
	'ep' => 0,
	'uf' => '',
	'db' => ''
);

my $usage = "\nusage: perl refineAlignment_left.pl [-option value]

options:  
-ia		input alignment sequence fasta file
-dt		sequence data type (nt or aa. default: $option{dt})
-ws		window size (default: $option{ws})
-ss		stride size (default: $option{ss})
-ms		match score (default: $option{ms})
-mm		mismatch score (default: $option{mm})
-gp		gap penalty (default: $option{gp})
-bp		bigining position for refine (default: $option{bp} that is the beginnig of alignment)
-ep		ending position for refine (default: the end of alignment)
-uf		unique sequence flag (default: false). indicating if the sequences in alignment are unique,
		if true, the duplicate info must be at the end of sequence name
-db		output detailed information in .log file for debugging (default: false)
		
";

GetOptions (\%option, 'ia=s', 'dt=s', 'ws=i', 'ss=i', 'ms=f', 'mm=f', 'gp=f', 'bp=i', 'ep=i', 'uf', 'db');

my $startTime = time();
my $inFile = $option{'ia'} or die $usage;
my $dataType = $option{'dt'} or die $usage;
my $window = $option{'ws'};
my $overlap = $option{'ss'};
my $match = $option{'ms'};
my $misMatch = $option{'mm'};
my $gapPenalty = $option{'gp'};
my $start = $option{'bp'};
my $end = $option{'ep'};
my $uniqFlag = $option{'uf'};
my $debug = $option{'db'};
my $gapMatch = 0;

my $refinedFile = $inFile;
$refinedFile =~ s/_collapsed_aligned_ordered/_collapse/;

my %score = (
	'M' => $match,
	'R' => $misMatch,
	'I' => $gapPenalty,
	'D' => $gapPenalty
);

my @chars;
if ($dataType eq 'nt') {
	@chars = qw(A C G T -);
}elsif ($dataType eq 'aa') {
	@chars = qw(I L V F M C A G P T S Y W Q N H E D K R - *);
}else {
	die "data type must be nt or aa\n";
}

my ($name, $seq, @seqNames, %nameSeq);
my $totalSeqs = 0;
my $fileLines = GetFileLines($inFile);
foreach my $line (@$fileLines) {
	if ($line =~ /^>(\S+)/) {
		$name = $1;
		push @seqNames, $name;
		$seq = '';
		my $duplicates = 1;
		if ($uniqFlag) {
			$duplicates = GetDuplicates ($name);
		}
		$totalSeqs += $duplicates;
	}else {
		$nameSeq{$name} = '' if !$nameSeq{$name};
		$nameSeq{$name} .= uc $line;
	}
}


# check length of alignment 
my $flag = my $alignLen = 0;
foreach my $name (keys %nameSeq) {
	if (!$flag) {
		$alignLen = length $nameSeq{$name};
		$flag++;
	}else {
		if (length $nameSeq{$name} != $alignLen) {
			die "lengths of alignment are not same\n";
		}
	}
}

my ($seqStartPosRef, $seqEndPosRef) = GetPos (\%nameSeq, $alignLen);
$end = $alignLen unless $end;
my $orig_start = $start;
my $orig_end = $end;
my $refineLen = $end - $start + 1;
my $statFile = $refinedFile.".log";
if ($debug) {
	open OUT, ">", $statFile or die "couldn't open $statFile: $!\n";
}

print "Total sequences: $totalSeqs, alignment length: $refineLen\n";
for (my $i = $start - 1; $i < $end - $window + $overlap; $i += $overlap) {
	my ($actualWindow, %nameFullwindowseq, %namePrefixwindowseq, %nameSuffixwindowseq, %nameMidwindowseq, %nameWindowseq);
	my $position = $i + 1;
	my $totalWindowseqCount = 0;
	if ($debug) {
		print OUT "======== Refining window from position $position ========\n";
		print OUT "*** Window sequences ***\n";
	}	
	# get window sequences
	foreach my $seqName (@seqNames) {
		my $seq = $nameSeq{$seqName};
		my $windowseq;
		if ($i > $end - $window) {
			$windowseq = substr ($seq, $i, $end - $i);
		}else {
			$windowseq = substr ($seq, $i, $window);
		}		
		$actualWindow = length $windowseq;
		$windowseq = reverse $windowseq;
		unless ($windowseq =~ /^\-+$/) {	# skip window with all gaps
			if ($windowseq =~ /^\-/) {
				$namePrefixwindowseq{$seqName} = $windowseq;
			}elsif ($windowseq =~ /\-$/) {
				$nameSuffixwindowseq{$seqName} = $windowseq;
			}else {
				$nameFullwindowseq{$seqName} = $windowseq;
			}
		
			my $duplicates = 1;
			if ($uniqFlag) {
				$duplicates = GetDuplicates ($seqName);
			}
			$totalWindowseqCount += $duplicates;
			if ($debug) {
				print OUT "$windowseq\n";
			}			
		}		
	}
	if ($debug) {
		print OUT "Total window sequences: $totalWindowseqCount\n";
	}
	
	my ($seqCount, $alignSeqs, $fullseqAlignseq, $seqAlignseq, $prefixseqAlignseq, $suffixseqAlignseq, %prefixStatus, %suffixStatus);
	
	# align full window sequences first
	my $type = 'full';
	my ($uniqSeqs, $fullseqCount, $totalFullseqs) = GetFullUniqSeqs (\%nameFullwindowseq, $uniqFlag);

	if (@$uniqSeqs) {
		if ($debug) {
			print OUT "Total full window sequences: $totalFullseqs\n";
			if (@$uniqSeqs) {
				foreach my $uniqSeq (@$uniqSeqs) {
					print OUT "$uniqSeq: $fullseqCount->{$uniqSeq}\n";
				}
			}
		}
		my $mostFreqSeq = my $scdFreqSeq = my $mostFreqAlignSeq = my $scdFreqAlignSeq = '';
		if (@$uniqSeqs == 1) {	# only one unique full window sequence
			$mostFreqAlignSeq = $mostFreqSeq = shift @$uniqSeqs;
			push @$alignSeqs, $mostFreqSeq;
			$seqAlignseq->{$type}->{$mostFreqSeq} = $mostFreqSeq;
			$seqCount->{$type}->{$mostFreqSeq} = $fullseqCount->{$mostFreqSeq};
		}else {	# at least 2 unique full window sequences, first do pairwise alignment for the two most frequency sequences
			$mostFreqSeq = shift @$uniqSeqs;
			$scdFreqSeq = shift @$uniqSeqs;
		
			my @mostFreqNasNoGaps = split //, $mostFreqSeq;
			my @scdFreqNasNoGaps = split //, $scdFreqSeq;
			my $mostFreqSeqLen = length $mostFreqSeq;
			my $scdFreqSeqLen = length $scdFreqSeq;
			my $value = CalculateSimilarity (\@mostFreqNasNoGaps, \@scdFreqNasNoGaps, $mostFreqSeqLen, $scdFreqSeqLen, $match, $misMatch, $gapPenalty);						
			my $editTranscripts = Traceback (\@mostFreqNasNoGaps, \@scdFreqNasNoGaps, $value, $mostFreqSeqLen, $scdFreqSeqLen, $match, $misMatch, $gapPenalty);
			
			my $scores;						
			foreach my $et (@$editTranscripts) {
				$scores += $score{$et};
				if ($et eq 'I') {
					$mostFreqAlignSeq .= '-';
					$scdFreqAlignSeq .= shift @scdFreqNasNoGaps;
				}elsif ($et eq 'D') {
					$mostFreqAlignSeq .= shift @mostFreqNasNoGaps;
					$scdFreqAlignSeq .= '-';
				}else {
					$mostFreqAlignSeq .= shift @mostFreqNasNoGaps;
					$scdFreqAlignSeq .= shift @scdFreqNasNoGaps;
				}
			}
			$seqAlignseq->{$type}->{$mostFreqSeq} = $mostFreqAlignSeq;
			$seqAlignseq->{$type}->{$scdFreqSeq} = $scdFreqAlignSeq;
			$seqCount->{$type}->{$mostFreqAlignSeq} = $fullseqCount->{$mostFreqSeq};
			$seqCount->{$type}->{$scdFreqAlignSeq} = $fullseqCount->{$scdFreqSeq};
			push @$alignSeqs, $mostFreqAlignSeq, $scdFreqAlignSeq;
			die "problem in mostFreqNasNoGaps" if (@mostFreqNasNoGaps);
			die "problem in scdFreqNasNoGaps" if (@scdFreqNasNoGaps);
		}
		# calculate initail profile that is only the most frequent read
		my $profAlignLen = length $mostFreqAlignSeq;
		my @mostFreqAlignNas = split //, $mostFreqAlignSeq;
		my $profileFreq;
		my $profilePosNaCount;
		for (my $i = 1; $i <= $profAlignLen; $i++) {
			foreach my $char (@chars) {
				if ($mostFreqAlignNas[$i-1] eq $char) {
					$profileFreq->{$i}->{$char} = 1;
				}else {
					$profileFreq->{$i}->{$char} = 0;
				}				
			}
			$profilePosNaCount->{$i} = $seqCount->{full}->{$mostFreqAlignSeq};
		}
		while (@$uniqSeqs) {	# at least 3 unique sequences, do profile alignment			
			my $readSeq = shift @$uniqSeqs;
			my $readLen = length $readSeq;
			my @readNasNoGaps = split //, $readSeq;
			my @profileSeqs = @$alignSeqs;						
			my $profAlignLen = length $profileSeqs[0];
			my $lastAlignSeq = $profileSeqs[$#profileSeqs];	# last added aligned sequence in previous step, used to calculate the current profile
			CalculateFreq($profileFreq, $lastAlignSeq, $profilePosNaCount, \@chars, $seqCount, $profAlignLen, \%prefixStatus, \%suffixStatus);			
			if ($debug) {
				print OUT "Full\n";
				foreach my $char (@chars) {
					print OUT "$char";
					for (my $i = 1; $i <= $profAlignLen; $i++) {
						my $freq = $profileFreq->{$i}->{$char};
						print OUT "\t";
						printf OUT ("%.8f", $freq);
					}
					print OUT "\n";
				}
				print OUT "\n";
			}			
			my $value = CalculateProfileSimilarity ($profileFreq, \@readNasNoGaps, \@chars, $profAlignLen, $readLen, $match, $misMatch, $gapPenalty, $gapMatch, $type);
			my $editTranscripts = ProfileTraceback ($profileFreq, \@readNasNoGaps, \@chars, $value, $readLen, $profAlignLen, $match, $misMatch, $gapPenalty, $gapMatch, $type);		
			my $readAlignSeq = GetReadalignseq ($editTranscripts, \@readNasNoGaps);
			my $len = length $readAlignSeq;
			$seqCount->{$type}->{$readAlignSeq} = $fullseqCount->{$readSeq};
			$seqAlignseq->{$type}->{$readSeq} = $readAlignSeq;			
			# re-align profile sequences and re-assign frequency if alignment length changes
			if (my $count = grep /D/, @$editTranscripts) {
				RealignProfile (\@profileSeqs, $len, $alignSeqs, $editTranscripts, $seqCount, $seqAlignseq, $profilePosNaCount, \%prefixStatus, \%suffixStatus);
				$profileFreq = ReassignProfFreq ($profileFreq, $editTranscripts, \@chars);				
			}
			push @$alignSeqs, $readAlignSeq;
		}	
		if ($debug) {						
			print OUT "*** Refined full sequences ***\n";
			print OUT join ("\n", @$alignSeqs), "\n";
		}				

		if (%namePrefixwindowseq) {
			my $type = 'prefix';
			# get unique prefix 
			my ($uniqSeqs, $prefixseqCount, $totalPrefixseqs) = GetUniqSeqs (\%namePrefixwindowseq, $uniqFlag);
			if ($debug) {
				print OUT "Total prefix window sequences: $totalPrefixseqs\n";
				if (@$uniqSeqs) {
					foreach my $uniqSeq (@$uniqSeqs) {
						print OUT "$uniqSeq: $prefixseqCount->{$uniqSeq}\n";
					}
				}
			}
			while (@$uniqSeqs) {	# do profile alignment for prefix
				my $readSeq = shift @$uniqSeqs;
				my $readLen = length $readSeq;
				my @readNasNoGaps = split //, $readSeq;
				my @profileSeqs = @$alignSeqs;	# array holds previously aligned sequences to calculate the frequency				
				my $profAlignLen = length $profileSeqs[0];
				unless (scalar @profileSeqs == 1) {	# at least two sequences, otherwise use the initial profile
					my $lastAlignSeq = $profileSeqs[$#profileSeqs];	# last added aligned sequence in previous step, used to calculate the current profile
					CalculateFreq($profileFreq, $lastAlignSeq, $profilePosNaCount, \@chars, $seqCount, $profAlignLen, \%prefixStatus, \%suffixStatus);
				}
				if ($debug) {
					print OUT "Prefix\n";
					foreach my $char (@chars) {
						print OUT "$char";
						for (my $i = 1; $i <= $profAlignLen; $i++) {
							my $freq = $profileFreq->{$i}->{$char};
							print OUT "\t";
							printf OUT ("%.8f", $freq);
						}
						print OUT "\n";
					}
					print OUT "\n";
				}
				my $value = CalculateProfileSimilarity ($profileFreq, \@readNasNoGaps, \@chars, $profAlignLen, $readLen, $match, $misMatch, $gapPenalty, $gapMatch, $type);
				my $editTranscripts = ProfileTraceback ($profileFreq, \@readNasNoGaps, \@chars, $value, $readLen, $profAlignLen, $match, $misMatch, $gapPenalty, $gapMatch, $type);
				my $readAlignSeq = GetReadalignseq ($editTranscripts, \@readNasNoGaps);				
				$seqCount->{prefix}->{$readAlignSeq} = $prefixseqCount->{$readSeq};
				$seqAlignseq->{prefix}->{$readSeq} = $readAlignSeq;
				my $len = length $readAlignSeq;
				for (my $i = 0; $i < $len; $i++) {	# count leading gaps for prefix sequence
					if (substr ($readAlignSeq, $i, 1) ne '-') {
						$prefixStatus{$readAlignSeq} = $i;
						last;
					}
				}
				# re-align profile sequences and re-assign frequency if alignment length changes
				if (my $count = grep /D/, @$editTranscripts) {
					RealignProfile (\@profileSeqs, $len, $alignSeqs, $editTranscripts, $seqCount, $seqAlignseq, $profilePosNaCount, \%prefixStatus, \%suffixStatus);
					$profileFreq = ReassignProfFreq ($profileFreq, $editTranscripts, \@chars);				
				}
				push @$alignSeqs, $readAlignSeq;	# put aligned read sequence into new alignment 
			}
			if ($debug) {
				print OUT "*** Refined full+prefix sequences ***\n";
				print OUT join ("\n", @$alignSeqs), "\n";	
			}	
		}
				
		if (%nameSuffixwindowseq) {
			my $type = 'suffix';
			my ($uniqSeqs, $suffixseqCount, $totalSuffixseqs) = GetUniqSeqs (\%nameSuffixwindowseq, $uniqFlag);
			if ($debug) {
				print OUT "Total suffix window sequences: $totalSuffixseqs\n";
				if (@$uniqSeqs) {
					foreach my $uniqSeq (@$uniqSeqs) {
						print OUT "$uniqSeq: $suffixseqCount->{$uniqSeq}\n";
					}
				}
			}			
			while (@$uniqSeqs) {	# do profile alignment for prefix				
				my $readSeq = shift @$uniqSeqs;
				my $readLen = length $readSeq;
				my @readNasNoGaps = split //, $readSeq;
				my @profileSeqs = @$alignSeqs;				
				my $profAlignLen = length $profileSeqs[0];
				unless (scalar @profileSeqs == 1) {	# at least two sequences, otherwise use the initial profile
					my $lastAlignSeq = $profileSeqs[$#profileSeqs];	# last added aligned sequence in previous step, used to calculate the current profile
					CalculateFreq($profileFreq, $lastAlignSeq, $profilePosNaCount, \@chars, $seqCount, $profAlignLen, \%prefixStatus, \%suffixStatus);
				}
				if ($debug) {
					print OUT "Suffix\n";
					foreach my $char (@chars) {
						print OUT "$char";
						for (my $i = 1; $i <= $profAlignLen; $i++) {
							my $freq = $profileFreq->{$i}->{$char};
							print OUT "\t";
							printf OUT ("%.8f", $freq);
						}
						print OUT "\n";
					}
					print OUT "\n";
				}
				# reverse the profile frequency
				my $reversedProfFreq = RevProfFreq ($profileFreq, $profAlignLen, \@chars);
				if ($debug) {
					print OUT "Reverse\n";
					foreach my $char (@chars) {
						print OUT "$char";
						for (my $i = 1; $i <= $profAlignLen; $i++) {
							my $freq = $reversedProfFreq->{$i}->{$char};
							print OUT "\t";
							printf OUT ("%.8f", $freq);
						}
						print OUT "\n";
					}
					print OUT "\n";
				}
				my @reverseReadNasNoGaps = reverse @readNasNoGaps;
				my $value = CalculateProfileSimilarity ($reversedProfFreq, \@reverseReadNasNoGaps, \@chars, $profAlignLen, $readLen, $match, $misMatch, $gapPenalty, $gapMatch, $type);
				my $editTranscripts = ProfileTraceback ($reversedProfFreq, \@reverseReadNasNoGaps, \@chars, $value, $readLen, $profAlignLen, $match, $misMatch, $gapPenalty, $gapMatch, $type);
				
				# reverse back
				@$editTranscripts = reverse @$editTranscripts;				
				my $readAlignSeq = GetReadalignseq ($editTranscripts, \@readNasNoGaps);	
				$seqCount->{suffix}->{$readAlignSeq} = $suffixseqCount->{$readSeq};
				$seqAlignseq->{suffix}->{$readSeq} = $readAlignSeq;
				my $len = length $readAlignSeq;
				for (my $i = $len - 1; $i >= 0; $i--) {	# count ending gaps for suffix sequence
					if (substr ($readAlignSeq, $i, 1) ne '-') {
						$suffixStatus{$readAlignSeq} = $len - $i - 1;
						last;
					}
				}				
				# re-align profile sequences and re-assign frequency if alignment length changes
				if (my $count = grep /D/, @$editTranscripts) {
					RealignProfile (\@profileSeqs, $len, $alignSeqs, $editTranscripts, $seqCount, $seqAlignseq, $profilePosNaCount, \%prefixStatus, \%suffixStatus);
					$profileFreq = ReassignProfFreq ($profileFreq, $editTranscripts, \@chars);				
				}
				push @$alignSeqs, $readAlignSeq;
			}
		}
		if ($debug) {
			print OUT "*** Refined all sequences in the window ***\n";
			print OUT join ("\n", @$alignSeqs), "\n";	
		}	
	
		my $newAlignLen = length $alignSeqs->[0];
		foreach my $seqName (@seqNames) {
			my $seq = $nameSeq{$seqName};
			my $newAlignSeq = '';		
			if ($nameFullwindowseq{$seqName}) {
				my $windowseq = my $pureWindowseq = $nameFullwindowseq{$seqName};
				$pureWindowseq =~ s/\-//g;
				$newAlignSeq = $seqAlignseq->{full}->{$pureWindowseq};
				$seqEndPosRef->{$seqName} = $seqEndPosRef->{$seqName} + $newAlignLen - $actualWindow;
			}elsif ($namePrefixwindowseq{$seqName}) {
				my $windowseq = my $pureWindowseq = $namePrefixwindowseq{$seqName};
				$pureWindowseq =~ s/\-//g;
				$newAlignSeq = $seqAlignseq->{prefix}->{$pureWindowseq};
				$seqStartPosRef->{$seqName} = $i + $prefixStatus{$newAlignSeq};
				$seqEndPosRef->{$seqName} = $seqEndPosRef->{$seqName} + $newAlignLen - $actualWindow;
			}elsif ($nameSuffixwindowseq{$seqName}) {
				my $windowseq = my $pureWindowseq = $nameSuffixwindowseq{$seqName};
				$pureWindowseq =~ s/\-//g;
				$newAlignSeq = $seqAlignseq->{suffix}->{$pureWindowseq};
				$seqEndPosRef->{$seqName} = $i + $newAlignLen - $suffixStatus{$newAlignSeq} - 1;
			}else {	# all gaps
				for (my $i = 0; $i < $newAlignLen; $i++) {
					$newAlignSeq .= '-';
				}
				if ($i < $seqStartPosRef->{$seqName} && $newAlignLen != $actualWindow) {	# gaps are before start position, will affect the newly aligned start and ending position
					$seqStartPosRef->{$seqName} = $seqStartPosRef->{$seqName} + $newAlignLen - $actualWindow;
					$seqEndPosRef->{$seqName} = $seqEndPosRef->{$seqName} + $newAlignLen - $actualWindow;
				}
			}
			$newAlignSeq = reverse $newAlignSeq;			
			substr ($seq, $i, $actualWindow, $newAlignSeq);
			$nameSeq{$seqName} = $seq;
		}
		$end = $end - $actualWindow + $newAlignLen;
	}	
}
if ($debug) {
	close OUT;
}
my $refineAlignLen = $end - $start + 1;
print "Alignment length after refinement: $refineAlignLen\n";

open REFINE, ">$refinedFile" or die "couldn't open $refinedFile: $!\n";
foreach my $seqName (@seqNames) {
	print REFINE ">$seqName\n";
	print REFINE "$nameSeq{$seqName}\n";
}
close REFINE;

my $finishTime = time();
my $duration = $finishTime - $startTime;


sub GetFileLines {
	my $file = shift;
	my $line = "";
	open IN, $file or die "couldn't open $file: $!\n";
	my @buffer = <IN>;
	close IN;
 	foreach my $element (@buffer) {
 		$line .= $element;
 	}
 	if ($line =~ /\r\n/) {
		$line =~ s/\r//g;
	}elsif ($line =~ /\r/) {
		$line =~ s/\r/\n/g;
	}
	my @fileLines = split /\n/, $line;
	return \@fileLines;
}

sub GetPos {
	my ($nameSeqRef, $len) = @_;
	my (%seqStartPos, %seqEndPos);
	foreach my $name (keys %$nameSeqRef) {
		my $seq = $nameSeqRef->{$name};
		for (my $i = 0; $i < $len; $i++) {
			if (substr ($seq, $i, 1) ne '-') {
				$seqStartPos{$name} = $i;
				last;
			}
		}
		for (my $i = $len - 1; $i >= 0; $i--) {
			if (substr ($seq, $i, 1) ne '-') {
				$seqEndPos{$name} = $i;
				last;
			}
		}
	}
	return (\%seqStartPos, \%seqEndPos);
}


sub GetDuplicates {
	my $seqName = shift;
	my $duplicates;
	if ($seqName =~ /(\d+)$/) {
		$duplicates = $1;
	}else {
		die "there is no duplicate information at the end of sequence name\n";
	}
	return $duplicates;
}


sub GetFullUniqSeqs {
	my $nameWindowseqRef = shift;
	my $uniqFlag = shift;
	my (%fullseqCount, %seqCount, @sortUniqSeqs);
	my $count = my $len = my $flag = 0;
	foreach my $name (keys %$nameWindowseqRef) {
		my $alignSeq = my $pureSeq = $nameWindowseqRef->{$name};
		if (!$count) {
			$len = length $alignSeq;
		}
		my $duplicates = 1;
		if ($uniqFlag) {
			$duplicates = GetDuplicates ($name);
		}
		$count += $duplicates;
		if ($alignSeq !~ /\-/) {
			if (!$fullseqCount{$alignSeq}) {
				$fullseqCount{$alignSeq} = 0;
			}
			$fullseqCount{$alignSeq} += $duplicates;
		}
		$pureSeq =~ s/\-//g;			
		if (!$seqCount{$pureSeq}) {
			$seqCount{$pureSeq} = 0;
		}
		$seqCount{$pureSeq} += $duplicates;		
	}
	my $firstseq = my $firstfullseq = "";
	my %pushed = ();
	foreach my $seq (sort{$seqCount{$b} <=> $seqCount{$a}} keys %seqCount) {
		if (length $seq < $len) {
			push @sortUniqSeqs, $seq;
			$pushed{$seq} = 1;
			$flag = 1;
		}
		last;
	}
	
	if ($flag) {
		foreach my $fullseq (sort{$fullseqCount{$b} <=> $fullseqCount{$a}} keys %fullseqCount) {
			push @sortUniqSeqs, $fullseq;
			$pushed{$fullseq} = 1;
			last;
		}
	}
	
	foreach my $seq (sort{$seqCount{$b} <=> $seqCount{$a}} keys %seqCount) {
		unless ($pushed{$seq}) {
			push @sortUniqSeqs, $seq;
		}		
	}
	return (\@sortUniqSeqs, \%seqCount, $count);
}


sub GetUniqSeqs {
	my $nameWindowseqRef = shift;
	my $uniqFlag = shift;
	my (%seqCount, @sortUniqSeqs);
	my $count = 0;
	foreach my $name (keys %$nameWindowseqRef) {
		my $alignSeq = my $pureSeq = $nameWindowseqRef->{$name};
		my $duplicates = 1;
		if ($uniqFlag) {
			$duplicates = GetDuplicates ($name);
		}
		$count += $duplicates;
		$pureSeq =~ s/\-//g;			
		if (!$seqCount{$pureSeq}) {
			$seqCount{$pureSeq} = 0;
		}
		$seqCount{$pureSeq} += $duplicates;
	}	
	foreach my $seq (sort{$seqCount{$b} <=> $seqCount{$a}} keys %seqCount) {
		push @sortUniqSeqs, $seq;
	}
	return (\@sortUniqSeqs, \%seqCount, $count);
}


sub CalculateSimilarity {
	my ($consNasNoGapsRef, $readNasNoGapsRef, $consLen, $readLen, $match, $misMatch, $gapPenalty) = @_;
	my $value;
	for (my $i = 0; $i <= $consLen; $i++) {
		for (my $j = 0; $j <= $readLen; $j++) {
			if ($i == 0) {
				$value->[$i]->[$j] = $j * $gapPenalty;
			}elsif ($j == 0) {
				$value->[$i]->[$j] = $i * $gapPenalty;
			}else {
				my $upper = $value->[$i-1]->[$j] + $gapPenalty;
				my $left = $value->[$i]->[$j-1] + $gapPenalty;
				my $diagonal;
				if ($consNasNoGapsRef->[$i-1] eq $readNasNoGapsRef->[$j-1]) {
					$diagonal = $value->[$i-1]->[$j-1] + $match;
				}else {
					$diagonal = $value->[$i-1]->[$j-1] + $misMatch;
				}
				$value->[$i]->[$j] = Max ($upper, $left, $diagonal);
			}
			#print "$value->[$i]->[$j]\t";
		}
		#print "\n";
	}
	return $value;
}


sub CalculateProfileSimilarity {
	my ($profileFreq, $readNasNoGapsRef, $charsRef, $profAlignLen, $readLen, $match, $misMatch, $gapPenalty, $gapMatch, $flag) = @_;
	my $value;
	for (my $i = 0; $i <= $readLen; $i++) {
		for (my $j = 0; $j <= $profAlignLen; $j++) {
			if ($i == 0 && $j == 0) {
				$value->[$i]->[$j] = 0;
			}elsif ($i == 0) {				
				if ($flag eq 'full') {	# full window sequence
					my $current = 0;
					foreach my $char (@$charsRef) {						
						if ($char eq '-') {
							$current += $gapMatch * $profileFreq->{$j}->{$char};
						}else {
							$current += $gapPenalty * $profileFreq->{$j}->{$char};
						}			
					}
					$value->[$i]->[$j] = $value->[$i]->[$j-1] + $current;
				}else {	# pre or suffix
					$value->[$i]->[$j] = 0;
				}				
			}elsif ($j == 0) {
				if ($flag eq 'full') {
					$value->[$i]->[$j] = $value->[$i-1]->[$j] + $gapPenalty;
				}else {
					$value->[$i]->[$j] = 0;
				}								
			}else {
				my $gapCurrent = my $charCurrent = 0;
				foreach my $char (@$charsRef) {
					if ($char eq '-') {
						$charCurrent += $gapPenalty * $profileFreq->{$j}->{$char};
						$gapCurrent += $gapMatch * $profileFreq->{$j}->{$char};
					}else {
						$gapCurrent += $gapPenalty * $profileFreq->{$j}->{$char};						
						if ($readNasNoGapsRef->[$i-1] eq $char) {
							$charCurrent += $match * $profileFreq->{$j}->{$char};
#							print "match, char: $char, readNasNoGapsRef->[$i-1]: $readNasNoGapsRef->[$i-1], charCurrent: $charCurrent\n";
						}else {
							$charCurrent += $misMatch * $profileFreq->{$j}->{$char};
#							print "mismatch, char: $char, readNasNoGapsRef->[$i-1]: $readNasNoGapsRef->[$i-1], charCurrent: $charCurrent\n";
						}
					}
#					print "char: $char, charCurrent: $charCurrent, gapCurrent: $gapCurrent\n";	
				}
				my $upper = $value->[$i-1]->[$j] + $gapPenalty;
				my $left = $value->[$i]->[$j-1] + $gapCurrent;
				my $diagonal = $value->[$i-1]->[$j-1] + $charCurrent;
				$value->[$i]->[$j] = Max ($upper, $left, $diagonal);
#				print "upper: $upper, left: $left, diagonal: $diagonal, value->[$i]->[$j]: $value->[$i]->[$j]\n";
			}
#			print OUT "$value->[$i]->[$j]\t";
		}
#		print OUT "\n";
	}
	return $value;
}



sub Max {
	my @values = @_;
	my $flag = 0;
	my $max;
	foreach my $value (@values) {
		if (!$flag) {
			$max = $value;
			$flag++;
		}elsif ($value > $max) {
			$max = $value;
		}
	}
	return $max;
}


sub GetReadalignseq {
	my ($editTranscripts, $readNasNoGapsRef) = @_;
	my $readAlignSeq = '';
	foreach my $et (@$editTranscripts) {
		if ($et eq 'I') {
			$readAlignSeq .= '-';
		}else {
			$readAlignSeq .= shift @$readNasNoGapsRef;
		}
	}
	return $readAlignSeq;
}



sub Traceback {
	my ($consNasNoGaps, $readNasNoGaps, $value, $i, $j, $match, $misMatch, $gapPenalty) = @_;
	my @editTranscripts;	
	while ($i != 0 || $j != 0) {
		#print "i: $i, j: $j\n";
		my ($t, $editTranscript);
		if ($i != 0 && $j != 0) {
			if (@$consNasNoGaps[$i-1] eq @$readNasNoGaps[$j-1]) {
				$t = $match;
				$editTranscript = 'M';
			}else {
				$t = $misMatch;
				$editTranscript = 'R';
			}
			if ($value->[$i]->[$j] == $value->[$i-1]->[$j-1] + $t) {
				unshift @editTranscripts, $editTranscript;
				$i--;
				$j--;
			}elsif ($value->[$i]->[$j] == $value->[$i]->[$j-1] + $gapPenalty) {
				$editTranscript = 'I';
				unshift @editTranscripts, $editTranscript;
				$j--;
			}elsif ($value->[$i]->[$j] == $value->[$i-1]->[$j] + $gapPenalty) {
				$editTranscript = 'D';
				unshift @editTranscripts, $editTranscript;
				$i--;
			}
		}elsif ($i == 0) {	
			$editTranscript = 'I';
			unshift @editTranscripts, $editTranscript;
			$j--;
		}elsif ($j == 0) {
			$editTranscript = 'D';
			unshift @editTranscripts, $editTranscript;
			$i--;
		}
	}
	return \@editTranscripts;
}


sub CalculateFreq {
	my ($sub_profileFreq, $sub_lastAlignSeq, $sub_profPosNaCount, $sub_chars, $sub_seqCount, $sub_profAlignLen, $sub_prefixStatus, $sub_suffixStatus) = @_;
	my $start = 0;
	my $end = $sub_profAlignLen;
	my @sub_lastAlignNas = split //, $sub_lastAlignSeq;
	if ($sub_prefixStatus->{$sub_lastAlignSeq}) {
		$start = $sub_prefixStatus->{$sub_lastAlignSeq};
	}elsif ($sub_suffixStatus->{$sub_lastAlignSeq}) {
		$end = $end - $sub_suffixStatus->{$sub_lastAlignSeq};
	}
	for (my $i = $start; $i < $end; $i++) {	# will not count leading and ending '-' for profile calculations
		my $posNaCount = 0;
		if (defined $sub_prefixStatus->{$sub_lastAlignSeq}) {	# last aligned read is prefix read
			$posNaCount = $sub_seqCount->{prefix}->{$sub_lastAlignSeq};
		}elsif (defined $sub_suffixStatus->{$sub_lastAlignSeq}) {	# last aligned read is suffix read
			$posNaCount = $sub_seqCount->{suffix}->{$sub_lastAlignSeq};
		}else {	# last aligned read is full window read
			$posNaCount = $sub_seqCount->{full}->{$sub_lastAlignSeq};
		}
		my $totalPosNaCount = $sub_profPosNaCount->{$i+1} + $posNaCount;
		foreach my $char (@$sub_chars) {
			my $charCount = $sub_profileFreq->{$i+1}->{$char} * $sub_profPosNaCount->{$i+1};
			if ($sub_lastAlignNas[$i] eq $char) {
				$charCount += $posNaCount;
			}
			$sub_profileFreq->{$i+1}->{$char} = $charCount / $totalPosNaCount;
		}
		$sub_profPosNaCount->{$i+1} = $totalPosNaCount;
	}
}


sub ProfileTraceback {
	my ($profileFreq, $readNasNoGapsRef, $charsRef, $value, $readLen, $profAlignLen, $match, $misMatch, $gapPanulty, $gapMatch, $type) = @_;
	my @editTranscripts;
	my $i = $readLen;
	my $j = $profAlignLen;
	while ($i != 0 || $j != 0) {
		#print "i: $i, j: $j\n";
		my $editTranscript;
		my $gapCurrent = my $charCurrent = 0;
		if ($i != 0 && $j != 0) {
			foreach my $char (@$charsRef) {
				if ($char eq '-') {
					$charCurrent += $gapPanulty * $profileFreq->{$j}->{$char};
					$gapCurrent += $gapMatch * $profileFreq->{$j}->{$char};
				}else {
					$gapCurrent += $gapPanulty * $profileFreq->{$j}->{$char};
					if ($readNasNoGapsRef->[$i-1] eq $char) {
						$charCurrent += $match * $profileFreq->{$j}->{$char};
					}else {
						$charCurrent += $misMatch * $profileFreq->{$j}->{$char};
					}
				}
			}
			if ($value->[$i]->[$j] == $value->[$i-1]->[$j-1] + $charCurrent) {
				$editTranscript = 'M';
				unshift @editTranscripts, $editTranscript;
				$i--;
				$j--;
			}elsif ($value->[$i]->[$j] == $value->[$i]->[$j-1] + $gapCurrent) {
				$editTranscript = 'I';
				unshift @editTranscripts, $editTranscript;
				$j--;
			}elsif ($value->[$i]->[$j] == $value->[$i-1]->[$j] + $gapPanulty) {
				$editTranscript = 'D';
				unshift @editTranscripts, $editTranscript;
				$i--;
			}
		}elsif ($i == 0) {	
			$editTranscript = 'I';
			unshift @editTranscripts, $editTranscript;
			$j--;
		}elsif ($j == 0) {
			$editTranscript = 'D';
			unshift @editTranscripts, $editTranscript;
			$i--;
		}
#		print "$editTranscript, value->[$i]->[$j]: $value->[$i]->[$j]\n";
	}
	return \@editTranscripts;
}



sub RealignProfile {
	my ($profileSeqsRef, $len, $alignSeqs, $editTranscripts, $seqCount, $seqAlignseq, $sub_profilePosNaCount, $prefixStatusRef, $suffixStatusRef) = @_;
				
	for (my $i = 0; $i < @$profileSeqsRef; $i++) {
		my $profileSeq = $profileSeqsRef->[$i];
		my @profileSeqNas = split //, $profileSeq;
		my $alignSeq = '';
		foreach my $et (@$editTranscripts) {
			if ($et eq 'D') {
				$alignSeq .= '-';
			}else {
				$alignSeq .= shift @profileSeqNas;
			}
		}
		$alignSeqs->[$i] = $alignSeq;
		my $pureSeq = $profileSeq;
		$pureSeq =~ s/\-//g;
		
		if ($seqCount->{prefix}->{$profileSeq}) {	# aligned prefix sequence					
			$seqCount->{prefix}->{$alignSeq} = $seqCount->{prefix}->{$profileSeq};
			for (my $i = 0; $i < $len; $i++) {	# re-count biginning gaps for newly aligned prefix sequence
				if (substr ($alignSeq, $i, 1) ne '-') {
					$prefixStatusRef->{$alignSeq} = $i;
					last;
				}
			}
			
			$seqAlignseq->{prefix}->{$pureSeq} = $alignSeq;
		}
		
		if ($seqCount->{suffix}->{$profileSeq}) {	# aligned suffix sequence, which may be same as the prefix above				
			$seqCount->{suffix}->{$alignSeq} = $seqCount->{suffix}->{$profileSeq};
			for (my $i = $len - 1; $i >= 0; $i--) {	# re-count ending gaps for newly aligned prefix sequence
				if (substr ($alignSeq, $i, 1) ne '-') {
					$suffixStatusRef->{$alignSeq} = $len - $i - 1;
					last;
				}
			}
			$seqAlignseq->{suffix}->{$pureSeq} = $alignSeq;
		}
		
		if ($seqCount->{full}->{$profileSeq}) {	# aligned full window sequence, which may be same as the prefix or suffix above
			$seqCount->{full}->{$alignSeq} = $seqCount->{full}->{$profileSeq};
			$seqAlignseq->{full}->{$pureSeq} = $alignSeq;
		}	
	}
	# calculate na count in each position for re-aligned profile	
	for (my $i = 0; $i < $len; $i++) {
		my $realignPosNaCount = 0;
		foreach my $alignSeq (@$alignSeqs) {
			my @alignSeqNas = split //, $alignSeq;
			if ($seqCount->{prefix}->{$alignSeq}) {
				unless ($i < $prefixStatusRef->{$alignSeq}) {	# not count leading gaps
					$realignPosNaCount += $seqCount->{prefix}->{$alignSeq};
				}
			}elsif ($seqCount->{suffix}->{$alignSeq}) {
				unless ($i > $len - $suffixStatusRef->{$alignSeq} - 1) {	# not count ending gaps
					$realignPosNaCount += $seqCount->{suffix}->{$alignSeq};
				}
			}else {
				$realignPosNaCount += $seqCount->{full}->{$alignSeq};
			}
		}
		$sub_profilePosNaCount->{$i+1} = $realignPosNaCount;
	}	
}

sub ReassignProfFreq {
	my ($sub_profileFreq, $sub_editTranscripts, $sub_chars) = @_;
	my $pos = 0;
	my $realignProfFreq;
	for (my $i = 0; $i < @$sub_editTranscripts; $i++) {
		$pos += 1;
		if ($sub_editTranscripts->[$i] eq 'D') {
			foreach my $char (@$sub_chars) {
				if ($char eq '-') {
					$realignProfFreq->{$i+1}->{$char} = 1;
				}else {
					$realignProfFreq->{$i+1}->{$char} = 0;
				}								
			}
			$pos--;
		}else {
			foreach my $char (@$sub_chars) {
				$realignProfFreq->{$i+1}->{$char} = $sub_profileFreq->{$pos}->{$char};							
			}
		}
	}
	return $realignProfFreq;
}

sub RevProfFreq {
	my ($sub_profileFreq, $sub_profAlignLen, $sub_chars) = @_;
	my $revProfFreq;	
	for (my $i = 1; $i <= $sub_profAlignLen; $i++) {
			foreach my $char (@$sub_chars) {
			$revProfFreq->{$i}->{$char} = $sub_profileFreq->{$sub_profAlignLen - $i + 1}->{$char};
		}
	}
	return $revProfFreq;
}