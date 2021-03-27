# AMP_alignment_process
A workflow and toolkit of perl scripts to manipulate multiple sequence alignments

### Usage

#### 1.  In a directory that contains nucleotide sequence fasta files, run 
```
perl scriptPath/run_AMP_alignment_process.pl –id directoryPath(default: .)
```
  - Renames file names beginning with “Hxxx” to “Vxxx” and renames sequence names within each the file with the standard nomenclature
  - Merges sequence files for the same sample ID, based on file names 
    - Creates subdirectory “GP/REN” to store merged sequence files and files produced by the following steps
  - Collapses identical sequences into unique sequences, adds the number of sequences that were collapsed to the sequence name
  - Reverses collapsed sequences
  - Muscle aligns reversed sequences
  - Reverses back aligned sequences
#### 2.  Manually inspect/edit refined alignments
#### 3.  In the directory that contains both the inspected alignments and HXB2 GP/REN sequence fasta file, run 
```
perl scriptPath/run_ref_sample_profile_alignment.pl –id directoryPath(default: .) 
–ref HXB2FileName
```
  - Merges HXB2 sequence to participant sequence alignment by profile alignment
#### 4.  In the directory that contains the inspected HXB2-included alignments, run 
```
perl scriptPath/run_gagenv_extraction_translation.pl –id directoryPath(default: .) 
–sp startPositionOfGagOrRevInHXB2(default: 1)
```
  - Extracts gene regions (Gag for GP; Env for REN) in alignment based on the start positions in HXB2 within alignment
    - Creates subdirectories of gene regions to store extracted alignments and files produced by following steps
  - Translates extracted nucleotide sequence alignment to a corresponding amino acid sequence alignment
  - Retrieves functional protein sequences by filtering out defective protein sequences (Sequence doesn’t start with “ATG”, premature stop codon locates before 95% of median protein sequence length, deletion is > 80% of protein sequence median length)
  - Collapse functional protein sequence alignment
#### 5.  Manually inspect/edit translated functional amino acid alignment
#### 6.  In the directory that contains reviewed functional protein sequence alignments, run
```
perl scriptPath/calcAAconsensus_0majority_ignore_gaps.pl inFileName
```
  - Calculate amino acid consensus sequence (0% majority and ignore gaps) for each subject's functional amino acid alignment
    - Creates two subdirectories. One contains a file of all subjects' consensus sequences. The other contains a file of all subjects' consensus sequences and individual sequence.
#### 7.  Align all subjects' functional amino acid consensus sequences via Muscle
#### 8.  Manually inspect/edit all subjects' amino acid consensus sequence alignment
#### 9.  Align reference sequence alignment and consensus sequence alignment via Muscle’s profile-profile alignment to make a master alignment
#### 10. Manually inspect/edit the master alignment
