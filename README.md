# AMP_alignment_process
A workflow and toolkit of perl scripts to manipulate multiple sequence alignments

### Usage

#### 1. In a directory that contains nucleotide sequence fasta files, run 
```
perl scriptPath/run_AMP_alignment_process.pl –id directoryPath(default: .)
```
  - Renames file names beginning with “Hxxx” to “Vxxx” and renames sequence names within each the file with the standard nomenclature
  - Merges sequence files for the same sample ID, based on file names 
    - Creates subdirectory “GP/REN” to store merged sequence files and files produced by the following steps
  - Collapses identical sequences into unique sequences, adds the number of sequences that were collapsed to the sequence name
  - Muscle aligns collapsed sequences
#### 2. Manually inspect/edit refined alignments
#### 3. In the directory that contains both the inspected alignments and a reference alignment, run 
```
perl scriptPath/run_ref_sample_profile_alignment.pl –id directoryPath(default: .) 
–ref referenceFileName
```
  - Merges reference alignment to participant sequence alignment by profile alignment
#### 4. Manually inspect/edit reference containing alignments
#### 5. In the directory that contains the inspected reference-included alignments, run 
```
perl scriptPath/run_extraction_translation.pl –id directoryPath(default: .) 
–sp startPositionOfGagOrRevInHXB2
```
  - Extracts gene regions (Gag, Pol, Prot for GP; Rev, Vpu, Env, Nef for REN) in alignment based on the start positions in HXB2 within alignment
    - Creates subdirectories of gene regions to store extracted alignments and files produced by following steps
  - Translates extracted nucleotide sequence alignment a corresponding amino acid sequence alignment
#### 6. Manually inspect/edit translated amino acid alignment and adjust corresponding nucleotide alignment as needed
#### 7. In the directory that contains alignments needed to be uncollapsed and the corresponding .name files, run 
```
perl scriptPath/uncollapse_seqs.pl inputFile nameFile
```
  - Uncollapses alignment into full alignment of individual sequences with sequence names
