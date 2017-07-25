# ACRDP Shellfish Water Transcriptome Project #
This pipeline was developed as part of the Molecular Genetics Lab at Pacific Biological Station (Nanaimo, BC) in the working group of Kristi Miller. The pipeline is developed for the purpose of analyzing data in the lab and comes with no warrantee or guarantees.   

## Requirements:
`fastqc` (currently: v0.11.5) https://www.bioinformatics.babraham.ac.uk/projects/fastqc/     
`multiqc`   http://multiqc.info/     
`cutadapt`  http://cutadapt.readthedocs.io/en/stable/index.html        
`SortMeRNA` http://bioinfo.lifl.fr/RNA/sortmerna/    

### 1. Quality check and trimming
#### Run FastQC on raw data   
Make new folder for results output    
`mkdir 02_raw_data/fastqc_results`    
Run fastqc and output into the subfolder.   
`fastqc 02_raw_data/S*.fastq.gz -o 02_raw_data/fastqc_results/`    

#### Summarize results with MultiQC
`multiqc -o ./02_raw_data/fastqc_results/ ./02_raw_data/fastqc_results/`   

#### Trim for quality and adapters, output interleaved files 
See details within the following script    
`./01_cutadapt.sh`

#### Run FastQC on trimmed data 
`mkdir 03_trimmed/fastqc_trimmed_results`
`fastqc 03_trimmed/*_ilvd_trimmed.fq.gz -o 03_trimmed/fastqc_trimmed_results/`

#### Summarize results with MultiQC
`multiqc -o 03_trimmed/fastqc_trimmed_results/ 03_trimmed/fastqc_trimmed_results/`     


### 2. Separate rRNA and metatranscriptomic data with SortMeRNA
#### Index the SortMeRNA fasta database
Note: Confirm the path to SortMeRNA installation in the following script:        
`./01_scripts/02a_index_SortMeRNA_db.sh`

#### Decompress interleaved fastq files
SortMeRNA requires decompressed, interleaved input.    
`for i in 03_trimmed/*.fq.gz ; do gunzip $i ; done`

#### Separate rRNA from metatranscriptome data
`./01_scripts/02b_SortMeRNA.sh`
This will put your data into `04_sorted`   
To see pertinent log information, run    
`tail -n 20 04_sorted/*.log | less`

### 3. Metatranscriptomic analysis
First assemble the reference transcriptome from the reads using `IDBA-Tran`    

Prepare IDBA-Tran by editing for longer reads (MiSeq)
"please modify the constant kMaxShortSequence in src/sequence/short_sequence.h to support longer read length"      
`static const uint32_t kMaxShortSequence = 128` to `static const uint32_t kMaxShortSequence = 350`    
More info: http://seqanswers.com/forums/showthread.php?t=29109

PE data needs to be in interleaved fasta format. So first format the output from SortMeRNA to fasta format from fastq.     
`for i in 04_sorted/*ilvd_trimmed.fq_non_rRNA.fq ; do fq2fa --paired $i ${i%.fq}.fa ; done`    

Combine all libraries into a single fasta file    
`cat 04_sorted/*ilvd_trimmed.fq_non_rRNA.fa > 04_sorted/all_samples_ilvd_trimmed_non_rRNA.fa`

Start the assembly with the interleaved .fa files.    
`idba_tran --read 04_sorted/all_samples_ilvd_trimmed_non_rRNA.fa --num_threads=8 --out 05_assembled`


### Trinity assembly
Next, try assembling with trinity
De-interleave the concatenated fasta file.   
`fastaq_deinterleave.1 `

First need to deinterleave the concatenated fasta    
`deinterleave_fasta.sh < 04_sorted/all_samples_ilvd_trimmed_non_rRNA.fa 04_sorted/all_samples_trimmed_non_rRNA_R1.fa 04_sorted/all_samples_trimmed_non_rRNA_R2.fa`

Then launch trinity
`Trinity --seqType fq --left ./04_sorted/all_samples_trimmed_non_rRNA_R1.fa --right ./04_sorted/all_samples_trimmed_non_rRNA_R2.fa --CPU 6 --max_memory 100G`
  


