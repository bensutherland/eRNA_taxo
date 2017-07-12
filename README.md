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

#### 
`./01_scripts/02b_SortMeRNA.sh`

