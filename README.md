# ACRDP Shellfish Water Transcriptome Project #
This pipeline was developed as part of the Molecular Genetics Lab at Pacific Biological Station (Nanaimo, BC) in the working group of Kristi Miller. The pipeline is developed for the purpose of analyzing data in the lab and comes with no warrantee or guarantees.   

## Requirements:
`fastqc` (currently: v0.11.5) https://www.bioinformatics.babraham.ac.uk/projects/fastqc/     
`multiqc`   http://multiqc.info/     
`cutadapt`  http://cutadapt.readthedocs.io/en/stable/index.html        
`SortMeRNA` http://bioinfo.lifl.fr/RNA/sortmerna/    
`idba-tran` https://github.com/loneknightpy/idba     
`bbmap`     https://sourceforge.net/projects/bbmap/    
`cd-hit-est`http://weizhongli-lab.org/cd-hit/     
`bowtie2`   http://bowtie-bio.sourceforge.net/bowtie2/index.shtml     
`express`   https://pachterlab.github.io/eXpress/index.html      


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

### 3. Metatranscriptome assembly
A) assemble the reference transcriptome from the reads using `IDBA-Tran`    

Prepare IDBA-Tran by editing for longer reads (MiSeq)
"please modify the constant kMaxShortSequence in src/sequence/short_sequence.h to support longer read length"      
`static const uint32_t kMaxShortSequence = 128` to `static const uint32_t kMaxShortSequence = 350`    
More info: http://seqanswers.com/forums/showthread.php?t=29109

IDBA-Tran requires interleaved fasta. Convert interleaved fastq to fasta    
`for i in 04_sorted/*ilvd_trimmed.fq_non_rRNA.fq ; do fq2fa --paired $i ${i%.fq}.fa ; done`    

Assemble each library individually to reduce computational load.   
`./01_scripts/03_assemble.sh`

Copy the assemblies with names over to next folder:   
`for i in 05_assembled/*/contig.fa ; do RENAME_ID=$(echo $i | awk -F/ '{ print $2"_contig.fa" }') ; echo "copying $RENAME_ID" ; cp $i 06_metatranscriptome/$RENAME_ID ; done`

Merge into one metatranscriptome file (with redundancy):    
`cat 06_metatranscriptome/*.fa > 06_metatranscriptome/assemblies_merged.fa`

Remove spaces in the fasta file name     
`sed 's/\ /\_/g' 06_metatranscriptome/assemblies_merged.fa > 06_metatranscriptome/assemblies_merged_for_dedupe.fa`

B) Merge assemblies    
Dedupe.sh from bbmap:   
`dedupe.sh in=06_metatranscriptome/assemblies_merged_for_dedupe.fa out=assemblies_merged_bbmap_reduced.fa threads=4 uniquenames=t`

*OR*   

cd-hit-est from cd-hit:    
`cd-hit-est -i ./assemblies_merged.fa -o assemblies_merged_cd-hit-est90_reduced.fa -c 0.9 -T 4 -M 9000`

### 4. Quantification
First index the transcriptome:   
`bowtie2-build -f 06_metatranscriptome/assemblies_merged_bbmap_reduced.fa 06_metatranscriptome/assemblies_merged_bbmap_reduced`    

Launch bowtie2 mapping script:   
`./01_scripts/04_bowtie2.sh`    

This will produce .bam files. 

Quantify bam files:
`./01_scripts/05_express.sh`    
This produces several xprs files that are used in the next stage.  

### 5. Normalization and Differential Expression
Work within R to import the named results.xprs files and collect all into a single matrix file.   
`01_scripts/prepare_gxlevels_matrix.R`   

Work within R using the script in the directory     
`01_scripts/06_edgeR_normalization.R`    


### 6. Annotate transcriptome
For this I suggest using Eric Normandeau's `go_enrichment` pipeline    
https://github.com/enormandeau/go_enrichment    
This pipeline annotates with uniprot's swissprot database and performs enrichment tests with goatools https://github.com/tanghaibao/goatools.    

# Optional: Annotate only expressed transcripts
Unwrap fasta file
Use fasta_unwrap.py from Eric Normandeau's scripts folder to unwrap the fasta file to obtain only the expressed contigs.

`for i in `cat expr_contigs.csv` ; do grep -A1 $i 06_met
atranscriptome/assemblies_merged_red0.95_bbmap_red_unwrapped.fa ; done > 06_metatranscriptome/selected_contigs_2017-09-09.fa`

