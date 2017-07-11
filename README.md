## ACRDP Shellfish Water Transcriptome Project ## 

### Requirements:
`java`    
`fastqc`    
`multiqc`    
`trimmomatic`    


### 1. Quality check and trimming
### Use fastqc to check quality
### put all files into one folder, and create a folder for the results
`mkdir 02_raw_data/fastqc_results`    
`fastqc 02_raw_data/S*.fastq.gz -o 02_raw_data/fastqc_results/`    

### Use multiqc to summarize fastqc results, and output to fastqc results folder
`multiqc -o ./02_raw_data/fastqc_results/ ./02_raw_data/fastqc_results/`

### Trim
`./01_trimming.sh`

### Use fastqc again to check trimmed quality
`mkdir 03_trimmed/fastqc_trimmed_results`
`fastqc 03_trimmed/*_paired_trimmed.fq.gz -o 03_trimmed/fastqc_trimmed_results/`


### Summarize with multiqc
 

