## ACRDP Shellfish Water Transcriptome Project ## 

### Requirements:
`java`    
`fastqc`    
`multiqc`    
`trimmomatic`    


### 1. Quality check and trimming
### Use fastqc to check quality
### put all files into one folder, and create a folder for the results
`mkdir fastqc_results`    
`fastqc S*.fastq.gz -o fastqc_results/`    

### Use multiqc to summarize fastqc results, and output to fastqc results folder
`multiqc -o ./02_raw_data/fastqc_results/ ./02_raw_data/fastqc_results/`

### Trim
`./01_trimming.sh`
