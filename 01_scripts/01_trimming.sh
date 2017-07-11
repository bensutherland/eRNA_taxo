#!/bin/bash
# Remove adapters and light quality trimming with Trimmomatic

# Global variables
RAW_FOLDER="02_raw_data"
TRIMMED_FOLDER="03_trimmed"
VECTORS="./00_archive/truseq_universal.fasta"
TRIMMOMATIC_PROGRAM="/home/ben/Programs/Trimmomatic-0.36/trimmomatic-0.36.jar"

# User set variable
NUM_CPU="10"

# Filtering and trimming data with trimmomatic
ls -1 $RAW_FOLDER/*.fastq.gz | 
    perl -pe 's/R[12]\_001\.fastq\.gz//' |
    sort -u |
    while read i
    do
        echo "Trimming $i"
        java -Xmx100G -jar $TRIMMOMATIC_PROGRAM PE \
            -threads $NUM_CPU \
            -phred33 \
            "$i"R1_001.fastq.gz \
            "$i"R2_001.fastq.gz \
            "$i"R1_paired_trimmed.fq.gz \
            "$i"R1_unpaired_trimmed.fq.gz \
            "$i"R2_paired_trimmed.fq.gz \
            "$i"R2_unpaired_trimmed.fq.gz \
            ILLUMINACLIP:$VECTORS:2:30:10 \
            SLIDINGWINDOW:20:2 \
            LEADING:2 \
            TRAILING:2 \
            MINLEN:80
    done

# Move trimmed files to trimmed folder
mv $RAW_FOLDER/*_trimmed.fq.gz $TRIMMED_FOLDER

