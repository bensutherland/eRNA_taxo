#!/bin/bash
# Remove adapters and light quality trimming with cutadapt 

# Global variables
RAW_FOLDER="02_raw_data"
TRIMMED_FOLDER="03_trimmed"
R1_ADAPTER="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
R2_ADAPTER="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"

# Where here R1_ADAPTER is the TruSeq Indexed Adapter
# and R2_ADAPTER is the revcomp of TruSeq Universal Adapter
# as suggested by cutadapt README

# Filtering and trimming data with trimmomatic
ls -1 $RAW_FOLDER/*.fastq.gz | 
    perl -pe 's/R[12]\_001\.fastq\.gz//' |
    sort -u |
    while read i
    do
        echo "Trimming $i"
        cutadapt -a $R1_ADAPTER -A $R2_ADAPTER \
        --interleaved \
        -o "$i"ilvd_trimmed.fq.gz
        #-o "$i"R1_trimmed.fq.gz -p "$i"R2_trimmed.fq.gz \
        "$i"R1_001.fastq.gz "$i"R2_001.fastq.gz \
        -q 20
    done


# Move trimmed files to trimmed folder
mv $RAW_FOLDER/*_trimmed.fq.gz $TRIMMED_FOLDER

