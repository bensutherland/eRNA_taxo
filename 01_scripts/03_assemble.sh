#!/bin/bash
# Assemble individual RNA-seq libraries in each individual assembly 

# Global variables
SORTED_FOLDER="04_sorted"
ASSEMBLY_FOLDER="05_assembled"

# User set variable
NUM_CPU="8"


# Running SortMeRNA
ls -1 $SORTED_FOLDER/*_non_rRNA.fa |
    sort -u |
    while read i
    do
        echo "Treating $i"
        FOLDER_NAME="${i%_L001_ilvd_trimmed.fq_non_rRNA.fa}"
        OUTPUT="./$ASSEMBLY_FOLDER/$( basename $FOLDER_NAME )"
        mkdir $OUTPUT
        idba_tran --read $i --num_threads=$NUM_CPU --out $OUTPUT

        # Clean up
        rm $OUTPUT/kmer
    done

