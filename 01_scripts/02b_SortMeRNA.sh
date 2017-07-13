#!/bin/bash
# Index the SortMeRNA databases in the main SortMeRNA directory 

# Global variables
TRIMMED_FOLDER="03_trimmed"
SORTMERNA_DB="/home/ben/Programs/sortmerna-2.1"

# User set variable
NUM_CPU="12"


# Running SortMeRNA
ls -1 $TRIMMED_FOLDER/*_ilvd_trimmed.fq |
    sort -u |
    while read i
    do
        echo "Treating $i"
        sortmerna --ref $SORTMERNA_DB/rRNA_databases/silva-bac-16s-id90.fasta,$SORTMERNA_DB/index/silva-bac-16s-db:$SORTMERNA_DB/rRNA_databases/silva-bac-23s-id98.fasta,$SORTMERNA_DB/index/silva-bac-23s-db:$SORTMERNA_DB/rRNA_databases/silva-arc-16s-id95.fasta,$SORTMERNA_DB/index/silva-arc-16s-db:$SORTMERNA_DB/rRNA_databases/silva-arc-23s-id98.fasta,$SORTMERNA_DB/index/silva-arc-23s-db:$SORTMERNA_DB/rRNA_databases/silva-euk-18s-id95.fasta,$SORTMERNA_DB/index/silva-euk-18s-db:$SORTMERNA_DB/rRNA_databases/silva-euk-28s-id98.fasta,$SORTMERNA_DB/index/silva-euk-28s:$SORTMERNA_DB/rRNA_databases/rfam-5s-database-id98.fasta,$SORTMERNA_DB/index/rfam-5s-db:$SORTMERNA_DB/rRNA_databases/rfam-5.8s-database-id98.fasta,$SORTMERNA_DB/index/rfam-5.8s-db \
        --reads $i --sam --num_alignments 1 --fastx --aligned "$i"_rRNA --other "$i"_non_rRNA --log -v -a $NUM_CPU -m 4096 --paired_in
    done
