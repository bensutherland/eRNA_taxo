#!/bin/bash
# Index the SortMeRNA databases in the main SortMeRNA directory 

# Global variables
SORTMERNA_DB="/home/ben/Programs/sortmerna-2.1"

# User set variable
NUM_CPU="10"

indexdb_rna --ref \
$SORTMERNA_DB/rRNA_databases/silva-bac-16s-id90.fasta,$SORTMERNA_DB/index/silva-bac-16s-db:\
$SORTMERNA_DB/rRNA_databases/silva-bac-23s-id98.fasta,$SORTMERNA_DB/index/silva-bac-23s-db:\
$SORTMERNA_DB/rRNA_databases/silva-arc-16s-id95.fasta,$SORTMERNA_DB/index/silva-arc-16s-db:\
$SORTMERNA_DB/rRNA_databases/silva-arc-23s-id98.fasta,$SORTMERNA_DB/index/silva-arc-23s-db:\
$SORTMERNA_DB/rRNA_databases/silva-euk-18s-id95.fasta,$SORTMERNA_DB/index/silva-euk-18s-db:\
$SORTMERNA_DB/rRNA_databases/silva-euk-28s-id98.fasta,$SORTMERNA_DB/index/silva-euk-28s:\
$SORTMERNA_DB/rRNA_databases/rfam-5s-database-id98.fasta,$SORTMERNA_DB/index/rfam-5s-db:\
$SORTMERNA_DB/rRNA_databases/rfam-5.8s-database-id98.fasta,$SORTMERNA_DB/index/rfam-5.8s-db

