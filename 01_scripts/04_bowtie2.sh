#!/bin/bash
# Map interleaved fasta to the reference transcriptome with bowtie2 
# Note: Requires that reference is indexed (see README.md) 

# Global variables
TRIMMED_FOLDER="04_sorted"
MAPPED_FOLDER="07_mapped"
REFERENCE="06_metatranscriptome/assemblies_merged_bbmap_reduced"

# User variables
NUM_THREADS="1"

# Map reads and add RG
ls -1 $TRIMMED_FOLDER/*non_rRNA.fq |
	sort -u |
	while read i
	do
	  echo $i
	  name=$(basename $i)
	  label=$(echo $name | perl -pe 's/\_L001\_ilvd\_trimmed\.fq\_non\_rRNA\.fq//')
	  ID="@RG\tID:${label}\tSM:${label}\tPL:Illumina"
	  #bowtie2 --end-to-end -k 40 --threads $NUM_THREADS --rg-id $ID -x $REFERENCE --interleaved $i -S $i.bowtie2.sam
	  samtools view -Sb $i.bowtie2.sam > $i.bowtie2.unsorted.bam
	  samtools sort -n $i.bowtie2.unsorted.bam $i.bowtie2.sorted
      #rm $i.bowtie2.sam
done

# clean up space
mv ./$TRIMMED_FOLDER/*.bam ./$MAPPED_FOLDER/
