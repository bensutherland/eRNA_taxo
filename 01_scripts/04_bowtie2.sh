#!/bin/bash
# Map interleaved fasta to the reference transcriptome with bowtie2 
# Note: Requires that reference is indexed (see README.md) 

# Global variables
TRIMMED_FOLDER="04_sorted"
MAPPED_FOLDER="07_mapped"

# Choose reference
#REFERENCE="06_metatranscriptome/assemblies_merged_bbmap_reduced"
#REFERENCE="06_metatranscriptome/two_libs_contig"
REFERENCE="06_metatranscriptome/assemblies_merged_red0.95_bbmap_red.fa"

# User variables
NUM_THREADS="2"

# Map reads and add RG
ls -1 $TRIMMED_FOLDER/*non_rRNA.fq |
	sort -u |
	while read i
	do
	  echo $i
	  name=$(basename $i)
	  label=$(echo $name | perl -pe 's/\_L001\_ilvd\_trimmed\.fq\_non\_rRNA\.fq//')
	  ID="@RG\tID:${label}\tSM:${label}\tPL:Illumina"
	  bowtie2 --end-to-end -k 40 --threads $NUM_THREADS --rg-id $ID -x $REFERENCE --interleaved $i -S $i.bowtie2.sam
          ##no longer needed?
	  # samtools view -Sb $i.bowtie2.sam > $i.bowtie2.unsorted.bam
          
          ## With sufficient filtering flags (-q 2 -F 4) #mapping qual > 2
          samtools view -Sb -F 4 -q 2 $i.bowtie2.sam > $i.bowtie2.unsorted.bam
          
          # Sort 
	  samtools sort -n $i.bowtie2.unsorted.bam -o $i.sorted.bam #samtools >= 1.5 

          ## In case using samtools <1.5
	  # samtools sort -n $i.bowtie2.unsorted.bam $i.bowtie2.sorted #samtools < 1.5

          rm $i.bowtie2.sam
done

# clean up space
mv ./$TRIMMED_FOLDER/*.bam ./$MAPPED_FOLDER/
