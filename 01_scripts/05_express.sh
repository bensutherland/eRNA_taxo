#!/bin/bash
# Obtain gene counts per individual using express

# Set environment variables
MAPPED_FOLDER="07_mapped"
COUNT_FOLDER="08_gx_levels"

# Choose reference
#REFERENCE="06_metatranscriptome/assemblies_merged_bbmap_reduced.fa"
#REFERENCE="06_metatranscriptome/assemblies_merged_red0.95_bbmap_red.fa"
REFERENCE="06_metatranscriptome/16_libs_contig.fa"

# Produce counts per individual with express 
ls -1 $MAPPED_FOLDER/*.sorted.bam |
    sort -u |
    while read i
    do
        echo "Counts for sample" $i
        name=$(basename $i)
        express $REFERENCE $i --max-read-len=300 --no-update-check
        mv results.xprs $COUNT_FOLDER/"$name"_results.xprs
        mv params.xprs $COUNT_FOLDER/"$name"_params.xprs 
    done

