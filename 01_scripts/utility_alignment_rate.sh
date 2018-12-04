#!/bin/bash
# Generate alignment rate table, execute from main repo directory 

grep -E 'overall alignment rate' 07_mapped/report_alignment_result.txt | awk -F"%" '{ print $1 }' - > 07_mapped/overall_alignment_rate.txt

grep -E '\) aligned concordantly exactly 1 time' 07_mapped/report_alignment_result.txt | awk '{ print $2 }' - | awk -F"(" '{ print $2 }' - | awk -F"%" '{ print $1 }' - > 07_mapped/align_concord_exact_once.txt

grep -E '\) aligned concordantly >1 times' 07_mapped/report_alignment_result.txt | awk '{ print $2 }' - | awk -F"(" '{ print $2 }' - | awk -F"%" '{ print $1 }' - > 07_mapped/align_concord_more_than_once.txt
