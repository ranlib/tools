#!/bin/bash

BAM=~/Analysis/varpipe4/tests//ERR552797/ERR552797_sdrcsm.bam
REF=~/Analysis/varpipe4/references/NC_000962.3/NC_000962.3.fasta 
BED=/home/dieterbest/Analysis/varpipe4/intervals/tbdb_intervals.bed

# ~/Analysis/tools/depth_of_coverage.py --bam $BAM \
# 				      --reference $REF \
# 				      --intervals $BED \
# 				      --lower_coverage 10 \
# 				      --output ERR552797 \
# 				      --csv ERR552797_depth_of_coverage.csv

# cross check with samtools
samtools depth -a -b $BED $BAM > samtools_depth.tsv

