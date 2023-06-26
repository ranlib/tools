#!/bin/bash

BAM=~/Analysis/varpipe3/data/Output_with_tbdb_reference/ERR552797/ERR552797_sdrcsm.bam
REF=~/Analysis/varpipe4/references/NC_000962.3/NC_000962.3.fasta 

~/Analysis/tools/depth_of_coverage.py --bam $BAM \
				      --reference $REF \
				      --intervals tbdb.bed \
				      --lower_coverage 10 \
				      --output ERR552797 \
				      --csv ERR552797_depth_of_coverage.csv
