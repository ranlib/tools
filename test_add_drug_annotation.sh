#!/bin/bash

BAM=~/Analysis/varpipe3/data/Output_with_tbdb_reference/ERR552797/ERR552797_sdrcsm.bam
REF=~/Analysis/varpipe4/references/NC_000962.3/NC_000962.3.fasta 

# rm -r ERR552797
# ~/Analysis/tools/depth_of_coverage.py --bam $BAM \
# 				      --reference $REF \
# 				      --intervals tbdb.bed \
# 				      --lower_coverage 10 \
# 				      --output ERR552797 \
# 				      --csv ERR552797_depth_of_coverage.csv

# ~/Analysis/tools/merge_bcftools_SnpSift.py --vcf CDC/ERR552797_DR_loci_raw_annotation.vcf --samplename ERR552797 --output ERR552797_cdc_snpsift.tsv

~/Analysis/tools/add_drug_annotation.py --tsv ERR552797_cdc_snpsift.tsv --minimum_coverage 10 --all_genes\
					--json tbdb.dr.json --bed tbdb.bed --coverage ERR552797_depth_of_coverage.csv \
					--output ERR552797_cdc_drugs.tsv --report ERR552797_cdc_interpretation.tsv

