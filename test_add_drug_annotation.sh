#!/bin/bash

#
# test TB mutations interpretation procedure
#

SAMPLE=ERR552797
BAM=./CDC//ERR552797_sdrcsm.bam
VCF=./CDC/ERR552797_DR_loci_raw_annotation.vcf
REF=~/Analysis/varpipe4/references/NC_000962.3/NC_000962.3.fasta 
BED=$HOME/Analysis/varpipe4/intervals/tbdb.bed
JSON=$HOME/Analysis/varpipe4/intervals/tbdb.dr.json

#
# determine depth of coverage in regions of interest
# rm -r ERR552797
# ~/Analysis/tools/depth_of_coverage.py --bam $BAM \
#  				      --reference $REF \
#  				      --intervals tbdb.bed \
#  				      --lower_coverage 10 \
#  				      --output ERR552797 \
#  				      --csv ERR552797_depth_of_coverage.tsv

# # vcf -> tsv
# ~/Analysis/tools/merge_bcftools_SnpSift.py --vcf $VCF --samplename ERR552797 --output ERR552797_cdc_snpsift.tsv

# this GATK command is basically like bcftools
# ie doesn't get the INFO field like SnpSift
#gatk VariantsToTable -V CDC/ERR552797_DR_loci_raw_annotation.vcf \
#     -F CHROM -F POS -F REF -F ALT -F FILTER -F GENE_NAME -F GENE_ID -F EFFECT -F IMPACT -F HGVS_C -F HGVS_P -F CDC_POS -GF AF -GF AD -GF DP \
#     -O ERR552797_cdc_snpsift.tsv

# run interpretation
# ~/Analysis/tools/add_drug_annotation.py --tsv ERR552797_cdc_snpsift.tsv \
# 					--coverage ERR552797_depth_of_coverage.csv \
# 					--minimum_coverage 10 \
# 					--filter_genes \
# 					--json tbdb.dr.json \
# 					--bed tbdb.bed \
# 					--output ERR552797_cdc_drugs.tsv \
# 					--report ERR552797_cdc_interpretation.tsv

# ~/Analysis/tools/vcf_to_pandas_dataframe.py \
#     --vcf all_variants.vcf \
#     --samplename ERR552797 \
#     --tsv all_variants.tsv \
#     --filter \
#     --verbose

# ~/Analysis/tools/vcf_to_pandas_dataframe.py \
#     --vcf CDC/ERR552797_DR_loci_raw_annotation.vcf \
#     --samplename ERR552797 \
#     --tsv CDC/ERR552797_DR_loci_raw_annotation.tsv \
#     --filter \
#     --verbose


# what if vcf has no mutations, just header
# ~/Analysis/tools/vcf_to_pandas_dataframe.py \
#     --vcf NC_000962.3_86_DR_loci_raw_annotation.vcf \
#     --samplename  NC_000962.3_86 \
#     --tsv NC_000962.3_86_DR_loci_raw_annotation.tsv \
#     --filter \
#     --verbose

#
# total rewrite of vcf_to_pandas_dataframe
# to account for all annotation for a variant
#
# ~/Analysis/tools/vcf_to_pandas_dataframe_2.py \
#     --vcf CDC/all_variants.vcf \
#     --samplename ERR552797 \
#     --tsv all_variants_2.tsv \
#     --verbose

~/Analysis/tools/vcf_to_pandas_dataframe_2.py \
    --vcf CDC/NC_000962.3_39/all_variants.vcf \
    --samplename NC_000962.3_39 \
    --tsv all_variants_NC_000962.3_2.tsv \
    --verbose

exit


#~/Analysis/tools/variant_interpretation.py --tsv ERR552797_cdc_snpsift.tsv \


#~/Analysis/tools/coverage.py --bam $BAM \
# 			     --bed $BED \
# 			     --minimum_coverage 10 \
# 			     --tsv ERR552797_depth_of_coverage_1.tsv

#
# CDC 
#
VCF=./CDC/ERR552797_DR_loci_raw_annotation.vcf.gz
#BED=amp.bed
#BED=cdc_regions.bed
# ~/Analysis/tools/variant_interpretation.py --vcf $VCF \
# 					   --bam $BAM \
# 					   --bed $BED \
# 					   --json $JSON \
# 					   --samplename "ERR552797" \
# 					   --minimum_coverage 10 \
# 					   --verbose \
# 					   --filter_genes \
# 					   --report ERR552797_tbdb_interpretation_1.tsv

# echo #########
# echo $SAMPLE
# echo #########

# # test with empty vcf file
# #SAMPLE=empty
# #VCF=./CDC/${SAMPLE}.vcf

# VCF=./CDC/all_variants.vcf # this file has also INS, INV, DUP, BND

# LAB_REPORT=./CDC/${SAMPLE}_lab_report_wDST_$(date +"%Y%m%d")_DB.csv

# #~/Analysis/tools/get_deletions_in_region.py --vcf $VCF --bed $BED --verbose


# #python3 -m scalene --outfile variant_interpretation.profile
# ####					   --filtered_vcf "filtered_vcf.vcf" \
# ~/Analysis/tools/variant_interpretation.py --vcf $VCF \
# 					   --bam $BAM \
# 					   --bed $BED \
# 					   --json $JSON \
# 					   --samplename $SAMPLE \
# 					   --minimum_coverage 10 \
#  					   --minimum_total_depth 10 \
#  					   --minimum_variant_depth 10 \
#  					   --minimum_allele_percentage 10 \
# 					   --filter_genes \
# 					   --report $LAB_REPORT


# #~/Analysis/tools/lineage_parser.py CDC/ERR552797_DR_loci_Final_annotation.txt CDC/lineage_markers.txt ERR552797_lineages.tsv ERR552797

# #~/Analysis/tools/get_lineage.py \
# #    --vcf CDC/ERR552797_DR_loci_raw_annotation.vcf \
# #    --lineages CDC/lineage_markers.txt \
# #    --report ERR552797_lineages.tsv \
# #    --samplename ERR552797

# OPERATOR=DB
# LIMS_REPORT=./CDC/${SAMPLE}_LIMSwDST_$(date +"%Y%m%d")_${OPERATOR}.csv
# ~/Analysis/tools/lims_report.py \
#     --lab $LAB_REPORT \
#     --lims $LIMS_REPORT \
#     --operator ${OPERATOR} \
#     --lineages ERR552797_lineages.tsv

# #csvtk -lt transpose lims.tsv 


# #
# # CDC, but try anothter test sample with a large deletion in katG
# #
# echo #########
# echo M21C00365
# echo #########

# SAMPLE=M21C00365
# VCF=./CDC/${SAMPLE}/all_variants.vcf
# BAM=./CDC/${SAMPLE}/M21C00365_s.bam
# LAB_REPORT=./CDC/${SAMPLE}/${SAMPLE}_lab_report_wDST_$(date +"%Y%m%d")_DB.csv

# ~/Analysis/tools/get_deletions_in_region.py --vcf $VCF --bed $BED --filter_SVs

# ~/Analysis/tools/variant_interpretation.py --vcf $VCF \
# 					   --bam $BAM \
# 					   --bed $BED \
# 					   --json $JSON \
# 					   --samplename $SAMPLE \
# 					   --minimum_coverage 10 \
#  					   --minimum_total_depth 10 \
#  					   --minimum_variant_depth 10 \
#  					   --minimum_allele_percentage 10 \
# 					   --filter_genes \
# 					   --filter_variants \
# 					   --debug \
# 					   --report $LAB_REPORT


# OPERATOR=DB
# LIMS_REPORT=./CDC/${SAMPLE}/${SAMPLE}_LIMSwDST_$(date +"%Y%m%d")_${OPERATOR}.csv
# ~/Analysis/tools/lims_report.py \
#     --lab $LAB_REPORT \
#     --lims $LIMS_REPORT \
#     --operator ${OPERATOR} \
#     --lineages ERR552797_lineages.tsv

# #
# # CDC, but try anothter test sample with a large deletion in katG
# #
# echo #########
# echo M21C00422
# echo #########
# SAMPLE=M21C00422
# VCF=./CDC/${SAMPLE}/all_variants.vcf
# BAM=./CDC/${SAMPLE}/M21C00422_s.bam
# LAB_REPORT=./CDC/${SAMPLE}/${SAMPLE}_lab_report_wDST_$(date +"%Y%m%d")_DB.csv

# ~/Analysis/tools/get_deletions_in_region.py --vcf $VCF --bed $BED --filter_SVs

# ~/Analysis/tools/variant_interpretation.py --vcf $VCF \
# 					   --bam $BAM \
# 					   --bed $BED \
# 					   --json $JSON \
# 					   --samplename $SAMPLE \
# 					   --minimum_coverage 10 \
#  					   --minimum_total_depth 10 \
#  					   --minimum_variant_depth 10 \
#  					   --minimum_allele_percentage 10 \
# 					   --filter_genes \
# 					   --filter_variants \
# 					   --debug \
# 					   --verbose \
# 					   --report $LAB_REPORT


# OPERATOR=DB
# LIMS_REPORT=./CDC/${SAMPLE}/${SAMPLE}_LIMSwDST_$(date +"%Y%m%d")_${OPERATOR}.csv
# ~/Analysis/tools/lims_report.py \
#     --lab $LAB_REPORT \
#     --lims $LIMS_REPORT \
#     --operator ${OPERATOR} \
#     --lineages ERR552797_lineages.tsv


#
# CDC, but try anothter test sample with a large deletion in katG
#
SAMPLE=NC_000962.3_39
echo #########
echo $SAMPLE
echo #########
VCF=./CDC/${SAMPLE}/all_variants.vcf
BAM=./CDC/${SAMPLE}/${SAMPLE}_s.bam
LAB_REPORT=./CDC/${SAMPLE}/${SAMPLE}_lab_report_wDST_$(date +"%Y%m%d")_DB.csv

#~/Analysis/tools/get_deletions_in_region.py --vcf $VCF --bed $BED --filter_SVs

~/Analysis/tools/variant_interpretation.py --vcf $VCF \
					   --bam $BAM \
					   --bed $BED \
					   --json $JSON \
					   --samplename $SAMPLE \
					   --minimum_coverage 10 \
 					   --minimum_total_depth 10 \
 					   --minimum_variant_depth 10 \
 					   --minimum_allele_percentage 10 \
					   --filter_genes \
					   --debug \
					   --verbose \
					   --report $LAB_REPORT


OPERATOR=DB
LIMS_REPORT=./CDC/${SAMPLE}/${SAMPLE}_LIMSwDST_$(date +"%Y%m%d")_${OPERATOR}.csv
~/Analysis/tools/lims_report.py \
    --lab $LAB_REPORT \
    --lims $LIMS_REPORT \
    --operator ${OPERATOR} \
    --lineages ERR552797_lineages.tsv

exit

#
# London
#
VCF=./London/ERR552797.targets.csq.vcf
BAM=./London/ERR552797.bam
BED=$HOME/Software/TBProfiler/db/tbdb.bed
LAB_REPORT=./London/${SAMPLE}_lab_report_wDST_$(date +"%Y%m%d")_DB.csv
$HOME/Analysis/tools/variant_interpretation.py --vcf $VCF \
 					       --bam $BAM \
 					       --json $JSON \
 					       --bed $BED \
 					       --samplename "ERR552797" \
 					       --minimum_coverage 10 \
 					       --minimum_total_depth 10 \
 					       --minimum_variant_depth 10 \
 					       --minimum_allele_percentage 10 \
					       --verbose \
					       --filter_genes \
					       --report $LAB_REPORT

