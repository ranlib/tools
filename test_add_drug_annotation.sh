#!/bin/bash

#
# test TB mutations interpretation procedure
#

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


#~/Analysis/tools/variant_interpretation.py --tsv ERR552797_cdc_snpsift.tsv \


#~/Analysis/tools/coverage.py --bam $BAM \
# 			     --bed $BED \
# 			     --minimum_coverage 10 \
# 			     --tsv ERR552797_depth_of_coverage_1.tsv


# CDC 
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

VCF=all_variants.vcf
#python3 -m scalene --outfile variant_interpretation.profile
####					   --filtered_vcf "filtered_vcf.vcf" \
~/Analysis/tools/variant_interpretation.py --vcf $VCF \
					   --bam $BAM \
					   --bed $BED \
					   --json $JSON \
					   --samplename "ERR552797" \
					   --minimum_coverage 10 \
 					   --minimum_total_depth 10 \
 					   --minimum_variant_depth 10 \
 					   --minimum_allele_percentage 10 \
					   --verbose \
					   --filter_genes \
					   --report ERR552797_tbdb_interpretation_2.tsv

# London
# VCF=./London/ERR552797.targets.csq.vcf
# BAM=./London/ERR552797.bam
# ~/Analysis/tools/variant_interpretation.py --vcf $VCF \
# 					   --bam $BAM \
# 					   --json $JSON \
# 					   --bed $BED \
# 					   --samplename "ERR552797" \
# 					   --minimum_coverage 10 \
# 					   --minimum_total_depth 10 \
# 					   --minimum_variant_depth 10 \
# 					   --minimum_allele_percentage 10 \
# 					   --report ERR552797_London_interpretation.tsv

#~/Analysis/tools/lineage_parser.py CDC/ERR552797_DR_loci_Final_annotation.txt CDC/lineage_markers.txt ERR552797_lineages.tsv ERR552797

#~/Analysis/tools/get_lineage.py \
#    --vcf CDC/ERR552797_DR_loci_raw_annotation.vcf \
#    --lineages CDC/lineage_markers.txt \
#    --report ERR552797_lineages.tsv \
#    --samplename ERR552797

~/Analysis/tools/lims_report.py \
    --lab ERR552797_tbdb_interpretation_2.tsv \
    --lims LIMSwDST_$(date +"%Y%m%d")_DB.csv \
    --operator DB \
    --lineages ERR552797_lineages.tsv

#~/Analysis/tools/lims_report.py \
#    --lab ERR552797_tbdb_interpretation_2.tsv \
#    --lims lims.tsv \
#    --operator DB \
#    --bed $BED

#csvtk -lt transpose lims.tsv 

				    
#~/Analysis/tools/get_deletions_in_region.py --vcf $VCF --bed $BED

