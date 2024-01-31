#!/bin/bash

#
# CDC, some simulated test sample with single mutation
#

# for drug annotation
JSON=$HOME/Analysis/varpipe4/intervals/tbdb.dr.json
BED=$HOME/Analysis/varpipe4/intervals/tbdb.bed

# for lims report
OPERATOR=DB

#
# test with empy vcf
#
SAMPLE=empty
echo #########
echo $SAMPLE
echo #########
VCF=./CDC/${SAMPLE}.vcf
LAB_REPORT=./CDC/${SAMPLE}_lab_report_wDST_$(date +"%Y%m%d")_DB.csv
BAM=./CDC/ERR552797_sdrcsm.bam
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
					   --report $LAB_REPORT
OPERATOR=DB
LIMS_REPORT=./CDC/${SAMPLE}_LIMSwDST_$(date +"%Y%m%d")_${OPERATOR}.csv
~/Analysis/tools/lims_report.py \
    --lab $LAB_REPORT \
    --lims $LIMS_REPORT \
    --operator ${OPERATOR} \
    --lineages ERR552797_lineages.tsv

#
# test with simulated data
#
for i in 4 5 6 39 ; do
   SAMPLE=NC_000962.3_${i}
   echo #########
   echo $SAMPLE
   echo #########
   VCF=./CDC/NC_000962.3_${i}/_LAST/out/vcf/all_variants.vcf
   BAM=./CDC/NC_000962.3_${i}/_LAST/out/bam/${SAMPLE}_s.bam
   LAB_REPORT=./CDC/${SAMPLE}/${SAMPLE}_lab_report_wDST_$(date +"%Y%m%d")_DB.csv
   #samtools index $BAM
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
done

exit 0
