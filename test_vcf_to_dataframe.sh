#!/bin/bash

#
# test vcf to dataframe conversion
#

OUTPUT=output/vcf_to_pandas_dataframe_all_annotations
mkdir -p $OUTPUT

BED=$HOME/Analysis/varpipe4/intervals/tbdb.bed
VCF=./CDC/ERR552797_DR_loci_raw_annotation.vcf

# ~/Analysis/tools/vcf_to_pandas_dataframe_all_annotations.py \
#     --vcf $VCF \
#     --samplename ERR552797 \
#     --bed $BED \
#     --tsv all_variants_all_annotations.tsv \
#     --tsv_filtered all_variants_NC_000962.3_filtered.tsv \
#     --filter \
#     --verbose

for i in 4 5 6 39 ; do
    VCF=./CDC/NC_000962.3_${i}/_LAST/out/vcf/all_variants.vcf
    ~/Analysis/tools/vcf_to_pandas_dataframe_all_annotations.py \
	--vcf $VCF \
	--samplename NC_000962.3_${i} \
	--bed $BED \
	--tsv $OUTPUT/all_variants_NC_000962.3_${i}.tsv \
	--tsv_filtered $OUTPUT/all_variants_NC_000962.3_${i}_filtered.tsv \
	--filter \
	--verbose
done

#
# for comparison the "old" vcf_to_pandas_dataframe
#
#
echo "###################################"
echo "# results of old vcf => dataframe #"
echo "###################################"
OUTPUT=output/vcf_to_pandas_dataframe
mkdir -p $OUTPUT

for i in 4 5 6 39 ; do 
    VCF=./CDC/NC_000962.3_${i}/_LAST/out/vcf/all_variants.vcf
    ~/Analysis/tools/vcf_to_pandas_dataframe.py \
	--vcf $VCF \
	--samplename NC_000962.3_${i} \
	--tsv $OUTPUT/all_variants_NC_000962.3_${i}.tsv \
	--filter \
	--verbose
done

exit 0
