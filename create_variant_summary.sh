#!/bin/bash
set -euxo pipefail

SAMPLE="ERR552797"

#
# London vcf file
# change chromosome name
# range filter
#
# massage London output vcf to show correct chromosome name
LONDON="London/ERR552797.targets.vcf"
cp ~/Analysis/TBprofiler/ERR552797.targets.vcf "$LONDON"
sed -i "s|^Chromosome|NC_000962.3|" "$LONDON"

# filter London file on CDC ranges
LONDON_filtered="${LONDON/.vcf/.range_filtered.vcf}"
vcf_filter_range.py -v "$LONDON" -b amp_bed.txt -o "$LONDON_filtered"

#
# CDC vcf file
#
CDC="CDC/ERR552797_DR_loci_raw_annotation.vcf"
cp ~/Analysis/varpipe3/data/Output_05_26_2023/ERR552797/ERR552797_DR_loci_raw_annotation.vcf "$CDC"

#
# parse vcf file into tsv file
#
#create_annotation.py -v "$CDC" -n "$SAMPLE"  > ${SAMPLE}_cdc.tsv
#create_annotation.py -v "$LONDON_filtered" -n "$SAMPLE"  > ${SAMPLE}_London.tsv
parse_annotation.py -v "$CDC" -l mutation_loci.txt -n "$SAMPLE" -o ${SAMPLE}_cdc.tsv
parse_annotation.py -v "$LONDON_filtered" -l mutation_loci.txt -n "$SAMPLE" -o ${SAMPLE}_London.tsv

#
# create report
#
echo "CDC"
create_variant_summary.py ${SAMPLE}_cdc.tsv
echo
echo "London"
create_variant_summary.py ${SAMPLE}_London.tsv

