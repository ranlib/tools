#!/bin/bash

## calculate Ts/Tv across all sites (includes AC=0)
## file must be b/gzipped
## ./tstv.sh file.vcf.gz

VCF=$1

# get count for transitions:
ag=$(zcat ${VCF} | awk '! /\#/' | awk '{if(length($4) == 1 && length($5) == 1) print}' | \
awk '($4 == "A" && $5 == "G" || $4 == "G" && $5 == "A")' | wc -l)
echo "AG:" $ag

ct=$(zcat ${VCF} | awk '! /\#/' | awk '{if(length($4) == 1 && length($5) == 1) print}' | \
awk '($4 == "C" && $5 == "T" || $4 == "T" && $5 == "C")' | wc -l)
echo "CT:" $ct

# transversions:
ac=$(zcat ${VCF} | awk '! /\#/' | awk '{if(length($4) == 1 && length($5) == 1) print}' | \
awk '($4 == "A" && $5 == "C" || $4 == "C" && $5 == "A")' | wc -l)
echo "AC:" $ac

at=$(zcat ${VCF} | awk '! /\#/' | awk '{if(length($4) == 1 && length($5) == 1) print}' | \
awk '($4 == "A" && $5 == "T" || $4 == "T" && $5 == "A")' | wc -l)
echo "AT:" $at

cg=$(zcat ${VCF} | awk '! /\#/' | awk '{if(length($4) == 1 && length($5) == 1) print}' | \
awk '($4 == "G" && $5 == "C" || $4 == "C" && $5 == "G")' | wc -l)
echo "CG:" $cg

gt=$(zcat ${VCF} | awk '! /\#/' | awk '{if(length($4) == 1 && length($5) == 1) print}' | \
awk '($4 == "G" && $5 == "T" || $4 == "T" && $5 == "G")' | wc -l)
echo "GT:" $gt


sum_ts=$((ag + ct))
sum_tv=$((ac + at + cg + gt))
ts_tv=$(awk "BEGIN {printf \"%.2f\",${sum_ts}/${sum_tv}}")
echo "Ts = $sum_ts"
echo "Tv = $sum_tv"
echo "Ts/Tv = $ts_tv"
