#!/bin/bash
R1=production/v1/WA0267320-WAPHL-M8782-231122/20240125_152712_wf_varpipe/call-wf_clockwork_decontamination/out/clean_reads_1/WA0267320-WAPHL-M8782-231122_clockwork_cleaned_1.fq.gz
R2=production/v1/WA0267320-WAPHL-M8782-231122/20240125_152712_wf_varpipe/call-wf_clockwork_decontamination/out/clean_reads_2/WA0267320-WAPHL-M8782-231122_clockwork_cleaned_2.fq.gz
N=WA0267320-WAPHL-M8782-231122

export KRAKEN2_DB_PATH=/data/kraken2

time docker run --rm -e KRAKEN2_DB_PATH=/data -v /data/kraken2:/data -v $PWD:/mnt -w /mnt staphb/kraken2 kraken2 \
     -db k2_standard_20230314 \
     --threads 16 \
     --report $N.kraken.report \
     --gzip-compressed \
     --unclassified-out $N.unclassified#.fq \
     --classified-out $N.classified#.fq \
     --paired $R1 $R2 > ${N}.Kraken.out

exit 0
