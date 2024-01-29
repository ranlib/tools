#!/usr/bin/bash

export CENTRIFUGE_INDEXES=/data/centrifuge

#
#
#
R1=production/v1/WA0267320-WAPHL-M8782-231122/20240125_152712_wf_varpipe/call-wf_clockwork_decontamination/out/clean_reads_1/WA0267320-WAPHL-M8782-231122_clockwork_cleaned_1.fq.gz
R2=production/v1/WA0267320-WAPHL-M8782-231122/20240125_152712_wf_varpipe/call-wf_clockwork_decontamination/out/clean_reads_2/WA0267320-WAPHL-M8782-231122_clockwork_cleaned_2.fq.gz
N=WA0267320-WAPHL-M8782-231122
#time centrifuge -x p_compressed+h+v -1 $R1 -2 $R2 --report-file ${N}.centrifuge.summary.report.tsv -S ${N}.centrifuge.classification.tsv
#time centrifuge -x p_compressed -1 $R1 -2 $R2 --report-file ${N}.centrifuge.summary.report.tsv -S ${N}.centrifuge.classification.tsv
#time centrifuge -x hpvc --threads 24 -1 $R1 -2 $R2 --report-file ${N}.centrifuge.summary.report.tsv -S ${N}.centrifuge.classification.tsv

# make kraken style report
#time centrifuge-kreport -x hpvc WA0267320-WAPHL-M8782-231122.centrifuge.classification.tsv > WA0267320-WAPHL-M8782-231122.centrifuge.classification.kraken_style.tsv 

#
#
#
R1=production/v2/WA0236590-WAPHL-M5916-231103/20240126_135124_wf_varpipe/call-wf_clockwork_decontamination/out/clean_reads_1/WA0236590-WAPHL-M5916-231103_clockwork_cleaned_1.fq.gz
R2=production/v2/WA0236590-WAPHL-M5916-231103/20240126_135124_wf_varpipe/call-wf_clockwork_decontamination/out/clean_reads_2/WA0236590-WAPHL-M5916-231103_clockwork_cleaned_2.fq.gz
N=WA0236590-WAPHL-M5916-231103
time centrifuge -x hpvc --threads 24 -1 $R1 -2 $R2 --report-file ${N}.centrifuge.summary.report.tsv -S ${N}.centrifuge.classification.tsv

exit 0

