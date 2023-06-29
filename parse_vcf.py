#!/usr/bin/env python
import vcf

vcf_file = 'CDC/ERR552797_DR_loci_raw_annotation.vcf'
vcf_reader = vcf.Reader(open(vcf_file, 'r'))
for record in vcf_reader:
    info = record.INFO
    for element,value  in info.items():
        print(element, value)
    format = record.FORMAT
    for element in format.split(":"):
        print(element)

    genotype = record.genotype
    print(genotype('ERR552797'))
