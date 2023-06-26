#!/usr/bin/env python
import os
import pandas
import argparse

parser = argparse.ArgumentParser(description="run SnpSift & bcftools for vcf -> tsv")
parser.add_argument("--vcf", "-v", type=str, help="vcf file", required=True)
parser.add_argument("--samplename", "-s", type=str, help="sample name", required=True)
parser.add_argument("--output", "-o", type=str, help="tsv output file", required=True)
args = parser.parse_args()

bcftools = os.path.basename(args.vcf) + ".bcftools"
command = f'bcftools query -f "[%CHROM\t%POS\t%AD\t%DP\t%AF]\n" {args.vcf} > {bcftools}'
os.system(command)

#sed -i 's|,|\t|g' ERR552797_DR_loci_raw_annotation_bcftools.tsv
#( echo -e "CHROM\tPOS\tAD_REF\tAD_ALT\tDP\tAF" ; cat ERR552797_DR_loci_raw_annotation_bcftools.tsv ) > temp && mv temp ERR552797_DR_loci_raw_annotation_bcftools.tsv

snpsift = os.path.basename(args.vcf) + ".SnpSift"
command = f'java -jar ~/Software/snpEff/SnpSift.jar extractFields {args.vcf} \
     "CHROM" \
     "POS" \
     "REF" \
     "ALT" \
     "FILTER" \
     "AF" \
     "AD" \
     "DP" \
     "ANN[0].GENE" \
     "ANN[0].GENEID" \
     "ANN[0].EFFECT" \
     "ANN[0].IMPACT" \
     "ANN[0].HGVS_C" \
     "ANN[0].HGVS_P" \
     "ANN[0].CDS_POS" > {snpsift}'

os.system(command)

d = pandas.read_csv(bcftools, header=None, sep="\t")
d.columns = ["CHROM", "POS", "AD",  "DP", "AF"]
d[["AD_REF", "AD_ALT"]] = d["AD"].str.split(",", expand=True)
d = d.drop("AD", axis=1)
for col in ["AD_ALT", "AD_REF", "DP"]:
    d[col] = d[col].astype('int')
d["AF"] = d["AD_ALT"]/d["DP"].fillna(0)

s = pandas.read_csv(snpsift, sep="\t")
s["AF"] = d["AF"]
#s["Sample ID"] = [args.samplename] * len(s.index)
s.insert(0, "Sample ID" , [args.samplename] * len(s.index) )
s.columns = [ e.replace("ANN[0].","") for e in s.columns ]

#Sample ID

#CHROM   POS     REF     ALT     Read Depth      Percent Alt Allele      Annotation      Variant Type    Nucleotide Change       Position within CDS   Amino acid Change       REF Amino acid  ALT Amino acid  Codon Position  Gene Name       Gene ID

#CHROM   POS     REF     ALT     FILTER  AF      AD      DP      GENE    GENEID  EFFECT  IMPACT  HGVS_C  HGVS_P  CDS_POS

s.rename( columns = { "DP": "Read Depth"}, inplace = True)
s.rename( columns = { "AF": "Percent Alt Allele"}, inplace = True)
s.rename( columns = { "GENE": "Gene Name"}, inplace = True)
s.rename( columns = { "GENEID": "Gene ID"}, inplace = True)
s.rename( columns = { "EFFECT": "Annotation"}, inplace = True)
s.rename( columns = { "HGVS_C": "Nucleotide Change"}, inplace = True)
s.rename( columns = { "HGVS_P": "Amino acid Change"}, inplace = True)
s.rename( columns = { "CDS_POS": "Position within CDS"}, inplace = True)
s.to_csv(args.output, index=False, sep="\t")
