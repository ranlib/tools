#!/usr/bin/env python
"""
read vcf file, loop over record, write tsv file
"""

import re
import vcf
import pandas
import argparse
from collections import OrderedDict

def vcf_to_pandas_dataframe(vcf_file: str, sample_name: str) -> pandas.DataFrame:
    """
    parse vcf file

    :param str vcf_file: filename of vcf file
    :param str sample_name: name of sample
    :return: pandas data frame
    """
    vcf_reader = vcf.Reader(filename=vcf_file)

    # get annotation terms from header
    # check there are mutations in input vcf file
    if "ANN" in vcf_reader.infos:
        annotation = vcf_reader.infos["ANN"].desc
        _, terms = re.split(r":[ ]+", annotation)
        ann_keys = terms.replace("'", "").split(" | ")
    else:
        print(f"<W> vcf_to_pandas_dataframe: no mutations in {vcf_file}, return empty dataframe")
        df = pandas.DataFrame()
        print(df)
        print(len(df.index))
        return df
        
    # loop over records
    ANNS = []
    for record in vcf_reader:
        vcf_item = OrderedDict()
        vcf_item["Sample ID"] = sample_name
        vcf_item["CHROM"] = record.CHROM
        vcf_item["POS"] = record.POS
        #vcf_item["ID"] = "." if not record.ID else record.ID
        vcf_item["REF"] = record.REF
        vcf_item["ALT"] = ','.join([ str(n) for n in record.ALT ])
        #vcf_item["QUAL"] = "." if not record.QUAL else record.QUAL
        vcf_item["FILTER"] = "." if not record.FILTER else record.FILTER
        info = record.INFO
        vcf_item["NANN"] = len(info["ANN"])
        item = info["ANN"][0]
        ann_item = dict(zip(ann_keys, item.split("|")))
        cDNA = ann_item["cDNA.pos / cDNA.length"].split("/")
        CDS = ann_item["CDS.pos / CDS.length"].split("/")
        AA = ann_item["AA.pos / AA.length"].split("/")

        ann_item["cDNA.pos"], ann_item["cDNA.length"] = [-1,-1]
        ann_item["CDS.pos"], ann_item["CDS.length"] = [-1,-1]
        ann_item["AA.pos"], ann_item["AA.length"] = [-1,-1]

        if len(cDNA) == 2:
            ann_item["cDNA.pos"], ann_item["cDNA.length"] = list( map(int, cDNA) )
        if len(CDS) == 2:
            ann_item["CDS.pos"], ann_item["CDS.length"] = list( map(int, CDS))
        if len(AA) == 2:
            ann_item["AA.pos"], ann_item["AA.length"] = list( map(int, AA))

        ann_item["Total Read Depth"] = record.genotype(sample_name)["DP"]
        ann_item["AD_REF"], ann_item["Variant Read Depth"] = record.genotype(sample_name)["AD"]
        #ann_item["Percent Alt Allele"] = record.genotype(sample_name)["AF"]
        ann_item["Percent Alt Allele"] = -1 if ann_item["Total Read Depth"] == 0 else ann_item["Variant Read Depth"] * 100.0 / ann_item["Total Read Depth"]

        # kick out some keys
        forget_this = [ "Allele",
                        "Feature_Type",
                        "Feature_ID",
                        "Transcript_BioType",
                        "Rank",
                        "cDNA.pos / cDNA.length",
                        "CDS.pos / CDS.length",
                        "AA.pos / AA.length",
                        "Distance",
                        "ERRORS / WARNINGS / INFO ",
                        "cDNA.pos",
                        "cDNA.length",
                        "CDS.length",
                        "AA.pos",
                        "AA.length",
                        "AD_REF"
                       ]
        for item in forget_this:
            if item in ann_item:
                ann_item.pop(item)

        # rename some keys
        ann_item['Gene Name'] = ann_item.pop("Gene_Name")
        ann_item['Gene ID'] = ann_item.pop("Gene_ID")
        ann_item['Nucleotide Change'] = ann_item.pop("HGVS.c")
        ann_item['Amino acid Change'] = ann_item.pop("HGVS.p")
        ann_item['Position within CDS'] = ann_item.pop("CDS.pos")

        ANNS.append(vcf_item | ann_item)

    df = pandas.DataFrame(ANNS)
    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="add drug annotation to tsv form of annotated vcf file", prog="add_drug_annotation", formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=80))
    parser.add_argument("--vcf", "-v", type=str, dest = "vcf_file", help="annotated vcf file", required=True)
    parser.add_argument("--sample_name", "-s", type=str, dest = "sample_name", help="sample name", required=True)
    parser.add_argument("--tsv", "-t", type=str, dest = "output_tsv", help="output tsv file", required=True)
    args = parser.parse_args()

    df = vcf_to_pandas_dataframe(args.vcf_file, args.sample_name)
    df.to_csv(args.output_tsv, sep="\t", index=False)
