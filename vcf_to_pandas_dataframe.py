#!/usr/bin/env python
"""
read vcf file, loop over record, write tsv file
"""

import re
import vcf
import pandas


def vcf_to_pandas_dataframe(vcf_file: str, sample_name: str) -> pandas.DataFrame:
    """
    parse vcf file

    :param str vcf_file: filename of vcf file
    :param str sample_name: name of sample
    :return: pandas data frame
    """

    # open vcf file
    with open(vcf_file, "r", encoding="ascii") as vcf_file_handle:
        vcf_reader = vcf.Reader(vcf_file_handle)

        # get annotation terms from header
        annotation = vcf_reader.infos["ANN"].desc
        _, terms = re.split(r":[ ]+", annotation)
        ann_keys = terms.replace("'", "").split(" | ")

        # loop over records
        ANNS = []
        for record in vcf_reader:
            info = record.INFO
            item = info["ANN"][0]
            ann_item = dict(zip(ann_keys, item.split("|")))
            cDNA = ann_item["cDNA.pos / cDNA.length"].split("/")
            CDS = ann_item["CDS.pos / CDS.length"].split("/")
            AA = ann_item["AA.pos / AA.length"].split("/")
            if len(cDNA) == 2:
                ann_item["cDNA.pos"], ann_item["cDNA.length"] = cDNA
            if len(CDS) == 2:
                ann_item["CDS.pos"], ann_item["CDS.length"] = CDS
            if len(AA) == 2:
                ann_item["AA.pos"], ann_item["AA.length"] = AA

            ann_item["DP"] = record.genotype(sample_name)["DP"]
            ann_item["AF"] = record.genotype(sample_name)["AF"]
            ann_item["AD_REF"], ann_item["AD_ALT"] = record.genotype(sample_name)["AD"]
            ANNS.append(ann_item)

        df = pandas.DataFrame(ANNS)
    return df


if __name__ == "__main__":
    vcf_file = "CDC/test.vcf"
    sample_name = "ERR552797"
    df = vcf_to_pandas_dataframe(vcf_file, sample_name)
    df.to_csv("test.tsv", sep="\t", index=False)
