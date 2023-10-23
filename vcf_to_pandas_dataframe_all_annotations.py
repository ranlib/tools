#!/usr/bin/env python
"""
read vcf file, loop over record, write tsv file
"""
import re
import pandas
import argparse
from collections import OrderedDict
import vcf
import hgvs.parser

def vcf_to_pandas_dataframe(vcf_file: str, samplename: str, filter_variants: bool, verbose: bool) -> pandas.DataFrame:
    """
    parse vcf file

    :param str vcf_file: filename of annotated vcf file
    :param str samplename: name of sample
    :param str bed_file: tab separated ascii file with 6 columns ["genome", "start", "stop", "locus", "gene", "chemical"]
    :return: pandas data frame
    """

    # prepare default output tsv file
    cols = [ "Sample ID","CHROM","POS","REF","ALT","FILTER","Allele","Annotation","Annotation_Impact","Gene_Name","Gene_ID","Feature_Type","Feature_ID","Transcript_BioType","Rank","HGVS.c","HGVS.p","cDNA.pos / cDNA.length","CDS.pos / CDS.length","AA.pos / AA.length","Distance","cDNA.pos","cDNA.length","CDS.pos","CDS.length","AA.pos","AA.length","Total Read Depth","AD_REF","Variant Read Depth","Percent Alt Allele",]
    df = pandas.DataFrame(columns=cols)

    # read in vcf file
    vcf_reader = vcf.Reader(filename=vcf_file)
    
    # get annotation terms from header
    # check there are mutations in input vcf file
    if "ANN" in vcf_reader.infos:
        annotation = vcf_reader.infos["ANN"].desc
        _, terms = re.split(r":[ ]+", annotation)
        annotation_keys = terms.replace("'", "").split(" | ")
    else:
        print(f"<W> vcf_to_pandas_dataframe: no mutations in {vcf_file}, return empty dataframe")
        return df

    # loop over records
    ANNS = []
    for record in vcf_reader:
        info = record.INFO
        genotype = record.genotype(samplename)
        
        if "ANN" in info:
            if len(info["ANN"]) > 0:
                for index, item in enumerate(info["ANN"]):
                    this_annotation_item = OrderedDict(zip(annotation_keys, item.split("|")))
                    
                    cDNA = this_annotation_item["cDNA.pos / cDNA.length"].split("/")
                    CDS = this_annotation_item["CDS.pos / CDS.length"].split("/")
                    AA = this_annotation_item["AA.pos / AA.length"].split("/")

                    for element in ["cDNA.pos", "cDNA.length", "CDS.pos", "CDS.length", "AA.pos", "AA.length"]:
                        this_annotation_item[element] = -1
                        
                    if len(cDNA) == 2:
                        this_annotation_item["cDNA.pos"], this_annotation_item["cDNA.length"] = list(map(int, cDNA))

                    if len(CDS) == 2:
                        this_annotation_item["CDS.pos"], this_annotation_item["CDS.length"] = list(map(int, CDS))

                    if len(AA) == 2:
                        this_annotation_item["AA.pos"], this_annotation_item["AA.length"] = list(map(int, AA))


                    if "SVTYPE" in info:
                        is_SV = True
                    else:
                        is_SV = False
                        
                    if "PRECISE" in info:
                        if is_SV:
                            is_precise_SV = True
                        else:
                            is_precise_SV = False
                    else:
                        is_precise_SV = None

                    if hasattr(genotype.data, "DP"):
                        this_annotation_item["Total Read Depth"] = genotype["DP"]
                    elif hasattr(genotype.data,"DR") and hasattr(genotype.data,"DV") and not is_precise_SV:
                        this_annotation_item["Total Read Depth"] = genotype["DR"] + genotype["DV"]
                    elif hasattr(genotype.data,"RR") and hasattr(genotype.data,"RV") and is_precise_SV:
                        this_annotation_item["Total Read Depth"] = genotype["RR"] + genotype["RV"]
                    else:
                        this_annotation_item["Total Read Depth"] = -1

                    if hasattr(genotype.data, "AD"):
                        this_annotation_item["AD_REF"], this_annotation_item["Variant Read Depth"] = genotype["AD"]
                    elif hasattr(genotype.data,"DV") and hasattr(genotype.data,"DV") and not is_precise_SV:
                        this_annotation_item["AD_REF"] = genotype["DR"]
                        this_annotation_item["Variant Read Depth"] = genotype["DV"]
                    elif hasattr(genotype.data,"RV") and hasattr(genotype.data,"RV") and is_precise_SV:
                        this_annotation_item["AD_REF"] = genotype["RR"]
                        this_annotation_item["Variant Read Depth"] = genotype["RV"]
                    else:
                        this_annotation_item["AD_REF"] = -1
                        this_annotation_item["Variant Read Depth"] = -1

                    this_annotation_item["Percent Alt Allele"] = -1 if this_annotation_item["Total Read Depth"] <= 0 else this_annotation_item["Variant Read Depth"] * 100.0 / this_annotation_item["Total Read Depth"]

                    vcf_item = OrderedDict()
                    vcf_item["Sample ID"] = samplename
                    vcf_item["CHROM"] = record.CHROM
                    vcf_item["POS"] = record.POS
                    vcf_item["ID"] = "." if not record.ID else record.ID
                    vcf_item["REF"] = record.REF
                    vcf_item["ALT"] = ",".join([str(n) for n in record.ALT])
                    vcf_item["QUAL"] = "." if not record.QUAL else record.QUAL
                    vcf_item["FILTER"] = "." if not record.FILTER else ";".join(record.FILTER)

                    hgvs_parser = hgvs.parser.Parser()
                    variant = hgvs_parser.parse_hgvs_variant(record.CHROM + ":" + this_annotation_item["HGVS.c"])

                    vcf_item["ref_concordant"] = variant.posedit.edit.ref == record.REF

                    ANNS.append(vcf_item | this_annotation_item)

                    # if variant not modifier, take first annotation in snpEff annotation list
                    if index == 0 and this_annotation_item["Annotation_Impact"] != "MODIFIER":
                        break

    # output dataframe
    df = pandas.DataFrame(ANNS)

    # some cleanup of data frame here
    if "Distance" in df:
        df['Distance'] = df['Distance'].replace(to_replace="",value="-1")
        df["Distance"] = df['Distance'].astype(int)
        
    # select distance > 0
    df = df[df["Distance"]>0]

    # select same ref
    df = df[df["ref_concordant"]]

    # take min distance
    if filter_variants:
        df = df[ df.Distance == df.Distance.min() ]
    
    df.reset_index(drop=True)

    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="vcf_to_pandas_dataframe", prog="vcf_to_pandas_dataframe", formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=80))
    parser.add_argument("--vcf", "-v", type=argparse.FileType("r"), dest="vcf_file", help="annotated vcf file", required=True)
    parser.add_argument("--samplename", "-s", type=str, dest="samplename", help="sample name", required=True)
    parser.add_argument("--filter", "-f", action="store_true", help="filter out variants which are not labeled PASS in vcf filter column")
    parser.add_argument("--verbose", action="store_true", help="turn on debugging output")
    parser.add_argument("--tsv", "-t", type=str, dest="output_tsv", help="output tsv file", required=True)
    args = parser.parse_args()

    df = vcf_to_pandas_dataframe(args.vcf_file.name, args.samplename, args.filter, args.verbose)
    df.to_csv(args.output_tsv, sep="\t", index=False)
