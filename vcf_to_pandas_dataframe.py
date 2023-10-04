#!/usr/bin/env python
"""
read vcf file, loop over record, write tsv file
"""
import re
import pandas
import argparse
from collections import OrderedDict
import vcf


def vcf_to_pandas_dataframe(vcf_file: str, samplename: str, filter_variants: bool, verbose: bool) -> pandas.DataFrame:
    """
    parse vcf file

    :param str vcf_file: filename of annotated vcf file
    :param str samplename: name of sample
    :param str bed_file: tab separated ascii file with 6 columns ["genome", "start", "stop", "locus", "gene", "chemical"]
    :return: pandas data frame
    """

    # loop over vcf file
    vcf_reader = vcf.Reader(filename=vcf_file)
    i = 0
    
    # get annotation terms from header
    # check there are mutations in input vcf file
    if "ANN" in vcf_reader.infos:
        annotation = vcf_reader.infos["ANN"].desc
        _, terms = re.split(r":[ ]+", annotation)
        annotation_keys = terms.replace("'", "").split(" | ")
    else:
        print(f"<W> vcf_to_pandas_dataframe: no mutations in {vcf_file}, return empty dataframe")
        df = pandas.DataFrame()
        return df

    # loop over records
    ANNS = []
    for record in vcf_reader:
        vcf_item = OrderedDict()
        vcf_item["Sample ID"] = samplename
        vcf_item["CHROM"] = record.CHROM
        vcf_item["POS"] = record.POS
        # vcf_item["ID"] = "." if not record.ID else record.ID
        vcf_item["REF"] = record.REF
        vcf_item["ALT"] = ",".join([str(n) for n in record.ALT])
        # vcf_item["QUAL"] = "." if not record.QUAL else record.QUAL
        vcf_item["FILTER"] = "." if not record.FILTER else ";".join(record.FILTER)
        #print(i, record, flush=True)
        i += 1
        
        if filter_variants and record.FILTER:
            print(f"<W> vcf_to_pandas_dataframe: filter out variant {record} with filter column entry {record.FILTER}")
            continue

        info = record.INFO
        if not (("ANN" in info) and (len(info["ANN"]) > 0)):
            annotation_item = dict(zip(annotation_keys, list(range(len(annotation_keys)))))
            annotation_item["cDNA.pos"], annotation_item["cDNA.length"] = [-1, -1]
            annotation_item["CDS.pos"], annotation_item["CDS.length"] = [-1, -1]
            annotation_item["AA.pos"], annotation_item["AA.length"] = [-1, -1]
            annotation_item["Total Read Depth"] = -1
            annotation_item["AD_REF"], annotation_item["Variant Read Depth"] = [-1, -1]
            annotation_item["Percent Alt Allele"] = -1 if annotation_item["Total Read Depth"] <= 0 else annotation_item["Variant Read Depth"] * 100.0 / annotation_item["Total Read Depth"]
        else:
            # check first entry in annotation list
            item = info["ANN"][0]
            annotation_item = dict(zip(annotation_keys, item.split("|")))

            # if first entry in annotation list is modifier then look for the closest annotation
            is_modifier = annotation_item["Annotation_Impact"] == "MODIFIER"

            is_SV = False
            if "SVTYPE" in info:
                is_SV = True

            #is_large_deletion = "<DEL>" in vcf_item["ALT"]
            #if is_modifier and not is_large_deletion:
            if is_modifier and not is_SV:
                # print(annotation_item)
                distances = []
                annotation_item_list = []
                for annotation in info["ANN"]:
                    annotation_item = dict(zip(annotation_keys, annotation.split("|")))
                    annotation_item_list.append(annotation_item)
                    distances.append(annotation_item["Distance"])

                distances.remove("")
                distances = list(map(int, distances))
                # closest_distance = min(distances)
                closest_index = distances.index(min(distances))
                annotation_item = annotation_item_list[closest_index]

            cDNA = annotation_item["cDNA.pos / cDNA.length"].split("/")
            CDS = annotation_item["CDS.pos / CDS.length"].split("/")
            AA = annotation_item["AA.pos / AA.length"].split("/")

            annotation_item["cDNA.pos"], annotation_item["cDNA.length"] = [-1, -1]
            annotation_item["CDS.pos"], annotation_item["CDS.length"] = [-1, -1]
            annotation_item["AA.pos"], annotation_item["AA.length"] = [-1, -1]

            if len(cDNA) == 2:
                annotation_item["cDNA.pos"], annotation_item["cDNA.length"] = list(map(int, cDNA))

            if len(CDS) == 2:
                annotation_item["CDS.pos"], annotation_item["CDS.length"] = list(map(int, CDS))

            if len(AA) == 2:
                annotation_item["AA.pos"], annotation_item["AA.length"] = list(map(int, AA))

            the_call = record.genotype(samplename)

            is_precise_SV = None
            if is_SV and "PRECISE" in info:
                is_precise_SV = True
            else:
                is_precise_SV = False

            if hasattr(the_call.data, "DP"):
                annotation_item["Total Read Depth"] = record.genotype(samplename)["DP"]
            elif hasattr(the_call.data,"DR") and hasattr(the_call.data,"DV") and not is_precise_SV:
                annotation_item["Total Read Depth"] = record.genotype(samplename)["DR"] + record.genotype(samplename)["DV"]
            elif hasattr(the_call.data,"RR") and hasattr(the_call.data,"RV") and is_precise_SV:
                annotation_item["Total Read Depth"] = record.genotype(samplename)["RR"] + record.genotype(samplename)["RV"]
            else:
                annotation_item["Total Read Depth"] = -1

            if hasattr(the_call.data, "AD"):
                annotation_item["AD_REF"], annotation_item["Variant Read Depth"] = record.genotype(samplename)["AD"]
            elif hasattr(the_call.data,"DV") and hasattr(the_call.data,"DV") and not is_precise_SV:
                annotation_item["AD_REF"] = record.genotype(samplename)["DR"]
                annotation_item["Variant Read Depth"] = record.genotype(samplename)["DV"]
            elif hasattr(the_call.data,"RV") and hasattr(the_call.data,"RV") and is_precise_SV:
                annotation_item["AD_REF"] = record.genotype(samplename)["RR"]
                annotation_item["Variant Read Depth"] = record.genotype(samplename)["RV"]
            else:
                annotation_item["AD_REF"], annotation_item["Variant Read Depth"] = [-1, -1]

            # annotation_item["Percent Alt Allele"] = record.genotype(samplename)["AF"]
            annotation_item["Percent Alt Allele"] = -1 if annotation_item["Total Read Depth"] <= 0 else annotation_item["Variant Read Depth"] * 100.0 / annotation_item["Total Read Depth"]
            ANNS.append(vcf_item | annotation_item)

    #print(i,flush=True)
    # output dataframe
    df = pandas.DataFrame(ANNS)
    # remove some columns
    forget_this = [
        "Allele",
        "Feature_Type",
        "Feature_ID",
        "Transcript_BioType",
        "Rank",
        "cDNA.pos / cDNA.length",
        "CDS.pos / CDS.length",
        "AA.pos / AA.length",
        "Distance",
        "cDNA.pos",
        "cDNA.length",
        "CDS.length",
        "AA.pos",
        "AA.length",
        "AD_REF",
    ]
    df = df.drop(forget_this, axis=1)
    df = df.drop(df.columns[df.columns.str.contains("ERRORS / WARNINGS / INFO")], axis=1)
    # rename some columns
    #new_names = {"Gene_Name": "Gene Name", "Gene_ID": "Gene ID", "HGVS.c": "Nucleotide Change", "HGVS.p": "Amino acid Change", "CDS.pos": "Position within CDS"}
    #df = df.rename(columns=new_names)
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
