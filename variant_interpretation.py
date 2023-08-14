#!/usr/bin/env python
"""
generate interpretation  
"""
import argparse
import json
import pandas
from sympy import Interval
from vcf_to_pandas_dataframe import vcf_to_pandas_dataframe
from coverage import calculate_average_depth
import io

def get_intervals(regions: pandas.DataFrame):
    """
    generate intervals
    """
    global rpoB_codon
    global mmpS5, mmpL5, rplC, pepQ, atpE, mmpR5
    global mmpR5_promoter, atpE_promoter, pepQ_promoter, rplC_promoter
    global rrl_rRNA_1, rrl_rRNA_1_complement

    for index, row in regions.iterrows():
        if row["gene"] == "mmpR5":
            mmpR5 = Interval(row["start"], row["stop"])
            mmpR5_promoter = Interval(row["start"] - 84, row["stop"] - 1)

        if row["gene"] == "atpE":
            atpE = Interval(row["start"], row["stop"])
            atpE_promoter = Interval(row["start"] - 48, row["stop"] - 1)

        if row["gene"] == "pepQ":
            pepQ = Interval(row["start"], row["stop"])
            pepQ_promoter = Interval(row["start"] - 33, row["stop"] - 1)

        if row["gene"] == "rplC":
            rplC = Interval(row["start"], row["stop"])
            rplC_promoter = Interval(row["start"] - 18, row["stop"] - 1)

        if row["gene"] == "mmpL5":
            mmpL5 = Interval(row["start"], row["stop"])

        if row["gene"] == "mmpS5":
            mmpS5 = Interval(row["start"], row["stop"])

        if row["gene"] == "rrl":
            rrl = Interval(row["start"], row["stop"])
            rrl_rRNA = Interval(1, row["stop"] - row["start"])
            # position in CDS
            rrl_rRNA_1 = Interval(2003, 2367).union(Interval(2449, 3056))
            # not in [2003, 2367] and not in [2449, 3056]
            rrl_rRNA_1_complement = rrl_rRNA_1.complement(rrl_rRNA)

        if row["gene"] == "rpoB":
            # position in CDS
            rpoB_codon = Interval(426, 452)


gene_list_1 = ["mmpR5", "atpE", "pepQ", "mmpL5", "mmpS5", "rrl", "rplC"]
gene_list_2 = ["katG", "pncA", "ethA", "gid", "rpoB"]
gene_list_3 = ["katG", "pncA", "ethA", "gid"]

looker_1_1 = {"": "no WHO confidence", "Assoc w R": "R", "Assoc w R - interim": "R - interim", "Not assoc w R": "S", "Not assoc w R - Interim": "S - interim", "Uncertain significance": "U"}

mdl_1_1 = {"": "no WHO confidence", "Assoc w R": "R", "Assoc w R - interim": "U", "Not assoc w R": "S", "Not assoc w R - Interim": "U", "Uncertain significance": "U"}

looker_2_2 = {"": "no WHO confidence", "Assoc w R": "R", "Assoc w R - interim": "R - interim", "Not assoc w R": "S", "Not assoc w R - Interim": "S - interim", "Uncertain significance": "U"}

mdl_2_2 = {"": "no WHO confidence", "Assoc w R": "R", "Assoc w R - interim": "U", "Not assoc w R": "S", "Not assoc w R - Interim": "S", "Uncertain significance": "S"}

looker_3_1 = {"": "no WHO confidence", "Assoc w R": "R", "Assoc w R - interim": "R - interim", "Not assoc w R": "S", "Not assoc w R - Interim": "S - interim", "Uncertain significance": "U"}

mdl_3_1 = {"": "no WHO confidence", "Assoc w R": "R", "Assoc w R - interim": "U", "Not assoc w R": "S", "Not assoc w R - Interim": "S", "Uncertain significance": "S"}


def get_interpretation_1_2(gene: str, genomic_position: int, cds_position: int, annotation: str) -> list[str]:
    """
    get interpretation for mmpR5, atpE, pepQ, mmpL5, mmpS5, rrl, rplC
    """
    is_synonymous = annotation == "synonymous_variant"
    is_nonsynonymous = annotation != "synonymous_variant"
    looker = mdl = ""

    if gene == "atpE":
        if atpE.contains(genomic_position) & is_synonymous:
            looker = mdl = "S"
        if atpE.contains(genomic_position) & is_nonsynonymous:
            looker = mdl = "U"
        if atpE_promoter.contains(genomic_position):
            looker = mdl = "U"

    if gene == "mmpL5":
        if mmpL5.contains(genomic_position) & is_synonymous:
            looker = mdl = "S"
        if mmpL5.contains(genomic_position) & is_nonsynonymous:
            looker = mdl = "U"

    if gene == "mmpR5":
        if mmpR5.contains(genomic_position) & is_synonymous:
            looker = mdl = "S"
        if mmpR5.contains(genomic_position) & is_nonsynonymous:
            looker = mdl = "U"
        if mmpR5_promoter.contains(genomic_position):
            looker = mdl = "U"

    if gene == "mmpS5":
        if mmpS5.contains(genomic_position) & is_synonymous:
            looker = mdl = "S"
        if mmpS5.contains(genomic_position) & is_nonsynonymous:
            looker = mdl = "U"

    if gene == "pepQ":
        if pepQ.contains(genomic_position) & is_synonymous:
            looker = mdl = "S"
        if pepQ.contains(genomic_position) & is_nonsynonymous:
            looker = mdl = "U"
        if pepQ_promoter.contains(genomic_position):
            looker = mdl = "U"

    if gene == "rplC":
        if rplC.contains(genomic_position) & is_synonymous:
            looker = mdl = "S"
        if rplC.contains(genomic_position) & is_nonsynonymous:
            looker = mdl = "U"
        if rplC_promoter.contains(genomic_position):
            looker = mdl = "U"

    if gene == "rrl":
        if rrl_rRNA_1.contains(cds_position):
            looker = mdl = "U"
        if rrl_rRNA_1_complement.contains(cds_position):
            looker = "S"
            mdl = "U"

    return [looker, mdl]


def get_interpretation_2_2_1(annotation: str) -> list[str]:
    """
    implementation of interpretation according to 2.2.1
    """
    is_synonymous = annotation == "synonymous_variant"
    is_nonsynonymous = annotation != "synonymous_variant"
    # ?
    effect_types = ["disruptive_inframe_insertion", "frameshift_variant", "stop_lost", "splice_region_variant", "missense_variant", "upstream_gene_variant", "disruptive_inframe_deletion", "nonsense_variant"]
    looker = mdl = ""
    if annotation in effect_types:
        looker = mdl = "U"
    else:
        if is_synonymous:
            looker = mdl = "S"
        if is_nonsynonymous:
            looker = "U"
            mdl = "S"
    return [looker, mdl]


def get_interpretation_2_2_2(cds_position: int, annotation: str) -> list[str]:
    """
    implementation of interpretation according to 2.2.2
    """
    is_synonymous = annotation == "synonymous_variant"
    is_nonsynonymous = annotation != "synonymous_variant"
    looker = mdl = ""
    if rpoB_codon.contains(cds_position):
        if is_synonymous:
            looker = mdl = "S"
        if is_nonsynonymous:
            looker = mdl = "R"
    else:
        if is_synonymous:
            looker = mdl = "S"
        if is_nonsynonymous:
            looker = "U"
            mdl = "S"

    return [looker, mdl]


def get_interpretation_3_2(annotation: str) -> list[str]:
    """
    implementation of interpretation according to 3.2
    """
    is_synonymous = annotation == "synonymous_variant"
    is_nonsynonymous = annotation != "synonymous_variant"
    looker = mdl = ""
    if is_synonymous:
        looker = mdl = "S"
    if is_nonsynonymous:
        looker = "U"
        mdl = "S"
    return [looker, mdl]


def add_drug_annotation(tsv: pandas.DataFrame, annotation: str) -> pandas.DataFrame:
    """
    add drug annotation to variants
    
    :param str inputfile: filename of vcf file in tsv format
    :param str annotatio: filename of json file wit drug annotation
    :return: pandas dataframe of input tsv file with drug annotation
    """
    with open(annotation, "r", encoding="utf-8") as jsonfile:
        json_annotation = json.load(jsonfile)

    #
    # add drug annotation from json file
    #

    # create placeholders
    tsv["antimicrobial"] = [""] * len(tsv.index)
    tsv["confidence"] = [""] * len(tsv.index)
    tsv["rationale"] = [""] * len(tsv.index)
    tsv["count"] = [1] * len(tsv.index)

    # create output dataframe
    tsv_out = pandas.DataFrame(columns=list(tsv.columns))

    # loop over input dataframe, fill output dataframe
    i = 0
    for index, row in tsv.iterrows():
        position = int(row["POS"])
        gene_id = row["Gene ID"]
        nucleotide_change = row["Nucleotide Change"]
        amino_acid_change = row["Amino acid Change"]
        drug_annotation = []
        if gene_id in json_annotation:
            if nucleotide_change in list(json_annotation[gene_id].keys()):
                if int(position) in json_annotation[gene_id][nucleotide_change]["genome_positions"]:
                    for index, entry in enumerate(json_annotation[gene_id][nucleotide_change]["annotations"]):
                        pair = ["", ""]
                        if "drug" in entry:
                            pair[0] = entry["drug"]
                        if "who_confidence" in entry:
                            pair[1] = entry["who_confidence"]
                        drug_annotation.append(pair)
            elif amino_acid_change in list(json_annotation[gene_id].keys()):
                if int(position) in json_annotation[gene_id][amino_acid_change]["genome_positions"]:
                    for index, entry in enumerate(json_annotation[gene_id][amino_acid_change]["annotations"]):
                        pair = ["", ""]
                        if "drug" in entry:
                            pair[0] = entry["drug"]
                        if "who_confidence" in entry:
                            pair[1] = entry["who_confidence"]
                        drug_annotation.append(pair)

        if len(drug_annotation) == 0:
            tsv_out.loc[i, tsv_out.columns] = row
            i = i + 1
        else:
            for idx, item in enumerate(drug_annotation):
                row["count"] = idx + 1
                row["antimicrobial"] = item[0]
                row["confidence"] = item[1]
                row["rationale"] = "WHO"
                tsv_out.loc[i, tsv_out.columns] = row
                i = i + 1
    return tsv_out


def run_interpretation(tsv: pandas.DataFrame, drug_info: {}, coverage_percentage: {}, coverage_average: {}, minimum_coverage: int, all_genes: bool):
    """
    interpretation
    """
    header = ["Sample ID", "Gene Name", "Gene ID", "POS", "Position within CDS", "Nucleotide Change", "Amino acid Change", "Annotation", "confidence", "antimicrobial", "Total Read Depth", "Variant Read Depth", "Percent Alt Allele", "rationale"]
    tsv_interpretation = tsv.loc[:, header]
    tsv_interpretation["average_coverage_in_region"] = [0] * len(tsv_interpretation.index)
    tsv_interpretation["percent_above_threshold"] = [0] * len(tsv_interpretation.index)
    tsv_interpretation["Percent Alt Allele"] = tsv_interpretation["Percent Alt Allele"]
    tsv_interpretation["Comment"] = [""] * len(tsv_interpretation.index)

    # If no antimicrobial information, take antimicrobial information for this range (gene)
    # Note: this could be a comma separated list of chemicals
    for index, row in tsv_interpretation.iterrows():
        if row["antimicrobial"] == "":
            tsv_interpretation.loc[index, "antimicrobial"] = drug_info[row["Gene Name"]] if row["Gene Name"] in drug_info else ""

    genes_count_dict = {}
    for index, row in tsv_interpretation.iterrows():
        if row["Gene Name"] in genes_count_dict:
            genes_count_dict[row["Gene Name"]] += 1
        else:
            genes_count_dict[row["Gene Name"]] = 1

        # 0. add region average coverage (depth) information
        tsv_interpretation.loc[index, "average_coverage_in_region"] = coverage_average[row["Gene Name"]] if row["Gene Name"] in coverage_average else -1
        tsv_interpretation.loc[index, "percent_above_threshold"] = coverage_percentage[row["Gene Name"]] if row["Gene Name"] in coverage_percentage else -1

        # 1.
        if row["Gene Name"] in gene_list_1:
            # 1.1
            if (row["confidence"] != "") & (row["antimicrobial"] != ""):
                tsv_interpretation.loc[index, "looker"] = looker_1_1[row["confidence"]]
                tsv_interpretation.loc[index, "mdl"] = mdl_1_1[row["confidence"]]

            # 1.2
            if (row["confidence"] == "") & (row["antimicrobial"] != ""):
                looker, mdl = get_interpretation_1_2(row["Gene Name"], row["POS"], row["Position within CDS"], row["Annotation"])
                tsv_interpretation.loc[index, "looker"] = looker
                tsv_interpretation.loc[index, "mdl"] = mdl
                tsv_interpretation.loc[index, "confidence"] = "no WHO annotation"
                tsv_interpretation.loc[index, "rationale"] = "expert rule 1.2"

        # 2.
        if row["Gene Name"] in gene_list_2:
            # 2.1
            if (row["confidence"] != "") & (row["antimicrobial"] != ""):
                tsv_interpretation.loc[index, "looker"] = looker_2_2[row["confidence"]]
                tsv_interpretation.loc[index, "mdl"] = mdl_2_2[row["confidence"]]

            # 2.2 "WHO expert rules"
            if (row["confidence"] == "") & (row["antimicrobial"] != ""):
                # 2.2.1 gene_list_2 without rpoB
                if row["Gene Name"] in gene_list_3:
                    looker, mdl = get_interpretation_2_2_1(row["Annotation"])
                    tsv_interpretation.loc[index, "looker"] = looker
                    tsv_interpretation.loc[index, "mdl"] = mdl
                    tsv_interpretation.loc[index, "confidence"] = "no WHO annotation"
                    tsv_interpretation.loc[index, "rationale"] = "expert rule 2.2.1"

                # 2.2.2 just rpoB
                if row["Gene Name"] == "rpoB":
                    looker, mdl = get_interpretation_2_2_2(row["Position within CDS"], row["Annotation"])
                    tsv_interpretation.loc[index, "looker"] = looker
                    tsv_interpretation.loc[index, "mdl"] = mdl
                    tsv_interpretation.loc[index, "confidence"] = "no WHO annotation"
                    tsv_interpretation.loc[index, "rationale"] = "expert rule 2.2.2"

        # 3.
        if row["Gene Name"] not in (gene_list_1 + gene_list_2):
            # 3.1
            if (row["confidence"] != "") & (row["antimicrobial"] != ""):
                tsv_interpretation.loc[index, "looker"] = looker_3_1[row["confidence"]]
                tsv_interpretation.loc[index, "mdl"] = mdl_3_1[row["confidence"]]

            # 3.2
            if (row["confidence"] == "") & (row["antimicrobial"] != ""):
                looker, mdl = get_interpretation_3_2(row["Annotation"])
                tsv_interpretation.loc[index, "looker"] = looker
                tsv_interpretation.loc[index, "mdl"] = mdl
                tsv_interpretation.loc[index, "confidence"] = "no WHO annotation"
                tsv_interpretation.loc[index, "rationale"] = "expert rule 3.2"

    # 4.
    # loop over genes of interest,
    # check if they have mutations,
    # is not in genes_count_dict ==> no mutations, i.e. WT
    # manufacture default row for output dataframe

    # keep only genes of interest
    if not all_genes:
        tsv_interpretation = tsv_interpretation[tsv_interpretation["Gene Name"].isin(gene_list_1 + gene_list_2)]

    # manufacture default rows for genes with no mutations and add them
    tsv_additional = pandas.DataFrame(columns=list(tsv_interpretation.columns))
    index = 0
    for sample in tsv_interpretation["Sample ID"].unique().tolist():
        for gene in gene_list_1 + gene_list_2:
            if gene not in genes_count_dict:
                average_coverage_in_region = coverage_average[gene] if gene in coverage_average else -1
                percent_above_threshold = coverage_percentage[gene] if gene in coverage_average else -1
                if average_coverage_in_region > minimum_coverage:
                    # 4.1 there is coverage
                    #for drug in drug_info[gene].split(","):
                    tsv_additional.loc[index, tsv_additional.columns] = ["N/A"] * len(tsv_additional.columns)
                    tsv_additional.loc[index, "Sample ID"] = sample
                    tsv_additional.loc[index, "Gene Name"] = gene
                    tsv_additional.loc[index, "rationale"] = "WT"
                    #tsv_additional.loc[index, "antimicrobial"] = drug
                    tsv_additional.loc[index, "antimicrobial"] =  drug_info[gene] if gene in drug_info else ""
                    tsv_additional.loc[index, "average_coverage_in_region"] = average_coverage_in_region
                    tsv_additional.loc[index, "percent_above_threshold"] = percent_above_threshold
                    tsv_additional.loc[index, "looker"] = "S"
                    tsv_additional.loc[index, "mdl"] = "WT"
                    index = index + 1
                else:
                    # 4.2 there is not enough coverage
                    #for drug in drug_info[gene].split(","):
                    tsv_additional.loc[index, tsv_additional.columns] = ["N/A"] * len(tsv_additional.columns)
                    tsv_additional.loc[index, "Sample ID"] = sample
                    tsv_additional.loc[index, "Gene Name"] = gene
                    tsv_additional.loc[index, "rationale"] = "Insufficient Coverage"
                    #tsv_additional.loc[index, "antimicrobial"] = drug
                    tsv_additional.loc[index, "antimicrobial"] = drug_info[gene] if gene in drug_info else ""
                    tsv_additional.loc[index, "average_coverage_in_region"] = average_coverage_in_region
                    tsv_additional.loc[index, "percent_above_threshold"] = percent_above_threshold
                    tsv_additional.loc[index, "looker"] = "Insufficient Coverage"
                    tsv_additional.loc[index, "mdl"] = "Insufficient Coverage"
                    index = index + 1

    tsv_final = pandas.concat([tsv_interpretation, tsv_additional])
    # break up rows with "," separated list of anti_microbials into several rows
    tsv_final = tsv_final.assign(antimicrobial=tsv_final["antimicrobial"].str.split(',')).explode("antimicrobial")

    return [tsv_final, genes_count_dict]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="variant interpretation",
                                     prog="variant_interpretation",
                                     formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=100))
    parser.add_argument("--vcf", "-v", type=str, help="annotated vcf file", required=True)
    parser.add_argument("--bam", "-b", type=str, help="annotated vcf file", required=True)
    parser.add_argument("--bed", "-i", type=str, help="bed file with regions of interest", required=True)
    parser.add_argument("--json", "-j", type=str, help="json file with drug annotation", required=True)
    parser.add_argument("--sample_name", "-s", type=str, dest = "sample_name", help="sample name", required=True)
    parser.add_argument("--minimum_coverage", type=int, help="minimum average coverage in region (default: %(default)s)", default=0)
    parser.add_argument("--minimum_total_depth", type=int, help="minimum total number of reads at variant location (default: %(default)s)", default=0)
    parser.add_argument("--minimum_variant_depth", type=int, help="minimum number of reads that support variant (default: %(default)s)", default=0)
    parser.add_argument("--report", "-r", type=str, help="another tsv output file", required=True)
    parser.add_argument("--all_genes", "-a", action="store_true", help="output results for all genes, not only genes of interest")
    args = parser.parse_args()

    # vcf -> tsv
    vcf_df = vcf_to_pandas_dataframe(args.vcf, args.sample_name)
    
    # get drug annotation
    tsv_out = add_drug_annotation(vcf_df, args.json)
    #tsv_out.to_csv(args.output, index=False, sep="\t")

    # get intervals
    # get drug information for region
    regions = pandas.read_csv(args.bed, header=None, sep="\t")
    regions.columns = ["genome", "start", "stop", "locus", "gene", "chemical"]
    print(regions)
    regions["chemical"] = regions["chemical"].astype("str")
    get_intervals(regions)
    drug_info = dict(zip(regions.gene, regions.chemical))

    # get coverage
    coverage = calculate_average_depth(args.bam, args.bed, args.minimum_coverage)
    gene_coverage = pandas.merge(coverage, regions, on="start", how="right")
    coverage_percentage = dict(zip(gene_coverage.gene, gene_coverage.percent_above_threshold))
    coverage_average = dict(zip(gene_coverage.gene, gene_coverage.average_coverage))
    print(coverage_average)
    
    # get interpretation
    tsv_final, genes_count_dict = run_interpretation(tsv_out, drug_info, coverage_percentage, coverage_average, args.minimum_coverage, args.all_genes)
    tsv_final.to_csv(args.report, index=False, sep="\t")
