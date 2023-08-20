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

    if gene == "mmpR5":
        if mmpR5.contains(genomic_position) & is_synonymous:
            looker = mdl = "S"
        if mmpR5.contains(genomic_position) & is_nonsynonymous:
            looker = mdl = "U"
        if mmpR5_promoter.contains(genomic_position):
            looker = mdl = "U"
        if annotation == "upstream_gene_variant" and not mmpR5_promoter.contains(genomic_position):
            looker = "U"
            mdl = "S"

    if gene == "atpE":
        if atpE.contains(genomic_position) & is_synonymous:
            looker = mdl = "S"
        if atpE.contains(genomic_position) & is_nonsynonymous:
            looker = mdl = "U"
        if atpE_promoter.contains(genomic_position):
            looker = mdl = "U"
        if annotation == "upstream_gene_variant" and not atpE_promoter.contains(genomic_position):
            looker = "U"
            mdl = "S"

    if gene == "pepQ":
        if pepQ.contains(genomic_position) & is_synonymous:
            looker = mdl = "S"
        if pepQ.contains(genomic_position) & is_nonsynonymous:
            looker = mdl = "U"
        if pepQ_promoter.contains(genomic_position):
            looker = mdl = "U"
        if annotation == "upstream_gene_variant" and not pepQ_promoter.contains(genomic_position):
            looker = "U"
            mdl = "S"

    if gene == "rplC":
        if rplC.contains(genomic_position) & is_synonymous:
            looker = mdl = "S"
        if rplC.contains(genomic_position) & is_nonsynonymous:
            looker = mdl = "U"
        if rplC_promoter.contains(genomic_position):
            looker = mdl = "U"
        if annotation == "upstream_gene_variant" and not rplC_promoter.contains(genomic_position):
            looker = "U"
            mdl = "S"

    if gene == "mmpL5":
        if mmpL5.contains(genomic_position) & is_synonymous:
            looker = mdl = "S"
        if mmpL5.contains(genomic_position) & is_nonsynonymous:
            looker = mdl = "U"
        if annotation == "upstream_gene_variant":
            looker = "U"
            mdl = "S"

    if gene == "mmpS5":
        if mmpS5.contains(genomic_position) & is_synonymous:
            looker = mdl = "S"
        if mmpS5.contains(genomic_position) & is_nonsynonymous:
            looker = mdl = "U"
        if annotation == "upstream_gene_variant":
            looker = "U"
            mdl = "S"

    if gene == "rrl":
        if rrl_rRNA_1.contains(cds_position):
            looker = mdl = "U"
        if rrl_rRNA_1_complement.contains(cds_position):
            #looker = "S"
            #mdl = "U"
            looker = "U"
            mdl = "S"

    return [looker, mdl]


# def get_interpretation_2_2_1(annotation: str) -> list[str]:
#     """
#     implementation of interpretation according to 2.2.1
#     """
#     is_synonymous = annotation == "synonymous_variant"
#     is_nonsynonymous = annotation != "synonymous_variant"
#     # ?
#     effect_types = ["disruptive_inframe_insertion", "frameshift_variant", "stop_lost", "splice_region_variant", "missense_variant", "upstream_gene_variant", "disruptive_inframe_deletion", "nonsense_variant"]
#     looker = mdl = ""
#     if annotation in effect_types:
#         looker = mdl = "U"
#     else:
#         if is_synonymous:
#             looker = mdl = "S"
#         if is_nonsynonymous:
#             looker = "U"
#             mdl = "S"
#     return [looker, mdl]

def get_interpretation_2_2_1(annotation: str) -> list[str]:
    """
    implementation of interpretation according to 2.2.1
    """
    is_synonymous = annotation == "synonymous_variant"
    is_nonsynonymous = annotation != "synonymous_variant"
    looker = mdl = ""

    if is_synonymous:
        looker = mdl = "S"

    if is_nonsynonymous or annotation == "upstream_gene_variant":
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
        if is_nonsynonymous or annotation == "upstream_gene_variant":
            looker = "U"
            mdl = "S"

    return [looker, mdl]


def get_interpretation_3_2(annotation: str) -> list[str]:
    """
    implementation of interpretation according to 3.2
    :param str annotation: filename of json file wit drug annotation
    :return: list with 2 strings
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
    :param str annotation: filename of json file wit drug annotation
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


def run_interpretation(tsv: pandas.DataFrame, drug_info: {}, coverage_percentage: {}, coverage_average: {}, minimum_coverage: int, filter_genes: bool, verbose: bool):
    """
    interpretation
    """
    header = ["Sample ID", "Gene Name", "Gene ID", "POS", "Position within CDS", "Nucleotide Change", "Amino acid Change", "Annotation", "confidence", "antimicrobial", "Total Read Depth", "Variant Read Depth", "Percent Alt Allele", "rationale"]
    tsv_mutations = tsv.loc[:, header]
    tsv_mutations["average_coverage_in_region"] = [0] * len(tsv_mutations.index)
    tsv_mutations["percent_above_threshold"] = [0] * len(tsv_mutations.index)
    tsv_mutations["Percent Alt Allele"] = tsv_mutations["Percent Alt Allele"]
    tsv_mutations["Comment"] = [""] * len(tsv_mutations.index)

    # If no antimicrobial information, take antimicrobial information for this range (gene)
    # Note: this could be a comma separated list of chemicals
    for index, row in tsv_mutations.iterrows():
        if row["antimicrobial"] == "":
            tsv_mutations.loc[index, "antimicrobial"] = drug_info[row["Gene Name"]] if row["Gene Name"] in drug_info else ""

    genes_with_mutations = {}
    for index, row in tsv_mutations.iterrows():
        if row["Gene Name"] in genes_with_mutations:
            genes_with_mutations[row["Gene Name"]] += 1
        else:
            genes_with_mutations[row["Gene Name"]] = 1

        # 0. add region average coverage (depth) information
        tsv_mutations.loc[index, "average_coverage_in_region"] = coverage_average[row["Gene Name"]] if row["Gene Name"] in coverage_average else -1
        tsv_mutations.loc[index, "percent_above_threshold"] = coverage_percentage[row["Gene Name"]] if row["Gene Name"] in coverage_percentage else -1

        # 1.
        if row["Gene Name"] in gene_list_1:
            # 1.1
            if (row["confidence"] != "") & (row["antimicrobial"] != ""):
                tsv_mutations.loc[index, "looker"] = looker_1_1[row["confidence"]]
                tsv_mutations.loc[index, "mdl"] = mdl_1_1[row["confidence"]]

            # 1.2
            if (row["confidence"] == "") & (row["antimicrobial"] != ""):
                looker, mdl = get_interpretation_1_2(row["Gene Name"], row["POS"], row["Position within CDS"], row["Annotation"])
                tsv_mutations.loc[index, "looker"] = looker
                tsv_mutations.loc[index, "mdl"] = mdl
                tsv_mutations.loc[index, "confidence"] = "no WHO annotation"
                tsv_mutations.loc[index, "rationale"] = "expert rule 1.2"

        # 2.
        if row["Gene Name"] in gene_list_2:
            # 2.1
            if (row["confidence"] != "") & (row["antimicrobial"] != ""):
                tsv_mutations.loc[index, "looker"] = looker_2_2[row["confidence"]]
                tsv_mutations.loc[index, "mdl"] = mdl_2_2[row["confidence"]]

            # 2.2 "WHO expert rules"
            if (row["confidence"] == "") & (row["antimicrobial"] != ""):
                # 2.2.1 gene_list_2 without rpoB
                if row["Gene Name"] in gene_list_3:
                    looker, mdl = get_interpretation_2_2_1(row["Annotation"])
                    tsv_mutations.loc[index, "looker"] = looker
                    tsv_mutations.loc[index, "mdl"] = mdl
                    tsv_mutations.loc[index, "confidence"] = "no WHO annotation"
                    tsv_mutations.loc[index, "rationale"] = "expert rule 2.2.1"

                # 2.2.2 just rpoB
                if row["Gene Name"] == "rpoB":
                    looker, mdl = get_interpretation_2_2_2(row["Position within CDS"], row["Annotation"])
                    tsv_mutations.loc[index, "looker"] = looker
                    tsv_mutations.loc[index, "mdl"] = mdl
                    tsv_mutations.loc[index, "confidence"] = "no WHO annotation"
                    tsv_mutations.loc[index, "rationale"] = "expert rule 2.2.2"

        # 3.
        if row["Gene Name"] not in (gene_list_1 + gene_list_2):
            # 3.1
            if (row["confidence"] != "") & (row["antimicrobial"] != ""):
                tsv_mutations.loc[index, "looker"] = looker_3_1[row["confidence"]]
                tsv_mutations.loc[index, "mdl"] = mdl_3_1[row["confidence"]]

            # 3.2
            if (row["confidence"] == "") & (row["antimicrobial"] != ""):
                looker, mdl = get_interpretation_3_2(row["Annotation"])
                tsv_mutations.loc[index, "looker"] = looker
                tsv_mutations.loc[index, "mdl"] = mdl
                tsv_mutations.loc[index, "confidence"] = "no WHO annotation"
                tsv_mutations.loc[index, "rationale"] = "expert rule 3.2"

    # 4.
    # loop over genes of interest,
    # check if they have mutations,
    # is not in genes_with_mutations ==> no mutations, i.e. WT
    # manufacture default row for output dataframe

    # keep only genes of interest
    #print(coverage_average.keys())
    if filter_genes:
        #tsv_mutations = tsv_mutations[tsv_mutations["Gene Name"].isin(gene_list_1 + gene_list_2)]
        number_of_entries = len(tsv_mutations.index)
        tsv_mutations = tsv_mutations[tsv_mutations["Gene Name"].isin(coverage_average.keys())]
        number_of_entries_after_filter = len(tsv_mutations.index)
        if verbose:
            print(f"<I> variant_interpretation: number of entries = {number_of_entries}, after filtering = {number_of_entries_after_filter}")

    # manufacture default rows for genes with no mutations and add them
    tsv_no_mutations = pandas.DataFrame(columns=list(tsv_mutations.columns))
    index = 0
    for sample in tsv_mutations["Sample ID"].unique().tolist():
        #for gene in gene_list_1 + gene_list_2:
        for gene in coverage_average:
            if gene not in genes_with_mutations:
                average_coverage_in_region = coverage_average[gene] if gene in coverage_average else -1
                percent_above_threshold = coverage_percentage[gene] if gene in coverage_average else -1
                if average_coverage_in_region > minimum_coverage:
                    # 4.1 there is coverage
                    #for drug in drug_info[gene].split(","):
                    tsv_no_mutations.loc[index, tsv_no_mutations.columns] = ["N/A"] * len(tsv_no_mutations.columns)
                    tsv_no_mutations.loc[index, "Sample ID"] = sample
                    tsv_no_mutations.loc[index, "Gene Name"] = gene
                    tsv_no_mutations.loc[index, "rationale"] = "WT"
                    #tsv_no_mutations.loc[index, "antimicrobial"] = drug
                    tsv_no_mutations.loc[index, "antimicrobial"] =  drug_info[gene] if gene in drug_info else ""
                    tsv_no_mutations.loc[index, "average_coverage_in_region"] = average_coverage_in_region
                    tsv_no_mutations.loc[index, "percent_above_threshold"] = percent_above_threshold
                    tsv_no_mutations.loc[index, "looker"] = "S"
                    tsv_no_mutations.loc[index, "mdl"] = "WT"
                    index = index + 1
                else:
                    # 4.2 there is not enough coverage
                    #for drug in drug_info[gene].split(","):
                    tsv_no_mutations.loc[index, tsv_no_mutations.columns] = ["N/A"] * len(tsv_no_mutations.columns)
                    tsv_no_mutations.loc[index, "Sample ID"] = sample
                    tsv_no_mutations.loc[index, "Gene Name"] = gene
                    tsv_no_mutations.loc[index, "rationale"] = "Insufficient Coverage"
                    #tsv_no_mutations.loc[index, "antimicrobial"] = drug
                    tsv_no_mutations.loc[index, "antimicrobial"] = drug_info[gene] if gene in drug_info else ""
                    tsv_no_mutations.loc[index, "average_coverage_in_region"] = average_coverage_in_region
                    tsv_no_mutations.loc[index, "percent_above_threshold"] = percent_above_threshold
                    tsv_no_mutations.loc[index, "looker"] = "Insufficient Coverage"
                    tsv_no_mutations.loc[index, "mdl"] = "Insufficient Coverage"
                    index = index + 1

    tsv_final = pandas.concat([tsv_mutations, tsv_no_mutations])
    # break up rows with "," separated list of anti_microbials into several rows
    tsv_final = tsv_final.assign(antimicrobial=tsv_final["antimicrobial"].str.split(',')).explode("antimicrobial")

    return [tsv_final, genes_with_mutations]


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
    parser.add_argument("--filter_genes", "-f", action="store_true", help="output only only genes of interest")
    parser.add_argument("--verbose", action="store_true", help="turn on debugging output")
    args = parser.parse_args()

    # vcf -> tsv
    vcf_df = vcf_to_pandas_dataframe(args.vcf, args.sample_name)
    if len(vcf_df.index) == 0:
        print(f"<W> variant_interpreation: no mutations in {args.vcf}, no interpretation report.")
    else:
        # get drug annotation
        tsv_out = add_drug_annotation(vcf_df, args.json)
        #tsv_out.to_csv(args.output, index=False, sep="\t")

        # get intervals
        # get drug information for region
        regions = pandas.read_csv(args.bed, header=None, sep="\t")
        regions.columns = ["genome", "start", "stop", "locus", "gene", "chemical"]
        #print(regions)
        regions["chemical"] = regions["chemical"].astype("str")
        get_intervals(regions)
        drug_info = dict(zip(regions.gene, regions.chemical))
    
        # get coverage
        coverage = calculate_average_depth(args.bam, args.bed, args.minimum_coverage)
        gene_coverage = pandas.merge(coverage, regions, on="start", how="right")
        coverage_percentage = dict(zip(gene_coverage.gene, gene_coverage.percent_above_threshold))
        coverage_average = dict(zip(gene_coverage.gene, gene_coverage.average_coverage))
        #print(coverage_average)
        #gene_coverage.to_csv("gene_coverage.tsv",index=False,sep="\t")
    
        # get interpretation
        tsv_final, genes_with_mutations = run_interpretation(tsv_out, drug_info, coverage_percentage, coverage_average, args.minimum_coverage, args.filter_genes, args.verbose)
        tsv_final.to_csv(args.report, index=False, sep="\t")
