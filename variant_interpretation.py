#!/usr/bin/env python
"""
generate severity code for variants
"""
import argparse
import json
import pandas
from sympy import Interval
from vcf_to_pandas_dataframe import vcf_to_pandas_dataframe
from coverage import calculate_average_depth
from get_deletions_in_region import get_deletions_in_region


def get_deletion_status(row: [], has_deletions: {}) -> bool:
    """
    get deletion status
    """
    return has_deletions[row["Gene_Name"]] if row["Gene_Name"] in has_deletions else False


def get_deletion_length(row: [], regions_list: []) -> int:
    """
    get deletion length
    """
    deletion_length = -1
    for item in regions_list:
        if item[1] == row["Gene_Name"]:
            deletion_length = item[4]
    return deletion_length


def get_gene_length(row: [], regions_list: []) -> int:
    """
    get gene length
    """
    gene_length = -1
    for item in regions_list:
        if item[1] == row["Gene_Name"]:
            gene_length = item[5]
    return gene_length


def region_coverage_qc(row: [], has_deletions: {}) -> str:
    """
    region coverage QC
    """
    condition1 = row["percent_above_threshold"] == 100.0
    condition2 = has_deletions[row["Gene_Name"]] if row["Gene_Name"] in has_deletions else False
    return "PASS" if condition1 or condition2 else "FAIL"


def variant_qc(row: [], minimum_allele_percentage: float, minimum_total_depth: int, minimum_variant_depth: int) -> str:
    """
    variant QC
    """
    if row["rationale"] == "WT":
        pass_all = True
    else:
        pass_allele_percentage = row["Percent Alt Allele"] > minimum_allele_percentage
        pass_total_depth = row["Total Read Depth"] >= minimum_total_depth
        pass_variant_depth = row["Variant Read Depth"] >= minimum_variant_depth
        pass_all = pass_allele_percentage and pass_total_depth and pass_variant_depth
    return "PASS" if pass_all else "FAIL"


def get_intervals(regions: pandas.DataFrame):
    """
    generate intervals
    """
    global rpoB_cds_region
    global mmpS5, mmpL5, rplC, pepQ, atpE, Rv0678
    global Rv0678_promoter, atpE_promoter, pepQ_promoter, rplC_promoter
    global rrl_rRNA_1, rrl_rRNA_1_complement

    for _, row in regions.iterrows():
        if row["gene"] == "Rv0678":
            Rv0678 = Interval(row["start"], row["stop"])
            Rv0678_promoter = Interval(row["start"] - 84, row["stop"] - 1)

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
            rpoB_cds_region = Interval(426, 452)


gene_list_1 = ["Rv0678", "atpE", "pepQ", "mmpL5", "mmpS5", "rrl", "rplC"]
gene_list_2 = ["katG", "pncA", "ethA", "gid", "rpoB"]
gene_list_3 = ["katG", "pncA", "ethA", "gid"]

looker_1_1 = {
    "": "no WHO confidence",
    "Assoc w R": "R",
    "Assoc w R - interim": "R - interim",
    "Not assoc w R": "S",
    "Not assoc w R - Interim": "S - interim",
    "Uncertain significance": "U",
}

mdl_1_1 = {
    "": "no WHO confidence",
    "Assoc w R": "R",
    "Assoc w R - interim": "R",
    "Not assoc w R": "S",
    "Not assoc w R - Interim": "S",
    "Uncertain significance": "U",
}

looker_2_1 = {
    "": "no WHO confidence",
    "Assoc w R": "R",
    "Assoc w R - interim": "R - interim",
    "Not assoc w R": "S",
    "Not assoc w R - Interim": "S - interim",
    "Uncertain significance": "U",
}

mdl_2_1 = {
    "": "no WHO confidence",
    "Assoc w R": "R",
    "Assoc w R - interim": "R",
    "Not assoc w R": "S",
    "Not assoc w R - Interim": "S",
    "Uncertain significance": "U",
}

looker_3_1 = {
    "": "no WHO confidence",
    "Assoc w R": "R",
    "Assoc w R - interim": "R - interim",
    "Not assoc w R": "S",
    "Not assoc w R - Interim": "S - interim",
    "Uncertain significance": "U",
}

mdl_3_1 = {
    "": "no WHO confidence",
    "Assoc w R": "R",
    "Assoc w R - interim": "R",
    "Not assoc w R": "S",
    "Not assoc w R - Interim": "S",
    "Uncertain significance": "U",
}


def get_interpretation_1_2(gene: str, genomic_position: int, cds_position: int, annotation: str) -> list[str]:
    """
    get interpretation for Rv0678, atpE, pepQ, mmpL5, mmpS5, rrl, rplC
    """
    is_synonymous = annotation == "synonymous_variant"
    is_nonsynonymous = annotation != "synonymous_variant"
    looker = mdl = ""

    if gene == "Rv0678":
        if Rv0678.contains(genomic_position) & is_synonymous:
            looker = mdl = "S"
        if Rv0678.contains(genomic_position) & is_nonsynonymous:
            looker = mdl = "U"
        if Rv0678_promoter.contains(genomic_position):
            looker = mdl = "U"
        if annotation == "upstream_gene_variant" and not Rv0678_promoter.contains(genomic_position):
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
            # looker = "S"
            # mdl = "U"
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

# def get_interpretation_2_2_1(annotation: str) -> list[str]:
#     """
#     implementation of interpretation according to 2.2.1
#     """
#     is_synonymous = annotation == "synonymous_variant"
#     is_nonsynonymous = annotation != "synonymous_variant"
#     looker = mdl = ""

#     if is_synonymous:
#         looker = mdl = "S"

#     if is_nonsynonymous or annotation == "upstream_gene_variant":
#         looker = "U"
#         mdl = "S"

#     return [looker, mdl]


def get_interpretation_2_2_1(annotation: str) -> list[str]:
    """
    implementation of interpretation according to 2.2.1
    """
    is_synonymous = annotation == "synonymous_variant"
    is_nonsynonymous = annotation != "synonymous_variant"
    # ?
    effect_types = ["feature_ablation", "disruptive_inframe_insertion", "frameshift_variant", "stop_lost", "splice_region_variant", "missense_variant", "upstream_gene_variant", "disruptive_inframe_deletion", "nonsense_variant"]
    looker = mdl = ""
    if annotation in effect_types:
        looker = mdl = "U"
    else:
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
    if rpoB_cds_region.contains(cds_position):
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


def get_interpretation_3_2_1(position: int) -> list[str]:
    """
    implementation of interpretation for rrs gene according to 3.2.1
    :param int position: genomic position in NC_000962.3
    :return: list with 2 strings
    """
    if position in [1473246, 1473247, 1473329]:
        looker = mdl = "U"
    else:
        looker = "U"
        mdl = "S"
    return [looker, mdl]


# def get_interpretation_3_2_2(annotation: str) -> list[str]:
#     """
#     implementation of interpretation according to 3.2.2
#     :param str annotation: drug annotation
#     :return: list with 2 strings
#     """
#     is_synonymous = annotation == "synonymous_variant"
#     is_nonsynonymous = annotation != "synonymous_variant"
#     looker = mdl = ""
#     if is_synonymous:
#         looker = mdl = "S"
#     if is_nonsynonymous:
#         looker = "U"
#         mdl = "S"
#     return [looker, mdl]


def get_interpretation_3_2_2(annotation: str, nucleotide_change: str) -> list[str]:
    """
    implementation of interpretation according to 3.2.2
    :param str annotation: drug annotation
    :param str nucleotide_change: nucleotide change
    :return: list of 2 strings
    """
    is_synonymous = annotation == "synonymous_variant"
    is_nonsynonymous = annotation != "synonymous_variant"
    is_short_deletion = is_large_deletion = False
    if "del" in nucleotide_change:
        if len(nucleotide_change.split("del")[1]) in pandas.Interval(0, 50, closed="neither"):
            is_short_deletion = True
        else:
            is_large_deletion = True
    looker = mdl = ""
    if is_synonymous:
        looker = mdl = "S"
    if is_nonsynonymous:
        if is_large_deletion:
            looker = mdl = "U"
        else:
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
    for _, row in tsv.iterrows():
        position = int(row["POS"])
        gene_id = row["Gene_ID"]
        nucleotide_change = row["HGVS.c"]
        amino_acid_change = row["HGVS.p"]
        drug_annotation = []
        if gene_id in json_annotation:
            if nucleotide_change in list(json_annotation[gene_id].keys()):
                if int(position) in json_annotation[gene_id][nucleotide_change]["genome_positions"]:
                    for _, entry in enumerate(json_annotation[gene_id][nucleotide_change]["annotations"]):
                        pair = ["", ""]
                        if "drug" in entry:
                            pair[0] = entry["drug"]
                        if "who_confidence" in entry:
                            pair[1] = entry["who_confidence"]
                        drug_annotation.append(pair)
            elif amino_acid_change in list(json_annotation[gene_id].keys()):
                if int(position) in json_annotation[gene_id][amino_acid_change]["genome_positions"]:
                    for _, entry in enumerate(json_annotation[gene_id][amino_acid_change]["annotations"]):
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
            # remove duplicate entries from json
            drug_annotation_uniq = list(set((x[0], x[1]) for x in drug_annotation))
            for idx, item in enumerate(drug_annotation_uniq):
                row["count"] = idx + 1
                row["antimicrobial"] = item[0]
                row["confidence"] = item[1]
                row["rationale"] = "WHO"
                tsv_out.loc[i, tsv_out.columns] = row
                i = i + 1
    return tsv_out


def run_interpretation(tsv: pandas.DataFrame, samplename: str, drug_info: {}, coverage_percentage: {}, coverage_average: {}, regions_list: [], minimum_allele_percentage: int, minimum_total_depth: int, minimum_variant_depth: int, filter_genes: bool, debug: bool = False, verbose: bool = True):
    """
    interpretation of mutations
    produces severities as output: looker, mdl_prelim
    """
    has_deletions = {item[1]: item[2] for item in regions_list}

    header = ["Sample ID", "Gene_Name", "Gene_ID", "POS", "CDS.pos", "HGVS.c", "HGVS.p", "Annotation", "confidence", "antimicrobial", "Total Read Depth", "Variant Read Depth", "Percent Alt Allele", "rationale"]
    tsv_mutations = tsv.loc[:, header]

    # new columns produced here
    tsv_mutations["average_coverage_in_region"] = [0] * len(tsv_mutations.index)
    tsv_mutations["percent_above_threshold"] = [0] * len(tsv_mutations.index)
    tsv_mutations["looker"] = [""] * len(tsv_mutations.index)
    tsv_mutations["mdl_prelim"] = [""] * len(tsv_mutations.index)

    # since "Gene_Name" can be a &-separated list use sets and calculate intersection
    # keep only genes of interest
    # if intersection empty keep or skip row depending in filter_genes flag
    tsv_mutations_filtered = pandas.DataFrame(columns=list(tsv_mutations.columns))
    genes_of_interest = set(coverage_average.keys())
    number_of_entries = len(tsv_mutations.index)
    idx = 0
    for _, row in tsv_mutations.iterrows():
        genes_affected_by_variant = set(row["Gene_Name"].split("&"))
        genes_in_intersection = genes_of_interest.intersection(genes_affected_by_variant)
        if genes_in_intersection:
            tsv_mutations_filtered.loc[idx, tsv_mutations_filtered.columns] = row
            tsv_mutations_filtered.loc[idx, "Gene_Name"] = "&".join(genes_in_intersection)
            idx += 1
        else:
            if not filter_genes:
                tsv_mutations_filtered.loc[idx, tsv_mutations_filtered.columns] = row
                idx += 1
            if verbose:
                print(f"<I> variant_interpretation:  no gene if interest in row: {row}")

    tsv_mutations_filtered = tsv_mutations_filtered.assign(Gene_Name=tsv_mutations_filtered["Gene_Name"].str.split("&")).explode("Gene_Name")
    tsv_mutations = tsv_mutations_filtered.reset_index(drop=True)
    number_of_entries_after_filter = len(tsv_mutations.index)
    if verbose:
        print(f"<I> variant_interpretation: number of entries = {number_of_entries}, after gene filtering = {number_of_entries_after_filter}")

    # If no antimicrobial information, take antimicrobial information for this range (gene)
    # Note: this could be a comma separated list of chemicals
    for index, row in tsv_mutations.iterrows():
        if row["antimicrobial"] == "":
            tsv_mutations.loc[index, "antimicrobial"] = drug_info[row["Gene_Name"]] if row["Gene_Name"] in drug_info else ""

    genes_with_mutations = {}
    for index, row in tsv_mutations.iterrows():
        if row["Gene_Name"] in genes_with_mutations:  # for large deletions Gene_Name could be a list
            genes_with_mutations[row["Gene_Name"]] += 1
        else:
            genes_with_mutations[row["Gene_Name"]] = 1

        # 0. add region average coverage (depth) information
        tsv_mutations.loc[index, "average_coverage_in_region"] = coverage_average[row["Gene_Name"]] if row["Gene_Name"] in coverage_average else -1
        tsv_mutations.loc[index, "percent_above_threshold"] = coverage_percentage[row["Gene_Name"]] if row["Gene_Name"] in coverage_percentage else -1

        # 1.
        if row["Gene_Name"] in gene_list_1:
            # 1.1
            if (row["confidence"] != "") and (row["antimicrobial"] != ""):
                tsv_mutations.loc[index, "looker"] = looker_1_1[row["confidence"]]
                tsv_mutations.loc[index, "mdl_prelim"] = mdl_1_1[row["confidence"]]

            # 1.2
            if (row["confidence"] == "") and (row["antimicrobial"] != ""):
                looker, mdl = get_interpretation_1_2(row["Gene_Name"], row["POS"], row["CDS.pos"], row["Annotation"])
                tsv_mutations.loc[index, "looker"] = looker
                tsv_mutations.loc[index, "mdl_prelim"] = mdl
                tsv_mutations.loc[index, "confidence"] = "no WHO annotation"
                tsv_mutations.loc[index, "rationale"] = "expert rule 1.2"

        # 2.
        if row["Gene_Name"] in gene_list_2:
            # 2.1
            if (row["confidence"] != "") and (row["antimicrobial"] != ""):
                tsv_mutations.loc[index, "looker"] = looker_2_1[row["confidence"]]
                tsv_mutations.loc[index, "mdl_prelim"] = mdl_2_1[row["confidence"]]

            # 2.2 "WHO expert rules"
            if (row["confidence"] == "") and (row["antimicrobial"] != ""):
                # 2.2.1 gene_list_2 without rpoB
                if row["Gene_Name"] in gene_list_3:
                    looker, mdl = get_interpretation_2_2_1(row["Annotation"])
                    tsv_mutations.loc[index, "looker"] = looker
                    tsv_mutations.loc[index, "mdl_prelim"] = mdl
                    tsv_mutations.loc[index, "confidence"] = "no WHO annotation"
                    tsv_mutations.loc[index, "rationale"] = "expert rule 2.2.1"

                # 2.2.2 just rpoB
                if row["Gene_Name"] == "rpoB":
                    looker, mdl = get_interpretation_2_2_2(row["CDS.pos"], row["Annotation"])
                    tsv_mutations.loc[index, "looker"] = looker
                    tsv_mutations.loc[index, "mdl_prelim"] = mdl
                    tsv_mutations.loc[index, "confidence"] = "no WHO annotation"
                    tsv_mutations.loc[index, "rationale"] = "expert rule 2.2.2"

        # 3.
        if row["Gene_Name"] not in (gene_list_1 + gene_list_2):
            # 3.1
            if (row["confidence"] != "") and (row["antimicrobial"] != ""):
                tsv_mutations.loc[index, "looker"] = looker_3_1[row["confidence"]]
                tsv_mutations.loc[index, "mdl_prelim"] = mdl_3_1[row["confidence"]]

            # 3.2
            if (row["confidence"] == "") and (row["antimicrobial"] != ""):
                # 3.2.1, just rrs
                if row["Gene_Name"] == "rrs":
                    looker, mdl = get_interpretation_3_2_1(row["POS"])
                    tsv_mutations.loc[index, "looker"] = looker
                    tsv_mutations.loc[index, "mdl_prelim"] = mdl
                    tsv_mutations.loc[index, "confidence"] = "no WHO annotation"
                    tsv_mutations.loc[index, "rationale"] = "expert rule 3.2.1"
                else:
                    # 3.2.2
                    looker, mdl = get_interpretation_3_2_2(row["Annotation"], row["HGVS.c"])
                    tsv_mutations.loc[index, "looker"] = looker
                    tsv_mutations.loc[index, "mdl_prelim"] = mdl
                    tsv_mutations.loc[index, "confidence"] = "no WHO annotation"
                    tsv_mutations.loc[index, "rationale"] = "expert rule 3.2.2"

    # 4.
    # loop over genes of interest,
    # check if they have mutations,
    # is not in genes_with_mutations ==> no mutations, i.e. WT
    # manufacture default row for output dataframe
    # manufacture default rows for genes with no mutations and add them
    tsv_no_mutations = pandas.DataFrame(columns=list(tsv_mutations.columns))
    index = 0
    for gene in coverage_average:
        if gene not in genes_with_mutations:
            average_coverage_in_region = coverage_average[gene] if gene in coverage_average else -1
            percent_above_threshold = coverage_percentage[gene] if gene in coverage_average else -1
            has_region_coverage = percent_above_threshold == 100.0
            if has_region_coverage or has_deletions[gene]:
                # 4.1 there is 100% region coverage or deletion
                # for drug in drug_info[gene].split(","):
                # tsv_no_mutations.loc[index, tsv_no_mutations.columns] = ["N/A"] * len(tsv_no_mutations.columns)
                tsv_no_mutations.loc[index, "Sample ID"] = samplename
                tsv_no_mutations.loc[index, "Gene_Name"] = gene
                tsv_no_mutations.loc[index, "Gene_ID"] = "N/A"
                tsv_no_mutations.loc[index, "POS"] = -1
                tsv_no_mutations.loc[index, "CDS.pos"] = -1
                tsv_no_mutations.loc[index, "HGVS.c"] = "N/A"
                tsv_no_mutations.loc[index, "HGVS.p"] = "N/A"
                tsv_no_mutations.loc[index, "Annotation"] = "N/A"
                tsv_no_mutations.loc[index, "confidence"] = "N/A"
                tsv_no_mutations.loc[index, "Total Read Depth"] = -1
                tsv_no_mutations.loc[index, "Variant Read Depth"] = -1
                tsv_no_mutations.loc[index, "Percent Alt Allele"] = -1.0
                tsv_no_mutations.loc[index, "rationale"] = "WT"
                # tsv_no_mutations.loc[index, "antimicrobial"] = drug
                tsv_no_mutations.loc[index, "antimicrobial"] = drug_info[gene] if gene in drug_info else "N/A"
                tsv_no_mutations.loc[index, "average_coverage_in_region"] = average_coverage_in_region
                tsv_no_mutations.loc[index, "percent_above_threshold"] = percent_above_threshold
                tsv_no_mutations.loc[index, "looker"] = "S"
                tsv_no_mutations.loc[index, "mdl_prelim"] = "WT"
                index = index + 1
            else:
                # 4.2 there is not enough region coverage
                # for drug in drug_info[gene].split(","):
                # tsv_no_mutations.loc[index, tsv_no_mutations.columns] = ["N/A"] * len(tsv_no_mutations.columns)
                tsv_no_mutations.loc[index, "Sample ID"] = samplename
                tsv_no_mutations.loc[index, "Gene_Name"] = gene
                tsv_no_mutations.loc[index, "Gene_ID"] = "N/A"
                tsv_no_mutations.loc[index, "POS"] = -1
                tsv_no_mutations.loc[index, "CDS.pos"] = -1
                tsv_no_mutations.loc[index, "HGVS.c"] = "N/A"
                tsv_no_mutations.loc[index, "HGVS.p"] = "N/A"
                tsv_no_mutations.loc[index, "Annotation"] = "N/A"
                tsv_no_mutations.loc[index, "confidence"] = "N/A"
                tsv_no_mutations.loc[index, "Total Read Depth"] = -1
                tsv_no_mutations.loc[index, "Variant Read Depth"] = -1
                tsv_no_mutations.loc[index, "Percent Alt Allele"] = -1.0
                tsv_no_mutations.loc[index, "rationale"] = "Insufficient Coverage"
                # tsv_no_mutations.loc[index, "antimicrobial"] = drug
                tsv_no_mutations.loc[index, "antimicrobial"] = drug_info[gene] if gene in drug_info else "N/A"
                tsv_no_mutations.loc[index, "average_coverage_in_region"] = average_coverage_in_region
                tsv_no_mutations.loc[index, "percent_above_threshold"] = percent_above_threshold
                tsv_no_mutations.loc[index, "looker"] = "Insufficient Coverage"
                tsv_no_mutations.loc[index, "mdl_prelim"] = "Insufficient Coverage"
                index = index + 1

    tsv_final = pandas.concat([tsv_mutations, tsv_no_mutations])

    # break up rows with "," separated list of anti_microbials into several rows
    tsv_final = tsv_final.assign(antimicrobial=tsv_final["antimicrobial"].str.split(",")).explode("antimicrobial")
    tsv_final = tsv_final.reset_index(drop=True)

    # create new column deletion status to check whether mutation in region with deletion
    if debug:
        # tsv_final.insert(tsv_final.columns.get_loc("looker"), "Deletion_status", tsv_final.apply(lambda row: get_deletion_status(row, has_deletions), axis=1) )
        tsv_final.insert(tsv_final.columns.get_loc("looker"), "Length of deletion [bp]", tsv_final.apply(lambda row: get_deletion_length(row, regions_list), axis=1))
        tsv_final.insert(tsv_final.columns.get_loc("looker"), "Length of gene [bp]", tsv_final.apply(lambda row: get_gene_length(row, regions_list), axis=1))

    # create new column "Breadth_of_coverage_QC" based on region coverage
    tsv_final.insert(tsv_final.columns.get_loc("looker"), "Breadth_of_coverage_QC", tsv_final.apply(lambda row: region_coverage_qc(row, has_deletions), axis=1))

    # create new column "Variant_QC"
    tsv_final.insert(tsv_final.columns.get_loc("looker"), "Variant_QC", tsv_final.apply(lambda row: variant_qc(row, minimum_allele_percentage, minimum_total_depth, minimum_variant_depth), axis=1))

    # update severity based on "Variant_QC" and "Breadth_of_coverage_QC"
    tsv_final["mdl_LIMSfinal"] = tsv_final["mdl_prelim"]
    tsv_final.loc[tsv_final["Variant_QC"] == "FAIL", "mdl_LIMSfinal"] = "WT"
    tsv_final.loc[(tsv_final["Breadth_of_coverage_QC"] == "FAIL") & (tsv_final["mdl_prelim"] != "R"), "mdl_LIMSfinal"] = "Insufficient Coverage"

    return [tsv_final, genes_with_mutations]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="variant interpretation", prog="variant_interpretation", formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=100))
    parser.add_argument("--vcf", "-v", type=argparse.FileType("r"), help="annotated vcf file", required=True)
    # parser.add_argument("--filtered_vcf", type=argparse.FileType("w"), help="vcf file filtered to keep only large deletions which overlap a region of interest", required=True)
    parser.add_argument("--bam", "-b", type=argparse.FileType("r"), help="annotated vcf file", required=True)
    parser.add_argument("--bed", "-i", type=argparse.FileType("r"), help="bed file with regions of interest", required=True)
    parser.add_argument("--json", "-j", type=argparse.FileType("r"), help="json file with drug annotation", required=True)
    parser.add_argument("--samplename", "-s", type=str, dest="samplename", help="sample name", required=True)
    parser.add_argument("--minimum_coverage", type=int, help="minimum average coverage in region (default: %(default)s)", default=0)
    parser.add_argument("--minimum_total_depth", type=int, help="minimum total number of reads at variant location (default: %(default)s)", default=0)
    parser.add_argument("--minimum_variant_depth", type=int, help="minimum number of reads that support variant (default: %(default)s)", default=0)
    parser.add_argument("--minimum_allele_percentage", type=float, help="minimum variant allele percentage  (default: %(default)s)", default=0)
    parser.add_argument("--report", "-r", type=argparse.FileType("w"), help="another tsv output file", required=True)
    parser.add_argument("--filter_genes", "-f", action="store_true", help="output only genes of interest")
    parser.add_argument("--filter_variants", action="store_true", help="take only variants with PASS in vcf filter column")
    parser.add_argument("--debug", action="store_true", help="add debugging information to output")
    parser.add_argument("--verbose", action="store_true", help="turn on debugging output")
    args = parser.parse_args()

    # filter vcf file, take only large deletions which overlap regions
    # intersect_vcf_bed(args.vcf, args.bed.name, args.filtered_vcf)

    # vcf -> tsv
    vcf_df = vcf_to_pandas_dataframe(args.vcf.name, args.samplename, args.filter_variants, args.verbose)
    # vcf_df = vcf_to_pandas_dataframe(args.filtered_vcf.name, args.samplename, args.filter_variants, args.verbose)
    # vcf_df.to_csv("vcf_df.tsv",index=False,sep="\t")

    if len(vcf_df.index) == 0:
        print(f"<W> variant_interpretation: no mutations in {args.vcf}, no interpretation report.")
    else:
        # get drug annotation
        tsv_out = add_drug_annotation(vcf_df, args.json.name)
        # tsv_out.to_csv("vcf_df_drugs.tsv", index=False, sep="\t")

        # get drug information for region
        regions = pandas.read_csv(args.bed.name, header=None, sep="\t")
        regions.columns = ["genome", "start", "stop", "locus", "gene", "chemical"]
        regions["chemical"] = regions["chemical"].astype("str")
        drug_info = dict(zip(regions.gene, regions.chemical))

        # get intervals
        get_intervals(regions)

        # get coverage
        coverage = calculate_average_depth(args.bam.name, args.bed.name, args.minimum_coverage)
        gene_coverage = pandas.merge(coverage, regions, on="start", how="right")
        coverage_percentage = dict(zip(gene_coverage.gene, gene_coverage.percent_above_threshold))
        coverage_average = dict(zip(gene_coverage.gene, gene_coverage.average_coverage))
        # gene_coverage.to_csv("gene_coverage.tsv",index=False,sep="\t")

        # check which regions contain deletions
        regions_list = get_deletions_in_region(args.vcf.name, args.bed.name, args.filter_variants)

        # get interpretation
        tsv_final, genes_with_mutations = run_interpretation(tsv_out, args.samplename, drug_info, coverage_percentage, coverage_average, regions_list, args.minimum_allele_percentage, args.minimum_total_depth, args.minimum_variant_depth, args.filter_genes, args.debug, args.verbose)
        new_names = {
            "Gene_Name": "Gene Name",
            "Gene_ID": "Gene ID",
            "HGVS.c": "Nucleotide Change",
            "HGVS.p": "Amino acid Change",
            "CDS.pos": "Position within CDS",
        }
        tsv_final = tsv_final.rename(columns=new_names)
        tsv_final.to_csv(args.report, index=False)
