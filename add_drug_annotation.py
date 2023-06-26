#!/usr/bin/env python

import pandas
import argparse
import json
from sympy import Interval

def get_intervals(regions: pandas.DataFrame):
    
    for index, row in regions.iterrows():
        if row["gene"] == "mmpR":
            mmpR = Interval(row["Start"], row["Stop"])        
            mmpR_promoter = Interval(row["Start"] - 84, row["Stop"] - 1)    

        if row["gene"] == "atpE":
            atpE = Interval(row["Start"], row["Stop"])        
            atpE_promoter = Interval(row["Start"] - 48, row["Stop"] - 1)

        if row["gene"] == "pepQ":
            pepQ = Interval(row["Start"], row["Stop"])        
            pepQ_promoter = Interval(row["Start"] - 33, row["Stop"] - 1)

        if row["gene"] == "rplC":
            replC = Interval(row["Start"], row["Stop"])        
            replC_promoter = Interval(row["Start"] - 18, row["Stop"] - 1)

        if row["gene"] == "mmpL5":
            mmpL5 = Interval(row["Start"], row["Stop"])        

        if row["gene"] == "mmpS5":
            mmpS5 = Interval(row["Start"], row["Stop"])        

        if row["gene"] == "rrl":
            rrl = Interval(row["Start"], row["Stop"])        
            rrl_rRNA = Interval(1, row["Stop"] - row["Start"])
            # position in CDS
            rrl_rRNA_1 = Interval(2003, 2367).union( Interval(2449, 3056) )
            #not in [2003, 2367] and not in [2449, 3056]
            rrl_rRNA_1_complement = rrl_rRNA_1.complement(rrl_rRNA)

        if row["gene"] == "rpoB":
            # position in CDS
            rpoB_codon = Interval(426, 452)

gene_list_1 = [ "mmpR", "atpE", "pepQ", "mmpL5", "mmpS5", "rrl", "rplC" ]
gene_list_2 = [ "katG", "pncA", "ethA", "gid", "rpoB" ]
gene_list_3 = [ "katG", "pncA", "ethA", "gid" ]

looker_1_1 = {
    "": "no WHO confidence",
    "Assoc w R": "R",
    "Assoc w R - interim": "R - interim",
    "Not assoc w R": "S",
    "Not assoc w R - Interim": "S - interim",
    "Uncertain significance": "U"
}

mdl_1_1 = {
    "": "no WHO confidence",
    "Assoc w R": "R",
    "Assoc w R - interim": "U",
    "Not assoc w R": "S",
    "Not assoc w R - Interim": "U",
    "Uncertain significance": "U"
}

looker_2_2 = {
    "": "no WHO confidence",
    "Assoc w R": "R",
    "Assoc w R - interim": "R - interim",
    "Not assoc w R": "S",
    "Not assoc w R - Interim": "S - interim",
    "Uncertain significance": "U"
}

mdl_2_2 = {
    "": "no WHO confidence",
    "Assoc w R": "R",
    "Assoc w R - interim": "U",
    "Not assoc w R": "S",
    "Not assoc w R - Interim": "S",
    "Uncertain significance": "S"
}

looker_3_1 = {
    "": "no WHO confidence",
    "Assoc w R": "R",
    "Assoc w R - interim": "R - interim",
    "Not assoc w R": "S",
    "Not assoc w R - Interim": "S - interim",
    "Uncertain significance": "U"
}

mdl_3_1 = {
    "": "no WHO confidence",
    "Assoc w R": "R",
    "Assoc w R - interim": "U",
    "Not assoc w R": "S",
    "Not assoc w R - Interim": "S",
    "Uncertain significance": "S"
}


def get_interpretation_1_2(gene: str, genomic_position: int, cds_position: int, annotation: str) -> list[str]:
    """
    get interpretation for mmpR, atpE, pepQ, mmpL5, mmpS5, rrl, rplC
    """
    is_synonymous = (annotation == "synonymous_variant")
    is_nonsynonymous = (annotation != "synonymous_variant")
    looker = mds = ""
    
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
            
    if gene == "mmpR":
        if mmpR.contains(genomic_position) & is_synonymous:
            looker = mdl = "S"
        if mmpR.contains(genomic_position) & is_nonsynonymous:
            looker = mdl = "U"
        if mmpR_promoter.contains(genomic_position):
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
        if rrl_rRNA_1.contains(cds_position) & is_synonymous:
            looker = mdl = "U"
        if rrl_rRNA_1_complement.contains(cds_position):
            looker = "S"
            mdl = "U"

    return [looker, mds]

def get_interpretation_2_2_1(annotation: str) -> list[str]:
    is_synonymous = (annotation == "synonymous_variant")
    is_nonsynonymous = (annotation != "synonymous_variant")
    # ?
    effect_types = [ "disruptive_inframe_insertion",
                     "frameshift_variant",
                     "stop_lost",
                     "splice_region_variant",
                     "missense_variant",
                     "upstream_gene_variant",
                     "disruptive_inframe_deletion",
                     "nonsense_variant" ]
    looker = mdl = ""
    if annotation in effect_types:
        looker = mdl = "S"
    else:
        if is_synonymous:
            looker = mdl = "S"
        if is_nonsynonymous:
            looker = "U"
            mdl = "S"
    return [looker, mdl]

def get_interpretation_2_2_2(cds_position: int, annotation: str) -> list[str]:
    is_synonymous = (annotation == "synonymous_variant")
    is_nonsynonymous = (annotation != "synonymous_variant")
    looker = mdl = ""
    if rpoB_codon.contains(cds_position ):
        if is_synonymous:
            looker = mdl = "S"
        if is_nonsynonymous:
            looker = "R"
            mdl = "R"
    else:
        if is_synonymous:
            looker = mdl = "S"
        if is_nonsynonymous:
            looker = "U"
            mdl = "S"
            
    return [looker, mdl]

def get_interpretation_3_2(annotation: str) -> list[str]:
    is_synonymous = (annotation == "synonymous_variant")
    is_nonsynonymous = (annotation != "synonymous_variant")
    looker = mdl = ""
    if is_synonymous:
        looker = mdl = "S"
    if is_nonsynonymous:
        looker = "U"
        mdl = "S"
    return [looker, mdl]


def get_coverage_for_gene(gene: str) -> int:
    regions = pandas.read_csv(args.bed, sep="\t", header=None)
    regions.columns = ["genome", "Start", "Stop", "locus", "gene", "chemical"]
    coverage = pandas.read_csv(args.coverage)
    coverage["Start"] = coverage["Start"] - 1
    gene_coverage = pandas.merge(coverage, regions, on="Start", how="right")
    filter = gene_coverage.gene == gene
    coverage_percentage = gene_coverage.loc[filter, "%_above_10"].values.tolist()
    coverage_average = gene_coverage.loc[filter, "average_coverage"].values.tolist()
    
    if len( coverage_percentage ) > 0:
        percentage_is_ok = coverage_percentage[0] > 0
    else:
        percentage_is_ok = False
        
    if len( coverage_average ) > 0:
        average_is_ok = coverage_average[0] > 0
    else:
        average_is_ok = False
        
    if len( coverage_average ) > 0:
        average_coverage_value = coverage_average[0]
    else:
        average_coverage_value = 0
        
    return average_coverage_value

def get_drug_information_for_gene(gene: str) -> list[str]:
    drug_info = pandas.read_csv(args.bed, sep="\t", header=None)
    drug_info.columns = ["genome", "Start", "Stop", "locus", "gene", "chemical"]
    filter = drug_info.gene == gene
    chemicals = drug_info.loc[filter, "chemical"]
    drugs = chemicals.values.tolist()
    if len(drugs) > 0:
        drugs = drugs[0].split(",")
    else:
        drugs = []
    return drugs
    
def add_drug_annotation(input: str, annotation: str):
    tsv = pandas.read_csv(input, sep="\t")

    with open(annotation,"r") as jsonfile:
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
    tsv_out = pandas.DataFrame(columns = list(tsv.columns))

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
                if int(position) in json_annotation[gene_id][nucleotide_change]['genome_positions']:
                    for index, entry in enumerate(json_annotation[gene_id][nucleotide_change]['annotations']):
                        pair = ['', '']
                        if 'drug' in entry:
                            pair[0] = entry['drug']
                        if 'who_confidence' in entry:
                            pair[1] = entry['who_confidence']
                        drug_annotation.append(pair)
            elif amino_acid_change in list(json_annotation[gene_id].keys()):
                if int(position) in json_annotation[gene_id][amino_acid_change]['genome_positions']:
                    for index, entry in enumerate(json_annotation[gene_id][amino_acid_change]['annotations']):
                        pair = ['', '']
                        if 'drug' in entry:
                            pair[0] = entry['drug']
                        if 'who_confidence' in entry:
                            pair[1] = entry['who_confidence']
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

    
def run_interpretation(tsv: pandas.DataFrame):
    #
    # interpretation
    #
    header = [ "Sample ID","Gene Name", "Gene ID", "POS", "Position within CDS", "Nucleotide Change", "Amino acid Change", "Annotation", "confidence", "antimicrobial", "Read Depth", "Percent Alt Allele", "rationale"]
    tsv_interpretation = tsv.loc[ :, header ]
    tsv_interpretation["depth"] = [0]*len(tsv_interpretation.index)
    
    genes_count_dict = {}
    
    for index, row in tsv_interpretation.iterrows():
        if row["Gene Name"] in genes_count_dict:
            genes_count_dict[row["Gene Name"]] = genes_count_dict[row["Gene Name"]] + 1
        else:
            genes_count_dict[row["Gene Name"]] = 1

        # 0. add region average coverage (depth) information
        tsv_interpretation.loc[index, "depth"] = get_coverage_for_gene(row["Gene Name"])
            
        # 1.
        if row["Gene Name"] in gene_list_1:
            # 1.1
            if (row["confidence"] != "") & (row["antimicrobial"] != ""):
                tsv_interpretation.loc[index, "looker"] = looker_1_1[ row["confidence"] ]
                tsv_interpretation.loc[index, "mdl"] = mdl_1_1[ row["confidence"] ]

            # 1.2
            if (row["confidence"] == "") & (row["antimicrobial"] != ""):
                looker_mdl = get_interpretation_1_2(row["Gene Name"], row["POS"], row["Position within CDS "], row["Annotation"])
                tsv_interpretation.loc[index, "looker"] = looker_mdl[0]
                tsv_interpretation.loc[index, "mdl"] = looker_mdl[1]
                row["confidence"] = "no WHO annotation"

        # 2.
        if row["Gene Name"] in gene_list_2:
            # 2.1
            if (row["confidence"] != "") & (row["antimicrobial"] != ""):
                tsv_interpretation.loc[index, "looker"] = looker_2_2[ row["confidence"] ]
                tsv_interpretation.loc[index, "mdl"] = mdl_2_2[ row["confidence"] ]
            
            # 2.2 "WHO expert rules"
            if (row["confidence"] == "") & (row["antimicrobial"] != ""):
                # 2.2.1 gene_list_2 without rpoB
                if row["Gene Name"] in gene_list_3:
                    looker_mdl = get_interpretation_2_2_1(row["Annotation"])
                    tsv_interpretation.loc[index, "looker"] = looker_mdl[0]
                    tsv_interpretation.loc[index, "mdl"] = looker_mdl[1]
                    row["confidence"] = "no WHO annotation"

                # 2.2.2 just rpoB
                if row["Gene Name"] == "rpoB":
                    looker_mdl = get_interpretation_2_2_2(row["Annotation"])
                    tsv_interpretation.loc[index, "looker"] = looker_mdl[0]
                    tsv_interpretation.loc[index, "mdl"] = looker_mdl[1]
                    row["confidence"] = "no WHO annotation"
                    
                
        # 3.
        if row["Gene Name"] not in (gene_list_1 + gene_list_2):
            # 3.1
            if (row["confidence"] != "") & (row["antimicrobial"] != ""):
                tsv_interpretation.loc[index, "looker"] = looker_3_1[ row["confidence"] ]
                tsv_interpretation.loc[index, "mdl"] = mdl_3_1[ row["confidence"] ]

            # 3.2
            if (row["confidence"] == "") & (row["antimicrobial"] != ""):
                looker_mdl = get_interpretation_3_2(row["Annotation"])
                tsv_interpretation.loc[index, "looker"] = looker_mdl[0]
                tsv_interpretation.loc[index, "mdl"] = looker_mdl[1]
                row["confidence"] = "no WHO annotation"

    # 4.
    # loop over genes of interest,
    # check if they have mutations,
    # is not in genes_count_dict ==> no mutations, i.e. WT
    # manufacture default row for output dataframe

    # keep only genes of interest
    tsv_interpretation = tsv_interpretation[ tsv_interpretation["Gene Name"].isin(gene_list_1 + gene_list_2) ]

    # fill additional rows
    tsv_additional = pandas.DataFrame(columns = list(tsv_interpretation.columns))
    index = 0
    for sample in tsv_interpretation["Sample ID"].unique().tolist():
        for gene in (gene_list_1 + gene_list_2):
            if gene not in genes_count_dict:
                average_coverage_in_region = get_coverage_for_gene(gene)
                if average_coverage_in_region > args.lower_coverage:
                    # 4.1 there is coverage
                    for drug in get_drug_information_for_gene(gene):
                        tsv_additional.loc[index, tsv_additional.columns] = ["N/A"] * len(tsv_additional.columns)
                        tsv_additional.loc[index, "Sample ID"] = sample
                        tsv_additional.loc[index, "Gene Name"] = gene
                        tsv_additional.loc[index, "rationale"] = "WT"
                        tsv_additional.loc[index, "antimicrobial"] = drug
                        tsv_additional.loc[index, "depth"] = average_coverage_in_region
                        tsv_additional.loc[index, "looker"] = "S"
                        tsv_additional.loc[index, "mdl"] = "S"
                        index = index + 1
                else:
                    # 4.2 there is not enough coverage
                    for drug in get_drug_information_for_gene(gene):
                        tsv_additional.loc[index, tsv_additional.columns] = ["N/A"] * len(tsv_additional.columns)
                        tsv_additional.loc[index, "Sample ID"] = sample
                        tsv_additional.loc[index, "Gene Name"] = gene
                        tsv_additional.loc[index, "rationale"] = "Insufficient Coverage"
                        tsv_additional.loc[index, "antimicrobial"] = drug
                        tsv_additional.loc[index, "depth"] = average_coverage_in_region
                        tsv_additional.loc[index, "looker"] = "U" # ?
                        tsv_additional.loc[index, "mdl"] = "U" # ?
                        index = index + 1
                
    tsv_final = pandas.concat([tsv_interpretation, tsv_additional])

    return [tsv_final, genes_count_dict]


    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="add drug annotation to tsv form of annotated vcf file")
    parser.add_argument("--tsv", "-t", type=str, help="tsv version of annotated vcf file", required=True)
    parser.add_argument("--json", "-j", type=str, help="json file with drug annotation", required=True)
    parser.add_argument("--bed", "-b", type=str, help="bed file with regions of interest", required=True)
    parser.add_argument("--coverage", "-c", type=str, help="csv file with coverage for regions of interest", required=True)
    parser.add_argument("--lower_coverage", "-l", type=int, help="lower coverage cutoff", required=True)
    parser.add_argument("--output", "-o", type=str, help="tsv output file", required=True)
    parser.add_argument("--report", "-r", type=str, help="another tsv output file", required=True)
    args = parser.parse_args()

    tsv_out = add_drug_annotation(args.tsv, args.json)
    tsv_out.to_csv(args.output, index=False, sep="\t")

    regions = pandas.read_csv(args.bed, header=None, sep="\t")
    regions.columns = ["genome", "Start", "Stop", "locus", "gene", "chemical"]
    get_intervals(regions)
    
    tsv_final, genes_count_dict = run_interpretation(tsv_out)
    tsv_final.to_csv(args.report, index=False, sep="\t")
    #print(genes_count_dict)    
