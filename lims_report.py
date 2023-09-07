#!/usr/bin/env python
"""
lab report ==> lims report
"""
import argparse
import datetime
import pandas

# interval of RRDR region for rpoB
RRDR = pandas.Interval(left=426, right=452, closed="both")

def get_chemical_evaluation(severity: str, antimicrobial: str, gene: str, annotation: str, amino_acid_change: str, position_within_cds: int) -> str:
    """
    get chemical evalution based on severity and other variables
    :param str severity
    :param str antimicrobial
    :param str gene
    :param str annotation
    :param str amino_acid_change
    :return chemical_evaluation
    """
    chemical_evaluation = ""
    
    if severity in ["R", "R-interim"]:
        if gene == "rpoB":
            if amino_acid_change in ["Leu430Pro", "Asp435Tyr", "His445Asn", "His445Ser", "His445Leu", "His445Cys", "Leu452Pro", "Ile491Phe"]:
                chemical_evaluation = "Predicted low-level resistance to rifampin; may test susceptible by phenotypic methods"
            else:
                chemical_evaluation = "Predicted resistance to rifampin"
        else:
            chemical_evaluation = f"Mutation(s) associated with resistance to {antimicrobial} detected"

    if severity == "U":
        chemical_evaluation = f"The detected mutation(s) have uncertain significance; resistance to {antimicrobial} cannot be ruled out"

    if severity in ["S", "S-interim"]:
        if antimicrobial != "rifampin":
            chemical_evaluation = f"No mutations associated with resistance to {antimicrobial} detected"
        else:  # rifampin
            if (annotation == "synonymous_variant") and (position_within_cds in RRDR):
                chemical_evaluation = "Predicted susceptibility to rifampin. The detected synonymous mutation(s) do not confer resistance"
            else:
                chemical_evaluation = "Predicted susceptibility to rifampin"

    if severity == "WT":
        chemical_evaluation = f"No mutations associated with resistance to {antimicrobial} detected"

    if severity == "Insufficient Coverage":
        chemical_evaluation = "Pending Retest"

    return chemical_evaluation


def get_gene_chemical_evaluation(chem_gene: pandas.DataFrame, gene: str) -> str:
    """
    get evaluation of chemical/gene combination
    """
    severities = chem_gene["mdl"].to_list()
    nt_changes = chem_gene["Nucleotide_Change"].to_list()
    aa_changes = chem_gene["Amino_acid_Change"].to_list()
    aa_changes_with_parens = [f"({x})" for x in aa_changes]
    annotations = chem_gene["Annotation"].to_list()
    positions_within_cds = chem_gene["Position_within_CDS"].to_list()

    rpoB_flag = False
    chem_gene_eval = []
    mut_counter = {"R": 0, "U": 0, "S": 0}
    for i, severity in enumerate(severities):
        if severity in ["R", "R-interim"]:
            chem_gene_eval.append(nt_changes[i] + " " + aa_changes_with_parens[i])
            mut_counter["R"] += 1

        if severity == "U":
            chem_gene_eval.append(nt_changes[i] + " " + aa_changes_with_parens[i])
            mut_counter["U"] += 1

        if severity in ["S", "S-interim"]:
            if (gene == "rpoB") and (annotations[i] == "synonymous_variant") and (positions_within_cds[i] in RRDR):
                chem_gene_eval.append(nt_changes[i] + " " + aa_changes_with_parens[i] + "[synonymous]")
            mut_counter["S"] += 1

        if severity == "WT":
            chem_gene_eval.append("No mutations detected")

        if severity == "Insufficient Coverage":
            chem_gene_eval.append("No sequence")

    if (mut_counter["R"] == 0) and (mut_counter["U"] == 0) and (mut_counter["S"] > 0):
        chem_gene_eval.append("No high confidence mutations detected")

    return "; ".join(chem_gene_eval)


def lims_report(lab_tsv: str, bed: str, lineage_name: str, operator: str, verbose: bool) -> pandas.DataFrame():
    """
    create lims report
    """
    # header of lims report
    header = [
        "MDL sample accession numbers",
        "M_DST_A01_ID",
        "M_DST_B01_INH",
        "M_DST_B02_katG",
        "M_DST_B03_fabG1",
        "M_DST_B04_inhA",
        "M_DST_C01_ETO",
        "M_DST_C02_ethA",
        "M_DST_C03_fabG1",
        "M_DST_C04_inhA",
        "M_DST_D01_RIF",
        "M_DST_D02_rpoB",
        "M_DST_E01_PZA",
        "M_DST_E02_pncA",
        "M_DST_F01_EMB",
        "M_DST_F02_embA",
        "M_DST_F03_embB",
        "M_DST_G01_AMK",
        "M_DST_G02_rrs",
        "M_DST_G03_eis",
        "M_DST_H01_KAN",
        "M_DST_H02_rrs",
        "M_DST_H03_eis",
        "M_DST_I01_CAP",
        "M_DST_I02_rrs",
        "M_DST_I03_tlyA",
        "M_DST_J01_MFX",
        "M_DST_J02_gyrA",
        "M_DST_J03_gyrB",
        "M_DST_K01_LFX",
        "M_DST_K02_gyrA",
        "M_DST_K03_gyrB",
        "M_DST_L01_BDQ",
        "M_DST_L02_Rv0678",
        "M_DST_L03_atpE",
        "M_DST_L04_pepQ",
        "M_DST_L05_mmpL5",
        "M_DST_L06_mmpS5",
        "M_DST_M01_CFZ",
        "M_DST_M02_Rv0678",
        "M_DST_M03_pepQ",
        "M_DST_M04_mmpL5",
        "M_DST_M05_mmpS5",
        "M_DST_N01_LZD",
        "M_DST_N02_rrl",
        "M_DST_N03_rplC",
        "Analysis date",
        "Operator",
    ]

    # map input colum for chemical to output lims report column
    chemical_to_column = {
        "isoniazid": "M_DST_B01_INH",
        "ethionamide": "M_DST_C01_ETO",
        "rifampin": "M_DST_D01_RIF",
        "pyrazinamide": "M_DST_E01_PZA",
        "ethambutol": "M_DST_F01_EMB",
        "amikacin": "M_DST_G01_AMK",
        "kanamycin": "M_DST_H01_KAN",
        "capreomycin": "M_DST_I01_CAP",
        "moxifloxacin": "M_DST_J01_MFX",
        "levofloxacin": "M_DST_K01_LFX",
        "bedaquiline": "M_DST_L01_BDQ",
        "clofazimine": "M_DST_M01_CFZ",
        "linezolid": "M_DST_N01_LZD",
    }

    # map input columns for chemical/gene combination to output lims report columns
    gene_chemical_to_column = {
        ("katG", "isoniazid"): "M_DST_B02_katG",
        ("fabG1", "isoniazid"): "M_DST_B03_fabG1",
        ("inhA", "isoniazid"): "M_DST_B04_inhA",
        ("ethA", "ethionamide"): "M_DST_C02_ethA",
        ("fabG1", "ethionamide"): "M_DST_C03_fabG1",
        ("inhA", "ethionamide"): "M_DST_C04_inhA",
        ("rpoB", "rifampin"): "M_DST_D02_rpoB",
        ("pncA", "pyrazinamide"): "M_DST_E02_pncA",
        ("embA", "ethambutol"): "M_DST_F02_embA",
        ("embB", "ethambutol"): "M_DST_F03_embB",
        ("rrs", "amikacin"): "M_DST_G02_rrs",
        ("eis", "amikacin"): "M_DST_G03_eis",
        ("rrs", "kanamycin"): "M_DST_H02_rrs",
        ("eis", "kanamycin"): "M_DST_H03_eis",
        ("rrs", "capreomycin"): "M_DST_I02_rrs",
        ("tlyA", "capreomycin"): "M_DST_I03_tlyA",
        ("gyrA", "moxifloxacin"): "M_DST_J02_gyrA",
        ("gyrB", "moxifloxacin"): "M_DST_J03_gyrB",
        ("gyrA", "levofloxacin"): "M_DST_K02_gyrA",
        ("gyrB", "levofloxacin"): "M_DST_K03_gyrB",
        ("Rv0678", "bedaquiline"): "M_DST_L02_Rv0678",
        ("atpE", "bedaquiline"): "M_DST_L03_atpE",
        ("pepQ", "bedaquiline"): "M_DST_L04_pepQ",
        ("mmpL5", "bedaquiline"): "M_DST_L05_mmpL5",
        ("mmpS5", "bedaquiline"): "M_DST_L06_mmpS5",
        ("Rv0678", "clofazimine"): "M_DST_M02_Rv0678",
        ("pepQ", "clofazimine"): "M_DST_M03_pepQ",
        ("mmpL5", "clofazimine"): "M_DST_M04_mmpL5",
        ("mmpS5", "clofazimine"): "M_DST_M05_mmpS5",
        ("rrl", "linezolid"): "M_DST_N02_rrl",
        ("rplC", "linezolid"): "M_DST_N03_rplC",
    }

    # input = lab report
    lab_cols = ["Sample ID", "Position within CDS", "Nucleotide Change", "Amino acid Change", "Annotation", "Gene Name", "antimicrobial", "Warning", "mdl"]
    lab = pandas.read_csv(lab_tsv, sep="\t", usecols=lab_cols)
    lab.columns = lab.columns.str.replace(" ", "_")
    
    # consider only variants which do not have a failed QC remark in Warning column of input lab report
    n_before = len(lab.index))
    # remove NAs
    lab = lab.fillna("")
    lab = lab[ ~lab["Warning"].str.contains("Mutation failed QC") ]
    n_after = len(lab.index)
    if verbose:
        print(f"<I> lims_report: remove variants with failed QC: number of variants before =  {n_before}, after = {n_after}")
        
    # mdl = category, determine sort order according to severity
    # sort by severity, keep only most severe
    lab["mdl"] = pandas.Categorical(lab["mdl"], ["R", "Insufficient Coverage", "R-interim", "U", "S-interim", "S", "WT"])
    lab_sorted = lab.sort_values(["antimicrobial", "mdl"], ascending=True)

    # output = lims report
    lims = pandas.DataFrame(columns=header)
    lims.loc[0, "MDL sample accession numbers"] = lab["Sample_ID"][0]
    lims.loc[0, "M_DST_A01_ID"] = lineage_name
    lims.loc[0, "Analysis date"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    lims.loc[0, "Operator"] = operator

    # get chemicals and genes of interest from bed file
    chemicals = pandas.read_csv(bed, sep="\t", header=None, usecols=[4, 5])
    chemicals.columns = ["gene", "chemicals"]
    antimicrobials = chemicals.assign(chemicals=chemicals["chemicals"].str.split(",")).explode("chemicals")
    chemicals = set(antimicrobials["chemicals"].to_list())
    genes = set(antimicrobials["gene"].to_list())

    # loop over chemicals and fill lims report
    for chemical in chemicals:  # chemical from bed file
        # print(chemical)
        if chemical not in chemical_to_column:  # hardcoded, make sure consistent with chemical from bed file
            if verbose:
                print(f"<W> lims_report: {chemical} not in list of chemicals, continuing.")
        else:
            filter_chem = "antimicrobial==@chemical"
            chem = lab_sorted.query(filter_chem)
            chem = chem.reset_index()
            if len(chem.index) == 0:
                if verbose:
                    print("<W> lims_report: no data for chemical {chemical}.")
            else:
                if verbose:
                    print(chem)
                # fill chemical header based on mutation with highest severity
                gene = chem.iloc[0]["Gene_Name"]
                severity = chem.iloc[0]["mdl"]
                annotation = chem.iloc[0]["Annotation"]
                amino_acid_change = chem.iloc[0]["Amino_acid_Change"]
                position_within_cds = chem.iloc[0]["Position_within_CDS"]
                lims.loc[0, chemical_to_column[chemical]] = get_chemical_evaluation(severity, chemical, gene, annotation, amino_acid_change, position_within_cds)

                for gene in genes:
                    key = (gene, chemical)
                    if key not in gene_chemical_to_column:
                        if verbose:
                            print(f"<W> lims_report: ({gene}, {chemical}) not in map.")
                    else:
                        # get information for gene/chemical information and fill column
                        filter_gene = "antimicrobial==@chemical & Gene_Name==@gene"
                        chem_gene = lab_sorted.query(filter_gene)
                        if verbose:
                            print(chem_gene)
                        lims.loc[0, gene_chemical_to_column[key]] = get_gene_chemical_evaluation(chem_gene, gene)

    return lims


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="lims report", prog="lims_report", formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=80))
    parser.add_argument("--lab", "-l", type=argparse.FileType('r'), dest="lab", help="input lab tsv report", required=True)
    parser.add_argument("--operator", "-o", type=str, dest="operator", help="input operator name", required=True)
    parser.add_argument("--lineages", "-s", type=argparse.FileType('r'), dest="lineages", help="input lineage tsv file", nargs="?")
    parser.add_argument("--bed", "-b", type=argparse.FileType('r'), dest="bed", help="input bed file", required=True)
    parser.add_argument("--lims", "-r", type=argparse.FileType('w'), dest="lims", help="output lims tsv report", required=True)
    parser.add_argument("--verbose", "-v", action="store_true", help="turn on debugging output")
    args = parser.parse_args()

    # get lineage from lineages file
    if args.lineages is not None:
        lineages = pandas.read_csv(args.lineages, sep="\t")
        lineage = lineages.loc[0, "Lineage Name"]
    else:
        lineage = "Lineage information not available"
        
    # create lims report
    df = lims_report(args.lab, args.bed, lineage, args.operator, args.verbose)
    df.to_csv(args.lims, sep="\t", index=False)
