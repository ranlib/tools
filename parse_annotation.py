#!/usr/bin/env python
"""
The script accepts a SnpEff annotated VCF file and the sample ID name (string) as input options
it parses files and creates a final annotation file
"""
import argparse
import csv
import re
import gzip
import io

parser = argparse.ArgumentParser(description="Parse annotated vcf file")
parser.add_argument("--name", "-n", help="Name of sample", required=True)
parser.add_argument("--vcf", "-v", help="Annotated vcf file", required=True)
parser.add_argument("--loci", "-l", help="tsv file with mutation loci", required=True)
parser.add_argument("--output", "-o", help="tsv output file", required=True)
args = parser.parse_args()

position = ""
reference = ""
alternate = ""
annotation = ""
variant = ""
read_depth = ""
perc_alt = ""
nucleotide_change = ""
nuc_change = ""
transcript_pos = ""
amino_acid_change = ""
orig_aacid = ""
new_aacid = ""
codon_pos = ""
gene_name = ""
gene_id = ""
transcript = ""

(genez, genezid, start, stop, gene_anot, strand) = ([], [], [], [], [], [])
dic = {"A": "T", "T": "A", "C": "G", "G": "C"}
ref_comp = ""
alt_comp = ""

with open(args.loci, "r") as L:
    for lines in L:
        if lines.startswith("H37Rv"):
            continue
        lined = lines.rstrip("\r\n").split("\t")
        genez.append(lined[0])
        genezid.append(lined[1])
        start.append(lined[2])
        stop.append(lined[3])
        gene_anot.append(lined[4])
        strand.append(lined[5])

header = ["Sample ID", "CHROM", "POS", "REF", "ALT", "Read Depth", "Percent Alt Allele", "Annotation", "Variant Type",
          "Nucleotide Change", "Position within CDS ", "Amino acid Change", "REF Amino acid", "ALT Amino acid",
          "Codon Position", "Gene Name", "Gene ID"]
with open(args.output, "w") as output:
    with io.TextIOWrapper(io.BufferedReader(gzip.open(args.vcf))) as vcf:
        writer = csv.DictWriter(output, fieldnames=header, delimiter="\t")
        writer.writeheader()

        for lines in vcf:
            # skip comment lines
            if lines.startswith("#"):
                continue

            fields = lines.rstrip("\r\n").split("\t")
            position = fields[1]
            reference = fields[3]
            alternate = fields[4]
            #qual = fields[5] # not used
            filterf = fields[6]
            info = fields[7]
            rformat = fields[8].split(":")  # format
            rarr = fields[9].split(":")  # format numbers
            res = dict(zip(rformat, rarr))

            # if not PASS info field not necessarily ; separated list of X=Y elements
            if filterf not in [ "PASS" ]:
                continue
            
            read_depth = res["DP"]
            if "AD" in res:
                (num_ref, num_alt) = res["AD"].split(",")
                perc_alt = round(float(num_alt) * 100.0 / float(read_depth), 2)
            else:
                perc_alt = 1000

            info_list = info.split(";")
            annot_dict = {}
            for element in info_list:
                element_list = element.split("=")
                if len(element_list) == 2:
                    annot_dict[ element_list[0] ] = element_list[1]
                    
            subannot = annot_dict["ANN"].split(",")
            smallannot = subannot[0].split("|")
            if smallannot[2] == "MODIFIER":
                for x in range(0, 3):
                    if (int(start[x]) - 1) < int(position) < (int(stop[x]) + 1):
                        annotation = gene_anot[x]
                        if genez[x] == "rrs":
                            nuc_change = str((int(position)) - (int(start[x]) - 1))
                            gene_id = "MTB000019"
                            nucleotide_change = "c." + nuc_change + reference + ">" + alternate
                        elif genez[x] == "rrl":
                            nuc_change = str((int(position)) - (int(start[x]) - 1))
                            gene_id = "MTB000020"
                            nucleotide_change = "c." + nuc_change + reference + ">" + alternate
                        elif genez[x] == "crfA":
                            nuc_change = str((int(position)) - (int(start[x]) - 1))
                            gene_id = "crfA"
                            nucleotide_change = "c." + nuc_change + reference + ">" + alternate
                        elif strand[x] == "forward":
                            gene_id = genezid[x]
                            nuc_change = str((int(position)) - (int(stop[x]) + 1))
                            nucleotide_change = "c." + nuc_change + reference + ">" + alternate
                        elif strand[x] == "reverse":
                            ref_comp = ""
                            alt_comp = ""
                            for char in reference:
                                ref_comp += dic[char]
                            for char in alternate:
                                alt_comp += dic[char]
                            gene_id = genezid[x]
                            nuc_change = str((int(start[x]) - 1) - int(position))
                            nucleotide_change = "c." + nuc_change + ref_comp + ">" + alt_comp
                        gene_name = genez[x]
                        amino_acid_change = "NA"
                        if len(fields[4]) > len(fields[3]):
                            if strand[x] == "forward":
                                nucleotide_change = "c." + nuc_change + "_" + str(
                                    int(nuc_change) + 1) + "ins" + alternate[len(reference):]
                                if genez[x] == "crfA":
                                    transcript_pos = nuc_change + "-" + str(int(nuc_change) + 1)
                            elif strand[x] == "reverse":
                                nucleotide_change = "c." + str(
                                    int(nuc_change) - 1) + "_" + nuc_change + "ins" + alt_comp[len(reference):][::-1]
                            variant = "Insertion"
                        elif len(fields[3]) > len(fields[4]):
                            if strand[x] == "forward":
                                if len(reference) - len(alternate) == 1:
                                    nucleotide_change = "c." + str(
                                        int(nuc_change) + len(alternate)) + "del" + reference[len(alternate):]
                                    if genez[x] == "crfA":
                                        transcript_pos = str(int(nuc_change) + len(alternate))
                                else:
                                    nucleotide_change = "c." + str(int(nuc_change) + len(alternate)) + "_" + str(
                                        int(nuc_change) + len(reference) - 1) + "del" + reference[len(alternate):]
                                    if genez[x] == "crfA":
                                        transcript_pos = str(int(nuc_change) + len(alternate)) + "-" + str(
                                            int(nuc_change) + len(reference) - 1)
                            elif strand[x] == "reverse":
                                if len(reference) - len(alternate) == 1:
                                    nucleotide_change = "c." + str(
                                        int(nuc_change) - len(reference) - 1) + "del" + ref_comp[len(alternate):][::-1]
                                else:
                                    nucleotide_change = "c." + str(int(nuc_change) - len(reference) - 1) + "_" + str(
                                        int(nuc_change) - len(alternate)) + "del" + ref_comp[len(alternate):][::-1]
                            variant = "Deletion"
                        elif len(fields[3]) > 1 and len(fields[3]) == len(fields[4]):
                            if strand[x] == "forward":
                                nucleotide_change = "c." + nuc_change + "_" + str(
                                    int(nuc_change) + 1) + "del" + reference + "ins" + alternate
                            elif strand[x] == "reverse":
                                nucleotide_change = "c." + str(
                                    int(nuc_change) - 1) + "_" + nuc_change + "del" + ref_comp[::-1] + "ins" + alt_comp[
                                                                                                               ::-1]
                            variant = "MNP"
                        else:
                            variant = "SNP"
                        transcript = "NA"
                        if genez[x] == "crfA":
                            transcript_pos = nuc_change
                        else:
                            transcript_pos = "NA"
                        orig_aacid = "NA"
                        new_aacid = "NA"
                        codon_pos = "NA"
                        break
                    else:
                        (nucposarray, newsubannot) = ([], [])
                        indx = len(subannot) - 1
                        for y in range(0, indx):
                            nucpos = subannot[y].split("|")[14]
                            if "downstream" not in subannot[y].split("|")[1] and len(nucpos) > 0:
                                nucposarray.append(int(nucpos))
                                newsubannot.append(subannot[y])
                        if len(newsubannot) > 0:
                            minnuc = str(min((nucposarray)))
                            ind = nucposarray.index(int(minnuc))
                            workannot = newsubannot[ind]
                            newsmallannot = workannot.split("|")
                            nucleotide_change = newsmallannot[9]
                            gene_name = newsmallannot[3] + " upstream"
                            gene_id = newsmallannot[4] + " upstream"
                        else:
                            nucleotide_change = smallannot[9]
                            gene_name = smallannot[3] + " downstream"
                            gene_id = smallannot[4] + " downstream"
                        annotation = "Non-Coding"
                        if len(fields[3]) > len(fields[4]):
                            variant = "Deletion"
                        elif len(fields[4]) > len(fields[3]):
                            variant = "Insertion"
                        elif len(fields[3]) > 1 and len(fields[3]) == len(fields[4]):
                            variant = "MNP"
                        else:
                            variant = "SNP"
                        amino_acid_change = "NA"
                        transcript = "NA"
                        transcript_pos = "NA"
                        orig_aacid = "NA"
                        new_aacid = "NA"
                        codon_pos = "NA"
            else:
                if len(smallannot[10]) < 13:
                    if smallannot[10][2:5] == smallannot[10][-3:]:
                        annotation = "Synonymous"
                    else:
                        annotation = "Non-synonymous"
                elif len(smallannot[10]) > 13:
                    pos11 = re.findall(r"\d+", smallannot[10])
                    ind11 = smallannot[10].index(pos11[0])
                    ind11len = len(pos11[0])
                    translen = ind11 + ind11len
                    if smallannot[10][2:ind11] == smallannot[10][translen:]:
                        annotation = "Synonymous"
                    else:
                        annotation = "Non-synonymous"
                nucleotide_change = smallannot[9]
                if len(smallannot[10]) < 1:
                    amino_acid_change = "NA"
                else:
                    amino_acid_change = smallannot[10]
                gene_name = smallannot[3]
                gene_id = smallannot[4]
                transcript = smallannot[6]
                if gene_name == "erm_37_":
                    gene_name = "erm(37)"
                if len(fields[3]) == len(fields[4]) and len(fields[3]) > 1:
                    variant = "MNP"
                elif "del" in nucleotide_change or "del" in amino_acid_change:
                    variant = "Deletion"
                elif "ins" in nucleotide_change or "ins" in amino_acid_change:
                    variant = "Insertion"
                elif "dup" in nucleotide_change or "dup" in amino_acid_change:
                    variant = "Insertion"
                else:
                    variant = "SNP"
                if variant in ("Insertion", "Deletion"):
                    new_aacid = "NA"
                    if "_" in smallannot[9]:
                        array1 = smallannot[9].split("_")
                        po1 = array1[0].split(".")
                        pos1 = po1[1]
                        pos2 = re.findall(r"\d+", array1[1])[0]
                        transcript_pos = pos1 + "-" + pos2
                    else:
                        transcript_pos = re.findall(r"\d+", smallannot[9])[0]
                    if "_" in smallannot[10]:
                        array2 = smallannot[10].split("_")
                        po11 = array2[0].split(".")
                        orig_aacid = po11[1][0:3]
                        pos11 = po11[1][3:]
                        pos12 = re.findall(r"\d+", array2[1])[0]
                        codon_pos = pos11 + "-" + pos12

                    else:
                        if len(smallannot[10]) > 0:
                            codon_pos = re.findall(r"\d+", smallannot[10])[0]
                            orig_aacid = smallannot[10][2:5]
                        else:
                            codon_pos = "NA"
                            orig_aacid = "NA"
                else:
                    if len(smallannot[10]) < 13:
                        orig_aacid = smallannot[10][2:5]
                    elif len(smallannot[10]) > 13:
                        pos11 = re.findall(r"\d+", smallannot[10])
                        ind11 = smallannot[10].index(pos11[0])
                        orig_aacid = smallannot[10][2:ind11]
                    if "*" in smallannot[10] or "?" in smallannot[10]:
                        new_aacid = "NA"
                    else:
                        if len(smallannot[10]) < 13:
                            new_aacid = smallannot[10][-3:]
                        elif len(smallannot[10]) > 13:
                            pos11 = re.findall(r"\d+", smallannot[10])
                            ind11 = smallannot[10].index(pos11[0])
                            transpos = len(pos11[0]) + ind11
                            new_aacid = smallannot[10][transpos:]
                    transcript_pos = re.findall(r"\d+", smallannot[9])[0]
                    if len(smallannot[10]) < 13:
                        codon_pos = re.findall(r"\d+", smallannot[10])[0]
                    elif len(smallannot[10]) > 13:
                        string_pos = re.findall(r"\d+", smallannot[10])[0]
                        string_pos_2 = int(string_pos) + 1
                        codon_pos = string_pos + "-" + str(string_pos_2)

                for x in range(0, 3):
                    if (int(start[x]) - 1) < int(position) < (int(stop[x]) + 1):
                        annotation = gene_anot[x]
                        if strand[x] == "forward":
                            gene_id = genezid[x]
                            nuc_change = str((int(position)) - (int(stop[x]) + 1))
                            nucleotide_change = "c." + nuc_change + reference + ">" + alternate
                        elif strand[x] == "reverse":
                            ref_comp = ""
                            alt_comp = ""
                            for char in reference:
                                ref_comp += dic[char]
                            for char in alternate:
                                alt_comp += dic[char]
                            gene_id = genezid[x]
                            nuc_change = str((int(start[x]) - 1) - int(position))
                            nucleotide_change = "c." + nuc_change + ref_comp + ">" + alt_comp
                        gene_name = genez[x]
                        amino_acid_change = "NA"
                        if len(fields[4]) > len(fields[3]):
                            if strand[x] == "forward":
                                nucleotide_change = "c." + nuc_change + "_" + str(
                                    int(nuc_change) + 1) + "ins" + alternate[len(reference):]
                            elif strand[x] == "reverse":
                                nucleotide_change = "c." + str(
                                    int(nuc_change) - 1) + "_" + nuc_change + "ins" + alt_comp[len(reference):][::-1]
                            variant = "Insertion"
                        elif len(fields[3]) > len(fields[4]):
                            if strand[x] == "forward":
                                if len(reference) - len(alternate) == 1:
                                    nucleotide_change = "c." + str(
                                        int(nuc_change) + len(alternate)) + "del" + reference[len(alternate):]
                                else:
                                    nucleotide_change = "c." + str(int(nuc_change) + len(alternate)) + "_" + str(
                                        int(nuc_change) + len(reference) - 1) + "del" + reference[len(alternate):]
                            elif strand[x] == "reverse":
                                if len(reference) - len(alternate) == 1:
                                    nucleotide_change = "c." + str(
                                        int(nuc_change) - len(reference) - 1) + "del" + ref_comp[len(alternate):][::-1]
                                else:
                                    nucleotide_change = "c." + str(int(nuc_change) - len(reference) - 1) + "_" + str(
                                        int(nuc_change) - len(alternate)) + "del" + ref_comp[len(alternate):][::-1]
                            variant = "Deletion"
                        elif len(fields[3]) > 1 and len(fields[3]) == len(fields[4]):
                            if strand[x] == "forward":
                                nucleotide_change = "c." + nuc_change + "_" + str(
                                    int(nuc_change) + 1) + "del" + reference + "ins" + alternate
                            elif strand[x] == "reverse":
                                nucleotide_change = "c." + str(
                                    int(nuc_change) - 1) + "_" + nuc_change + "del" + ref_comp[::-1] + "ins" + alt_comp[
                                                                                                               ::-1]
                            variant = "MNP"
                        else:
                            variant = "SNP"
                        transcript = "NA"
                        transcript_pos = "NA"
                        orig_aacid = "NA"
                        new_aacid = "NA"
                        codon_pos = "NA"
                        break

            record = [args.name, fields[0], position, reference, alternate, read_depth, perc_alt, annotation, variant,
                      nucleotide_change, transcript_pos, amino_acid_change, orig_aacid, new_aacid, codon_pos, gene_name,
                      gene_id]
            record_dict = dict(zip(header, record))
            writer.writerow(record_dict)
