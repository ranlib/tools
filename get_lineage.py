#!/usr/bin/env python
""" 
lineage caller
based on SNPs
"""
import sys
import argparse
import vcf

parser = argparse.ArgumentParser(description="get lineage",
                                 prog="get_lineage",
                                 formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=100))
parser.add_argument("--vcf", "-v", type=argparse.FileType("r"), help="vcf file", required=True)
parser.add_argument("--lineages", "-l", type=argparse.FileType("r"), help="lineage information", required=True)
parser.add_argument("--samplename", "-s", type=str, dest = "samplename", help="sample name", required=True)
parser.add_argument("--report", "-r", type=argparse.FileType("w"), help="tsv output file", required=True)
parser.add_argument("--verbose", action="store_true", help="turn on debugging output")
args = parser.parse_args()

sublinn = ""
prevlin = []
prevsub = []

tribes = ["lineages", "Indo-Oceanic", "East-Asian", "East-African-Indian", "Euro-American", "West-Africa 1", "West-Africa 2", "Ethiopian"]

(concord , discord, concord1, discord1) = (0, 0, 0, 0)

discordance = False
sublinneage = False
animals = False

linfour = ""
hrv37 = ""
CAP = ""
BCG = ""
BOV = ""
ORYX = ""

# read in mutation informtion
(lineage, position, ref, alt) = ([], [], [], [])

vcf_reader = vcf.Reader(filename=args.vcf.name)
for record in vcf_reader:
    lineage.append(record.CHROM)
    position.append(record.POS)
    ref.append(record.REF)
    alt.append(",".join([str(n) for n in record.ALT]))

# read in lineage information
with open(args.lineages.name, "r", encoding="utf-8") as fh2:
    count = 0
    for lines in fh2:
        count += 1
        fields = lines.rstrip("\r\n").split("\t")
        if fields[2] == "931123":
            linfour = fields[2]
        if fields[2] == "1759252":
            hrv37 = fields[2]
        if fields[2] == "2831482":
            CAP = fields[2]
            animals = True
            print("SNP" + " " + fields[2] + " " + "suggests M. caprae ")
        if fields[2] == "2289073":
            BOV = fields[2]
            animals = True
            print("SNP" + " " + fields[2] + " " + "suggests M. bovis ")
        if fields[2] == "8624":
            BCG = fields[2]
            animals = True
            print("SNP" + " " + fields[2] + " " + "suggests M. bovis-BCG ")
        if fields[2] == "2726378":
            ORYX = fields[2]
            animals = True
            print("SNP" + " " + fields[2] + " " + "suggests M. orygis ")
        if fields[2] in position:
            ind = position.index(fields[2])
            if alt[ind] == fields[4]:
                if len(lineage[ind]) > 1:
                    sublin = lineage[ind]
                    prevsub.append(sublin)
                    sublinn = prevsub[0]
                    print("SNP" + " " + position[ind] + " " + "suggests sub-lineage: " + lineage[ind])
                    for item in prevsub:
                        if len(sublinn) < len(item):
                            sublinn = item
                else:
                    lin = lineage[ind]
                    prevlin.append(lin)
                    print("SNP" + " " + position[ind] + " " + "suggests lineage: " + lineage[ind])

# write output report
with open(args.report.name, "w", encoding="utf-8") as fh3:

    print("Sample ID" + "\t" + "Lineage" + "\t" + "Lineage Name", file=fh3)

    split_first = ["NA"]
    if len(prevsub) > 0:
        split_first = sublinn.split(".")
        sublinneage = True

    if len(prevlin) == 0:
        for item in prevsub:
            split_lin = item.split(".")
            if split_lin[0] != split_first[0]:
                discordance = True
            if split_lin[1] != split_first[1]:
                discordance = True
                
        if discordance:
            print("no precise lineage inferred")
            print(args.samplename + "\t" + "mixed lineage(s)" + "\t" + "mixed lineage(s)", file=fh3)
            sys.exit(1)
        else:
            if len(split_first) > 1 and animals is False:
                print("Lineage: " + split_first[0] + " :  " + tribes[int(split_first[0])])
                print(args.samplename + "\t" + split_first[0] + "\t" + tribes[int(split_first[0])], file=fh3)
            elif len(split_first) > 1 and animals is True:
                print("no precise lineage inferred")
                print(args.samplename + "\t" + "mixed lineage(s)" + "\t" + "mixed lineage(s)", file=fh3)
                sys.exit(1)
            elif len(linfour) < 2 and animals is False:
                print("Absence of SNP 931123 suggests lineage 4")
                print("Lineage: " + "4" + " :  " + "Euro-American")
                if len(hrv37) > 2:
                    print(args.samplename + "\t" + "4" + "\t" + "Euro American", file=fh3)
                elif len(hrv37) < 2:
                    print("Absence of SNP 1759252 suggests sublineage 4.9")
                    print(args.samplename + "\t" + "4" + "\t" + "Euro American", file=fh3)
            elif len(linfour) < 2 and animals is True:
                print("no precise lineage inferred")
                print(args.samplename + "\t" + "mixed lineage(s)" + "\t" + "mixed lineage(s)", file=fh3)
                sys.exit(1)
            elif len(BCG) > 0:
                print("Lineage: " + "Bovis-BCG")
                print(args.samplename + "\t" + "Bovis-BCG" + "\t" + "M.bovis-BCG", file=fh3)
            elif len(BOV) > 0:
                print("Lineage: " + "Bovis")
                print(args.samplename + "\t" + "Bovis" + "\t" + "M.bovis", file=fh3)
            elif len(ORYX) > 0:
                print("Lineage: " + "Oryx")
                print(args.samplename + "\t" + "Oryx" + "\t" + "M.orygis", file=fh3)
            elif len(CAP) > 0:
                print("Lineage: " + "Caprae")
                print(args.samplename + "\t" + "Caprae" + "\t" + "M.caprae", file=fh3)
            else:
                print("No Informative SNPs detected")
                print(args.samplename + "\t" + "No Informative SNPs detected" + "\t" + "No Informative SNPs detected", file=fh3)
    else:
        if len(prevlin) > 1:
            for j in range(0, len(prevlin)):
                if prevlin[0] != prevlin[j]:
                    discordance = True
            if discordance:
                print("no concordance between predicted lineage and sublineage(s)")
                print(args.samplename + "\t" + "mixed lineage(s)" + "\t" + "mixed lineage(s)" + "\t" + "NA", file=fh3)
                sys.exit(1)
        else:
            if len(sublinn) < 1:
                print("Lineage: " + prevlin[0] + " " + tribes[int(prevlin[0])])
                print(args.samplename + "\t" + prevlin[0] + "\t" + tribes[int(prevlin[0])], file=fh3)
            elif len(sublinn) > 1:
                for item in prevsub:
                    split_lin = item.split(".")
                    if split_lin[0] != prevlin[0]:
                        discordance = True
                    if split_lin[0] != split_first[0]:
                        discordance = True
                if discordance:
                    print("no precise lineage inferred")
                    print(args.samplename + "\t" + "mixed lineage(s)" + "\t" + "mixed lineage(s)", file=fh3)
                    sys.exit(1)
                elif animals is False:
                    print("Lineage: " + prevlin[0] + " " + tribes[int(prevlin[0])])
                    print(args.samplename + "\t" + prevlin[0] + "\t" + tribes[int(prevlin[0])], file=fh3)
                else:
                    print("no precise lineage inferred")
                    print(args.samplename + "\t" + "mixed lineage(s)" + "\t" + "mixed lineage(s)", file=fh3)
                    sys.exit(1)
