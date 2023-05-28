#!/usr/bin/env python
""" 
The scripts accepts the stats, target coverage and final annotation file
It parses and merge those to the summary output file
"""
import sys

input1 = sys.argv[1]

print("Variant Summary:")

fh3 = open(input1, "r")
print("POS" + "\t" + "Gene Name" + "\t" + "Nucleotide Change" + "\t" + "Amino acid Change" + "\t" + "Read Depth" + "\t" + "Percent Alt Allele" + "\t" + "Annotation")

for lines in fh3:
    if lines.startswith("Sample ID"):
        continue

    lined = lines.rstrip("\r\n").split("\t")
    linedd = lined[2] + "\t" + lined[15] + "\t" + lined[9] + "\t" + lined[11] + "\t" + lined[5] + "\t" + lined[6] + "\t" + lined[7]
    
    if lined[15] == "rrl" or lined[15] == "ahpC" or lined[15] == "ahpC upstream":
        print(linedd)
    elif lined[15] == "atpE" and lined[7] == "Non-synonymous":
        print(linedd)
    elif lined[15] == "pepQ" and lined[7] == "Non-synonymous":
        print(linedd)
    elif lined[15] == "mmpR" and lined[7] == "Non-synonymous":
        print(linedd)
    elif lined[15] == "inhA" and lined[7] == "Non-synonymous":
        print(linedd)
    elif "rplC" in lined[15] and lined[7] == "Non-synonymous":
        print(linedd)
    elif lined[15] == "tlyA" and lined[7] == "Non-synonymous":
        print(linedd)
    elif lined[15] == "embB" and lined[7] == "Non-synonymous":
        if int(lined[2]) < 4246524:
            print(linedd)
        elif 4246586 < int(lined[2]) < 4248314:
            print(linedd)
        elif 4248329 < int(lined[2]) < 4249653:
            print(linedd)
        elif 4249692 < int(lined[2]):
            print(linedd)
    elif lined[15] == "gyrA" and lined[7] == "Non-synonymous":
        if "-" in lined[14]:
            if 87 < int(lined[14].split("-")[1]) < 95:
                print(linedd)
        elif 87 < int(lined[14]) < 95:
            print(linedd)
    elif lined[15] == "gyrB" and lined[7] == "Non-synonymous":
        if "-" in lined[14]:
            if 445 < int(lined[14].split("-")[1]) < 508:
                print(linedd)
        elif 445 < int(lined[14]) < 508:
            print(linedd)
    elif lined[15] == "ethA":
        if lined[7] == "Non-synonymous":
            print(linedd)
        elif "1" in lined[14] and len(lined[14]) == 1 and lined[7] == "Synonymous":
            print(linedd)
    elif lined[15] == "katG":
        if lined[7] == "Non-synonymous":
            print(linedd)
        elif "1" in lined[14] and len(lined[14]) == 1 and lined[7] == "Synonymous":
            print(linedd)
    elif lined[15] == "eis" or lined[15] == "eis upstream":
        if lined[7] == "Non-synonymous" or lined[7] == "Non-Coding":
            print(linedd)
        elif "1" in lined[14] and len(lined[14]) == 1 and lined[7] == "Synonymous":
            print(linedd)
    elif lined[15] == "pncA" or lined[15] == "pncA upstream":
        if lined[7] == "Non-synonymous" or lined[7] == "Non-Coding":
            print(linedd)
        elif "1" in lined[14] and len(lined[14]) == 1 and lined[7] == "Synonymous":
            print(linedd)
    elif lined[15] == "rrs":
        if "1401" in lined[9] or "1402" in lined[9] or "1484" in lined[9]:
            print(linedd)
    elif lined[15] == "fabG1" or lined[15] == "fabG1 upstream":
        if "c.-17" in lined[9] or "c.-15" in lined[9] or "c.-8" in lined[9] or lined[14] == "203":
            print(linedd)
    elif lined[15] == "rpoB":
        if "170" in lined[14] and lined[7] == "Non-synonymous":
            print(linedd)
        elif 453 > int(lined[14].split("-")[0]) > 425:
            print(linedd)
        elif "491" in lined[14] and lined[7] == "Non-synonymous":
            print(linedd)
        elif "-" in lined[14]:
            if int(lined[14].split("-")[0]) < 170 < int(lined[14].split("-")[1]) and lined[7] == "Non-synonymous":
                print(linedd)
            elif 453 > int(lined[14].split("-")[1]) > 426:
                print(linedd)
            elif int(lined[14].split("-")[0]) < 426 and int(lined[14].split("-")[1]) > 452:
                print(linedd)
            elif int(lined[14].split("-")[0]) < 491 < int(lined[14].split("-")[1]) and lined[7] == "Non-synonymous":
                print(linedd)

fh3.close()
