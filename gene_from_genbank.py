#!/usr/bin/env python
"""
get nucleotide sequence for gene
given genbank file
Example:
gene_from_genbank.py -g ~/Analysis/varpipe4/references/NC_000962.3/NC_000962.3.gb -f katg.fa -s katG
"""
import argparse
from Bio import SeqIO

# Set up argument parser
parser = argparse.ArgumentParser(description="get gene nucleotide sequence from genbank file")
parser.add_argument("--genbank", "-g", type=str, help="Name of input genbank file", required = True)
parser.add_argument("--fasta", "-f", type=str, help="Name of output fasta file", required  = True)
parser.add_argument("--gene", "-s", type=str, help="Name of gene", required  = True)
args = parser.parse_args()

def main():
    """
    main
    """
    with open(args.fasta, 'w', encoding="ascii") as output_fasta:
        for rec in SeqIO.parse(args.genbank, "genbank"):
            if rec.features:
                for feature in rec.features:
                    if feature.type == "CDS":
                        if 'gene' in feature.qualifiers:
                            if feature.qualifiers['gene'][0] == args.gene:
                                output_fasta.write(">%s|%s from %s\n%s\n" % (feature.qualifiers['gene'][0],feature.qualifiers['product'][0],rec.name,feature.location.extract(rec).seq))

if __name__ == "__main__":
    main()
