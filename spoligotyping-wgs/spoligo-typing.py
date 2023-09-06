#!/usr/bin/env python
import argparse
import numpy as np
import pandas as pd
import pathlib
from snipgenie import tools


def get_spoligotype(filename, reads_limit=500000, threshold=2, spacers="dr_spacers.fa"):
    """Get spoligotype from reads"""
    # convert reads to fasta
    fasta = pathlib.Path(filename).with_suffix(".fa")
    tools.fastq_to_fasta(filename,fasta, reads_limit)
    # make blast db from reads
    tools.make_blast_database(fasta)
    # blast spacers to db
    bl = tools.blast_fasta(fasta, spacers, evalue=0.1, maxseqs=100000, show_cmd=False)
    # filter hits
    bl = bl[(bl.qcovs > 95) & (bl.mismatch < 2)]
    # group resulting table to get hits per spacer sequence
    x = bl.groupby("qseqid").agg({"pident": np.size}).reset_index()
    # filter
    x = x[x.pident >= threshold]
    found = list(x.qseqid)
    # convert hits to binary code
    s = "".join([ str(int(item in found)) for item in list(range(1,44)) ])
    # s = []
    # for i in range(1, 44):
    #     if i in found:
    #         s.append("1")
    #     else:
    #         s.append("0")
    # s = "".join(s)
    return s


def get_sb_number(binary_str, Mbovis="Mbovis.org-database-overview-2023-08-31_0019.tsv"):
    """Get SB number from binary pattern using database reference"""
    df = pd.read_csv(Mbovis, sep="\t")
    x = df[df["Binary"] == binary_str]
    if x.empty:
        return
    else:
        return x.iloc[0].SBNumber

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="spoligo-typing", prog="spoligo-typing", formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=100))
    parser.add_argument("--fastq", "-f", type=argparse.FileType("r"), help="fastq file", required=True)
    parser.add_argument("--spacers", "-s", type=argparse.FileType("r"), help="input file with spacer sequences", required=True)
    parser.add_argument("--database", "-d", type=argparse.FileType("r"), help="tsv file of Mbovis.org database", required=True)
    args = parser.parse_args()
    b = get_spoligotype(args.fastq.name, reads_limit=500000, threshold=2, spacers=args.spacers.name)
    print(b)

    # test
    #b="1001111111111111111111111111000010110001111"

    sb_number = get_sb_number(b, Mbovis=args.database.name)
    print(sb_number)
    
