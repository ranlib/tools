#!/usr/bin/env python
import argparse
import numpy as np
import pandas as pd
from snipgenie import tools

def get_spoligotype(filename, reads_limit=500000, threshold=2):
    """Get spoligotype from reads"""

    ref = 'dr_spacers.fa'
    #convert reads to fasta
    tools.fastq_to_fasta(filename, 'temp.fa', reads_limit)
    #make blast db from reads
    tools.make_blast_database('temp.fa')
    #blast spacers to db
    bl = tools.blast_fasta('temp.fa', ref, evalue=0.1,
                           maxseqs=100000, show_cmd=False)
    #filter hits
    bl=bl[(bl.qcovs>95) & (bl.mismatch<2)]
    #group resulting table to get hits per spacer sequence
    x = bl.groupby('qseqid').agg({'pident':np.size}).reset_index()
    #filter
    x = x[x.pident>=threshold]
    found = list(x.qseqid)
    s=[]
    #convert hits to binary code
    for i in range(1,44):
        if i in found:
            s.append('1')
        else:
            s.append('0')
    s =''.join(s)  
    return s

def get_sb_number(binary_str):
    """Get SB number from binary pattern using database reference"""

    df = pd.read_csv('Mbovis.org_db.csv')
    x = df[df['binary'] == binary_str]
    if len(x) == 0:
        return
    else:
        return x.iloc[0].SB



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="spoligo-typing",
                                     prog="spoligo-typing",
                                     formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=100))
    parser.add_argument("--fastq", "-f", type=argparse.FileType("r"), help="fastq file", required=True)
    parser.add_argument("--report", "-r", type=argparse.FileType("w"), help="another tsv output file", required=True)
    parser.add_argument("--verbose", action="store_true", help="turn on debugging output")
    args = parser.parse_args()
    b = get_spoligotype(args.fastq.name)
    print(b)
    print(get_sb_number(b))
    

    
