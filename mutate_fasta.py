#!/usr/bin/env python
import os
import argparse
import pathlib
import pandas
from Bio import SeqIO
from Bio.Seq import Seq

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="mutate_fasta",
                                     prog="mutate_fasta",
                                     formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=100))
    parser.add_argument("--fasta", "-f", type=str, help="fasta file", required=True)
    parser.add_argument("--excel", "-e", type=str, help="excel file with list of mutations", required=True)
    parser.add_argument("--output", "-o", type=str, help="output: compare true reference with reference nucleotide in input", required=True)
    args = parser.parse_args()

    d = pandas.read_excel(args.excel)
    cols_input = [ "final_annotation.Reference",
                   "final_annotation.Position",
                   "final_annotation.ReferenceNucleotide",
                   "final_annotation.AlternativeNucleotide",
                   "final_annotation.Gene",
                   "final_annotation.Effect"]
    di = d[cols_input].sort_values("final_annotation.Position")
    #print(di)

    with open(args.fasta, encoding="ascii") as handle:
        genome = SeqIO.read(handle, "fasta")

    cols_output = ["position",
                   "true_ref",
                   "ref",
                   "alt",
                   "len_ref",
                   "len_alt",
                   "len_seq",
                   "len_new_seq",
                   "agree"]
    do = pandas.DataFrame(columns = cols_output)
    
    counters = {"SNP": 0 , "MNP": 0, "Insertion": 0, "Deletion": 0}
    for index, row in di.iterrows():
        sequence = str(genome.seq)
        pos = row["final_annotation.Position"]
        ref = row["final_annotation.ReferenceNucleotide"].upper()
        alt = row["final_annotation.AlternativeNucleotide"].upper()

        if len(ref) > 1 and len(alt) > 1:
            print(pos,ref,alt)
            
            # remove common suffix
            p = [ref[::-1], alt[::-1]]
            commonprefix = os.path.commonprefix(p)
            [new_ref, new_alt] = [ x[len(commonprefix):] for x in p ]
            ref = new_ref[::-1]
            alt = new_alt[::-1]
            print(pos,ref,alt)
                    
            # remove common prefix
            p = [ref, alt]
            commonprefix = os.path.commonprefix(p)
            [new_ref, new_alt] = [ x[len(commonprefix)-1:] for x in p ]
            ref = new_ref
            alt = new_alt
            pos = pos + len(commonprefix) - 1
            print(pos,ref,alt)
            
        #true_ref = sequence[pos-1]
        true_ref = sequence[pos-1: pos+len(ref)-1]
        #print(pos, seq[pos], ref, alt)
        do.loc[index, "position"] = pos
        do.loc[index, "true_ref"] = true_ref
        #do.loc[index, "ref"] = ref
        #do.loc[index, "alt"] = alt
        do.loc[index, "ref"] = ref
        do.loc[index, "alt"] = alt
        do.loc[index, "len_ref"] = len(ref)
        do.loc[index, "len_alt"] = len(alt)
        #do.loc[index, "agree"] = true_ref == list(ref)[0]
        do.loc[index, "agree"] = true_ref == ref
        if len(ref) == len(alt):
            if len(ref) == 1:
                counters["SNP"] += 1
                new_seq = sequence[:pos-1] + alt + sequence[pos:]
            else:
                counters["MNP"] += 1
                new_seq = sequence[:pos-1] + alt + sequence[pos+len(alt)-1:]
        elif len(ref) < len(alt):
            counters["Insertion"] += 1
            new_seq = sequence[:pos-1] + alt + sequence[pos+len(ref)-1:]
        elif len(ref) > len(alt):
            counters["Deletion"] += 1
            new_seq = sequence[:pos-1] + alt + sequence[pos+len(ref)-1:]

        do.loc[index, "len_seq"] = len(sequence)
        do.loc[index, "len_new_seq"] = len(new_seq)
                
        new_genome_id = genome.id + "_" + str(index+1)
        new_genome = SeqIO.SeqRecord(seq=Seq(new_seq), id=new_genome_id, description=new_genome_id)
        new_filename = pathlib.Path(args.fasta).name.replace(".fasta", "_" +str(index+1) + ".fasta")
        do.loc[index, "filename"] = new_filename
        
        #r = SeqIO.write(new_genome, new_filename , 'fasta')
        #if r!=1:
        #    print('Error while writing sequence:  ' + new_genome.id)

    # temprarily remove some columns
    cols_to_drop = ["true_ref", "len_ref", "len_alt", "len_seq", "len_new_seq", "agree"]
    do.drop(cols_to_drop, axis=1, inplace=True)
    
    do.to_csv(args.output, sep="\t", index=False)
