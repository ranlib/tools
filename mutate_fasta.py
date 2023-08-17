#!/usr/bin/env python
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
    cols = [ "final_annotation.Reference",
             "final_annotation.Position",
             "final_annotation.ReferenceNucleotide",
             "final_annotation.AlternativeNucleotide",
             "final_annotation.Gene",
             "final_annotation.Effect"]
    di = d[cols].sort_values("final_annotation.Position")
    #print(di)

    with open(args.fasta, encoding="ascii") as handle:
        record = SeqIO.read(handle, "fasta")

    cols = ["position",
            "true_ref",
            "ref",
            "alt",
            "len_ref",
            "len_alt",
            "len_seq",
            "len_new_seq",
            "agree"]
    do = pandas.DataFrame(columns = cols)
    counters = {"SNP": 0 , "MNP": 0, "Insertion": 0, "Deletion": 0}
    for index, row in di.iterrows():
        sequence = str(record.seq)
        pos = row["final_annotation.Position"]
        ref = row["final_annotation.ReferenceNucleotide"].upper()
        alt = row["final_annotation.AlternativeNucleotide"].upper()
        true_ref = sequence[pos-1]
        #print(pos, seq[pos], ref, alt)
        do.loc[index, "position"] = pos
        do.loc[index, "true_ref"] = true_ref
        do.loc[index, "ref"] = ref
        do.loc[index, "alt"] = alt
        do.loc[index, "len_ref"] = len(ref)
        do.loc[index, "len_alt"] = len(alt)
        do.loc[index, "agree"] = true_ref == list(ref)[0]
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
                
        new_record_id = record.id + "_" + str(index+1)
        new_record = SeqIO.SeqRecord(seq=Seq(new_seq), id=new_record_id, description=new_record_id)
        new_filename = pathlib.Path(args.fasta).name.replace(".fasta", "_" +str(index+1) + ".fasta")
        r = SeqIO.write(new_record, new_filename , 'fasta')
        if r!=1:
            print('Error while writing sequence:  ' + new_record.id)

            
    do.to_csv(args.output, sep="\t")
