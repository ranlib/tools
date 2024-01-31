#!/usr/bin/env python
"""
turn input tbdb json into tsv file
"""
import argparse
import json
import pandas

def tbdb_json_to_tsv(annotation: str) -> pandas.DataFrame:
    """
    turn input tbdb tsv into tsv file
    """
    with open(annotation, "r", encoding="utf-8") as jsonfile:
        json_annotation = json.load(jsonfile)

    tsv_out = pandas.DataFrame(columns=["gene_id", "mutation", "drug", "who_confidence", "start", "stop"])

    i = 0 
    for gene_id in json_annotation:
        for mutation in json_annotation[gene_id]:
            positions = json_annotation[gene_id][mutation]["genome_positions"]
            if not isinstance(positions, list):
                #print(gene_id, mutation, -1, -1)
                # skip this entry for now, just 10 or so
                continue
        
            annotations = json_annotation[gene_id][mutation]["annotations"]
            if not isinstance(annotations, list):
                #print(annotations)
                # skip this entry for now
                continue
            
            for annotation in annotations:
                # take only annotations which have drug information
                if "drug" in annotation and annotation["type"] == "who_confidence":
                    tsv_out.loc[i, "gene_id"] = gene_id
                    tsv_out.loc[i, "mutation"] = mutation
                    tsv_out.loc[i, "start"] = positions[0]
                    tsv_out.loc[i, "stop"] = positions[-1]
                    tsv_out.loc[i, "drug"] = annotation["drug"]
                    tsv_out.loc[i, "who_confidence"] = annotation["who_confidence"]
                i = i + 1

    return tsv_out

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="variant interpretation", prog="variant_interpretation", formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=100))
    parser.add_argument("--json", "-j", type=argparse.FileType("r"), help="json file with drug annotation", required=True)
    parser.add_argument("--tsv", "-t", type=argparse.FileType("w"), help="tsv output file", required=True)
    parser.add_argument("--verbose", action="store_true", help="turn on debugging output")
    args = parser.parse_args()
    tsv = tbdb_json_to_tsv(args.json.name)
    tsv.to_csv(args.tsv.name, sep="\t", index=False)
