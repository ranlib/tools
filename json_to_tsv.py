#!/usr/bin/env python
"""
take json input file for wdl workflow
make a tsv out of it for import to terra.bio
"""

import os
import argparse
import pathlib
import re
import pandas

parser = argparse.ArgumentParser(prog="json_to_tsv", description="json -> tsv file", formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=80))
parser.add_argument("-j", "--json", required=True, type=pathlib.Path, help="Input json file")
parser.add_argument("-t", "--tsv", required=True, type=pathlib.Path, help="Output tsv file")
parser.add_argument("-s", "--sample_name", required=True, type=str, help="Sample name")
parser.add_argument("-n", "--dataset_name", required=True, type=str, help="Dataset name")
args = parser.parse_args()

j = pandas.read_json(args.json)
j = j.drop(columns=j.columns[j.columns.str.contains("task_")])
j.columns = [re.sub("^.*\.", "", s) for s in j.columns]
# j.to_csv("wf_varpipe.tsv",sep="\t", index=False)

jj = j.copy(deep=True)
jj.drop(1, axis=0, inplace=True)
jj.drop(2, axis=0, inplace=True)

i = 1
for read1, read2 in zip(j["read1"], j["read2"]):
    col_1 = f"lane{i}_read1"
    col_2 = f"lane{i}_read2"
    jj.loc[0, col_1] = [os.path.basename(read1)]
    jj.loc[0, col_2] = [os.path.basename(read2)]
    i += 1

jj.drop("read1", axis=1, inplace=True)
jj.drop("read2", axis=1, inplace=True)

for col in ["bed", "config", "json", "reference"]:
    if col in jj.columns:
        jj[col] = [os.path.basename(item) for item in jj[col]]

data_set_name = args.dataset_name  # "ERR552797_split_fastq"
sample_name = args.sample_name  # "ERR552797"

data_set = f"entity:{data_set_name}_id"
jj.insert(0, data_set, sample_name)

jj.to_csv(args.tsv, sep="\t", index=False)
