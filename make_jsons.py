#!/usr/bin/env python
"""
Purpose:
Make varpipe.wdl jason input file for list of fastq file

Input: 1) ABSOLUTE directory path to fastq files 2) json file template

Output: for each sample json file with sample name

Note: sample name created from fastq file name with the assumption
that fastq file name has the structure samplename_XXXXX.fastq.gz

"""
import sys
import os
import json
import re
import pathlib
import argparse

parser = argparse.ArgumentParser(prog="make_json", description="make input json files for varpipe for fastq files in input directory")
parser.add_argument("--dir", "-d", required=True, type=str, help="Input directory with fastq files")
parser.add_argument("--json_template", "-t", required=True, type=str, help="Output json file")
parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
args = parser.parse_args()

files = []
for filepath in pathlib.Path(args.dir).glob("**/*.fastq.gz"):
    files.append(str(filepath.absolute()))

if args.verbose:
    print(files)

with open(args.json_template, "r", encoding="ascii") as stream:
    try:
        json_template = json.load(stream)
    except json.JSONDecodeError as exc:
        print(exc)

R1s = [r for r in files if "_R1" in r]
# print(R1s)
if len(R1s) == 0:
    print(f"<E> length of fastq list = {len(R1s)}")
    sys.exit(1)
    
json_list = []
for r1 in R1s:
    forward_reads = [r1]
    reverse_reads = [r1.replace("_R1", "_R2")]
    json_template["wf_varpipe.read1"] = forward_reads
    json_template["wf_varpipe.read2"] = reverse_reads
    json_template["wf_varpipe.samplename"] = re.sub("_.*$", "", os.path.basename(r1))
    json_template["wf_varpipe.outdir"] = json_template["wf_varpipe.samplename"]
    json_list.append(json_template.copy())

for item in json_list:
    samplename = item["wf_varpipe.samplename"]
    output_json = pathlib.Path(samplename + ".json")
    # print(output_json)
    with open(output_json, "w", encoding="ascii") as f:
        json.dump(item, f, indent=4)
