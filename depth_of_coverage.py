#!/usr/bin/env python
"""
create a coverage report for genomic intervals
"""
import sys
import os
import shutil
import pandas
import argparse
from pathlib import Path

parser = argparse.ArgumentParser(prog="depth of coverage")
parser.add_argument("--bam",help="bam file",required=True)
parser.add_argument("--reference",help="reference genome",required=True)
parser.add_argument("--intervals",help="intervals",required=True)
parser.add_argument("--lower_coverage",help="lower coverage",required=True)
parser.add_argument("--output",help="output base name",required=True)
parser.add_argument("--csv",help="output csv file",required=True)
args = parser.parse_args()

# create output directory
if os.path.isdir(args.output):
    print(f"<E> directory {args.output} already exists.")
    sys.exit(1)
else:
    os.makedirs(args.output)

# create faidx file
c = f"samtools faidx {args.reference}"
os.system(c)

# create dict file
dict_file = Path(args.reference)
dict_file = dict_file.with_suffix(".dict")
if not dict_file.exists():
    c = f"gatk CreateSequenceDictionary --REFERENCE {args.reference}"
    os.system(c)

# index bam file
c = f"gatk BuildBamIndex --INPUT {args.bam}"
os.system(c)

# coverage with lower coverage cut off at args.lower_coverage
output = os.path.join(args.output, args.output)
c = f"gatk DepthOfCoverage --input {args.bam} --reference {args.reference} --intervals {args.intervals} --output {output} --summary-coverage-threshold {args.lower_coverage}"
os.system(c)

# coverage with lower coverage cut off at 1
output1 = os.path.join(args.output, args.output + "_1")
c = f"gatk DepthOfCoverage --input {args.bam} --reference {args.reference} --intervals {args.intervals} --output {output1} --summary-coverage-threshold 1"
os.system(c)

# write csv file
f = output1 + ".sample_interval_summary"
d = pandas.read_csv(f)
d = d.drop(columns=d.columns[d.columns.str.contains('granular', regex=True)])
d = d.drop(columns=d.columns[d.columns.str.contains('cvg', regex=True)])
header = [ e.replace(args.output+"_", "") for e in d.columns ]
d.columns = header
d[['Target', 'Start-Stop']] = d['Target'].str.split(':', 1, expand=True)
d[['Start','Stop']] = d['Start-Stop'].str.split('-', 1, expand=True)
d = d.drop(["Start-Stop"], axis=1)
d = d.reindex(columns=['Target', 'Start', 'Stop', 'total_coverage', 'average_coverage','%_above_1'])

g = output + ".sample_interval_summary"
df = pandas.read_csv(g)
d["%_above_10"] = df.iloc[:, -1]

# write to file
d.to_csv(args.csv, index=False)

# clean up
shutil.rmtree(args.output)
