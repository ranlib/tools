#!/usr/bin/env python
"""
compare keys of 2 json files
"""
import json
import argparse
from pprint import pprint

parser = argparse.ArgumentParser(prog="diff_json_files", description="compare keys in 2 json files")
parser.add_argument("--json1", "-j1", required=True, type=str, help="json file 1")
parser.add_argument("--json2", "-j2", required=True, type=str, help="json file 2")
args = parser.parse_args()

with open(args.json1, 'r') as json1:
  data1 = json.load(json1)
  
with open(args.json2, 'r') as json2:
  data2 = json.load(json2)

#print(data1.keys())
#print(data2.keys())

set1 = set(data1.keys())
set2 = set(data2.keys())

# elements of set1 that are not in set2
not_in_set2 = set1.difference(set2)
pprint(not_in_set2)

# elements of set2 that are not in set1
not_in_set1 = set2.difference(set1)
pprint(not_in_set1)
