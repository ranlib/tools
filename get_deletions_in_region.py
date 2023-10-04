#!/usr/bin/env python
"""
determine which regions have at least one large deletion
"""
import argparse
import pandas
import sympy
import vcf
import json
from pprint import pprint

def get_deletions_in_region(vcf_file: str, bed: str, filter_SVs: bool = True) -> {}:
    """
    determine which regions have at least one large deletion
    """

    # read list or regions
    # make list of intervals
    # attach gene and has_deletion information
    regions = pandas.read_csv(bed, header = None, sep="\t")
    regions.columns = ["genome", "start", "stop", "locus", "gene", "chemical"]
    regions_list: list[ [pandas.Interval, str, bool, int, int, int] ] = []
    for index, row in regions.iterrows():
        regions_list.append([ pandas.Interval(left=row["start"], right=row["stop"], closed="both"), row["gene"], False, -1, -1, -1] )

    regions_array = pandas.arrays.IntervalArray([item[0] for item in regions_list])

    # loop over vcf
    # check DEL objects
    # check if in regions_array
    vcf_reader = vcf.Reader(filename=vcf_file)
    for record in vcf_reader:
        for item in record.ALT:
            if "SVTYPE" in record.INFO:
                is_SV = True
            else:
                is_SV = False
                
            if is_SV and record.INFO["SVTYPE"] == "DEL" and ( not filter_SVs or not record.FILTER):
                deletion = pandas.Interval(record.POS, record.INFO["END"])
                contained_in_region = regions_array.overlaps(deletion)
                if any(contained_in_region):
                    index_list = [ idx for idx, value  in enumerate(contained_in_region) if value ]
                    for index in index_list:
                        regions_list[index][2] = True
                        regions_list[index][3] = deletion.length
                        s_deletion = sympy.Interval(deletion.left, deletion.right)
                        s_region = sympy.Interval(regions_list[index][0].left, regions_list[index][0].right)
                        overlap = s_deletion.intersect(s_region)
                        regions_list[index][4] = overlap.measure
                        regions_list[index][5] = s_region.measure
                            
    return regions_list

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="get_deletions_in_region",
                                     prog="get_deletions_in_region",
                                     formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=100))
    parser.add_argument("--vcf", "-v", type=argparse.FileType("r"), help="vcf file", required=True)
    parser.add_argument("--bed", "-b", type=argparse.FileType("r"), help="bed file with regions of interest", required=True)
    parser.add_argument("--filter_SVs", "-f", action="store_true", help="filter out SVs that do not have a PASS in filter column")
    parser.add_argument("--verbose", action="store_true", help="turn on debugging output")
    args = parser.parse_args()
    regions_list = get_deletions_in_region(args.vcf.name, args.bed.name, args.filter_SVs)
    has_deletions = { item[1]:item[2] for item in regions_list }
    if args.verbose:
        #print(json.dumps(has_deletions,sort_keys=True, indent=4))
        pprint(regions_list)
