#!/usr/bin/env python
"""
determine which regions have at least one large deletion
"""
import argparse
import pandas
import vcf

def get_deletions_in_region(vcf_file: str, bed: str) -> {}:
    """
    determine which regions have at least one large deletion
    """

    # read list or regions
    regions = pandas.read_csv(bed, sep="\t")
    regions.columns = ["genome", "start", "stop", "locus", "gene", "chemical"]
    # make list of intervals
    regions_list = []
    genes_list = []
    for index, row in regions.iterrows():
        regions_list.append(pandas.Interval(left=row["start"], right=row["stop"], closed="both"))
        genes_list.append(row["gene"])

    regions_array = pandas.arrays.IntervalArray(regions_list)

    # loop over vcf
    # check DEL objects
    # check if in regions_array
    vcf_reader = vcf.Reader(filename=vcf_file)
    index_list = [ False for i in genes_list ]
    for record in vcf_reader:
        #print(record)
        for item in record.ALT:
            if str(item) == "<DEL>":
                contained_in_region = regions_array.contains(record.POS)
                if any(contained_in_region):
                    indices = [ idx for idx, value  in enumerate(contained_in_region) if value ]
                    #print(f"deletion at {record.POS} in region {indices}")
                    for index in indices:
                        index_list[index] = True

    return dict(zip(genes_list,index_list))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="get_deletions_in_region",
                                     prog="get_deletions_in_region",
                                     formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=100))
    parser.add_argument("--vcf", "-v", type=argparse.FileType("r"), help="vcf file", required=True)
    parser.add_argument("--bed", "-b", type=argparse.FileType("r"), help="bed file with regions of interest", required=True)
    parser.add_argument("--verbose", action="store_true", help="turn on debugging output")
    args = parser.parse_args()
    has_deletions = get_deletions_in_region(args.vcf.name, args.bed.name)
    if args.verbose:
        print(has_deletions)
