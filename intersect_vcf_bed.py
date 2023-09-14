#!/usr/bin/env python

import pandas
import argparse
import vcf

def intersect_vcf_bed(vcf_file, bed_file, output_vcf):
    # turn regions into pandas intervals array
    regions = pandas.read_csv(bed_file, header=None, sep="\t")
    regions.columns = ["genome", "start", "stop", "gene_id", "gene", "chemical"]
    intervals = []
    for index, row in regions.iterrows():
        intervals.append( pandas.Interval(row["start"], row["stop"], closed="both") )
    interval_array = pandas.arrays.IntervalArray(intervals)
    # loop over vcf and filter out large deletions not overlapping any region
    vcf_reader = vcf.Reader(filename=vcf_file.name)
    output = vcf.Writer(output_vcf, vcf_reader)
    del_counter = 0
    del_counter_overlapping = 0
    for record in vcf_reader:
        alt_alleles = [str(allele) for allele in record.ALT]
        keep_deletion = False
        if '<DEL>' in alt_alleles: # use only DEL with PASS in filter column
            del_counter += 1
            deletion = pandas.Interval(record.POS, record.INFO["END"], closed="both")
            for interval in interval_array:
                if interval.overlaps(deletion):
                    keep_deletion = True
            if keep_deletion and not record.FILTER:
                del_counter_overlapping += 1
                output.write_record(record)
        else:
            output.write_record(record)
    print(f"number of deletions = {del_counter}, passed {del_counter_overlapping}")
    return 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="intersect vcf with bed file", prog="intersect_vcf_bed", formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=100))
    parser.add_argument("--vcf", "-v", type=argparse.FileType("r"), help="annotated vcf file", required=True)
    parser.add_argument("--bed", "-b", type=argparse.FileType("r"), help="bed file with regions of interest", required=True)
    parser.add_argument("--out", "-o", type=argparse.FileType("w"), help="filtered output vcf file", required=True)
    args = parser.parse_args()
    intersect_vcf_bed(args.vcf, args.bed, args.out)
