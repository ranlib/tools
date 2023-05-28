#!/usr/bin/env python
"""
filter vcf file on intervals in bedfile
"""
import vcfpy
import argparse

def read_bed_file(file_path):
    intervals = []
    with open(file_path, "r") as file:
        for line in file:
            if line.startswith("#"):  # Skip comment lines if any
                continue
            fields = line.strip().split("\t")
            if len(fields) >= 3:
                start = int(fields[1])
                end = int(fields[2])
                intervals.append((start, end))
    return intervals


def is_position_in_intervals(position, intervals):
    for interval in intervals:
        start, end = interval
        if start <= position <= end:
            return True
    return False


if __name__ == "__main__":
    # Create argument parser
    parser = argparse.ArgumentParser(description="Intersect vcf file with bed file.")
    parser.add_argument("--vcf", "-v", help="vcf file from London", required=True)
    parser.add_argument("--bed", "-b", help="bed file from CDC", required=True)
    parser.add_argument("--out", "-o", help="output filtered vcf file", required=True)
    args = parser.parse_args()

    # Read BED file and extract intervals
    intervals = read_bed_file(args.bed)

    # Read in vcf file
    vcf_reader = vcfpy.Reader.from_path(args.vcf)

    # Filter the VCF file for the positions of the genes
    filtered_vcf_records = []
    for record in vcf_reader:
        if is_position_in_intervals(record.POS, intervals):
            filtered_vcf_records.append(record)

    # Write the filtered VCF records to a new VCF file
    vcf_writer = vcfpy.Writer.from_path(args.out, vcf_reader.header)
    for record in filtered_vcf_records:
        vcf_writer.write_record(record)

    # Close the VCF writer
    vcf_writer.close()
