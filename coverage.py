#!/usr/bin/env python
"""
get average coverage per region
get % region covered above threshold
"""
import pysam
import pandas
import pathlib
import argparse

def calculate_average_depth(bam_file: str, bed_file: str , minimal_coverage: int) -> pandas.DataFrame:
    """
    get average coverage per region
    get % region covered above threshold
    """
    
    bed_regions = pandas.read_csv(bed_file, sep="\t", header=None)
    bed_regions.columns = ["genome", "start", "stop", "gene_id", "gene", "chemical"]
    d = bed_regions.loc[:, ["genome", "start", "stop"]]
    
    with pysam.AlignmentFile(bam_file, "rb")  as samfile:

        total_depths = []
        average_depths = []
        coverage_percentages = []
        coverage_percentages_1 = []

        for _, region in bed_regions.iterrows():
            chrom, start, end, gene_id,  gene, chemical = region
            total_depth = 0
            total_bases = 0
            covered_bases = 0
            covered_bases_1 = 0

            for pileupcolumn in samfile.pileup(chrom, start, end, truncate=True, stepper="samtools"):
                depth = pileupcolumn.nsegments
                total_depth += depth
                total_bases += 1
                #print(f"{chrom}\t{start}\t{end}\t{pileupcolumn.pos}\t{depth}")
                
                if depth >= minimal_coverage:
                    covered_bases += 1

                if depth >= 1:
                    covered_bases_1 += 1

            total_depths.append(total_depth)
                
            average_depth = total_depth / total_bases if total_bases > 0 else 0
            average_depths.append(average_depth)

            coverage_percentage = (covered_bases / total_bases) * 100 if total_bases > 0 else 0
            coverage_percentages.append(coverage_percentage)

            coverage_percentage_1 = (covered_bases_1 / total_bases) * 100 if total_bases > 0 else 0
            coverage_percentages_1.append(coverage_percentage_1)

        d["total_coverage"] = total_depths
        d["average_coverage"] = [ round(i,2) for i in average_depths ]
        d["%_above_1"] = [  round(i,2) for i in coverage_percentages_1 ]
        d["percent_above_threshold"] = [ round(i,2) for i in coverage_percentages ]

    return d

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="calculate depth of coverage", prog="coverage", formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=80))
    parser.add_argument("--bam", "-b", type=pathlib.Path, help="bam file", required=True)
    parser.add_argument("--bed", "-i", type=pathlib.Path, help="intervals bed file", required=True)
    parser.add_argument("--minimum_coverage", "-m", type=int, help="minimum coverage", required=True)
    parser.add_argument("--tsv", "-t", type=str, help="output tsv file", required=True)
    args = parser.parse_args()

    d = calculate_average_depth(args.bam, args.bed, args.minimum_coverage)
    d.to_csv(args.tsv, index=False, sep="\t")
