#!/usr/bin/env python
"""
get average coverage per region
get % region covered above threshold
"""
import pysam
import pandas

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
                
            average_depth = total_depth / total_bases
            average_depths.append(average_depth)

            coverage_percentage = (covered_bases / total_bases) * 100
            coverage_percentages.append(coverage_percentage)

            coverage_percentage_1 = (covered_bases_1 / total_bases) * 100
            coverage_percentages_1.append(coverage_percentage_1)

        d["total_coverage"] = total_depths
        d["average_coverage"] = average_depths
        d["percentage >= 1"] = coverage_percentages_1
        d["percentage >= threshold"] = coverage_percentages

    return d

if __name__ == "__main__":
    # Usage example
    bam_file_path = "/home/dieterbest/Analysis/varpipe4/tests/ERR552797/ERR552797_sdrcsm.bam"
    bed_file_path = "/home/dieterbest/Analysis/varpipe4/intervals/tbdb.bed"
    minimal_coverage = 10 
    d = calculate_average_depth(bam_file_path, bed_file_path, minimal_coverage)
    d.to_csv("coverage.tsv", index=False, sep="\t")
