import argparse
import subprocess

def run_freebayes(bam_file, reference_fasta, output_file, ploidy):
    # Build FreeBayes command
    command = f'freebayes -f {reference_fasta} {bam_file} --vcf {output_file} --ploidy {ploidy}'

    # Run FreeBayes using subprocess
    process = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    _, stderr = process.communicate()

    # Check for errors
    if process.returncode != 0:
        print(f"FreeBayes returned an error:\n{stderr.decode()}")
    else:
        print(f"Variant calling completed successfully. VCF output written to {output_file}")

if __name__ == '__main__':
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Run FreeBayes on an alignment BAM file')
    parser.add_argument('--bam_file', type=str, required=True, help='Path to the alignment BAM file')
    parser.add_argument('--reference_fasta', type=str, required=True, help='Path to the reference genome FASTA file')
    parser.add_argument('--output_file', type=str, required=True, help='Path to the output VCF file')
    parser.add_argument('--ploidy', type=int, default=2, help='Ploidy level for variant calling')
    args = parser.parse_args()

    # Call the FreeBayes function with the provided arguments
    run_freebayes(args.bam_file, args.reference_fasta, args.output_file, args.ploidy)
