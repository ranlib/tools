import argparse
import vcfpy

def modify_chromosome_name(record, new_chromosome):
    # Modify the chromosome name in the record
    record.CHROM = new_chromosome

def read_vcf_file(input_file, output_file=None, new_chromosome=None):
    reader = vcfpy.Reader.from_path(input_file)
    
    writer = None
    if output_file:
        writer = vcfpy.Writer.from_path(output_file, header=reader.header)

    for record in reader:
        # Modify the chromosome name if specified
        if new_chromosome:
            modify_chromosome_name(record, new_chromosome)
        
        # Process or modify the variant information if needed
        
        # Write the record to the output VCF file if provided
        if writer:
            writer.write_record(record)

    reader.close()

    if writer:
        writer.close()

def main():
    parser = argparse.ArgumentParser(description='Read and optionally write VCF file.')
    parser.add_argument('--input_file', help='Path to the input VCF file')
    parser.add_argument('--output_file', help='Path to the output VCF file')
    parser.add_argument('--new_chromosome', help='New chromosome name')

    args = parser.parse_args()
    input_file_path = args.input_file
    output_file_path = args.output_file
    new_chromosome = args.new_chromosome

    if not input_file_path:
        parser.error('The --input_file argument is required.')

    read_vcf_file(input_file_path, output_file_path, new_chromosome)

if __name__ == '__main__':
    main()
