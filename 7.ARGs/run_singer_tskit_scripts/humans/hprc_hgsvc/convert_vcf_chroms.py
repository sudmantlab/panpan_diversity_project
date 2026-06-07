import sys
import gzip
import pysam
import re

def read_alias_file(alias_file):
    """Read the alias file and create a mapping from RefSeq to UCSC names"""
    refseq_to_ucsc = {}
    with open(alias_file, 'r') as f:
        header = f.readline().strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            refseq = parts[1].strip()
            ucsc = parts[4].strip()
            if refseq and ucsc:
                refseq_to_ucsc[refseq] = ucsc
    return refseq_to_ucsc

def convert_vcf(input_vcf, output_vcf, refseq_to_ucsc):
    """Convert chromosome names in VCF file with proper BGZF compression"""
    # Open input file (handles both compressed/uncompressed)
    if input_vcf.endswith('.gz'):
        infile = gzip.open(input_vcf, 'rt')
    else:
        infile = open(input_vcf, 'r')
    
    # Create BGZF-compressed output
    with pysam.BGZFile(output_vcf, 'w') as bgzf:
        for line in infile:
            if line.startswith('##contig=<'):
                # Process contig header lines with robust parsing
                line = line.strip()
                content = line[10:-1]  # Remove '##contig=<' and '>'
                attributes = content.split(',')
                new_attrs = []
                for attr in attributes:
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        if key == 'ID' and value in refseq_to_ucsc:
                            new_attrs.append(f'ID={refseq_to_ucsc[value]}')
                        else:
                            new_attrs.append(attr)
                    else:
                        new_attrs.append(attr)
                new_line = f"##contig=<{','.join(new_attrs)}>\n"
                bgzf.write(new_line.encode())
            
            elif line.startswith('#'):
                # Preserve other header lines
                bgzf.write(line.encode())
            
            else:
                # Process data lines
                fields = line.split('\t')
                chrom = fields[0]
                if chrom in refseq_to_ucsc:
                    fields[0] = refseq_to_ucsc[chrom]
                bgzf.write('\t'.join(fields).encode())
    
    infile.close()
    # Create tabix index
    pysam.tabix_index(output_vcf, preset="vcf", force=True)
    print(f"Successfully created BGZF-compressed VCF and index: {output_vcf}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python convert_vcf_chroms.py <alias_file> <input_vcf> <output_vcf.gz>")
        sys.exit(1)
    
    alias_file = sys.argv[1]
    input_vcf = sys.argv[2]
    output_vcf = sys.argv[3]
    
    if not output_vcf.endswith('.gz'):
        print("Error: Output file must end with .gz for BGZF compression")
        sys.exit(1)
    
    refseq_to_ucsc = read_alias_file(alias_file)
    convert_vcf(input_vcf, output_vcf, refseq_to_ucsc)
