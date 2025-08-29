#!/usr/bin/env python3
"""
Simple processor for Volkova VCF files without external dependencies.
"""

import os
import glob
import re
from collections import defaultdict
import csv

# ID83 categories based on COSMIC classification
ID83_CATEGORIES = [
    'DEL.C.1.1', 'DEL.C.1.2', 'DEL.C.1.3', 'DEL.C.1.4', 'DEL.C.1.5', 'DEL.C.1.6+',
    'DEL.T.1.1', 'DEL.T.1.2', 'DEL.T.1.3', 'DEL.T.1.4', 'DEL.T.1.5', 'DEL.T.1.6+',
    'INS.C.1.0', 'INS.C.1.1', 'INS.C.1.2', 'INS.C.1.3', 'INS.C.1.4', 'INS.C.1.5+',
    'INS.T.1.0', 'INS.T.1.1', 'INS.T.1.2', 'INS.T.1.3', 'INS.T.1.4', 'INS.T.1.5+',
    'DEL.repeats.2.1', 'DEL.repeats.2.2', 'DEL.repeats.2.3', 'DEL.repeats.2.4', 'DEL.repeats.2.5', 'DEL.repeats.2.6+',
    'DEL.repeats.3.1', 'DEL.repeats.3.2', 'DEL.repeats.3.3', 'DEL.repeats.3.4', 'DEL.repeats.3.5', 'DEL.repeats.3.6+',
    'DEL.repeats.4.1', 'DEL.repeats.4.2', 'DEL.repeats.4.3', 'DEL.repeats.4.4', 'DEL.repeats.4.5', 'DEL.repeats.4.6+',
    'DEL.repeats.5+.1', 'DEL.repeats.5+.2', 'DEL.repeats.5+.3', 'DEL.repeats.5+.4', 'DEL.repeats.5+.5', 'DEL.repeats.5+.6+',
    'INS.repeats.2.0', 'INS.repeats.2.1', 'INS.repeats.2.2', 'INS.repeats.2.3', 'INS.repeats.2.4', 'INS.repeats.2.5+',
    'INS.repeats.3.0', 'INS.repeats.3.1', 'INS.repeats.3.2', 'INS.repeats.3.3', 'INS.repeats.3.4', 'INS.repeats.3.5+',
    'INS.repeats.4.0', 'INS.repeats.4.1', 'INS.repeats.4.2', 'INS.repeats.4.3', 'INS.repeats.4.4', 'INS.repeats.4.5+',
    'INS.repeats.5+.0', 'INS.repeats.5+.1', 'INS.repeats.5+.2', 'INS.repeats.5+.3', 'INS.repeats.5+.4', 'INS.repeats.5+.5+',
    'DEL.MH.2.1', 'DEL.MH.3.1', 'DEL.MH.3.2', 'DEL.MH.4.1', 'DEL.MH.4.2', 'DEL.MH.4.3',
    'DEL.MH.5+.1', 'DEL.MH.5+.2', 'DEL.MH.5+.3', 'DEL.MH.5+.4', 'DEL.MH.5+.5+'
]


def parse_vcf_file(vcf_path):
    """Parse VCF file and extract basic indel information."""
    indels = []
    
    with open(vcf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            
            chrom = parts[0]
            pos = int(parts[1])
            ref = parts[3]
            alt = parts[4]
            
            # Skip if not an indel
            if len(ref) == len(alt):
                continue
            
            indels.append({
                'chrom': chrom,
                'pos': pos,
                'ref': ref,
                'alt': alt,
                'type': 'DEL' if len(ref) > len(alt) else 'INS',
                'length': abs(len(ref) - len(alt))
            })
    
    return indels


def main():
    # First, let's analyze what types of indels we have
    vcf_files = glob.glob('IND/*.vcf')
    print(f"Found {len(vcf_files)} VCF files")
    
    # Count basic statistics
    total_indels = 0
    deletions = 0
    insertions = 0
    length_dist = defaultdict(int)
    
    # Sample a few files first
    sample_files = vcf_files[:5]
    
    for vcf_file in sample_files:
        print(f"\nProcessing {os.path.basename(vcf_file)}...")
        indels = parse_vcf_file(vcf_file)
        
        print(f"  Found {len(indels)} indels")
        
        for indel in indels:
            total_indels += 1
            if indel['type'] == 'DEL':
                deletions += 1
            else:
                insertions += 1
            
            length_dist[indel['length']] += 1
            
            # Print first few examples
            if total_indels <= 10:
                print(f"  Example: {indel['type']} length={indel['length']}, "
                      f"REF={indel['ref'][:20]}{'...' if len(indel['ref']) > 20 else ''}, "
                      f"ALT={indel['alt'][:20]}{'...' if len(indel['alt']) > 20 else ''}")
    
    print(f"\nSummary from {len(sample_files)} files:")
    print(f"Total indels: {total_indels}")
    print(f"Deletions: {deletions}")
    print(f"Insertions: {insertions}")
    print(f"\nLength distribution:")
    for length in sorted(length_dist.keys())[:10]:
        print(f"  Length {length}: {length_dist[length]} indels")
    
    # Now let's look at the actual sequences in more detail
    print("\n\nDetailed analysis of first file:")
    first_file = vcf_files[0]
    indels = parse_vcf_file(first_file)
    
    # Save detailed info to CSV
    with open('indel_examples.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Type', 'Length', 'REF', 'ALT', 'Chrom', 'Pos'])
        
        for i, indel in enumerate(indels[:20]):
            writer.writerow([
                indel['type'],
                indel['length'],
                indel['ref'],
                indel['alt'],
                indel['chrom'],
                indel['pos']
            ])
    
    print(f"Saved first 20 indels to indel_examples.csv for inspection")


if __name__ == "__main__":
    main()