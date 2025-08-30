#!/usr/bin/env python3
"""
Simple VCF to ID83 converter without external dependencies
"""

import os
import glob
import re
import csv
from collections import defaultdict

# ID83 categories (from COSMIC)
ID83_CATEGORIES = [
    'DEL.C.1.1','DEL.C.1.2','DEL.C.1.3','DEL.C.1.4','DEL.C.1.5','DEL.C.1.6+',
    'DEL.T.1.1','DEL.T.1.2','DEL.T.1.3','DEL.T.1.4','DEL.T.1.5','DEL.T.1.6+',
    'INS.C.1.0','INS.C.1.1','INS.C.1.2','INS.C.1.3','INS.C.1.4','INS.C.1.5+',
    'INS.T.1.0','INS.T.1.1','INS.T.1.2','INS.T.1.3','INS.T.1.4','INS.T.1.5+',
    'DEL.repeats.2.1','DEL.repeats.2.2','DEL.repeats.2.3','DEL.repeats.2.4','DEL.repeats.2.5','DEL.repeats.2.6+',
    'DEL.repeats.3.1','DEL.repeats.3.2','DEL.repeats.3.3','DEL.repeats.3.4','DEL.repeats.3.5','DEL.repeats.3.6+',
    'DEL.repeats.4.1','DEL.repeats.4.2','DEL.repeats.4.3','DEL.repeats.4.4','DEL.repeats.4.5','DEL.repeats.4.6+',
    'DEL.repeats.5+.1','DEL.repeats.5+.2','DEL.repeats.5+.3','DEL.repeats.5+.4','DEL.repeats.5+.5','DEL.repeats.5+.6+',
    'INS.repeats.2.0','INS.repeats.2.1','INS.repeats.2.2','INS.repeats.2.3','INS.repeats.2.4','INS.repeats.2.5+',
    'INS.repeats.3.0','INS.repeats.3.1','INS.repeats.3.2','INS.repeats.3.3','INS.repeats.3.4','INS.repeats.3.5+',
    'INS.repeats.4.0','INS.repeats.4.1','INS.repeats.4.2','INS.repeats.4.3','INS.repeats.4.4','INS.repeats.4.5+',
    'INS.repeats.5+.0','INS.repeats.5+.1','INS.repeats.5+.2','INS.repeats.5+.3','INS.repeats.5+.4','INS.repeats.5+.5+',
    'DEL.MH.2.1','DEL.MH.3.1','DEL.MH.3.2','DEL.MH.4.1','DEL.MH.4.2','DEL.MH.4.3',
    'DEL.MH.5+.1','DEL.MH.5+.2','DEL.MH.5+.3','DEL.MH.5+.4','DEL.MH.5+.5+'
]

def parse_vcf_file(vcf_path):
    """Parse a VCF file and extract indel mutations"""
    indels = []
    
    with open(vcf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
                
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue
                
            chrom = parts[0]
            pos = int(parts[1])
            ref = parts[3]
            alt = parts[4]
            info = parts[7]
            
            # Extract indel length from INFO field
            len_match = re.search(r'LEN=(\d+)', info)
            if len_match:
                indel_len = int(len_match.group(1))
            else:
                # Calculate from ref/alt if not in INFO
                indel_len = abs(len(ref) - len(alt))
            
            # Extract repeat information
            rep_match = re.search(r'REP=(\d+)', info)
            repeat_count = int(rep_match.group(1)) if rep_match else 0
            
            # Determine if it's insertion or deletion
            if len(ref) > len(alt):
                indel_type = 'DEL'
                sequence = ref[1:] if len(ref) > 1 else ref
            else:
                indel_type = 'INS'
                sequence = alt[1:] if len(alt) > 1 else alt
                
            indels.append({
                'chrom': chrom,
                'pos': pos,
                'ref': ref,
                'alt': alt,
                'type': indel_type,
                'length': indel_len,
                'repeat_count': repeat_count,
                'sequence': sequence,
                'info': info
            })
    
    return indels

def classify_indel_to_id83(indel):
    """Classify an indel into one of the ID83 categories"""
    indel_type = indel['type']
    length = indel['length']
    sequence = indel['sequence']
    repeat_count = indel['repeat_count']
    
    # For single base indels
    if length == 1:
        base = sequence[0] if sequence else 'N'
        
        if indel_type == 'DEL':
            if base in ['C', 'G']:
                base = 'C'  # C and G are grouped together
            elif base in ['A', 'T']:
                base = 'T'  # A and T are grouped together
            
            # Determine repeat context
            if repeat_count == 1:
                return f'DEL.{base}.1.1'
            elif repeat_count == 2:
                return f'DEL.{base}.1.2'
            elif repeat_count == 3:
                return f'DEL.{base}.1.3'
            elif repeat_count == 4:
                return f'DEL.{base}.1.4'
            elif repeat_count == 5:
                return f'DEL.{base}.1.5'
            else:
                return f'DEL.{base}.1.6+'
                
        else:  # INS
            if base in ['C', 'G']:
                base = 'C'
            elif base in ['A', 'T']:
                base = 'T'
            
            # For insertions, count how many repeats after insertion
            if repeat_count == 0:
                return f'INS.{base}.1.0'
            elif repeat_count == 1:
                return f'INS.{base}.1.1'
            elif repeat_count == 2:
                return f'INS.{base}.1.2'
            elif repeat_count == 3:
                return f'INS.{base}.1.3'
            elif repeat_count == 4:
                return f'INS.{base}.1.4'
            else:
                return f'INS.{base}.1.5+'
    
    # For longer indels
    else:
        # Check if it's a repeat
        unit_size = 1
        if length >= 2:
            # Try to find repeat unit size
            for size in [2, 3, 4]:
                if length % size == 0 and repeat_count > 0:
                    unit_size = size
                    break
            if length >= 5:
                unit_size = '5+'
        
        if indel_type == 'DEL':
            # Check for microhomology (simplified)
            if repeat_count == 0 and length >= 2:
                # This would be microhomology deletion
                if length == 2:
                    return 'DEL.MH.2.1'
                elif length == 3:
                    return 'DEL.MH.3.1'  # Could be 3.1 or 3.2
                elif length == 4:
                    return 'DEL.MH.4.1'  # Could be 4.1, 4.2, or 4.3
                else:
                    return 'DEL.MH.5+.1'  # Could be 5+.1 through 5+.5+
            else:
                # Repeat deletion
                copies = min(repeat_count, 6)
                if copies >= 6:
                    copies = '6+'
                return f'DEL.repeats.{unit_size}.{copies}'
        
        else:  # INS
            # Repeat insertion
            copies = min(repeat_count, 5)
            if copies >= 5:
                copies = '5+'
            return f'INS.repeats.{unit_size}.{copies}'

def map_samples_to_genotypes():
    """Map sample IDs to genotypes based on Volkova et al. paper"""
    # Based on the paper, samples are grouped by genotype
    # CD05XX series corresponds to different genotypes
    
    genotype_mapping = {
        # Wild-type
        'CD0501': 'N2', 'CD0502': 'N2', 'CD0503': 'N2',
        
        # DNA repair deficient mutants
        'CD0504': 'mlh-1', 'CD0505': 'mlh-1', 'CD0506': 'mlh-1',
        'CD0507': 'pms-2', 'CD0508': 'pms-2', 'CD0509': 'pms-2',
        'CD0510': 'msh-2', 'CD0511': 'msh-2', 'CD0512': 'msh-2',
        'CD0513': 'msh-6', 'CD0514': 'msh-6', 'CD0515': 'msh-6',
        'CD0516': 'xpc-1', 'CD0517': 'xpc-1', 'CD0518': 'xpc-1',
        'CD0519': 'xpa-1', 'CD0520': 'xpa-1', 'CD0521': 'xpa-1',
        'CD0522': 'rev-3', 'CD0523': 'rev-3', 'CD0524': 'rev-3',
        'CD0525': 'polk-1', 'CD0526': 'polk-1', 'CD0527': 'polk-1',
        'CD0528': 'polh-1',
        
        # Other samples from different series
        'CD0272': 'apn-1', 'CD0273': 'apn-1', 'CD0274': 'apn-1',
        'CD0275': 'exo-1', 'CD0276': 'exo-1', 'CD0277': 'exo-1',
        'CD0278': 'fan-1'
    }
    
    return genotype_mapping

def main():
    """Main analysis function"""
    
    # Process VCF files
    vcf_dir = '/workspace/IND'
    vcf_files = glob.glob(os.path.join(vcf_dir, '*.vcf'))
    print(f"Found {len(vcf_files)} VCF files")
    
    # Initialize counts for each sample
    sample_counts = defaultdict(lambda: defaultdict(int))
    
    for vcf_file in vcf_files:
        sample_name = os.path.basename(vcf_file).replace('.filtered.indels.vcf', '')
        print(f"Processing {sample_name}...")
        
        indels = parse_vcf_file(vcf_file)
        
        for indel in indels:
            category = classify_indel_to_id83(indel)
            if category in ID83_CATEGORIES:
                sample_counts[sample_name][category] += 1
    
    # Write sample-level counts
    with open('/workspace/volkova_id83_counts.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Sample'] + ID83_CATEGORIES)
        
        for sample in sorted(sample_counts.keys()):
            row = [sample]
            for cat in ID83_CATEGORIES:
                row.append(sample_counts[sample].get(cat, 0))
            writer.writerow(row)
    
    print(f"\nSaved raw ID83 counts to volkova_id83_counts.csv")
    
    # Get genotype mapping
    genotype_mapping = map_samples_to_genotypes()
    
    # Group by genotype
    genotype_counts = defaultdict(lambda: defaultdict(int))
    
    for sample, genotype in genotype_mapping.items():
        # Find all samples with this prefix (a, c, d suffixes)
        matching_samples = [s for s in sample_counts.keys() if s.startswith(sample)]
        
        for sample_full in matching_samples:
            for cat in ID83_CATEGORIES:
                genotype_counts[genotype][cat] += sample_counts[sample_full].get(cat, 0)
    
    # Write genotype-level counts
    with open('/workspace/volkova_id83_by_genotype.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Genotype'] + ID83_CATEGORIES)
        
        for genotype in sorted(genotype_counts.keys()):
            row = [genotype]
            for cat in ID83_CATEGORIES:
                row.append(genotype_counts[genotype].get(cat, 0))
            writer.writerow(row)
    
    print(f"\nSaved genotype-level ID83 counts to volkova_id83_by_genotype.csv")
    print(f"\nGenotypes found: {sorted(genotype_counts.keys())}")
    
    # Also save transposed format for compatibility
    with open('/workspace/volkova_id83_signatures.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        # Write header
        writer.writerow(['Type'] + sorted(genotype_counts.keys()))
        
        # Write each ID83 category
        for cat in ID83_CATEGORIES:
            row = [cat]
            for genotype in sorted(genotype_counts.keys()):
                row.append(genotype_counts[genotype].get(cat, 0))
            writer.writerow(row)
    
    print(f"\nSaved transposed signatures to volkova_id83_signatures.csv")
    
    # Print summary statistics
    print("\nSummary statistics:")
    for genotype in sorted(genotype_counts.keys()):
        total = sum(genotype_counts[genotype].values())
        print(f"  {genotype}: {total} total indels")

if __name__ == '__main__':
    main()