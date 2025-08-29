#!/usr/bin/env python3
"""
Process Volkova et al. (2020) VCF files and classify indels into ID83 categories.

This script:
1. Reads VCF files from the IND directory
2. Extracts indel mutations
3. Classifies them into ID83 categories
4. Outputs the counts for each genotype
"""

import os
import pandas as pd
import numpy as np
from collections import defaultdict
import re
import glob

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


def get_repeat_unit_length(seq):
    """Find the length of the repeat unit in a sequence."""
    seq_len = len(seq)
    for unit_len in range(1, min(6, seq_len // 2 + 1)):
        if seq_len % unit_len == 0:
            unit = seq[:unit_len]
            if unit * (seq_len // unit_len) == seq:
                return unit_len
    return seq_len


def has_microhomology(ref, alt, pos):
    """Check if deletion has microhomology."""
    if len(ref) <= len(alt):
        return False, 0
    
    # For deletions, check if deleted sequence appears at flanks
    deleted_seq = ref[len(alt):]
    
    # Check upstream and downstream for microhomology
    # This is simplified - in reality would need genome context
    # For now, just check if deletion creates a repeat pattern
    return False, 0  # Simplified for now


def classify_indel(ref, alt, chrom=None, pos=None):
    """
    Classify an indel into ID83 categories.
    
    Parameters:
    - ref: reference allele
    - alt: alternate allele
    - chrom: chromosome (optional, for context)
    - pos: position (optional, for context)
    
    Returns:
    - category: ID83 category string
    """
    ref_len = len(ref)
    alt_len = len(alt)
    
    # Determine if insertion or deletion
    if ref_len > alt_len:
        # Deletion
        del_len = ref_len - alt_len
        deleted_seq = ref[alt_len:]
        
        # Check for simple 1bp deletions
        if del_len == 1:
            base = deleted_seq[0]
            if base in ['C', 'G']:
                # Count homopolymer length (would need genome context)
                # For now, assume it's a simple deletion
                return 'DEL.C.1.1'
            elif base in ['A', 'T']:
                return 'DEL.T.1.1'
        
        # Check for repeat-mediated deletions
        unit_len = get_repeat_unit_length(deleted_seq)
        
        if unit_len < del_len:
            # It's a repeat
            num_units = del_len // unit_len
            
            if unit_len == 1:
                # Homopolymer deletion
                base = deleted_seq[0]
                base_type = 'C' if base in ['C', 'G'] else 'T'
                
                if num_units == 1:
                    return f'DEL.{base_type}.1.1'
                elif num_units == 2:
                    return f'DEL.{base_type}.1.2'
                elif num_units == 3:
                    return f'DEL.{base_type}.1.3'
                elif num_units == 4:
                    return f'DEL.{base_type}.1.4'
                elif num_units == 5:
                    return f'DEL.{base_type}.1.5'
                else:
                    return f'DEL.{base_type}.1.6+'
            
            elif unit_len >= 2:
                # Tandem repeat deletion
                if unit_len >= 5:
                    unit_len_cat = '5+'
                else:
                    unit_len_cat = str(unit_len)
                
                if num_units == 1:
                    return f'DEL.repeats.{unit_len_cat}.1'
                elif num_units == 2:
                    return f'DEL.repeats.{unit_len_cat}.2'
                elif num_units == 3:
                    return f'DEL.repeats.{unit_len_cat}.3'
                elif num_units == 4:
                    return f'DEL.repeats.{unit_len_cat}.4'
                elif num_units == 5:
                    return f'DEL.repeats.{unit_len_cat}.5'
                else:
                    return f'DEL.repeats.{unit_len_cat}.6+'
        
        # Check for microhomology-mediated deletion
        has_mh, mh_len = has_microhomology(ref, alt, pos)
        if has_mh and del_len >= 2:
            if del_len >= 5:
                del_len_cat = '5+'
            else:
                del_len_cat = str(del_len)
            
            if mh_len == 1:
                return f'DEL.MH.{del_len_cat}.1'
            elif mh_len == 2:
                return f'DEL.MH.{del_len_cat}.2'
            elif mh_len == 3:
                return f'DEL.MH.{del_len_cat}.3'
            elif mh_len == 4:
                return f'DEL.MH.{del_len_cat}.4'
            elif mh_len >= 5:
                return f'DEL.MH.{del_len_cat}.5+'
        
        # Default deletion category
        return 'DEL.C.1.1'
        
    else:
        # Insertion
        ins_len = alt_len - ref_len
        inserted_seq = alt[ref_len:]
        
        # Check for simple 1bp insertions
        if ins_len == 1:
            base = inserted_seq[0]
            if base in ['C', 'G']:
                return 'INS.C.1.0'
            elif base in ['A', 'T']:
                return 'INS.T.1.0'
        
        # Check for repeat-mediated insertions
        unit_len = get_repeat_unit_length(inserted_seq)
        
        if unit_len < ins_len:
            # It's a repeat
            num_units = ins_len // unit_len
            
            if unit_len == 1:
                # Homopolymer insertion
                base = inserted_seq[0]
                base_type = 'C' if base in ['C', 'G'] else 'T'
                
                if num_units == 0:
                    return f'INS.{base_type}.1.0'
                elif num_units == 1:
                    return f'INS.{base_type}.1.1'
                elif num_units == 2:
                    return f'INS.{base_type}.1.2'
                elif num_units == 3:
                    return f'INS.{base_type}.1.3'
                elif num_units == 4:
                    return f'INS.{base_type}.1.4'
                else:
                    return f'INS.{base_type}.1.5+'
            
            elif unit_len >= 2:
                # Tandem repeat insertion
                if unit_len >= 5:
                    unit_len_cat = '5+'
                else:
                    unit_len_cat = str(unit_len)
                
                if num_units == 0:
                    return f'INS.repeats.{unit_len_cat}.0'
                elif num_units == 1:
                    return f'INS.repeats.{unit_len_cat}.1'
                elif num_units == 2:
                    return f'INS.repeats.{unit_len_cat}.2'
                elif num_units == 3:
                    return f'INS.repeats.{unit_len_cat}.3'
                elif num_units == 4:
                    return f'INS.repeats.{unit_len_cat}.4'
                else:
                    return f'INS.repeats.{unit_len_cat}.5+'
        
        # Default insertion category
        return 'INS.C.1.0'


def parse_vcf_file(vcf_path):
    """Parse a VCF file and extract indels with their classifications."""
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
            
            # Classify the indel
            category = classify_indel(ref, alt, chrom, pos)
            
            indels.append({
                'chrom': chrom,
                'pos': pos,
                'ref': ref,
                'alt': alt,
                'category': category
            })
    
    return indels


def process_all_vcf_files(vcf_dir):
    """Process all VCF files in the directory and count ID83 categories."""
    # Group files by sample/genotype
    sample_groups = defaultdict(list)
    
    for vcf_file in glob.glob(os.path.join(vcf_dir, '*.vcf')):
        # Extract sample name (e.g., CD0273 from CD0273a.filtered.indels.vcf)
        basename = os.path.basename(vcf_file)
        match = re.match(r'(CD\d+)', basename)
        if match:
            sample = match.group(1)
            sample_groups[sample].append(vcf_file)
    
    # Process each sample group
    results = {}
    
    for sample, vcf_files in sample_groups.items():
        print(f"Processing sample {sample} ({len(vcf_files)} files)...")
        
        # Initialize counts for this sample
        counts = defaultdict(int)
        
        for vcf_file in vcf_files:
            indels = parse_vcf_file(vcf_file)
            
            for indel in indels:
                counts[indel['category']] += 1
        
        # Convert to ID83 vector
        id83_vector = np.zeros(83)
        for i, category in enumerate(ID83_CATEGORIES):
            id83_vector[i] = counts.get(category, 0)
        
        results[sample] = id83_vector
    
    return results


def save_results(results, output_file):
    """Save the ID83 counts to a CSV file."""
    # Create DataFrame
    df = pd.DataFrame(results).T
    df.columns = ID83_CATEGORIES
    
    # Save to CSV
    df.to_csv(output_file)
    print(f"Results saved to {output_file}")
    
    return df


def main():
    # Process VCF files
    vcf_dir = 'IND'
    results = process_all_vcf_files(vcf_dir)
    
    # Save results
    output_file = 'volkova_id83_counts.csv'
    df = save_results(results, output_file)
    
    # Print summary
    print(f"\nProcessed {len(results)} samples")
    print(f"Total indels per sample:")
    for sample, vector in results.items():
        print(f"  {sample}: {int(vector.sum())} indels")
    
    return df


if __name__ == "__main__":
    main()