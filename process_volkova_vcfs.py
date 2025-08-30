#!/usr/bin/env python3
"""
Process Volkova et al. (2020) VCF files to classify indels into ID83 categories.

This script reads the VCF files from the IND/ directory and classifies each indel
mutation according to the ID83 classification system used in COSMIC.
"""

import os
import glob
import pandas as pd
import numpy as np
from collections import defaultdict
import re

# ID83 classification categories (based on your notebook)
indel_types_83_str = [
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

def classify_indel(ref, alt, info_dict):
    """
    Classify an indel mutation into one of the 83 ID83 categories.
    
    Args:
        ref (str): Reference sequence
        alt (str): Alternate sequence
        info_dict (dict): INFO field parsed into dictionary
        
    Returns:
        str: ID83 classification category
    """
    # Get basic information
    ref_len = len(ref)
    alt_len = len(alt)
    indel_len = abs(ref_len - alt_len)
    
    # Get repeat information from INFO
    repeat_count = info_dict.get('REP', 0)
    
    # Determine mutation type
    if ref_len > alt_len:
        # Deletion
        if indel_len == 1:
            # 1bp deletion
            if ref[0] in ['C', 'G']:
                # C/G context
                if indel_len == 1:
                    if indel_len <= 6:
                        return f'DEL.C.1.{indel_len}'
                    else:
                        return 'DEL.C.1.6+'
            else:
                # T/A context
                if indel_len <= 6:
                    return f'DEL.T.1.{indel_len}'
                else:
                    return 'DEL.T.1.6+'
        else:
            # >1bp deletion
            if repeat_count > 0:
                # Repeat context
                if indel_len <= 6:
                    return f'DEL.repeats.{indel_len}.{repeat_count}'
                else:
                    return f'DEL.repeats.5+.{repeat_count}'
            else:
                # Microhomology context
                if indel_len <= 6:
                    return f'DEL.MH.{indel_len}.{repeat_count}'
                else:
                    return f'DEL.MH.5+.{repeat_count}'
    else:
        # Insertion
        if indel_len == 1:
            # 1bp insertion
            if alt[0] in ['C', 'G']:
                # C/G context
                if indel_len <= 5:
                    return f'INS.C.1.{indel_len-1}'
                else:
                    return 'INS.C.1.5+'
            else:
                # T/A context
                if indel_len <= 5:
                    return f'INS.T.1.{indel_len-1}'
                else:
                    return 'INS.T.1.5+'
        else:
            # >1bp insertion
            if repeat_count > 0:
                # Repeat context
                if indel_len <= 5:
                    return f'INS.repeats.{indel_len}.{repeat_count-1}'
                else:
                    return f'INS.repeats.5+.{repeat_count-1}'
            else:
                # Microhomology context
                if indel_len <= 5:
                    return f'INS.MH.{indel_len}.{repeat_count}'
                else:
                    return f'INS.MH.5+.{repeat_count}'

def parse_info_field(info_str):
    """Parse INFO field string into dictionary."""
    info_dict = {}
    if info_str == '.':
        return info_dict
    
    for item in info_str.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info_dict[key] = value
        else:
            info_dict[item] = True
    
    return info_dict

def process_vcf_file(vcf_path):
    """
    Process a single VCF file and return indel counts by ID83 category.
    
    Args:
        vcf_path (str): Path to VCF file
        
    Returns:
        dict: Counts of indels by ID83 category
    """
    indel_counts = defaultdict(int)
    
    with open(vcf_path, 'r') as f:
        for line in f:
            line = line.strip()
            
            # Skip header lines
            if line.startswith('#'):
                continue
            
            # Parse VCF line
            fields = line.split('\t')
            if len(fields) < 8:
                continue
            
            chrom, pos, var_id, ref, alt, qual, filter_val, info = fields[:8]
            
            # Skip if not PASS filter
            if filter_val != 'PASS':
                continue
            
            # Parse INFO field
            info_dict = parse_info_field(info)
            
            # Classify indel
            try:
                indel_class = classify_indel(ref, alt, info_dict)
                indel_counts[indel_class] += 1
            except Exception as e:
                print(f"Warning: Could not classify indel {ref}->{alt} in {vcf_path}: {e}")
                continue
    
    return dict(indel_counts)

def process_all_vcfs(ind_dir="IND"):
    """
    Process all VCF files in the IND directory.
    
    Args:
        ind_dir (str): Directory containing VCF files
        
    Returns:
        pd.DataFrame: Matrix of indel counts (samples x ID83 categories)
    """
    vcf_files = glob.glob(os.path.join(ind_dir, "*.vcf"))
    
    if not vcf_files:
        raise FileNotFoundError(f"No VCF files found in {ind_dir}")
    
    print(f"Found {len(vcf_files)} VCF files to process...")
    
    # Process each VCF file
    all_counts = {}
    
    for vcf_file in vcf_files:
        sample_name = os.path.basename(vcf_file).replace('.filtered.indels.vcf', '')
        print(f"Processing {sample_name}...")
        
        try:
            counts = process_vcf_file(vcf_file)
            all_counts[sample_name] = counts
        except Exception as e:
            print(f"Error processing {vcf_file}: {e}")
            continue
    
    # Convert to DataFrame
    df = pd.DataFrame(all_counts).T
    
    # Fill missing categories with 0
    for category in indel_types_83_str:
        if category not in df.columns:
            df[category] = 0
    
    # Reorder columns to match ID83 order
    df = df.reindex(columns=indel_types_83_str)
    
    # Fill NaN values with 0
    df = df.fillna(0)
    
    return df

def convert_to_id83(df):
    """
    Convert ID94 counts to ID83 by merging microhomology categories.
    
    Args:
        df (pd.DataFrame): DataFrame with ID94 categories
        
    Returns:
        pd.DataFrame: DataFrame with ID83 categories
    """
    # Create a copy to avoid modifying original
    df_id83 = df.copy()
    
    # Merge microhomology categories into repeat categories
    if 'INS.MH.2.1' in df_id83.columns and 'INS.repeats.2.0' in df_id83.columns:
        df_id83['INS.repeats.2.0'] += df_id83['INS.MH.2.1']
    
    if 'INS.MH.3.1' in df_id83.columns and 'INS.repeats.3.0' in df_id83.columns:
        df_id83['INS.repeats.3.0'] += df_id83['INS.MH.3.1']
        if 'INS.MH.3.2' in df_id83.columns:
            df_id83['INS.repeats.3.0'] += df_id83['INS.MH.3.2']
    
    if 'INS.MH.4.1' in df_id83.columns and 'INS.repeats.4.0' in df_id83.columns:
        df_id83['INS.repeats.4.0'] += df_id83['INS.MH.4.1']
        if 'INS.MH.4.2' in df_id83.columns:
            df_id83['INS.repeats.4.0'] += df_id83['INS.MH.4.2']
        if 'INS.MH.4.3' in df_id83.columns:
            df_id83['INS.repeats.4.0'] += df_id83['INS.MH.4.3']
    
    if 'INS.MH.5+.1' in df_id83.columns and 'INS.repeats.5+.0' in df_id83.columns:
        df_id83['INS.repeats.5+.0'] += df_id83['INS.MH.5+.1']
        if 'INS.MH.5+.2' in df_id83.columns:
            df_id83['INS.repeats.5+.0'] += df_id83['INS.MH.5+.2']
        if 'INS.MH.5+.3' in df_id83.columns:
            df_id83['INS.repeats.5+.0'] += df_id83['INS.MH.5+.3']
        if 'INS.MH.5+.4' in df_id83.columns:
            df_id83['INS.repeats.5+.0'] += df_id83['INS.MH.5+.4']
        if 'INS.MH.5+.5+' in df_id83.columns:
            df_id83['INS.repeats.5+.0'] += df_id83['INS.MH.5+.5+']
    
    # Keep only the first 83 columns (ID83 categories)
    df_id83 = df_id83.iloc[:, :83]
    
    return df_id83

def main():
    """Main function to process all VCF files."""
    print("Processing Volkova et al. (2020) VCF files...")
    
    try:
        # Process all VCF files
        df_counts = process_all_vcfs()
        
        print(f"\nProcessed {len(df_counts)} samples")
        print(f"Total indel categories: {len(df_counts.columns)}")
        
        # Save raw counts
        output_file = "volkova_indel_counts_id83.tsv"
        df_counts.to_csv(output_file, sep='\t')
        print(f"\nSaved indel counts to {output_file}")
        
        # Show summary statistics
        print("\nSummary of indel counts:")
        print(f"Total indels across all samples: {df_counts.sum().sum()}")
        print(f"Mean indels per sample: {df_counts.sum(axis=1).mean():.1f}")
        print(f"Range of indels per sample: {df_counts.sum(axis=1).min()}-{df_counts.sum(axis=1).max()}")
        
        # Show top categories
        total_by_category = df_counts.sum()
        top_categories = total_by_category.nlargest(10)
        print(f"\nTop 10 indel categories:")
        for category, count in top_categories.items():
            print(f"  {category}: {count}")
        
        # Convert to ID83 if needed
        if df_counts.shape[1] > 83:
            print("\nConverting to ID83 categories...")
            df_id83 = convert_to_id83(df_counts)
            output_id83 = "volkova_indel_counts_id83_converted.tsv"
            df_id83.to_csv(output_id83, sep='\t')
            print(f"Saved ID83 counts to {output_id83}")
        
        return df_counts
        
    except Exception as e:
        print(f"Error: {e}")
        return None

if __name__ == "__main__":
    df_result = main()