#!/usr/bin/env python3
"""
Simplified script to process Volkova et al. (2020) VCF files and classify indels into ID83 categories.

This version only uses Python standard library modules.
"""

import os
import glob
from collections import defaultdict
import re

# ID83 classification categories
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

def classify_indel_simple(ref, alt, info_dict):
    """
    Simplified indel classification based on the VCF structure.
    
    Args:
        ref (str): Reference sequence
        alt (str): Alternate sequence  
        info_dict (dict): Parsed INFO field
        
    Returns:
        str: ID83 classification category
    """
    ref_len = len(ref)
    alt_len = len(alt)
    indel_len = abs(ref_len - alt_len)
    
    # Get repeat information from INFO
    repeat_count = info_dict.get('REP', 0)
    try:
        repeat_count = int(repeat_count)
    except (ValueError, TypeError):
        repeat_count = 0
    
    # Get context base (first base of ref for deletions, first base of alt for insertions)
    if ref_len > alt_len:  # Deletion
        context_base = ref[0]
    else:  # Insertion
        context_base = alt[0]
    
    if ref_len > alt_len:
        # Deletion
        if indel_len == 1:
            # 1bp deletion
            if context_base in ['C', 'G']:
                return 'DEL.C.1.1'
            else:
                return 'DEL.T.1.1'
        else:
            # >1bp deletion
            if repeat_count > 0:
                # Repeat context
                if indel_len <= 6:
                    if indel_len == 2:
                        return f'DEL.repeats.2.{min(repeat_count, 6)}'
                    elif indel_len == 3:
                        return f'DEL.repeats.3.{min(repeat_count, 6)}'
                    elif indel_len == 4:
                        return f'DEL.repeats.4.{min(repeat_count, 6)}'
                    else:
                        return f'DEL.repeats.5+.{min(repeat_count, 6)}'
                else:
                    return f'DEL.repeats.5+.{min(repeat_count, 6)}'
            else:
                # Microhomology context
                if indel_len <= 6:
                    if indel_len == 2:
                        return 'DEL.MH.2.1'
                    elif indel_len == 3:
                        return f'DEL.MH.3.{min(repeat_count, 2)}'
                    elif indel_len == 4:
                        return f'DEL.MH.4.{min(repeat_count, 3)}'
                    else:
                        return f'DEL.MH.5+.{min(repeat_count, 5)}'
                else:
                    return f'DEL.MH.5+.{min(repeat_count, 5)}'
    else:
        # Insertion
        if indel_len == 1:
            # 1bp insertion
            if context_base in ['C', 'G']:
                return 'INS.C.1.0'
            else:
                return 'INS.T.1.0'
        else:
            # >1bp insertion
            if repeat_count > 0:
                # Repeat context
                if indel_len <= 5:
                    if indel_len == 2:
                        return f'INS.repeats.2.{max(0, repeat_count-1)}'
                    elif indel_len == 3:
                        return f'INS.repeats.3.{max(0, repeat_count-1)}'
                    elif indel_len == 4:
                        return f'INS.repeats.4.{max(0, repeat_count-1)}'
                    else:
                        return f'INS.repeats.5+.{max(0, repeat_count-1)}'
                else:
                    return f'INS.repeats.5+.{max(0, repeat_count-1)}'
            else:
                # Microhomology context
                if indel_len <= 5:
                    if indel_len == 2:
                        return 'INS.MH.2.1'
                    elif indel_len == 3:
                        return f'INS.MH.3.{min(repeat_count, 2)}'
                    elif indel_len == 4:
                        return f'INS.MH.4.{min(repeat_count, 3)}'
                    else:
                        return f'INS.MH.5+.{min(repeat_count, 5)}'
                else:
                    return f'INS.MH.5+.{min(repeat_count, 5)}'

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

def process_vcf_file_simple(vcf_path):
    """
    Process a single VCF file with simplified indel classification.
    
    Args:
        vcf_path (str): Path to VCF file
        
    Returns:
        dict: Counts of indels by ID83 category
    """
    indel_counts = defaultdict(int)
    total_variants = 0
    classified_variants = 0
    
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
            total_variants += 1
            
            # Skip if not PASS filter
            if filter_val != 'PASS':
                continue
            
            # Parse INFO field
            info_dict = parse_info_field(info)
            
            # Classify indel
            try:
                indel_class = classify_indel_simple(ref, alt, info_dict)
                indel_counts[indel_class] += 1
                classified_variants += 1
            except Exception as e:
                print(f"Warning: Could not classify indel {ref}->{alt} in {vcf_path}: {e}")
                continue
    
    print(f"  Total variants: {total_variants}, Classified: {classified_variants}")
    return dict(indel_counts)

def process_all_vcfs_simple(ind_dir="IND"):
    """
    Process all VCF files in the IND directory with simplified classification.
    
    Args:
        ind_dir (str): Directory containing VCF files
        
    Returns:
        dict: Dictionary of indel counts by sample
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
            counts = process_vcf_file_simple(vcf_file)
            all_counts[sample_name] = counts
        except Exception as e:
            print(f"Error processing {vcf_file}: {e}")
            continue
    
    return all_counts

def save_results_simple(all_counts, output_file="volkova_indel_counts_simple.tsv"):
    """Save results in TSV format."""
    
    # Get all unique categories
    all_categories = set()
    for counts in all_counts.values():
        all_categories.update(counts.keys())
    
    # Sort categories to match ID83 order
    sorted_categories = []
    for cat in indel_types_83_str:
        if cat in all_categories:
            sorted_categories.append(cat)
    
    # Add any missing categories
    for cat in sorted(all_categories):
        if cat not in sorted_categories:
            sorted_categories.append(cat)
    
    # Write TSV file
    with open(output_file, 'w') as f:
        # Header
        f.write("Sample\t" + "\t".join(sorted_categories) + "\n")
        
        # Data rows
        for sample, counts in all_counts.items():
            row = [sample]
            for category in sorted_categories:
                row.append(str(counts.get(category, 0)))
            f.write("\t".join(row) + "\n")
    
    print(f"Saved results to {output_file}")

def create_summary_report_simple(all_counts):
    """Create a comprehensive summary report of the indel analysis."""
    print("\n" + "="*80)
    print("VOLKOVA ET AL. (2020) INDEL ANALYSIS SUMMARY")
    print("="*80)
    
    print(f"\nSamples processed: {len(all_counts)}")
    
    # Total counts
    total_indels = sum(sum(counts.values()) for counts in all_counts.values())
    print(f"\nTotal indels across all samples: {total_indels:,}")
    
    # Per-sample statistics
    indels_per_sample = {sample: sum(counts.values()) for sample, counts in all_counts.items()}
    values = list(indels_per_sample.values())
    print(f"Mean indels per sample: {sum(values)/len(values):.1f}")
    print(f"Range: {min(values):.0f} - {max(values):.0f}")
    
    # Top categories across all samples
    category_totals = defaultdict(int)
    for counts in all_counts.values():
        for category, count in counts.items():
            category_totals[category] += count
    
    # Sort by total count
    sorted_categories = sorted(category_totals.items(), key=lambda x: x[1], reverse=True)
    
    print(f"\nTop 15 indel categories:")
    for i, (category, count) in enumerate(sorted_categories[:15], 1):
        percentage = (count / total_indels) * 100
        print(f"  {i:2d}. {category:<20} {count:>6,} ({percentage:5.1f}%)")
    
    # Sample distribution
    print(f"\nSample indel counts:")
    sorted_samples = sorted(indels_per_sample.items(), key=lambda x: x[1], reverse=True)
    for sample, count in sorted_samples[:10]:
        print(f"  {sample}: {count:>4,}")
    
    return {
        'total_indels': total_indels,
        'samples_processed': len(all_counts),
        'top_categories': sorted_categories[:15],
        'indels_per_sample': indels_per_sample
    }

def main():
    """Main function to process all VCF files."""
    print("Processing Volkova et al. (2020) VCF files with simplified classification...")
    
    try:
        # Process all VCF files
        all_counts = process_all_vcfs_simple()
        
        # Create summary report
        summary = create_summary_report_simple(all_counts)
        
        # Save results
        save_results_simple(all_counts)
        
        print(f"\nAnalysis complete!")
        print(f"Processed {summary['samples_processed']} samples")
        print(f"Total indels: {summary['total_indels']:,}")
        
        return all_counts, summary
        
    except Exception as e:
        print(f"Error: {e}")
        return None, None

if __name__ == "__main__":
    all_counts, summary_result = main()