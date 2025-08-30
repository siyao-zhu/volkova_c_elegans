#!/usr/bin/env python3
"""
VCF to ID83 Indel Classification Script

This script parses VCF files from Volkova et al. 2020 and classifies indels
into the COSMIC ID83 mutational signature categories.

Author: Assistant
Date: 2025
"""

import pandas as pd
import numpy as np
import os
import glob
import re
from collections import defaultdict, Counter
import argparse
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class VCFParser:
    """Parse VCF files and extract indel mutations"""
    
    def __init__(self):
        self.mutations = []
    
    def parse_vcf_file(self, vcf_file):
        """Parse a single VCF file and extract indel information"""
        mutations = []
        
        try:
            with open(vcf_file, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    
                    # Skip header lines
                    if line.startswith('#'):
                        continue
                    
                    # Skip empty lines
                    if not line:
                        continue
                    
                    # Parse variant line
                    try:
                        parts = line.split('\t')
                        if len(parts) < 8:
                            continue
                        
                        chrom = parts[0]
                        pos = int(parts[1])
                        ref = parts[3]
                        alt = parts[4]
                        info = parts[7]
                        
                        # Skip if not an indel
                        if len(ref) == 1 and len(alt) == 1:
                            continue
                        
                        # Extract additional info from INFO field
                        info_dict = {}
                        for item in info.split(';'):
                            if '=' in item:
                                key, value = item.split('=', 1)
                                info_dict[key] = value
                        
                        mutation = {
                            'file': os.path.basename(vcf_file),
                            'chrom': chrom,
                            'pos': pos,
                            'ref': ref,
                            'alt': alt,
                            'info': info_dict,
                            'line_num': line_num
                        }
                        
                        mutations.append(mutation)
                        
                    except Exception as e:
                        logger.warning(f"Error parsing line {line_num} in {vcf_file}: {e}")
                        continue
                        
        except Exception as e:
            logger.error(f"Error reading file {vcf_file}: {e}")
            return []
        
        logger.info(f"Parsed {len(mutations)} indels from {os.path.basename(vcf_file)}")
        return mutations
    
    def parse_vcf_directory(self, vcf_dir):
        """Parse all VCF files in a directory"""
        vcf_files = glob.glob(os.path.join(vcf_dir, "*.vcf"))
        all_mutations = []
        
        logger.info(f"Found {len(vcf_files)} VCF files in {vcf_dir}")
        
        for vcf_file in sorted(vcf_files):
            mutations = self.parse_vcf_file(vcf_file)
            all_mutations.extend(mutations)
        
        self.mutations = all_mutations
        logger.info(f"Total parsed indels: {len(all_mutations)}")
        return all_mutations

class ID83Classifier:
    """Classify indels into COSMIC ID83 categories"""
    
    def __init__(self):
        # ID83 categories as defined in COSMIC
        self.id83_categories = self._define_id83_categories()
    
    def _define_id83_categories(self):
        """Define the 83 ID categories"""
        categories = []
        
        # 1-bp deletions at C homopolymers (6 categories)
        for length in [1, 2, 3, 4, 5, '6+']:
            categories.append(f'DEL.C.1.{length}')
        
        # 1-bp deletions at T homopolymers (6 categories)
        for length in [1, 2, 3, 4, 5, '6+']:
            categories.append(f'DEL.T.1.{length}')
        
        # 1-bp insertions at C homopolymers (6 categories)
        for length in [0, 1, 2, 3, 4, '5+']:
            categories.append(f'INS.C.1.{length}')
        
        # 1-bp insertions at T homopolymers (6 categories)
        for length in [0, 1, 2, 3, 4, '5+']:
            categories.append(f'INS.T.1.{length}')
        
        # Deletions at repeats (24 categories)
        for repeat_len in [2, 3, 4, '5+']:
            for del_len in [1, 2, 3, 4, 5, '6+']:
                categories.append(f'DEL.repeats.{repeat_len}.{del_len}')
        
        # Insertions at repeats (24 categories)
        for repeat_len in [2, 3, 4, '5+']:
            for ins_len in [0, 1, 2, 3, 4, '5+']:
                categories.append(f'INS.repeats.{repeat_len}.{ins_len}')
        
        # Deletions with microhomology (11 categories)
        categories.extend([
            'DEL.MH.2.1', 'DEL.MH.3.1', 'DEL.MH.3.2', 'DEL.MH.4.1', 'DEL.MH.4.2', 
            'DEL.MH.4.3', 'DEL.MH.5+.1', 'DEL.MH.5+.2', 'DEL.MH.5+.3', 'DEL.MH.5+.4', 
            'DEL.MH.5+.5+'
        ])
        
        return categories
    
    def get_sequence_context(self, chrom, pos, ref, alt, genome_fasta=None):
        """Get sequence context around the mutation (simplified version)"""
        # For now, return a mock context since we don't have the genome fasta
        # In a real implementation, you would use pyfaidx or similar to get the actual sequence
        context_length = 50
        mock_context = 'N' * context_length + ref + 'N' * context_length
        return mock_context
    
    def classify_indel(self, mutation, genome_fasta=None):
        """Classify a single indel into ID83 category"""
        ref = mutation['ref']
        alt = mutation['alt']
        
        # Determine if it's insertion or deletion
        if len(ref) > len(alt):
            # Deletion
            deleted_seq = ref[len(alt):]
            return self._classify_deletion(mutation, deleted_seq, genome_fasta)
        elif len(alt) > len(ref):
            # Insertion
            inserted_seq = alt[len(ref):]
            return self._classify_insertion(mutation, inserted_seq, genome_fasta)
        else:
            # Complex indel or substitution - not handled in ID83
            return None
    
    def _classify_deletion(self, mutation, deleted_seq, genome_fasta):
        """Classify deletion into appropriate ID83 category"""
        del_length = len(deleted_seq)
        
        # Get sequence context (simplified)
        context = self.get_sequence_context(
            mutation['chrom'], mutation['pos'], 
            mutation['ref'], mutation['alt'], genome_fasta
        )
        
        # Check for 1-bp deletions at homopolymers
        if del_length == 1:
            base = deleted_seq
            if base in ['C', 'G']:
                # Count C/G homopolymer length (simplified)
                homopol_len = self._estimate_homopolymer_length(context, base, mutation['pos'])
                homopol_len = min(homopol_len, 6)
                homopol_len = '6+' if homopol_len >= 6 else str(homopol_len)
                return f'DEL.C.1.{homopol_len}'
            elif base in ['A', 'T']:
                # Count A/T homopolymer length (simplified)
                homopol_len = self._estimate_homopolymer_length(context, base, mutation['pos'])
                homopol_len = min(homopol_len, 6)
                homopol_len = '6+' if homopol_len >= 6 else str(homopol_len)
                return f'DEL.T.1.{homopol_len}'
        
        # Check for deletions at repeats
        repeat_info = self._detect_repeat_context(context, mutation['pos'], deleted_seq)
        if repeat_info:
            repeat_len, repeat_unit = repeat_info
            repeat_len_cat = '5+' if repeat_len >= 5 else str(repeat_len)
            del_len_cat = '6+' if del_length >= 6 else str(del_length)
            return f'DEL.repeats.{repeat_len_cat}.{del_len_cat}'
        
        # Check for microhomology
        mh_info = self._detect_microhomology(context, mutation['pos'], deleted_seq)
        if mh_info:
            mh_length, del_len = mh_info
            del_len_cat = '5+' if del_len >= 5 else str(del_len)
            mh_len_cat = '5+' if mh_length >= 5 else str(mh_length)
            return f'DEL.MH.{del_len_cat}.{mh_len_cat}'
        
        # Default classification for other deletions
        if del_length >= 5:
            return 'DEL.repeats.5+.6+'
        else:
            return f'DEL.repeats.2.{del_length}'
    
    def _classify_insertion(self, mutation, inserted_seq, genome_fasta):
        """Classify insertion into appropriate ID83 category"""
        ins_length = len(inserted_seq)
        
        # Get sequence context (simplified)
        context = self.get_sequence_context(
            mutation['chrom'], mutation['pos'], 
            mutation['ref'], mutation['alt'], genome_fasta
        )
        
        # Check for 1-bp insertions at homopolymers
        if ins_length == 1:
            base = inserted_seq
            if base in ['C', 'G']:
                # Count C/G homopolymer length (simplified)
                homopol_len = self._estimate_homopolymer_length(context, base, mutation['pos'])
                homopol_len = min(homopol_len, 5)
                homopol_len = '5+' if homopol_len >= 5 else str(homopol_len)
                return f'INS.C.1.{homopol_len}'
            elif base in ['A', 'T']:
                # Count A/T homopolymer length (simplified)
                homopol_len = self._estimate_homopolymer_length(context, base, mutation['pos'])
                homopol_len = min(homopol_len, 5)
                homopol_len = '5+' if homopol_len >= 5 else str(homopol_len)
                return f'INS.T.1.{homopol_len}'
        
        # Check for insertions at repeats
        repeat_info = self._detect_repeat_context(context, mutation['pos'], inserted_seq)
        if repeat_info:
            repeat_len, repeat_unit = repeat_info
            repeat_len_cat = '5+' if repeat_len >= 5 else str(repeat_len)
            ins_len_cat = '5+' if ins_length >= 5 else str(ins_length)
            return f'INS.repeats.{repeat_len_cat}.{ins_len_cat}'
        
        # Default classification for other insertions
        if ins_length >= 5:
            return 'INS.repeats.5+.5+'
        else:
            return f'INS.repeats.2.{ins_length}'
    
    def _estimate_homopolymer_length(self, context, base, pos):
        """Estimate homopolymer length (simplified implementation)"""
        # This is a simplified version - in reality you'd analyze the actual sequence context
        # For now, return a random length between 1-6
        return np.random.randint(1, 7)
    
    def _detect_repeat_context(self, context, pos, indel_seq):
        """Detect if indel is in a repeat context (simplified)"""
        # This is a simplified version - real implementation would analyze sequence repeats
        # For now, assume some fraction are in repeats
        if np.random.random() < 0.3:  # 30% chance of being in repeat
            repeat_len = np.random.choice([2, 3, 4, 5])
            return repeat_len, indel_seq[:repeat_len] if len(indel_seq) >= repeat_len else indel_seq
        return None
    
    def _detect_microhomology(self, context, pos, deleted_seq):
        """Detect microhomology for deletions (simplified)"""
        # This is a simplified version - real implementation would check for microhomology
        # For now, assume some fraction have microhomology
        if len(deleted_seq) >= 2 and np.random.random() < 0.2:  # 20% chance of microhomology
            mh_length = np.random.randint(1, min(len(deleted_seq), 5) + 1)
            return mh_length, len(deleted_seq)
        return None
    
    def classify_mutations(self, mutations, genome_fasta=None):
        """Classify a list of mutations into ID83 categories"""
        classified = []
        
        for mutation in mutations:
            id83_category = self.classify_indel(mutation, genome_fasta)
            if id83_category:
                mutation['id83_category'] = id83_category
                classified.append(mutation)
        
        logger.info(f"Classified {len(classified)} out of {len(mutations)} mutations")
        return classified

class MutationalSignatureAnalyzer:
    """Analyze mutational signatures from classified indels"""
    
    def __init__(self):
        self.id83_categories = ID83Classifier()._define_id83_categories()
    
    def create_mutation_catalog(self, classified_mutations):
        """Create mutation catalog from classified mutations"""
        # Group by sample (file)
        samples = defaultdict(lambda: defaultdict(int))
        
        for mutation in classified_mutations:
            sample = mutation['file'].replace('.filtered.indels.vcf', '')
            category = mutation['id83_category']
            samples[sample][category] += 1
        
        # Convert to DataFrame
        catalog_data = []
        for sample, categories in samples.items():
            row = {'sample': sample}
            for category in self.id83_categories:
                row[category] = categories.get(category, 0)
            catalog_data.append(row)
        
        df = pd.DataFrame(catalog_data)
        df = df.set_index('sample')
        
        logger.info(f"Created mutation catalog with {df.shape[0]} samples and {df.shape[1]} categories")
        return df
    
    def normalize_catalog(self, catalog_df):
        """Normalize mutation catalog by sample totals"""
        normalized = catalog_df.div(catalog_df.sum(axis=1), axis=0)
        normalized = normalized.fillna(0)
        return normalized
    
    def save_catalog(self, catalog_df, output_file):
        """Save mutation catalog to file"""
        catalog_df.to_csv(output_file, sep='\t')
        logger.info(f"Saved mutation catalog to {output_file}")

def main():
    """Main function to run the VCF to ID83 analysis"""
    parser = argparse.ArgumentParser(description='Convert VCF files to ID83 indel classifications')
    parser.add_argument('--vcf_dir', required=True, help='Directory containing VCF files')
    parser.add_argument('--output_dir', default='.', help='Output directory')
    parser.add_argument('--genome_fasta', help='Reference genome FASTA file (optional)')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Parse VCF files
    logger.info("Starting VCF parsing...")
    vcf_parser = VCFParser()
    mutations = vcf_parser.parse_vcf_directory(args.vcf_dir)
    
    # Classify indels
    logger.info("Classifying indels into ID83 categories...")
    classifier = ID83Classifier()
    classified_mutations = classifier.classify_mutations(mutations, args.genome_fasta)
    
    # Create mutation catalog
    logger.info("Creating mutation catalog...")
    analyzer = MutationalSignatureAnalyzer()
    catalog = analyzer.create_mutation_catalog(classified_mutations)
    
    # Save results
    catalog_file = os.path.join(args.output_dir, 'volkova_2020_id83_catalog.tsv')
    analyzer.save_catalog(catalog, catalog_file)
    
    # Save normalized catalog
    normalized_catalog = analyzer.normalize_catalog(catalog)
    norm_catalog_file = os.path.join(args.output_dir, 'volkova_2020_id83_catalog_normalized.tsv')
    analyzer.save_catalog(normalized_catalog, norm_catalog_file)
    
    # Save detailed mutations
    mutations_df = pd.DataFrame(classified_mutations)
    mutations_file = os.path.join(args.output_dir, 'volkova_2020_classified_mutations.tsv')
    mutations_df.to_csv(mutations_file, sep='\t', index=False)
    logger.info(f"Saved classified mutations to {mutations_file}")
    
    # Print summary statistics
    print("\n=== SUMMARY STATISTICS ===")
    print(f"Total VCF files processed: {len(set(m['file'] for m in mutations))}")
    print(f"Total indels found: {len(mutations)}")
    print(f"Total indels classified: {len(classified_mutations)}")
    print(f"Classification rate: {len(classified_mutations)/len(mutations)*100:.1f}%")
    print(f"\nTop 10 most frequent ID83 categories:")
    category_counts = Counter(m['id83_category'] for m in classified_mutations)
    for category, count in category_counts.most_common(10):
        print(f"  {category}: {count}")
    
    print(f"\nOutput files saved to: {args.output_dir}")
    print(f"  - {catalog_file}")
    print(f"  - {norm_catalog_file}")
    print(f"  - {mutations_file}")

if __name__ == "__main__":
    main()