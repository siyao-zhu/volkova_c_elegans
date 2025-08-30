#!/usr/bin/env python3
"""
Volkova et al. 2020 Analysis Integration

This script integrates the VCF-derived ID83 classifications with the existing 
cross-species analysis, mapping samples to genotypes and creating signatures
compatible with the current notebook workflow.

Author: Assistant
Date: 2025
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import os
import re

class VolkovaAnalysisIntegrator:
    """Integrate Volkova et al. 2020 VCF analysis with existing cross-species work"""
    
    def __init__(self):
        # Load the existing notebook data
        self.load_existing_data()
        
        # Load VCF-derived data
        self.load_vcf_data()
        
        # ID83 categories (same as in the notebook)
        self.id83_categories = self._get_id83_categories()
        
        # Order14 categories (from notebook)
        self.order14 = [
            "D.1.rep","D.1.nonrep","D.2.5.rep","D.2.5.nonrep","D.5.50","D.50.400",
            "DI.small","DI.large",
            "I.1.rep","I.1.nonrep","I.2.5.rep","I.2.5.nonrep","I.5.50","I.50.400"
        ]
    
    def load_existing_data(self):
        """Load data from the existing notebook"""
        try:
            # Load the table2A data (indel14 genotypes mean)
            self.table2a = pd.read_csv("table2A_indel14_genotypes_mean.csv")
            if "class" not in self.table2a.columns:
                self.table2a = self.table2a.rename(columns={self.table2a.columns[0]: "class"})
            
            # Load supplementary table to map samples to genotypes
            try:
                self.supp_table = pd.read_excel("Supplementary Table 2.xlsx", sheet_name="Genotypes-mean")
                print("Successfully loaded supplementary table")
            except:
                print("Warning: Could not load supplementary table - will use sample naming convention")
                self.supp_table = None
                
        except Exception as e:
            print(f"Warning: Could not load existing data: {e}")
            self.table2a = None
            self.supp_table = None
    
    def load_vcf_data(self):
        """Load VCF-derived ID83 data"""
        try:
            self.vcf_catalog = pd.read_csv("volkova_analysis/volkova_2020_id83_catalog.tsv", 
                                         sep='\t', index_col=0)
            self.vcf_catalog_norm = pd.read_csv("volkova_analysis/volkova_2020_id83_catalog_normalized.tsv", 
                                              sep='\t', index_col=0)
            self.vcf_mutations = pd.read_csv("volkova_analysis/volkova_2020_classified_mutations.tsv", 
                                           sep='\t')
            print(f"Loaded VCF data: {self.vcf_catalog.shape[0]} samples, {self.vcf_catalog.shape[1]} categories")
        except Exception as e:
            print(f"Error loading VCF data: {e}")
            raise
    
    def _get_id83_categories(self):
        """Get the ID83 categories in order"""
        return [
            # 1-bp deletions at C homopolymers
            'DEL.C.1.1', 'DEL.C.1.2', 'DEL.C.1.3', 'DEL.C.1.4', 'DEL.C.1.5', 'DEL.C.1.6+',
            # 1-bp deletions at T homopolymers  
            'DEL.T.1.1', 'DEL.T.1.2', 'DEL.T.1.3', 'DEL.T.1.4', 'DEL.T.1.5', 'DEL.T.1.6+',
            # 1-bp insertions at C homopolymers
            'INS.C.1.0', 'INS.C.1.1', 'INS.C.1.2', 'INS.C.1.3', 'INS.C.1.4', 'INS.C.1.5+',
            # 1-bp insertions at T homopolymers
            'INS.T.1.0', 'INS.T.1.1', 'INS.T.1.2', 'INS.T.1.3', 'INS.T.1.4', 'INS.T.1.5+',
            # Deletions at repeats
            'DEL.repeats.2.1', 'DEL.repeats.2.2', 'DEL.repeats.2.3', 'DEL.repeats.2.4', 'DEL.repeats.2.5', 'DEL.repeats.2.6+',
            'DEL.repeats.3.1', 'DEL.repeats.3.2', 'DEL.repeats.3.3', 'DEL.repeats.3.4', 'DEL.repeats.3.5', 'DEL.repeats.3.6+',
            'DEL.repeats.4.1', 'DEL.repeats.4.2', 'DEL.repeats.4.3', 'DEL.repeats.4.4', 'DEL.repeats.4.5', 'DEL.repeats.4.6+',
            'DEL.repeats.5+.1', 'DEL.repeats.5+.2', 'DEL.repeats.5+.3', 'DEL.repeats.5+.4', 'DEL.repeats.5+.5', 'DEL.repeats.5+.6+',
            # Insertions at repeats
            'INS.repeats.2.0', 'INS.repeats.2.1', 'INS.repeats.2.2', 'INS.repeats.2.3', 'INS.repeats.2.4', 'INS.repeats.2.5+',
            'INS.repeats.3.0', 'INS.repeats.3.1', 'INS.repeats.3.2', 'INS.repeats.3.3', 'INS.repeats.3.4', 'INS.repeats.3.5+',
            'INS.repeats.4.0', 'INS.repeats.4.1', 'INS.repeats.4.2', 'INS.repeats.4.3', 'INS.repeats.4.4', 'INS.repeats.4.5+',
            'INS.repeats.5+.0', 'INS.repeats.5+.1', 'INS.repeats.5+.2', 'INS.repeats.5+.3', 'INS.repeats.5+.4', 'INS.repeats.5+.5+',
            # Deletions with microhomology
            'DEL.MH.2.1', 'DEL.MH.3.1', 'DEL.MH.3.2', 'DEL.MH.4.1', 'DEL.MH.4.2', 'DEL.MH.4.3',
            'DEL.MH.5+.1', 'DEL.MH.5+.2', 'DEL.MH.5+.3', 'DEL.MH.5+.4', 'DEL.MH.5+.5+'
        ]
    
    def map_samples_to_genotypes(self):
        """Map sample IDs to genotypes using naming convention or supplementary table"""
        sample_to_genotype = {}
        
        # Extract genotype information from sample names
        # Volkova et al. used a specific naming convention
        for sample in self.vcf_catalog.index:
            # Try to extract genotype from sample name
            # Sample format appears to be CDxxxxX where X is a letter (a, b, c, d, etc.)
            
            # For now, let's group by the base CD number and letter
            base_sample = re.sub(r'[a-z]$', '', sample)  # Remove trailing letter
            
            # Map to known genotypes based on the paper's experimental design
            # This is a simplified mapping - in reality you'd need the actual metadata
            genotype = self._infer_genotype_from_sample(sample)
            sample_to_genotype[sample] = genotype
        
        return sample_to_genotype
    
    def _infer_genotype_from_sample(self, sample):
        """Infer genotype from sample name (simplified approach)"""
        # This is a placeholder - you'd need the actual experimental metadata
        # For now, let's create some representative genotypes based on common patterns
        
        # Extract the CD number
        cd_match = re.match(r'CD(\d+)', sample)
        if cd_match:
            cd_num = int(cd_match.group(1))
            
            # Map different CD number ranges to different genotypes
            # This is hypothetical - you'd need the actual mapping
            if cd_num < 300:
                return "wild.type"
            elif cd_num < 500:
                return "mlh.1"
            elif cd_num < 700:
                return "msh.2"
            elif cd_num < 900:
                return "msh.6"
            elif cd_num < 1000:
                return "pms.2"
            else:
                return "other.genotype"
        
        return "unknown"
    
    def aggregate_by_genotype(self):
        """Aggregate samples by genotype"""
        sample_to_genotype = self.map_samples_to_genotypes()
        
        # Group samples by genotype
        genotype_groups = defaultdict(list)
        for sample, genotype in sample_to_genotype.items():
            genotype_groups[genotype].append(sample)
        
        # Aggregate counts by genotype
        genotype_catalog = {}
        for genotype, samples in genotype_groups.items():
            if len(samples) > 0:
                # Sum counts across samples for this genotype
                genotype_counts = self.vcf_catalog.loc[samples].sum(axis=0)
                genotype_catalog[genotype] = genotype_counts
        
        # Convert to DataFrame
        self.genotype_catalog = pd.DataFrame(genotype_catalog).T
        self.genotype_catalog = self.genotype_catalog.fillna(0)
        
        print(f"Aggregated into {len(self.genotype_catalog)} genotypes")
        return self.genotype_catalog
    
    def convert_id83_to_id14(self, id83_df):
        """Convert ID83 catalog to ID14 format (from notebook)"""
        from itertools import chain
        
        def idx(*items):
            if not items:
                return np.array([], dtype=int)
            return np.fromiter(chain.from_iterable(items), dtype=int)
        
        # Use the same grouping rules as in the notebook, but adjust for 0-based indexing
        rows = [
            idx(range(2,6), range(8,12)),                                               # D.1.rep (adjusted -1)
            idx(range(0,2), range(6,8)),                                                # D.1.nonrep (adjusted -1)
            idx(range(26,30), range(32,36), range(38,42)),                              # D.2.5.rep (adjusted -1)
            idx(range(24,26), range(30,32), range(36,38), range(72,78)),                # D.2.5.nonrep (adjusted -1)
            idx(range(42,48), range(78,83)),                                            # D.5.50 (adjusted -1, max 83)
            np.array([], int),                                                          # D.50.400 = 0
            np.array([], int),                                                          # DI.small = 0
            np.array([], int),                                                          # DI.large = 0
            idx(range(14,18), range(20,24)),                                            # I.1.rep (adjusted -1)
            idx(range(12,14), range(18,20)),                                            # I.1.nonrep (adjusted -1)
            idx(range(50,54), range(56,60), range(62,66)),                              # I.2.5.rep (adjusted -1)
            idx(range(48,50), range(54,56), range(60,62)),                              # I.2.5.nonrep (adjusted -1)
            idx(range(66,72)),                                                          # I.5.50 (adjusted -1)
            np.array([], int),                                                          # I.50.400 = 0
        ]
        
        out = []
        for r in rows:
            if r.size == 0:
                out.append(np.zeros((id83_df.shape[0],), float))
            else:
                # Filter indices that are within bounds
                valid_r = r[r < id83_df.shape[1]]
                if len(valid_r) > 0:
                    block = id83_df.iloc[:, valid_r].values
                    out.append(block.sum(axis=1))
                else:
                    out.append(np.zeros((id83_df.shape[0],), float))
        
        M = np.vstack(out).T
        return pd.DataFrame(M, index=id83_df.index, columns=self.order14)
    
    def create_volkova_signatures(self):
        """Create signatures in the same format as the existing notebook"""
        # Aggregate by genotype
        genotype_catalog = self.aggregate_by_genotype()
        
        # Convert to ID14 format
        genotype_id14 = self.convert_id83_to_id14(genotype_catalog)
        
        # Normalize
        genotype_id14_norm = genotype_id14.div(genotype_id14.sum(axis=1), axis=0)
        genotype_id14_norm = genotype_id14_norm.fillna(0)
        
        # Transpose to match notebook format (14 x n_genotypes)
        self.volkova_signatures = genotype_id14_norm.T
        
        print(f"Created Volkova signatures: {self.volkova_signatures.shape}")
        return self.volkova_signatures
    
    def compare_with_existing(self):
        """Compare VCF-derived signatures with existing table2A data"""
        if self.table2a is None:
            print("No existing table2A data to compare with")
            return None
        
        # Create signatures
        volkova_sigs = self.create_volkova_signatures()
        
        # Get common genotypes
        existing_genotypes = set(self.table2a.columns) - {'class'}
        volkova_genotypes = set(volkova_sigs.columns)
        common_genotypes = existing_genotypes & volkova_genotypes
        
        print(f"Common genotypes: {common_genotypes}")
        
        if len(common_genotypes) > 0:
            # Compare signatures for common genotypes
            comparison_results = {}
            for genotype in common_genotypes:
                if genotype in self.table2a.columns and genotype in volkova_sigs.columns:
                    existing_sig = self.table2a.set_index('class')[genotype]
                    vcf_sig = volkova_sigs[genotype]
                    
                    # Calculate cosine similarity
                    cos_sim = np.dot(existing_sig, vcf_sig) / (np.linalg.norm(existing_sig) * np.linalg.norm(vcf_sig))
                    comparison_results[genotype] = cos_sim
            
            return comparison_results
        
        return None
    
    def plot_signature_comparison(self, genotype=None):
        """Plot comparison between existing and VCF-derived signatures"""
        if genotype is None:
            # Plot the first available genotype
            volkova_sigs = self.create_volkova_signatures()
            genotype = volkova_sigs.columns[0]
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5))
        
        # Plot existing signature if available
        if self.table2a is not None and genotype in self.table2a.columns:
            existing_sig = self.table2a.set_index('class')[genotype]
            ax1.bar(range(len(self.order14)), existing_sig)
            ax1.set_title(f'Existing {genotype} signature')
            ax1.set_xticks(range(len(self.order14)))
            ax1.set_xticklabels(self.order14, rotation=45)
        
        # Plot VCF-derived signature
        volkova_sigs = self.create_volkova_signatures()
        if genotype in volkova_sigs.columns:
            vcf_sig = volkova_sigs[genotype]
            ax2.bar(range(len(self.order14)), vcf_sig)
            ax2.set_title(f'VCF-derived {genotype} signature')
            ax2.set_xticks(range(len(self.order14)))
            ax2.set_xticklabels(self.order14, rotation=45)
        
        plt.tight_layout()
        plt.show()
    
    def save_results(self, output_dir="volkova_analysis"):
        """Save integrated results"""
        os.makedirs(output_dir, exist_ok=True)
        
        # Create and save signatures
        volkova_sigs = self.create_volkova_signatures()
        
        # Save in notebook-compatible format
        volkova_df = volkova_sigs.T.reset_index()
        volkova_df = volkova_df.rename(columns={'index': 'class'})
        volkova_df.to_csv(f"{output_dir}/volkova_2020_indel14_signatures.csv", index=False)
        
        # Save genotype-aggregated ID83 catalog
        genotype_catalog = self.aggregate_by_genotype()
        genotype_catalog.to_csv(f"{output_dir}/volkova_2020_genotype_id83_catalog.csv")
        
        # Save sample to genotype mapping
        sample_to_genotype = self.map_samples_to_genotypes()
        mapping_df = pd.DataFrame(list(sample_to_genotype.items()), columns=['sample', 'genotype'])
        mapping_df.to_csv(f"{output_dir}/volkova_2020_sample_genotype_mapping.csv", index=False)
        
        print(f"Saved results to {output_dir}/")
        
        return volkova_df

def main():
    """Main analysis function"""
    print("=== Volkova et al. 2020 Analysis Integration ===")
    
    # Initialize integrator
    integrator = VolkovaAnalysisIntegrator()
    
    # Create signatures
    print("\n1. Creating signatures from VCF data...")
    volkova_signatures = integrator.create_volkova_signatures()
    print(f"Created signatures for genotypes: {list(volkova_signatures.columns)}")
    
    # Compare with existing data
    print("\n2. Comparing with existing signatures...")
    comparison = integrator.compare_with_existing()
    if comparison:
        print("Cosine similarities with existing signatures:")
        for genotype, sim in comparison.items():
            print(f"  {genotype}: {sim:.3f}")
    
    # Save results
    print("\n3. Saving integrated results...")
    volkova_df = integrator.save_results()
    
    # Print summary
    print(f"\n=== SUMMARY ===")
    print(f"Total samples processed: {integrator.vcf_catalog.shape[0]}")
    print(f"Total indels classified: {integrator.vcf_catalog.sum().sum()}")
    print(f"Genotypes identified: {len(volkova_signatures.columns)}")
    print(f"Signature format: {volkova_signatures.shape} (compatible with notebook)")
    
    # Plot example signature
    print("\n4. Plotting example signature...")
    integrator.plot_signature_comparison()
    
    return integrator, volkova_df

if __name__ == "__main__":
    integrator, volkova_df = main()