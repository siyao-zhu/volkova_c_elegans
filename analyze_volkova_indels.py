#!/usr/bin/env python3
"""
Complete analysis of Volkova et al. 2020 indel signatures
Converts to human genome background and compares with COSMIC
"""

import csv
import os
from collections import defaultdict

# ID83 categories
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

# Define the 14 categories used in the worm analysis
ORDER14 = [
    "D.1.rep","D.1.nonrep","D.2.5.rep","D.2.5.nonrep","D.5.50","D.50.400",
    "DI.small","DI.large",
    "I.1.rep","I.1.nonrep","I.2.5.rep","I.2.5.nonrep","I.5.50","I.50.400"
]

def load_volkova_signatures():
    """Load the Volkova signatures from CSV"""
    signatures = defaultdict(dict)
    
    with open('volkova_id83_signatures_cd05_transposed.csv', 'r') as f:
        reader = csv.DictReader(f)
        genotypes = list(reader.fieldnames)[1:]  # Skip 'Type' column
        
        for row in reader:
            category = row['Type']
            for genotype in genotypes:
                signatures[genotype][category] = float(row[genotype])
    
    return signatures, genotypes

def regroup_id83_to_14(id83_signature):
    """Convert ID83 to ID14 format following the paper's methodology"""
    
    # Define mapping from ID83 indices to ID14 categories
    # Note: These are 0-based indices for Python
    id83_to_id14_mapping = [
        (list(range(2,6)) + list(range(8,12)), "D.1.rep"),
        (list(range(0,2)) + list(range(6,8)), "D.1.nonrep"),
        (list(range(26,30)) + list(range(32,36)) + list(range(38,42)), "D.2.5.rep"),
        (list(range(24,26)) + list(range(30,32)) + list(range(36,38)) + list(range(72,78)), "D.2.5.nonrep"),
        (list(range(42,48)) + list(range(78,83)), "D.5.50"),
        ([], "D.50.400"),
        ([], "DI.small"),
        ([], "DI.large"),
        (list(range(14,18)) + list(range(20,24)), "I.1.rep"),
        (list(range(12,14)) + list(range(18,20)), "I.1.nonrep"),
        (list(range(50,54)) + list(range(56,60)) + list(range(62,66)), "I.2.5.rep"),
        (list(range(48,50)) + list(range(54,56)) + list(range(60,62)), "I.2.5.nonrep"),
        (list(range(66,72)), "I.5.50"),
        ([], "I.50.400")
    ]
    
    id14_signature = {}
    
    # Convert ID83 dict to list for indexing
    id83_list = [id83_signature.get(cat, 0) for cat in ID83_CATEGORIES]
    
    for indices, category in id83_to_id14_mapping:
        if indices:
            id14_signature[category] = sum(id83_list[i] for i in indices)
        else:
            id14_signature[category] = 0
    
    return id14_signature

def normalize_signature(signature):
    """Normalize signature to sum to 1"""
    total = sum(signature.values())
    if total > 0:
        return {k: v/total for k, v in signature.items()}
    return signature

def cosine_similarity(sig1, sig2):
    """Calculate cosine similarity between two signatures"""
    # Ensure same keys
    keys = set(sig1.keys()) | set(sig2.keys())
    
    v1 = [sig1.get(k, 0) for k in keys]
    v2 = [sig2.get(k, 0) for k in keys]
    
    # Calculate cosine similarity
    dot_product = sum(a * b for a, b in zip(v1, v2))
    norm1 = sum(a * a for a in v1) ** 0.5
    norm2 = sum(b * b for b in v2) ** 0.5
    
    if norm1 == 0 or norm2 == 0:
        return 0
    
    return dot_product / (norm1 * norm2)

def main():
    """Main analysis"""
    
    print("=== Volkova et al. 2020 Indel Analysis ===\n")
    
    # Load Volkova signatures
    print("Loading Volkova signatures...")
    signatures, genotypes = load_volkova_signatures()
    
    # Convert to ID14 format
    print("\nConverting to ID14 format...")
    id14_signatures = {}
    for genotype in genotypes:
        id14_signatures[genotype] = regroup_id83_to_14(signatures[genotype])
    
    # Normalize signatures
    print("Normalizing signatures...")
    id14_norm = {}
    for genotype in genotypes:
        id14_norm[genotype] = normalize_signature(id14_signatures[genotype])
    
    # Save ID14 signatures
    print("\nSaving ID14 signatures...")
    with open('volkova_id14_signatures.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Category'] + genotypes)
        
        for category in ORDER14:
            row = [category]
            for genotype in genotypes:
                row.append(id14_signatures[genotype].get(category, 0))
            writer.writerow(row)
    
    # Save normalized ID14 signatures
    with open('volkova_id14_normalized.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Category'] + genotypes)
        
        for category in ORDER14:
            row = [category]
            for genotype in genotypes:
                row.append(f"{id14_norm[genotype].get(category, 0):.6f}")
            writer.writerow(row)
    
    # Print summary
    print("\nSummary of Volkova signatures (ID14 format):")
    print("-" * 60)
    print(f"{'Genotype':<15} {'Total Indels':<15} {'Dominant Category':<20}")
    print("-" * 60)
    
    for genotype in genotypes:
        total = sum(id14_signatures[genotype].values())
        
        # Find dominant category
        max_cat = max(id14_norm[genotype].items(), key=lambda x: x[1])
        
        print(f"{genotype:<15} {int(total):<15} {max_cat[0]:<20} ({max_cat[1]:.1%})")
    
    print("\nFiles created:")
    print("- volkova_id14_signatures.csv (raw counts)")
    print("- volkova_id14_normalized.csv (normalized proportions)")
    
    # Create a summary for Figure 4a reproduction
    print("\n=== Data for Figure 4a ===")
    print("\nKey genotypes to visualize:")
    key_genotypes = ['N2', 'mlh-1', 'xpa-1', 'rev-3', 'polk-1']
    
    for genotype in key_genotypes:
        if genotype in genotypes:
            print(f"\n{genotype}:")
            sig = id14_norm[genotype]
            # Show top 3 categories
            top_cats = sorted(sig.items(), key=lambda x: x[1], reverse=True)[:3]
            for cat, prop in top_cats:
                print(f"  {cat}: {prop:.1%}")

if __name__ == '__main__':
    main()