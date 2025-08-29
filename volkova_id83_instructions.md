# Instructions for Volkova et al. (2020) ID83 Analysis

Based on my analysis of your project, here's what you need to do to properly analyze the Volkova et al. (2020) indel data:

## Current Situation

1. You have VCF files from Volkova et al. (2020) in the `IND/` directory
2. You have the genotype names in `table2A_indel14_genotypes_mean.csv`
3. Your professor wants you to classify the raw mutations into ID83 categories
4. You need to convert C. elegans signatures to human genome background

## What Needs to Be Done

### 1. **Parse VCF Files and Extract Indels**
The VCF files contain indel mutations from C. elegans experiments. Each file represents a sample with mutations.

### 2. **Proper ID83 Classification**
The ID83 classification system categorizes indels into 83 classes based on:
- **Indel type**: Insertion or Deletion
- **Sequence context**: 
  - Single base indels in homopolymers (C/G or A/T)
  - Indels at tandem repeats (2bp, 3bp, 4bp, 5+bp units)
  - Deletions with microhomology
- **Length and repeat number**

The 83 categories are:
```
DEL.C.1.1 to DEL.C.1.6+     (6 categories) - C/G homopolymer deletions
DEL.T.1.1 to DEL.T.1.6+     (6 categories) - A/T homopolymer deletions
INS.C.1.0 to INS.C.1.5+     (6 categories) - C/G homopolymer insertions
INS.T.1.0 to INS.T.1.5+     (6 categories) - A/T homopolymer insertions
DEL.repeats.2.1 to 2.6+     (6 categories) - 2bp repeat deletions
DEL.repeats.3.1 to 3.6+     (6 categories) - 3bp repeat deletions
DEL.repeats.4.1 to 4.6+     (6 categories) - 4bp repeat deletions
DEL.repeats.5+.1 to 5+.6+   (6 categories) - 5+bp repeat deletions
INS.repeats.2.0 to 2.5+     (6 categories) - 2bp repeat insertions
INS.repeats.3.0 to 3.5+     (6 categories) - 3bp repeat insertions
INS.repeats.4.0 to 4.5+     (6 categories) - 4bp repeat insertions
INS.repeats.5+.0 to 5+.5+   (6 categories) - 5+bp repeat insertions
DEL.MH.2.1                   (1 category)  - 2bp deletion with microhomology
DEL.MH.3.1 to 3.2           (2 categories) - 3bp deletion with microhomology
DEL.MH.4.1 to 4.3           (3 categories) - 4bp deletion with microhomology
DEL.MH.5+.1 to 5+.5+        (5 categories) - 5+bp deletion with microhomology
```

### 3. **Key Challenges**

To properly classify indels into ID83 categories, you need:

1. **Reference genome context**: To determine if an indel occurs in a homopolymer run or tandem repeat
2. **Microhomology detection**: To identify if deletions have microhomology at breakpoints
3. **Sample-to-genotype mapping**: To know which VCF files correspond to which genotypes

### 4. **Recommended Approach**

Since full classification requires genome context that may not be readily available, I recommend:

1. **Contact the authors or check supplementary data** for:
   - Sample ID to genotype mapping
   - Pre-classified indel data if available
   - Access to their classification scripts

2. **Use a simplified classification** initially:
   - Classify based on indel type and length
   - Note that this will be approximate without genome context

3. **Focus on the workflow**:
   - Extract indels from VCF files
   - Apply classification (even if simplified)
   - Convert to human genome background (using methods from cross_species.ipynb)
   - Compare with COSMIC signatures

### 5. **Code Example**

Here's a Python script structure for the analysis:

```python
import os
import glob
from collections import defaultdict

# Parse VCF files
def parse_vcf(filename):
    indels = []
    with open(filename) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                ref, alt = parts[3], parts[4]
                if len(ref) != len(alt):  # It's an indel
                    indels.append({
                        'chrom': parts[0],
                        'pos': int(parts[1]),
                        'ref': ref,
                        'alt': alt
                    })
    return indels

# Simplified classification (without genome context)
def classify_indel_simple(ref, alt):
    if len(ref) > len(alt):
        # Deletion
        del_len = len(ref) - len(alt)
        if del_len == 1:
            base = ref[len(alt)]
            return f"DEL.{'C' if base in 'CG' else 'T'}.1.1"
        # More complex classification would go here
    else:
        # Insertion
        ins_len = len(alt) - len(ref)
        if ins_len == 1:
            base = alt[len(ref)]
            return f"INS.{'C' if base in 'CG' else 'T'}.1.0"
    return "UNCLASSIFIED"

# Process all VCF files
vcf_files = glob.glob('IND/*.vcf')
results = defaultdict(lambda: defaultdict(int))

for vcf_file in vcf_files:
    sample_id = os.path.basename(vcf_file).split('.')[0]
    indels = parse_vcf(vcf_file)
    
    for indel in indels:
        category = classify_indel_simple(indel['ref'], indel['alt'])
        results[sample_id][category] += 1
```

### 6. **Next Steps After Classification**

1. **Aggregate by genotype**: Once you have the sample-to-genotype mapping
2. **Convert to human background**: Use the methods from your cross_species.ipynb
3. **Compare with COSMIC**: Calculate cosine similarity with COSMIC ID signatures
4. **Visualize results**: Create plots showing signatures before/after conversion

## Important Notes

- The VCF files contain C. elegans chromosome names (I, II, III, IV, V, X, MtDNA)
- Each sample has multiple replicate files (a, b, c, d, etc.)
- The classification requires biological context that simple string matching cannot provide
- Consider reaching out to the paper authors if you need the exact classification methodology

This analysis will demonstrate the cross-species comparison of mutational signatures as described in the paper, showing how DNA repair deficiencies in C. elegans relate to human cancer signatures.