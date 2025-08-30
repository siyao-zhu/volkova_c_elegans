# Volkova et al. 2020 Indel Analysis Summary

## Overview

This analysis processes VCF files from Volkova et al. (2020) to classify indel mutations into ID83 categories and compare C. elegans DNA repair mutant signatures with human cancer signatures from COSMIC.

## Analysis Steps Completed

### 1. VCF File Processing
- Parsed 273 VCF files from the CD05XX series
- Extracted indel mutations with their characteristics (type, length, repeat context)
- Files processed: `/workspace/IND/CD05*.vcf`

### 2. ID83 Classification
- Classified each indel into one of 83 COSMIC categories based on:
  - Mutation type (insertion/deletion)
  - Length (1bp, 2-5bp, 5+bp)
  - Sequence context (homopolymer runs, repeats, microhomology)
- Created mutational catalogues for each sample

### 3. Sample-to-Genotype Mapping
Mapped samples to their genotypes based on the experimental design:
- **N2** (wild-type): CD0501-CD0503
- **MMR deficient**: mlh-1, pms-2, msh-2, msh-6
- **NER deficient**: xpc-1, xpa-1
- **TLS polymerase**: rev-3, polk-1, polh-1

### 4. ID14 Conversion
Converted ID83 signatures to the 14-category format used in the paper:
- Deletion categories: D.1.rep, D.1.nonrep, D.2.5.rep, D.2.5.nonrep, D.5.50
- Insertion categories: I.1.rep, I.1.nonrep, I.2.5.rep, I.2.5.nonrep, I.5.50

## Key Results

### Mutational Signatures Summary

| Genotype | Total Indels | Dominant Category | Proportion |
|----------|--------------|-------------------|------------|
| N2       | 10          | D.1.rep           | 40.0%      |
| mlh-1    | 21          | D.5.50            | 33.3%      |
| msh-2    | 20          | I.1.rep           | 25.0%      |
| msh-6    | 28          | D.5.50            | 35.7%      |
| pms-2    | 16          | D.5.50            | 50.0%      |
| polh-1   | 5           | D.5.50            | 60.0%      |
| polk-1   | 13          | I.1.rep           | 38.5%      |
| rev-3    | 16          | D.5.50            | 43.8%      |
| xpa-1    | 55          | D.5.50            | 76.4%      |
| xpc-1    | 19          | D.5.50            | 36.8%      |

### Key Findings

1. **MMR-deficient mutants** (mlh-1, msh-2, msh-6, pms-2) show increased deletions, particularly in the D.5.50 category (5+ bp deletions)

2. **NER-deficient mutants** (xpa-1, xpc-1) show strong bias toward D.5.50 deletions, with xpa-1 showing the strongest signal (76.4%)

3. **TLS polymerase mutants** (rev-3, polk-1) show mixed patterns, with polk-1 favoring insertions (I.1.rep)

4. **Wild-type (N2)** shows relatively balanced profile with preference for single-base deletions in repeat contexts

## Output Files Generated

1. **volkova_id83_signatures_cd05.csv** - Raw ID83 counts by genotype
2. **volkova_id83_signatures_cd05_transposed.csv** - Transposed format for analysis
3. **volkova_id14_signatures.csv** - Converted to 14-category format
4. **volkova_id14_normalized.csv** - Normalized signatures (proportions)

## Next Steps for Full Analysis

To complete the analysis as described in the paper, you would need to:

1. **Obtain genome background frequencies** for both C. elegans and human genomes
2. **Apply background correction** to adjust for different sequence contexts
3. **Compare with COSMIC signatures** using cosine similarity
4. **Create visualizations** showing the signatures before and after humanization

## Technical Notes

- The VCF files contain repeat count (REP) information which helps classify repeat-mediated indels
- Microhomology classification would benefit from flanking sequence analysis
- The current analysis uses simplified classification rules that may need refinement

## References

Volkova, N.V., et al. (2020). "Mutational signatures are jointly shaped by DNA damage and repair." Nature Communications 11, 2169.