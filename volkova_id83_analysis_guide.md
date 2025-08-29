# Volkova et al. (2020) ID83 Analysis Guide

## Overview
Your professor has requested that you process the raw VCF files from Volkova et al. (2020) and classify indels into ID83 categories instead of using the pre-processed ID14 data.

## Key Changes Required

### Previous Approach (ID14)
- Used pre-processed data from `table2A_indel14_genotypes_mean.csv`
- 14 indel categories
- Already aggregated counts

### New Approach (ID83)
- Process raw VCF files directly
- 83 indel categories (more detailed)
- Better resolution for understanding mutational processes

## Data Sources

### Option 1: Figshare (Recommended)
- URL: https://figshare.com/articles/dataset/Strain_and_sample_description_/16754007
- Contains filtered VCF files for each C. elegans strain
- Easier to use as files are already processed

### Option 2: ENA
- Accession numbers: ERP000975 and ERP004086
- Raw sequencing data (would need to call variants first)

## Implementation Steps

1. **Download VCF Files**
   ```bash
   # Create directory for VCF files
   mkdir volkova_vcf_files
   cd volkova_vcf_files
   
   # Download from figshare (manual download required)
   # Navigate to the figshare link and download all VCF files
   ```

2. **Install Required Python Packages**
   ```bash
   pip install pyvcf3 pysam pandas numpy matplotlib seaborn
   ```

3. **Run the Analysis**
   - The new sections in `cross_species.ipynb` provide all necessary functions
   - Key functions:
     - `process_vcf_to_id83()`: Process individual VCF files
     - `create_id83_matrix()`: Create count matrix for all samples
     - `normalize_to_human_background()`: Essential for cross-species comparison
     - `id83_to_id14()`: Convert to ID14 for comparison with previous analysis

4. **Update File Paths**
   ```python
   # In the notebook, update this line:
   results = analyze_volkova_vcfs('/path/to/your/vcf/directory/')
   ```

## ID83 Classification System

The ID83 system classifies indels based on:

### 1bp Indels
- Base type (C/G vs A/T)
- Homopolymer run length (0-6+)
- Total: 24 categories (12 DEL + 12 INS)

### Multi-bp Indels
- Length (2bp, 3bp, 4bp, 5+bp)
- Context:
  - Tandem repeats (with repeat count)
  - Microhomology (deletions only)
- Total: 59 categories

## Important Considerations

1. **Reference Genome**: The classification requires sequence context. The simplified version in the notebook uses placeholder values. For accurate classification, you'll need:
   - C. elegans reference genome (ce11 or WBcel235)
   - Proper context extraction around each indel

2. **Normalization**: Critical for cross-species comparison
   - Accounts for different sequence compositions between C. elegans and human
   - Current implementation uses simplified factors
   - Consider using actual genome-wide context frequencies

3. **Validation**: Compare ID14 results from both approaches
   - Convert ID83 → ID14 using provided function
   - Results should be similar but not identical due to:
     - Different classification methods
     - Potential differences in variant filtering

## Expected Outputs

1. `volkova_id83_raw_counts.csv`: Raw counts in 83 categories
2. `volkova_id83_normalized.csv`: Normalized to human background
3. `volkova_id14_normalized.csv`: Converted to ID14 for comparison

## Troubleshooting

- **VCF parsing errors**: Check VCF format, ensure files are not corrupted
- **Missing categories**: Some ID83 categories may have zero counts
- **Normalization issues**: Verify genome context frequencies are reasonable

## Next Steps

After running the analysis:
1. Compare ID83 results with COSMIC signatures
2. Visualize signatures (similar to Fig4a in your goal)
3. Compare with your previous ID14-based results
4. Document any interesting differences or insights