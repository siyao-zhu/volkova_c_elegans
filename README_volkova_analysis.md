# Volkova et al. (2020) Indel Analysis - Integration Guide

This guide explains how to process the raw VCF files from Volkova et al. (2020) and integrate them with your existing cross-species indel signature analysis.

## Overview

The Volkova study used C. elegans as a model organism to experimentally determine mutational signatures of DNA repair and damage. Your professor's feedback emphasizes that **renormalizing to the human background is essential** for interpreting these signatures in the context of human cancer genomics.

## What These Scripts Do

### 1. `process_volkova_vcfs_improved.py`
- **Input**: Raw VCF files from the `IND/` directory
- **Output**: Indel counts classified into ID83 COSMIC categories
- **Key Features**:
  - Parses VCF files and extracts indel mutations
  - Classifies each indel into the appropriate ID83 category
  - Handles edge cases and provides detailed error reporting
  - Creates comprehensive summary statistics

### 2. `volkova_analysis_notebook.py`
- **Input**: Processed indel counts from step 1
- **Output**: Integrated analysis ready for your existing pipeline
- **Key Features**:
  - Converts ID83 to ID14 categories (matching your existing analysis)
  - Applies background normalization (C. elegans → Human)
  - Creates visualizations comparing before/after signatures
  - Saves results in formats compatible with your existing code

## Quick Start

### Step 1: Process VCF Files
```bash
python process_volkova_vcfs_improved.py
```

This will:
- Process all VCF files in the `IND/` directory
- Classify indels into ID83 categories
- Save results to `volkova_indel_counts_id83_improved.tsv`
- Generate a summary report

### Step 2: Run Analysis Pipeline
```bash
python volkova_analysis_notebook.py
```

This will:
- Load the processed data
- Convert to ID14 categories
- Apply background normalization
- Create comparison plots
- Save results for integration

## Output Files

After running both scripts, you'll have:

1. **`volkova_indel_counts_id83_improved.tsv`** - Raw indel counts in ID83 categories
2. **`volkova_indel_signatures_id14.tsv`** - Signatures in ID14 categories (C. elegans background)
3. **`volkova_indel_signatures_humanized.tsv`** - Signatures in ID14 categories (Human background)
4. **`volkova_mean_signatures_comparison.tsv`** - Mean signatures comparison
5. **`volkova_indel_analysis_summary.txt`** - Detailed analysis summary

## Integration with Your Existing Analysis

### Loading the Data
```python
# Load the processed Volkova data
volkova_ce = pd.read_csv("volkova_indel_signatures_id14.tsv", sep="\t", index_col=0)
volkova_human = pd.read_csv("volkova_indel_signatures_humanized.tsv", sep="\t", index_col=0)

# These are now compatible with your existing functions
```

### Using with Your Existing Functions
The processed data follows the same format as your existing analysis:

- **Rows**: ID14 categories (D.1.rep, D.1.nonrep, etc.)
- **Columns**: Samples (CD0268c, CD0268d, etc.)
- **Values**: Normalized indel frequencies

### Background Normalization
The key insight from your professor is implemented as:

```python
# S_humanized[c] = S_worm[c] / BG_worm[c] * BG_human[c]
conv = (M / bw) * bh
```

Where:
- `M` = Original C. elegans signatures
- `bw` = C. elegans background frequencies
- `bh` = Human background frequencies
- `conv` = Humanized signatures

## Understanding the ID83 Classification

The ID83 system classifies indels based on:

1. **Type**: Deletion (DEL) or Insertion (INS)
2. **Context**: C/G vs T/A bases
3. **Length**: 1bp, 2-5bp, 5-50bp, 50-400bp
4. **Mechanism**: Repeat-mediated vs Microhomology-mediated
5. **Repeat count**: Number of repeat units affected

### Example Classifications
- `DEL.C.1.1` = 1bp deletion in C/G context
- `INS.repeats.2.0` = 2bp insertion in repeat context
- `DEL.MH.3.1` = 3bp deletion with microhomology

## Troubleshooting

### Common Issues

1. **VCF files not found**
   - Ensure VCF files are in the `IND/` directory
   - Check file permissions

2. **Classification errors**
   - Review the VCF file format
   - Check INFO field parsing

3. **Background data missing**
   - The scripts use placeholder values for demonstration
   - Replace with your actual background data

### Debug Mode
Add debug prints to see detailed classification:
```python
# In process_volkova_vcfs_improved.py
print(f"Classifying: {ref}->{alt}, REP={repeat_count}")
print(f"Result: {indel_class}")
```

## Next Steps for Your Analysis

1. **Run the scripts** to process your VCF files
2. **Load the results** into your existing notebook
3. **Compare signatures** before/after human background conversion
4. **Match with COSMIC** using your existing cosine similarity functions
5. **Create Figure 4a** showing the transformation
6. **Analyze the impact** of background normalization on signature matching

## Key Insights from Your Professor

> "BTW I was reading your updates from the past weeks. For data in this paper, you need to obtain the raw mutations from their experiments and annotate them into ID83 classes. The vcf files should be available according to https://www.nature.com/articles/s41467-020-15912-7#data-availability."

This means:
1. **Use raw VCF files** (not pre-processed data)
2. **Classify into ID83** (COSMIC standard)
3. **Apply background normalization** (essential for interpretation)

## Contact and Support

If you encounter issues:
1. Check the error messages in the terminal output
2. Verify VCF file format and content
3. Ensure all required Python packages are installed
4. Review the classification logic for edge cases

The scripts are designed to be robust and provide detailed feedback about any issues encountered during processing.