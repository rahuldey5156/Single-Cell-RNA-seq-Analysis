# Seurat Implementation: Single-Cell Analysis of Mouse ESCs

This directory contains the R-based re-implementation of the mESC analysis using the **Seurat v5** framework. This workflow replicates the biological findings of the Bioconductor/scran pipeline, specifically focusing on transcriptional heterogeneity in **Serum, 2i, and a2i** conditions and the identification of the **2C-like totipotent state**.

## Technical Environment
The analysis is performed using Seurat’s object-oriented structure, utilizing the `v5` assay version for optimized data handling and integration.

### Dependencies
The following R packages are required to execute the pipeline:
- **Core Analysis**: `Seurat`, `SeuratObject`
- **Data Manipulation**: `dplyr`, `stringr`, `AnnotationDbi`, `org.Mm.eg.db`
- **Visualization**: `ggplot2`, `pheatmap`, `RColorBrewer`, `gridExtra`

---

## Bioinformatics Pipeline

### 1. Data Ingestion & Metadata Integration
Raw counts from `counts_matrix.csv` (869 cells) were rounded to integers to satisfy Seurat's expectations for count-based models. Metadata from `E-MTAB-2600.targets.txt` was mapped to cell IDs to define experimental groups (Serum, 2i, a2i).

### 2. Stringent Quality Control
Data-driven filtering was applied using a **3 Median Absolute Deviation (MAD)** approach:
- **Library Size & Feature Counts**: Outliers were identified and removed based on log-transformed total counts and detected genes.
- **Mitochondrial Interference**: Percent mitochondrial loading was calculated for every cell; apoptotic cells exceeding the 3-MAD threshold were excluded.
- **ERCC Spikes**: Technical noise was monitored via 92 ERCC spike-in transcripts.

### 3. Normalization & Variance Modeling
- **Normalization**: Counts were normalized using a global-scaling method (`LogNormalize`) with a scale factor of 10,000.
- **HVG Selection**: 2,000 Highly Variable Genes (HVGs) were identified using the **Variance Stabilizing Transformation (vst)** method, which models the mean-variance relationship.

### 4. Dimensionality Reduction & Manifold Learning
- **PCA**: Linear dimensionality reduction was performed on scaled HVG expression data to capture principal components of variance.
- **t-SNE**: Non-linear embedding was optimized with a **perplexity of 10** to resolve the distinct transcriptional architectures of the three culture conditions.

### 5. Totipotency Signature (2C-like State)
The rare totipotent subpopulation was identified through correlation analysis:
- **Reference**: *Zscan4a* expression was used as the anchor for the 2C-state.
- **Correlation**: Spearman’s Rank Correlation ($\rho$) was calculated for all genes against *Zscan4a*.
- **Identification**: Genes with $\rho > 0.4$ were extracted as the core **2C-like signature**.

---

## Directory Contents & Outputs

### Documentation & Source
- **`seurat_analysis.R`**: The complete analysis script.
- **`tx2gene.csv`**: Mapping file for Ensembl-to-Symbol conversion.

### Data Exports
- **`hsc_hvg.tsv`**: List of top Highly Variable Genes identified by Seurat.
- **`hsc_cor.tsv`**: Spearman correlation coefficients for the Zscan4a signature.

### Visualizations
- **QC & Filtering**: `qc_before_filtering.png`, `qc_violin_postfilter.png`, `gene_filtering.png`.
- **Dimensionality Reduction**: `pca_condition.png`, `tsne_condition.png`, `tsne_perplexity_condition.png`.
- **Biological Analysis**: 
  - `heatmap_zscan4a_correlated.png`: The definitive 2C-like signature heatmap.
  - `cell_cycle_scores.png` & `tsne_cell_cycle.png`: Assessment of cell cycle phase distribution.
  - `heatmap_hvg100.png`: Top 100 HVG expression across conditions.

---

## Execution
To run the analysis and regenerate all plots and tables:
```bash
Rscript seurat_analysis.R
