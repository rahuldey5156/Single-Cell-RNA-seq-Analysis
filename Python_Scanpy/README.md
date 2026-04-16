# Scanpy Implementation: Single-Cell Analysis of Mouse ESCs (E-MTAB-2600)

This implementation provides a high-fidelity Python-based reproduction of the Bioconductor workflow used to analyze mouse Embryonic Stem Cells (mESCs). The analysis focuses on characterizing transcriptional heterogeneity across **Serum, 2i, and a2i** conditions and identifying the transient **2C-like totipotent state**.

## Technical Configuration & Dependencies

The environment is strictly controlled to ensure reproducibility of the Seurat-v3 variance modeling and to resolve package version conflicts between Scanpy's plotting backends and PyArrow.

- **Framework**: Scanpy (v1.10.0+)
- **Numerical Ops**: NumPy, SciPy, Pandas
- **Visualization**: Matplotlib (Agg backend), Seaborn
- **Gene Annotation**: MyGene API integration
- **Critical Requirements**: 
    - `scikit-misc`: Required for LOESS (Locally Estimated Scatterplot Smoothing).
    - `pyarrow < 22`: Required for Streamlit-compatible data serialization.

---

##  Bioinformatics Pipeline

### 1. Pre-processing & Stringent Quality Control
Data was ingested as an AnnData object from a raw counts matrix of 32,723 genes. QC was performed using a data-driven **3 Median Absolute Deviation (MAD)** approach:
- **Library Size**: Low-count outliers (potential debris) and high-count outliers (potential doublets) were excluded.
- **Transcriptional Complexity**: Cells with low gene counts were removed to ensure sufficient depth for correlation analysis.
- **Mitochondrial Interference**: A 37-gene mitochondrial signature was used to filter apoptotic cells.
- **Technical Controls**: 92 ERCC spike-ins were monitored to evaluate technical noise floors.

### 2. Normalization & Variance Stabilization
To normalize for varying sequencing depth, counts were scaled to a median-normalized target sum of $10^4$ and transformed using $log(1+x)$. 
- **Cell Cycle Correction**: Cells were scored for S and G2M phases using the standard mouse cell cycle marker set.
- **HVG Selection**: 190 Highly Variable Genes (HVGs) were identified using the `seurat_v3` flavor, which uses a variance-stabilizing transformation to identify genes with high biological signal relative to their mean expression.

### 3. Dimensionality Reduction & Manifold Learning
- **Principal Component Analysis (PCA)**: Performed on HVGs to reduce dimensionality while preserving 50 axes of variance.
- **t-SNE Embedding**: A perplexity scan was performed to optimize the visualization of the global data structure. **Perplexity 10** was selected to clearly resolve the separation between culture conditions while preserving local cluster density.

### 4. Identification of the 2C-like Signature
The rare totipotent-like subpopulation was identified through a correlation-based approach:
- **Reference Gene**: *Zscan4a* (a definitive marker of the 2C-state).
- **Correlation Metric**: Spearman’s Rank Correlation ($\rho$) was calculated for all HVGs against *Zscan4a*.
- **Signature Extraction**: Genes with $\rho > 0.4$ were identified as the core 2C-like signature (e.g., *Zscan4c*, *Zscan4f*, *Pramel20*, *Obox4*).

---

## Summary of Results
- **Post-QC Population**: 480 high-quality single cells.
- **Condition Distribution**: 2i (231 cells), Serum (170 cells), a2i (79 cells).
- **Cell Cycle Profile**: Predominantly G2M and S phase, characteristic of rapidly cycling pluripotent cells.
- **Biological Insight**: The signature genes show a clear, localized expression pattern, successfully capturing the rare 2C-like subpopulation as previously described in the original Bioconductor study.

## Usage
```bash
# Install dependencies
pip install -r requirements.txt

# Run analysis (generates outputs in figures/ directory)
python single_cell_analysis.py
