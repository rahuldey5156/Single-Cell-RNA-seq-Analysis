# Single-Cell RNA-seq Analysis of Mouse Embryonic Stem Cells

A fully reproducible R-based analysis of single-cell RNA-sequencing data from mouse embryonic stem cells (mESCs) grown under three culture conditions. The analysis identifies transcriptional differences between conditions at single-cell resolution and characterises a rare 2C-like totipotent cell subpopulation using Zscan4 gene family expression and correlation-based gene discovery.

This project reproduces and extends key findings from:

> Kolodziejczyk AA et al. (2015). *Single Cell RNA-Sequencing of Pluripotent States Unlocks Modular Transcriptional Variation.* Cell Stem Cell, 17(4):471–85. [DOI: 10.1016/j.stem.2015.09.011](https://doi.org/10.1016/j.stem.2015.09.011)

---

## Background

Embryonic stem cells were cultured under three conditions — **serum**, **2i+LIF**, and **a2i** — and profiled by single-cell RNA-seq (scRNA-seq). The dataset (ArrayExpress accession [E-MTAB-2600](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-2600/)) comprises 869 individual cells with paired-end sequencing libraries.

A key biological question is whether a rare **2C-like cell state** (previously identified in bulk RNA-seq) can be detected at single-cell resolution. 2C-like cells express a signature enriched for Zscan4-family genes and are thought to represent a transient totipotent-like state that ES cells cycle through in culture.

---

## Repository Structure

```
.
├── single_cell_es.Rmd          # Main R Markdown analysis notebook
├── single_cell_functions2.R    # Global plotting variables and helper functions
├── E-MTAB-2600_targets.txt     # Sample annotation file (from ArrayExpress SDRF)
├── tx2gene.csv                 # Transcript-to-gene mapping (Ensembl/Biomart)
├── hsc_hvg.tsv                 # Output: highly variable genes table
├── hsc_cor.tsv                 # Output: correlated gene pairs table
├── single_cell_es.html         # Rendered HTML report
├── single_cell_es.pdf          # Rendered PDF report
├── .gitignore
└── README.md
```

> **Note:** The pre-processed Kallisto quantification object (`single_es_tximport.rds`, ~31 MB) is not included due to file size constraints. See the **Data** section below for instructions on how to reproduce it.

---

## Analysis Workflow

### 1. Upstream Processing (performed prior to R analysis)

Raw FASTQ files were quality-checked with **FastQC** and **MultiQC**. Adapter trimming was performed with **Trimmomatic** (settings: `ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36`). Quantification against the mouse transcriptome (Ensembl) concatenated with ERCC spike sequences was performed using **Kallisto** (pseudoalignment, `--bias` flag). Parallel execution of the ~869 alignment jobs was managed using **GNU parallel**.

### 2. Data Loading and QC

- Kallisto output was imported using **tximport**, collapsing transcript-level counts to gene level via `tx2gene.csv`
- ERCC spike-ins (92 sequences) and mitochondrial genes (13) were identified for QC metrics
- Cells were filtered using a **3 MAD (Median Absolute Deviation)** threshold across four metrics:
  - Library size (filtered from below — low counts indicate damaged or empty cells)
  - Number of detected genes (filtered from below)
  - Mitochondrial read proportion (filtered from above — high values indicate cytoplasmic RNA loss)
  - ERCC spike proportion (filtered from above — high values indicate failed library preparation)
- 354 of 869 cells were removed; **515 cells** passed QC

### 3. Gene Filtering and Normalisation

- Genes with mean count < 1 across all cells were removed (32,631 → 20,730 genes retained)
- ERCC spikes were moved to an `altExp` slot within the `SingleCellExperiment` object to avoid interfering with normalisation
- Library size normalisation was performed using **scran** `computeSumFactors()`, with cells pooled by culture condition group size to avoid introducing condition-level bias into clustering

### 4. Gene Annotation

- Ensembl IDs were converted to gene symbols using **org.Mm.eg.db** (Bioconductor mouse annotation)
- Ambiguous mappings were individually reviewed; the Ensembl ID is retained in the `ENSEMBL` slot for traceability

### 5. Cell Cycle Assignment

- Cell cycle phase (G1, S, G2M) was assigned using **Cyclone** (scran package) with built-in mouse cell cycle gene sets
- Cells were assigned to a phase only if they scored above 0.5 for that phase and below 0.5 for both others, ensuring mutually exclusive assignments
- Phase was stored in the `CellCycle` slot of the SCE object
- No strong partitioning of cell cycle phase by culture condition was observed

### 6. Highly Variable Gene (HVG) Detection

- Variance was decomposed into biological and technical components using a loess trend fitted to the mean-variance relationship
- HVG selection criteria: FDR ≤ 0.2 and biological variance component > 0.5
- **190 HVGs** were identified
- Correlated gene pairs among HVGs were identified using `correlatePairs()` (scran); pairs with FDR ≤ 0.01 formed the **chosen feature set** for dimensionality reduction

### 7. Dimensionality Reduction

- **PCA** and **tSNE** were applied using the chosen HVG feature set (scater)
- tSNE perplexity was tuned across values of 5, 10, 15, and 20 using a reusable helper function (`run_tsne_perplexity_scan()` defined in `single_cell_functions2.R`); perplexity = 10 was selected as results were stable across the range
- `set.seed(100)` was applied throughout for reproducibility
- Three well-separated clusters corresponding to culture conditions (serum, 2i, a2i) were clearly resolved

### 8. 2C-like Cell Identification

- All **Zscan4-family genes** were extracted by name pattern matching and visualised in a heatmap
- Spearman correlation of all genes against the **Zscan4a** expression profile identified 45 co-expressed genes (Spearman r > 0.4)
  - Spearman correlation was used over Pearson as it is more robust to the zero-inflated, non-normally distributed counts typical in scRNA-seq
  - A threshold of 0.4 was chosen as a deliberately inclusive cutoff to avoid missing genuine 2C-like co-expressed genes that show noisy single-cell expression
- The correlated gene set includes known 2C markers (Eif1ad family, Pramel genes, Obox pseudogenes) predominantly expressed in a small cluster of **2i cells**, consistent with published findings

---

## Key Results

| Metric | Value |
|---|---|
| Starting cells | 869 |
| Cells after QC | 515 (510 after bulk removal) |
| Cells per condition (post-QC) | 2i: 243 · a2i: 89 · serum: 178 |
| Genes after filtering | 20,730 |
| HVGs identified | 190 |
| Chosen feature set (correlated HVGs) | 190 genes |
| Variance explained by PCA1+2 | ~31% |
| Zscan4a-correlated 2C-like genes | 45 |

Culture conditions separate cleanly by tSNE and PCA. A small subpopulation of 2i cells expresses the Zscan4-family 2C-like signature, consistent with the known ability of ES cells to cycle through a transient totipotent-like state.

---

## Code Refactoring and Version Control

This repository demonstrates a version-controlled development workflow with two stages of commits:

- **Initial commit** — working analysis code rendered from the original tutorial
- **Refactoring commits** — a series of improvements made to the code:
  - Removed duplicate and defunct package calls (`scde` was listed but discontinued)
  - Replaced repetitive tSNE perplexity scan blocks with a reusable `run_tsne_perplexity_scan()` function
  - Replaced an iterative `while` loop for cell cycle assignment with vectorised `dplyr::case_when()`
  - Replaced deprecated `multiplot()` with `gridExtra::grid.arrange()`
  - Added biological reasoning comments explaining key analytical thresholds (MAD=3, correlation > 0.4)
  - Added `sessionInfo()` for full reproducibility documentation

---

## Requirements

### R version

R ≥ 4.0 recommended.

### R Packages

Install Bioconductor packages:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "SingleCellExperiment",
  "scater",
  "scran",
  "scuttle",
  "tximport",
  "org.Mm.eg.db",
  "annotate",
  "limma",
  "edgeR"
))
```

Install CRAN packages:

```r
install.packages(c(
  "ggplot2",
  "plotly",
  "Rtsne",
  "gplots",
  "pheatmap",
  "gridExtra",
  "RColorBrewer",
  "Cairo",
  "pkgconfig",
  "dplyr"
))
```

### Mac (Apple Silicon) users

If you encounter compilation errors for `irlba` or `scater`, install gfortran via Homebrew and configure R to find it:

```bash
brew install gcc
mkdir -p ~/.R
cat > ~/.R/Makevars << 'EOF'
FC = /opt/homebrew/bin/gfortran
F77 = /opt/homebrew/bin/gfortran
FLIBS = -L/opt/homebrew/lib/gcc/14/gcc/aarch64-apple-darwin23/14 -L/opt/homebrew/lib/gcc/14 -lgfortran -lquadmath -lm
EOF
```

Then retry the package installation.

---

## Data

The raw data are publicly available from ArrayExpress under accession **[E-MTAB-2600](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-2600/)**.

To reproduce the upstream Kallisto quantification from raw FASTQ files:

1. Download FASTQ files using links in `E-MTAB-2600.sdrf.txt`
2. Trim adapters with Trimmomatic
3. Build a Kallisto index from the Ensembl mouse transcript FASTA + ERCC spike FASTA (concatenated with `cat`)
4. Run Kallisto quantification (paired-end, `--bias` flag) for each sample using GNU parallel
5. Import with `tximport` using `tx2gene.csv` to collapse to gene level and save as `single_es_tximport.rds`

The `tx2gene.csv` mapping file was generated from Ensembl via Biomart (mouse GRCm38, Ensembl gene ID to transcript ID mapping).

---

## Reproducing the Analysis

```bash
git clone https://github.com/rahuldey5156/Single-Cell-RNA-seq-Analysis.git
cd Single-Cell-RNA-seq-Analysis
```

Place `single_es_tximport.rds` in the working directory (see Data section above), then render from the terminal:

```bash
Rscript -e "rmarkdown::render('single_cell_es.Rmd', output_format='html_document')"
```

All random seeds are fixed via `set.seed(100)` (defined in `single_cell_functions2.R`).

---

## References

1. Kolodziejczyk AA et al. (2015). Single Cell RNA-Sequencing of Pluripotent States Unlocks Modular Transcriptional Variation. *Cell Stem Cell* 17(4):471–85.
2. Kim JK, Kolodziejczyk AA et al. (2015). Characterizing noise structure in single-cell RNA-seq distinguishes genuine from technical stochastic allelic expression. *Nat Commun* 6:8687.
3. Lun AT, McCarthy DJ, Marioni JC (2016). A step-by-step workflow for low-level analysis of single-cell RNA-seq data with Bioconductor. *F1000Res* 5:2122.

---

## License

This analysis code is shared for educational and reproducibility purposes. The underlying data are subject to the terms of [E-MTAB-2600](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-2600/).
