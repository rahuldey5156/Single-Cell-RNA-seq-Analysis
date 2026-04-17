# mESC Single-Cell RNA-seq Analysis (E-MTAB-2600)

This repository contains a comparative study of mouse Embryonic Stem Cells (mESCs) across three culture conditions: **Serum, 2i, and a2i**. 

The primary goal is the identification of the rare **2C-like totipotent state** through two independent bioinformatics pipelines.

##  Project Structure

- [**Python / Scanpy Workflow**](./Python_Scanpy/): A modern, scalable implementation using the Scanpy ecosystem.
- [**R / Bioconductor Workflow**](./R_Bioconductor/): The classic genomic analysis pipeline using `scran` and `scater`.
- [**R / Seurat Workflow**](./R_Seurat/): A high-performance, object-oriented implementation using the Seurat v5 framework.

##  Repository Structure & Organization

The repository is organized into independent workflows. All the 3 pipelines utilize the same raw data but employ language-specific libraries and statistical methodologies.

```
Single-Cell-RNA-seq-Analysis/
├── README.md                   # Root documentation (this file)
│
├── Python_Scanpy/              # Python Implementation
│   ├── counts_matrix.csv       # Raw gene expression data
│   ├── README.md               # Python technical documentation
│   ├── single_cell_analysis.py # Python analysis source code
│   ├── requirements.txt        # Python dependency list
│   └── figures/                # Generated QC and analysis plots
│
├── R_Bioconductor/             # R Implementation (scran/scater)
│   ├── README.md               # Bioconductor technical documentation
│   ├── E-MTAB-2600.targets.txt # Metadata (Conditions: 2i, a2i, Serum)
│   ├── single_cell_es.Rmd      # R analysis source code
│   ├── single_cell_es.html     # Rendered HTML report
│   ├── single_cell_es.pdf      # Final scientific report
│   ├── single_cell_functions2.R # Modularized R functions
│   ├── hsc_cor.tsv             # Spearman correlation results
│   ├── hsc_hvg.tsv             # Highly Variable Genes list
│   └── tx2gene.csv             # Transcript-to-gene mapping file
│
└── R_Seurat/                   # R Implementation (Seurat v5)
    ├── README.md               # Seurat technical documentation
    ├── counts_matrix.csv       # Local copy of expression data
    ├── E-MTAB-2600.targets.txt # Local copy of metadata
    ├── tx2gene.csv             # Local copy of mapping file
    ├── seurat_analysis.R       # Seurat analysis source code
    ├── hsc_cor.tsv             # Spearman correlation results (Seurat)
    ├── hsc_hvg.tsv             # Highly Variable Genes list (Seurat)
    └── [Visualizations]        # Comprehensive PNG outputs (QC, PCA, t-SNE, Heatmaps)
```

## Biological Summary
We analyzed 869 cells to identify a signature subpopulation expressing the *Zscan4* gene family. These cells represent a transient peak of totipotency within pluripotent cultures.

## Unified Pipeline Objectives
- Quality Control: Data-driven filtering (3 MADs) for library size, feature complexity, and mitochondrial content.

- Normalization: Correcting for sequencing depth to ensure comparable gene expression profiles.

- Feature Selection: Identifying top Highly Variable Genes (HVGs) via language-specific variance modeling (e.g., VST in Seurat vs. Seurat_v3 in Scanpy).

- Dimensionality Reduction: Utilizing PCA and t-SNE to visualize the separation of culture conditions and the emergence of the 2C-like cluster.

## Implementation Mapping
| Feature | Bioconductor | Scanpy | Seurat |
|--------|-------------|--------|--------|
| HVG Selection | `modelGeneVar` | `highly_variable_genes` (Seurat_v3) | `FindVariableFeatures` (VST) |
| Annotation | biomaRt | MyGene API | org.Mm.eg.db |
| Data Structure | SingleCellExperiment | AnnData | Seurat Object |

## Methodological Comparison

### Highly Variable Genes (HVGs)
- Bioconductor (`modelGeneVar`) uses a variance decomposition approach.
- Seurat (VST) stabilizes variance using a transformation.
- Scanpy (Seurat_v3 flavor) mimics Seurat but may differ slightly due to implementation.

Observation: The number of HVGs selected differed slightly across methods, which affected downstream clustering.

---

### Clustering
- Seurat uses Louvain/Leiden clustering via `FindClusters`.
- Scanpy typically uses Leiden clustering.
- Parameter choices (resolution) significantly influenced cluster granularity.

Observation: Scanpy produced slightly more clusters at the same resolution compared to Seurat.

---

### Normalization
- Seurat: `LogNormalize`
- Scanpy: `normalize_total` + `log1p`
- Bioconductor: `scran` normalization

Observation: Differences in normalization led to minor variation in gene expression scaling but did not drastically change biological interpretation.

## Data Source
- **Dataset:** ArrayExpress [E-MTAB-2600](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-2600)
- **Input:** Ensembl-indexed counts matrix and metadata targets file.
