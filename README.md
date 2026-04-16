# mESC Single-Cell RNA-seq Analysis (E-MTAB-2600)

This repository contains a comparative study of mouse Embryonic Stem Cells (mESCs) across three culture conditions: **Serum, 2i, and a2i**. 

The primary goal is the identification of the rare **2C-like totipotent state** through two independent bioinformatics pipelines.

##  Project Structure

- [**Python / Scanpy Workflow**](./Python_Scanpy/): A modern, scalable implementation using the Scanpy ecosystem.
- [**R / Bioconductor Workflow**](./R_Bioconductor/): The classic genomic analysis pipeline using `scran` and `scater`.

##  Repository Structure & Organization

The repository is organized into independent workflows. Both pipelines utilize the same raw data but employ language-specific libraries and statistical methodologies.

```
Single-Cell-RNA-seq-Analysis/
├── README.md                   # Root documentation (this file)
├── Python_Scanpy/              # Python Implementation
│   ├── counts_matrix.csv       # Raw gene expression data
│   ├── README.md               # Python technical documentation
│   ├── single_cell_analysis.py # Python analysis source code
│   ├── requirements.txt        # Python dependency list
│   └── figures/                # Generated QC and analysis plots
└── R_Bioconductor/             # R Implementation
    ├── README.md               # R technical documentation
    ├── E-MTAB-2600.targets.txt # Metadata (Conditions: 2i, a2i, Serum)
    ├── single_cell_es.Rmd      # R analysis source code
    ├── single_cell_es.html     # Rendered HTML report
    ├── single_cell_es.pdf      # Final scientific report
    ├── single_cell_functions2.R # Modularized R functions
    ├── hsc_cor.tsv             # Spearman correlation results (R)
    ├── hsc_hvg.tsv             # Highly Variable Genes list (R)
    └── tx2gene.csv             # Transcript-to-gene mapping file
```

## Biological Summary
We analyzed 869 cells to identify a signature subpopulation expressing the *Zscan4* gene family. These cells represent a transient peak of totipotency within pluripotent cultures.

## Data Source
- **Dataset:** ArrayExpress [E-MTAB-2600](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-2600)
- **Input:** Ensembl-indexed counts matrix and metadata targets file.
