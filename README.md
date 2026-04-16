# mESC Single-Cell RNA-seq Analysis (E-MTAB-2600)

This repository contains a comparative study of mouse Embryonic Stem Cells (mESCs) across three culture conditions: **Serum, 2i, and a2i**. 

The primary goal is the identification of the rare **2C-like totipotent state** through two independent bioinformatics pipelines.

##  Project Structure

- [**Python / Scanpy Workflow**](./python_scanpy/): A modern, scalable implementation using the Scanpy ecosystem.
- [**R / Bioconductor Workflow**](./r_bioconductor/): The classic genomic analysis pipeline using `scran` and `scater`.

##  Repository Structure & Organization

The repository is organized into independent workflows. Both pipelines utilize the same raw data but employ language-specific libraries and statistical methodologies.

```
Single-Cell-RNA-seq-Analysis/
├── README.md                   # Root documentation (this file)
├── counts_matrix.csv           # Shared raw gene expression data (Ensembl IDs)
├── E-MTAB-2600.targets.txt     # Shared metadata (Conditions: 2i, a2i, Serum)
│
├── r_bioconductor/             # R Workflow Directory
│   ├── README.md               # Technical details for R implementation
│   ├── single_cell_es.Rmd      # Analysis source code
│   ├── single_cell_es.pdf      # Formatted scientific report & findings
│   ├── single_cell_functions2.R # Modularized R helper functions
│   └── tx2gene.csv             # Transcript-to-gene mapping file
│
└── python_scanpy/              # Python Workflow Directory
    ├── README.md               # Technical details for Python implementation
    ├── single_cell_analysis.py  # Analysis source code
    └── requirements.txt        # Minimalist dependency list (Scanpy/Scikit-misc)
```

## Biological Summary
We analyzed 869 cells to identify a signature subpopulation expressing the *Zscan4* gene family. These cells represent a transient peak of totipotency within pluripotent cultures.

## Data Source
- **Dataset:** ArrayExpress [E-MTAB-2600](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-2600)
- **Input:** Ensembl-indexed counts matrix and metadata targets file.
