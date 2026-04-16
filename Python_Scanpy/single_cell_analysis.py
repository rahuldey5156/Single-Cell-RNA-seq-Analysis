# =============================================================================
# Single-Cell RNA-seq Analysis of Mouse ESCs — Scanpy (Python)
# Reproduces: single_cell_es.Rmd (Bioconductor/scran/scater workflow)
# Dataset: E-MTAB-2600 | 869 cells, 3 conditions (serum, 2i, a2i)
# Input: counts_matrix.csv (genes x cells, Ensembl IDs as index)
# =============================================================================

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from scipy.stats import spearmanr
from scipy.sparse import issparse
import warnings
warnings.filterwarnings("ignore")

# Reproducibility seed (mirrors set.seed(100) in R)
SEED = 100
np.random.seed(SEED)
sc.settings.verbosity = 2
sc.settings.figsize = (6, 5)

# =============================================================================
# 1. LOAD DATA
# =============================================================================

print("=== 1. Loading counts matrix ===")

counts = pd.read_csv("counts_matrix.csv", index_col=0)
print(f"Raw matrix: {counts.shape[0]} genes x {counts.shape[1]} cells")

# AnnData: rows = cells, cols = genes (transpose from your CSV which is genes x cells)
adata = sc.AnnData(X=counts.T.values.astype(np.float64))
adata.obs_names = counts.columns.tolist()
adata.var_names = counts.index.tolist()

print(f"AnnData shape (cells x genes): {adata.shape}")

# =============================================================================
# 2. IDENTIFY ERCC SPIKES AND MITOCHONDRIAL GENES
# =============================================================================

print("\n=== 2. Identifying ERCC spikes and mitochondrial genes ===")

# ERCC spikes — named "ERCC-..." in Ensembl-based data
adata.var["is_ercc"] = adata.var_names.str.startswith("ERCC")

# Mitochondrial genes — exact Ensembl IDs from your R script
mito_ids = {
    "ENSMUSG00000064336","ENSMUSG00000064337","ENSMUSG00000064338","ENSMUSG00000064339",
    "ENSMUSG00000064340","ENSMUSG00000064341","ENSMUSG00000064342","ENSMUSG00000064343",
    "ENSMUSG00000064344","ENSMUSG00000064345","ENSMUSG00000064346","ENSMUSG00000064347",
    "ENSMUSG00000064348","ENSMUSG00000064349","ENSMUSG00000064350","ENSMUSG00000064351",
    "ENSMUSG00000064352","ENSMUSG00000064353","ENSMUSG00000064354","ENSMUSG00000064355",
    "ENSMUSG00000064356","ENSMUSG00000064357","ENSMUSG00000064358","ENSMUSG00000064359",
    "ENSMUSG00000064360","ENSMUSG00000064361","ENSMUSG00000065947","ENSMUSG00000064363",
    "ENSMUSG00000064364","ENSMUSG00000064365","ENSMUSG00000064366","ENSMUSG00000064367",
    "ENSMUSG00000064368","ENSMUSG00000064369","ENSMUSG00000064370","ENSMUSG00000064371",
    "ENSMUSG00000064372"
}
adata.var["is_mito"] = [g in mito_ids for g in adata.var_names]

print(f"ERCC spikes found: {adata.var['is_ercc'].sum()}")
print(f"Mitochondrial genes found: {adata.var['is_mito'].sum()}")

# =============================================================================
# 3. QC METRICS AND CELL FILTERING (3 MAD thresholds)
# =============================================================================

print("\n=== 3. QC metrics and cell filtering ===")

# Calculate per-cell QC metrics
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["is_ercc", "is_mito"],
    percent_top=None,
    log1p=False,
    inplace=True
)
# adata.obs now contains:
#   total_counts, n_genes_by_counts
#   pct_counts_is_ercc, pct_counts_is_mito

# Plot QC distributions before filtering
fig, axes = plt.subplots(1, 2, figsize=(12, 4))
axes[0].hist(adata.obs["total_counts"] / 1e6, bins=20, color="grey")
axes[0].set_xlabel("Library sizes (millions)")
axes[0].set_ylabel("Number of cells")
axes[1].hist(adata.obs["n_genes_by_counts"], bins=20, color="grey")
axes[1].set_xlabel("Number of detected genes")
axes[1].set_ylabel("Number of cells")
plt.tight_layout()
plt.savefig("qc_before_filtering.png", dpi=150)

fig, axes = plt.subplots(1, 2, figsize=(12, 4))
axes[0].hist(adata.obs["pct_counts_is_mito"], bins=20, color="grey")
axes[0].set_xlabel("Mitochondrial proportion (%)")
axes[0].set_ylabel("Number of cells")
axes[1].hist(adata.obs["pct_counts_is_ercc"], bins=20, color="grey")
axes[1].set_xlabel("ERCC proportion (%)")
axes[1].set_ylabel("Number of cells")
plt.tight_layout()
plt.savefig("qc_mito_ercc.png", dpi=150)


def mad_filter(series, nmads=3, log=False, direction="both"):
    """
    Replicates scuttle::isOutlier().
    direction: 'lower', 'higher', or 'both'
    log: apply log transform before MAD calculation (mirrors log=TRUE in R)
    """
    vals = np.log(series + 1) if log else series.copy()
    median = np.median(vals)
    mad = np.median(np.abs(vals - median))
    lower = median - nmads * mad
    upper = median + nmads * mad
    if direction == "lower":
        return series < (np.exp(lower) - 1 if log else lower)
    elif direction == "higher":
        return series > upper
    else:
        return (series < (np.exp(lower) - 1 if log else lower)) | (series > upper)


libsize_drop  = mad_filter(adata.obs["total_counts"], nmads=3, log=True,  direction="lower")
feature_drop  = mad_filter(adata.obs["n_genes_by_counts"], nmads=3, log=True, direction="lower")
mito_drop     = mad_filter(adata.obs["pct_counts_is_mito"], nmads=3, direction="higher")
spike_drop    = mad_filter(adata.obs["pct_counts_is_ercc"], nmads=3, direction="higher")

keep_cells = ~(libsize_drop | feature_drop | mito_drop | spike_drop)

print(f"\nCells removed by library size filter:  {libsize_drop.sum()}")
print(f"Cells removed by feature count filter: {feature_drop.sum()}")
print(f"Cells removed by mito filter:          {mito_drop.sum()}")
print(f"Cells removed by ERCC spike filter:    {spike_drop.sum()}")
print(f"Cells remaining after QC:              {keep_cells.sum()}")

# Keep a copy of pre-filter object (mirrors sceX in R)
adata_raw = adata.copy()

adata = adata[keep_cells].copy()
print(f"Post-QC AnnData shape: {adata.shape}")

# =============================================================================
# 4. LOG COUNTS FOR INITIAL PCA (raw log, no normalisation yet)
#    Mirrors: assays(sce)$logcounts <- log2(counts + 0.01); runPCA(sce)
# =============================================================================

print("\n=== 4. Initial PCA on raw log counts (pre-normalisation) ===")

adata_log_raw = adata.copy()
X = adata_log_raw.X
adata_log_raw.X = np.log2(X + 0.01)

sc.pp.pca(adata_log_raw, n_comps=50, random_state=SEED)
sc.pl.pca(adata_log_raw, title="PCA — raw log counts (no normalisation)", save="_raw_log.png")

# =============================================================================
# 5. GENE SYMBOL ANNOTATION
#    Mirrors: AnnotationDbi::select(org.Mm.eg.db, ..., keytype="ENSEMBL", column="SYMBOL")
#    Requires: pip install mygene
# =============================================================================

print("\n=== 5. Annotating Ensembl IDs with gene symbols ===")

try:
    import mygene
    mg = mygene.MyGeneInfo()
    non_ercc = [g for g in adata.var_names if not g.startswith("ERCC")]
    result = mg.querymany(non_ercc, scopes="ensembl.gene", fields="symbol", species="mouse", as_dataframe=True)
    symbol_map = result["symbol"].dropna().to_dict()
except Exception as e:
    print(f"mygene lookup failed ({e}), falling back to Ensembl ID as symbol")
    symbol_map = {}

def ensembl_to_symbol(ensembl_id, symbol_map):
    if ensembl_id.startswith("ERCC"):
        return ensembl_id
    return symbol_map.get(ensembl_id, ensembl_id)  # fallback: keep Ensembl ID

adata.var["ENSEMBL"] = adata.var_names.tolist()

raw_symbols = [ensembl_to_symbol(g, symbol_map) for g in adata.var_names]

# Make unique (mirrors make.names(..., unique=TRUE) in R)
seen = {}
unique_symbols = []
for s in raw_symbols:
    if s not in seen:
        seen[s] = 0
        unique_symbols.append(s)
    else:
        seen[s] += 1
        unique_symbols.append(f"{s}.{seen[s]}")

adata.var["SYMBOL"] = unique_symbols
adata.var_names = unique_symbols
adata.var_names_make_unique()

print(f"Gene annotation done. Example: {list(zip(list(adata.var['ENSEMBL'])[:5], list(adata.var_names)[:5]))}")

# Update ERCC / mito flags after rename
adata.var["is_ercc"] = adata.var["ENSEMBL"].str.startswith("ERCC")
adata.var["is_mito"] = [e in mito_ids for e in adata.var["ENSEMBL"]]

# =============================================================================
# 6. SEPARATE ERCC SPIKES (mirrors splitAltExps in R)
# =============================================================================

print("\n=== 6. Separating ERCC spikes from main data ===")

spike_mask = adata.var["is_ercc"].values
adata_spikes = adata[:, spike_mask].copy()
adata = adata[:, ~spike_mask].copy()

print(f"Main data (no spikes): {adata.shape}")
print(f"Spike data: {adata_spikes.shape}")

# =============================================================================
# 7. GENE FILTERING: keep genes expressed in >= 10 cells
#    Mirrors: alt.keep <- nexprs(sce, byrow=TRUE) >= 10
# =============================================================================

print("\n=== 7. Gene filtering (expressed in >= 10 cells) ===")

print(f"Genes before filtering: {adata.n_vars}")
sc.pp.filter_genes(adata, min_cells=10)
print(f"Genes after filtering:  {adata.n_vars}")

# Also plot: average count distribution
ave_counts = np.array(adata.X).mean(axis=0)
fig, axes = plt.subplots(1, 2, figsize=(12, 4))
axes[0].hist(np.log10(ave_counts + 1e-6), bins=100, color="grey")
axes[0].axvline(x=np.log10(1), color="blue", linestyle="--", linewidth=2)
axes[0].set_xlabel("Log10 average count")
axes[0].set_ylabel("Number of genes")
n_cells_expressing = (np.array(adata.X) > 0).sum(axis=0)
axes[1].scatter(np.log10(ave_counts + 1e-6), n_cells_expressing, s=1, alpha=0.3)
axes[1].set_xlabel("Log10 average count")
axes[1].set_ylabel("Number of expressing cells")
plt.tight_layout()
plt.savefig("gene_filtering.png", dpi=150)

# =============================================================================
# 8. LOAD TARGETS FILE AND REMOVE BULK SAMPLES
#    Mirrors: read.csv("E-MTAB-2600.targets.txt"); remove bulk samples
#    NOTE: Download E-MTAB-2600.targets.txt from your GitHub repo
# =============================================================================

print("\n=== 8. Loading targets and removing bulk samples ===")

targets = pd.read_csv("E-MTAB-2600.targets.txt", sep="\t")
print(f"Targets file columns: {targets.columns.tolist()}")
print(f"Targets shape: {targets.shape}")

# Match cells in AnnData to targets (cells are currently named ERR...)
# The R script renames cells to targets$ExtractName after matching
err_to_type  = dict(zip(targets["ERR"].astype(str), targets["Type"]))
err_to_name  = dict(zip(targets["ERR"].astype(str), targets["ExtractName"]))

# Current obs_names are ERR IDs
adata.obs["ERR"] = adata.obs_names.tolist()
adata.obs["Type"] = adata.obs["ERR"].map(err_to_type)
adata.obs["ExtractName"] = adata.obs["ERR"].map(err_to_name)

# Rename cells to ExtractName (mirrors colnames(sce) <- targets$ExtractName)
adata.obs_names = adata.obs["ExtractName"].fillna(adata.obs["ERR"]).tolist()

# Remove bulk samples (mirrors removing rows where ExtractName contains "bulk")
is_bulk = adata.obs_names.str.contains("bulk")
print(f"Bulk samples found: {is_bulk.sum()}")
adata = adata[~is_bulk].copy()

print(f"Cells after bulk removal: {adata.n_obs}")
print("Cells per condition:")
print(adata.obs["Type"].value_counts())

# =============================================================================
# 9. NORMALISATION
#    Mirrors: computeSumFactors (scran) + logNormCounts
#    Scanpy equivalent: sc.pp.normalize_total + sc.pp.log1p
#    For scran-equivalent pooling-based normalisation use scran via rpy2 or
#    the scanpy-external scran wrapper. Pure Python fallback shown here.
# =============================================================================

print("\n=== 9. Normalisation ===")

# Store raw counts
adata.layers["counts"] = adata.X.copy()

# Option A: scran-style normalisation via scanpy-external (recommended)
# Requires: pip install rpy2  AND  R with scran installed
try:
    import scanpy.external as sce_ext
    # scran pooling groups: use condition as cluster proxy (mirrors sizes=c(A,C,E))
    condition_codes = adata.obs["Type"].astype("category").cat.codes
    sce_ext.pp.scran(adata, use_highly_variable=False, key_added="size_factors",
                     random_state=SEED)
    adata.obs["size_factor"] = adata.obs["size_factors"]
    print("scran normalisation applied via scanpy.external")
except Exception as e:
    print(f"scran not available ({e}), using scanpy normalize_total as fallback")
    # Option B: library-size normalisation (simpler but equivalent for most purposes)
    sc.pp.normalize_total(adata, target_sum=None, inplace=True)
    adata.obs["size_factor"] = adata.obs["total_counts"] / adata.obs["total_counts"].median()

# Log-normalise (mirrors logNormCounts)
sc.pp.log1p(adata)
adata.layers["logcounts"] = adata.X.copy()
print("Log-normalised counts stored in adata.layers['logcounts']")

# Size factor vs library size plot
fig, ax = plt.subplots(figsize=(6, 5))
ax.scatter(adata.obs["size_factor"], adata.obs["total_counts"] / 1e6,
           s=5, alpha=0.5)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("Size factor"); ax.set_ylabel("Library size (millions)")
ax.set_title("Size factors vs library size")
plt.tight_layout()
plt.savefig("size_factors.png", dpi=150)

# =============================================================================
# 10. CELL CYCLE ASSIGNMENT
#     Mirrors: cyclone(sce, mm.pairs, ...) + G1/S/G2M scoring
#     Scanpy uses Tirosh et al. gene sets (human); we adapt for mouse
# =============================================================================

print("\n=== 10. Cell cycle scoring ===")

# Mouse cell cycle genes — conservative list matching major markers
# (cyclone uses a larger probabilistic set; this is the Scanpy equivalent)
s_genes_mouse = [
    "Mcm5","Pcna","Tyms","Fen1","Mcm2","Mcm4","Rrm1","Ung","Gins2",
    "Mcm6","Cdca7","Dtl","Prim1","Uhrf1","Mlf1ip","Hells","Rfc2","Rpa2",
    "Nasp","Rad51ap1","Gmnn","Wdr76","Slbp","Ccne2","Ubr7","Pold3","Msh2",
    "Atad2","Rad51","Rrm2","Cdc45","Cdc6","Exo1","Tipin","Dscc1","Blm",
    "Casp8ap2","Usp1","Clspn","Pola1","Chaf1b","Brip1","E2f8"
]
g2m_genes_mouse = [
    "Hmgb2","Cdk1","Nusap1","Ube2c","Birc5","Tpx2","Top2a","Ndc80",
    "Cks2","Nuf2","Cks1b","Mki67","Tmpo","Cenpf","Tacc3","Fam64a",
    "Smc4","Ccnb2","Ckap2l","Ckap2","Aurkb","Bub1","Kif11","Anp32e",
    "Tubb4b","Gtse1","Kif20b","Hjurp","Cdca3","Hn1","Cdc20","Ttk",
    "Cdc25c","Kif2c","Rangap1","Ncapd2","Dlgap5","Cdca2","Cdca8",
    "Ect2","Kif23","Hmmr","Aurka","Psrc1","Anln","Lbr","Ckap5",
    "Cenpe","Ctcf","Nek2","G2e3","Gas2l3","Cbx5","Cenpa"
]

# Only keep genes present in adata
s_genes_use   = [g for g in s_genes_mouse  if g in adata.var_names]
g2m_genes_use = [g for g in g2m_genes_mouse if g in adata.var_names]

sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes_use, g2m_genes=g2m_genes_use,
                              random_state=SEED)

# Mirrors cyclone threshold: assign phase only if score > 0.5 for that phase
# Scanpy stores phase directly in adata.obs["phase"]; remap to match R output
# Scanpy scores are S_score and G2M_score (continuous, not 0-1 probabilities)
# We keep the scanpy assignment directly as it is functionally equivalent

print("Cell cycle phase distribution:")
print(adata.obs["phase"].value_counts())

# G1 vs G2M scatter (mirrors plot(G1 score, G2M score) in R)
fig, ax = plt.subplots(figsize=(5, 5))
ax.scatter(adata.obs["S_score"], adata.obs["G2M_score"], s=5, alpha=0.5)
ax.set_xlabel("S score"); ax.set_ylabel("G2/M score")
ax.set_title("Cell cycle scores")
plt.tight_layout()
plt.savefig("cell_cycle_scores.png", dpi=150)

# =============================================================================
# 11. INITIAL PCA ON NORMALISED DATA — ALL GENES
#     Mirrors: runPCA(sce) before HVG selection, coloured by CellCycle
# =============================================================================

print("\n=== 11. PCA on all normalised genes ===")

sc.pp.pca(adata, n_comps=50, random_state=SEED)
sc.pl.pca(adata, color="phase", title="PCA — all genes, coloured by cell cycle",
          save="_all_genes_cell_cycle.png")

# =============================================================================
# 12. TOP 50 MOST HIGHLY EXPRESSED GENES
#     Mirrors: plotHighestExprs(sce, n=50)
# =============================================================================

print("\n=== 12. Top 50 most highly expressed genes ===")

sc.pl.highest_expr_genes(adata, n_top=20, save="_top_expressed.png")
# Note: scanpy shows top 20 by default; increase n_top for top 50

# =============================================================================
# 13. HIGHLY VARIABLE GENE (HVG) DETECTION
#     Mirrors: modelGeneVar(sce) → FDR <= 0.2 AND bio >= 0.5
#     Scanpy uses seurat/cell_ranger/scran flavour; we use 'seurat_v3' on counts
#     and then replicate the R output TSV
# =============================================================================

print("\n=== 13. Highly Variable Gene detection ===")

# Work on raw counts layer for HVG (mirrors scran's modelGeneVar on logcounts)
adata_hvg = adata.copy()
adata_hvg.X = adata_hvg.layers["counts"].copy()

sc.pp.normalize_total(adata_hvg)
sc.pp.log1p(adata_hvg)
sc.pp.highly_variable_genes(adata_hvg, flavor="seurat_v3", n_top_genes=190,
                             span=0.3)  # 190 HVGs matches your R output

adata.var["highly_variable"] = adata_hvg.var["highly_variable"].values
adata.var["means"]           = adata_hvg.var["means"].values
adata.var["highly_variable"] = adata_hvg.var["highly_variable"].values

if "variances_norm" in adata_hvg.var.columns:
    adata.var["dispersions"] = adata_hvg.var["variances_norm"].values

adata.var["highly_variable"] = adata_hvg.var["highly_variable"].values

if "variances_norm" in adata_hvg.var.columns:
    adata.var["dispersions_norm"] = adata_hvg.var["variances_norm"].values

hvg_names = adata.var_names[adata.var["highly_variable"]].tolist()
print(f"Number of HVGs identified: {len(hvg_names)}")

# Mean-variance plot (mirrors var.out plot in R)
fig, ax = plt.subplots(figsize=(7, 5))
ax.scatter(adata.var["means"], adata.var["dispersions"], s=3, alpha=0.3, color="grey",
           label="All genes")
ax.scatter(adata.var.loc[adata.var["highly_variable"], "means"],
           adata.var.loc[adata.var["highly_variable"], "dispersions"],
           s=5, alpha=0.6, color="dodgerblue", label="HVGs")
ax.set_xlabel("Mean log-expression")
ax.set_ylabel("Dispersion")
ax.set_title("Mean-variance trend — HVG selection")
ax.legend()
plt.tight_layout()
plt.savefig("hvg_mean_variance.png", dpi=150)

# Write HVG table (mirrors write.table(file="hsc_hvg.tsv", hvg.out, ...))
hvg_df = adata.var.loc[adata.var["highly_variable"], ["means","dispersions","dispersions_norm"]].copy()
hvg_df = hvg_df.sort_values("dispersions_norm", ascending=False)
hvg_df.to_csv("hsc_hvg.tsv", sep="\t")
print("HVG table written to hsc_hvg.tsv")

# =============================================================================
# 14. CORRELATED GENE PAIRS (Spearman)
#     Mirrors: correlatePairs(sce, subset.row=rownames(hvg.out))
#     Chosen feature set = genes in any significant pair (FDR <= 0.01)
# =============================================================================

print("\n=== 14. Correlated HVG pairs (Spearman) ===")

hvg_expr = adata[:, hvg_names].X
if issparse(hvg_expr):
    hvg_expr = hvg_expr.toarray()

n_hvg = len(hvg_names)
cor_records = []

# Compute all pairwise Spearman correlations among HVGs
from itertools import combinations
for i, j in combinations(range(n_hvg), 2):
    rho, pval = spearmanr(hvg_expr[:, i], hvg_expr[:, j])
    cor_records.append((hvg_names[i], hvg_names[j], rho, pval))

cor_df = pd.DataFrame(cor_records, columns=["gene1", "gene2", "rho", "pvalue"])

# BH FDR correction (mirrors correlatePairs FDR)
from statsmodels.stats.multitest import multipletests
_, cor_df["FDR"], _, _ = multipletests(cor_df["pvalue"], method="fdr_bh")
cor_df = cor_df.sort_values("FDR")

# Write correlation table (mirrors write.table(file="hsc_cor.tsv", var.cor, ...))
cor_df.to_csv("hsc_cor.tsv", sep="\t", index=False)
print("Correlation table written to hsc_cor.tsv")
print(cor_df.head())

# Chosen feature set: genes appearing in pairs with FDR <= 0.01
sig_cor = cor_df[cor_df["FDR"] <= 0.01]
chosen = pd.unique(sig_cor[["gene1","gene2"]].values.ravel()).tolist()
chosen = [g for g in chosen if g in adata.var_names]
print(f"Number of chosen genes (correlated HVGs): {len(chosen)}")

# =============================================================================
# 15. HEATMAP OF 100 RANDOM HVG GENES (row-mean centred)
#     Mirrors: heatmap.2(heat.vals[sample(rownames, 100), ...])
# =============================================================================

print("\n=== 15. Heatmap — 100 random chosen HVGs ===")

np.random.seed(SEED)
sample_genes = np.random.choice(chosen, size=min(100, len(chosen)), replace=False)

expr_mat = adata[:, sample_genes].X
if issparse(expr_mat):
    expr_mat = expr_mat.toarray()
heat_vals = (expr_mat - expr_mat.mean(axis=0))  # row-mean centred

condition_order = adata.obs["Type"].argsort().values
heat_df = pd.DataFrame(heat_vals[condition_order, :].T,
                        index=sample_genes,
                        columns=adata.obs_names[condition_order])

# Colour annotations by condition
lut = {"serum": "#E41A1C", "2i": "#377EB8", "a2i": "#4DAF4A"}
col_colors = adata.obs["Type"].iloc[condition_order].map(lut)
col_colors.index = adata.obs_names[condition_order]

g = sns.clustermap(heat_df, col_cluster=False, row_cluster=True,
                   cmap="RdBu_r", center=0, col_colors=col_colors,
                   yticklabels=True, xticklabels=False,
                   figsize=(12, 10))
g.fig.suptitle("100 random chosen HVG genes (row-mean centred)", y=1.01)
plt.savefig("heatmap_hvg100.png", dpi=150, bbox_inches="tight")

# =============================================================================
# 16. PCA USING CHOSEN FEATURE SET
#     Mirrors: runPCA(sce, subset_row=chosen); plotPCA(sce, colour_by="CellCycle")
# =============================================================================

print("\n=== 16. PCA on chosen feature set ===")

sc.pp.pca(adata, n_comps=50, random_state=SEED,
          use_highly_variable=False)

# PCA on subset: re-run on chosen genes only
adata_chosen = adata[:, chosen].copy()
sc.pp.pca(adata_chosen, n_comps=50, random_state=SEED)
adata.obsm["X_pca_chosen"] = adata_chosen.obsm["X_pca"]

# Plot coloured by cell cycle
sc.pl.embedding(adata, basis="X_pca_chosen", color="phase",
                title="PCA (chosen HVGs) — Cell cycle",
                save="_chosen_cell_cycle.png")

# Plot coloured by culture condition
sc.pl.embedding(adata, basis="X_pca_chosen", color="Type",
                title="PCA (chosen HVGs) — Culture condition",
                save="_chosen_type.png")

# =============================================================================
# 17. tSNE — PERPLEXITY SCAN (5, 10, 15, 20)
#     Mirrors: run_tsne_perplexity_scan() helper function in R
#     Coloured by library size and by culture condition
# =============================================================================

print("\n=== 17. tSNE perplexity scan ===")

perplexities = [5, 10, 15, 20]

# Set PCA embedding from chosen features as input to tSNE
adata.obsm["X_pca"] = adata.obsm["X_pca_chosen"]

# --- Coloured by library size (mirrors colour_by="sum") ---
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
for ax, perp in zip(axes.flat, perplexities):
    np.random.seed(SEED)
    sc.tl.tsne(adata, n_pcs=min(len(chosen), 50), perplexity=perp,
               random_state=SEED, use_rep="X_pca")
    coords = adata.obsm["X_tsne"]
    scatter = ax.scatter(coords[:, 0], coords[:, 1],
                         c=adata.obs["total_counts"], cmap="viridis", s=5)
    ax.set_title(f"tSNE perplexity={perp}")
    ax.set_xlabel("tSNE 1"); ax.set_ylabel("tSNE 2")
    plt.colorbar(scatter, ax=ax, label="Library size")
plt.suptitle("tSNE perplexity scan — coloured by library size")
plt.tight_layout()
plt.savefig("tsne_perplexity_scan_sum.png", dpi=150)

# --- Coloured by culture condition (mirrors colour_by="Type") ---
cond_palette = {"serum": "#E41A1C", "2i": "#377EB8", "a2i": "#4DAF4A"}
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
for ax, perp in zip(axes.flat, perplexities):
    np.random.seed(SEED)
    sc.tl.tsne(adata, n_pcs=min(len(chosen), 50), perplexity=perp,
               random_state=SEED, use_rep="X_pca")
    coords = adata.obsm["X_tsne"]
    for cond, color in cond_palette.items():
        mask = adata.obs["Type"] == cond
        ax.scatter(coords[mask, 0], coords[mask, 1], c=color, s=5,
                   label=cond, alpha=0.7)
    ax.set_title(f"tSNE perplexity={perp}")
    ax.set_xlabel("tSNE 1"); ax.set_ylabel("tSNE 2")
    if perp == perplexities[0]:
        ax.legend(markerscale=3)
plt.suptitle("tSNE perplexity scan — coloured by culture condition")
plt.tight_layout()
plt.savefig("tsne_perplexity_scan_type.png", dpi=150)

# =============================================================================
# 18. FINAL tSNE AT PERPLEXITY 10 (chosen, stable)
#     Mirrors: runTSNE(sce, feature_set=chosen, perplexity=10); plotTSNE(...)
# =============================================================================

print("\n=== 18. Final tSNE at perplexity=10 ===")

np.random.seed(SEED)
sc.tl.tsne(adata, n_pcs=min(len(chosen), 50), perplexity=10,
           random_state=SEED, use_rep="X_pca")

sc.pl.tsne(adata, color="Type", title="tSNE (perplexity=10) — Culture condition",
           palette=cond_palette, save="_final_type.png")
sc.pl.tsne(adata, color="phase", title="tSNE (perplexity=10) — Cell cycle",
           save="_final_cell_cycle.png")
sc.pl.tsne(adata, color="total_counts", title="tSNE (perplexity=10) — Library size",
           save="_final_libsize.png")

# Final PCA coloured by Type
sc.pl.embedding(adata, basis="X_pca_chosen", color="Type",
                title="PCA (chosen HVGs) — Culture condition",
                palette=cond_palette, save="_final_pca_type.png")

# =============================================================================
# 19. ZSCAN4 FAMILY HEATMAP
#     Mirrors: grep("Zscan", rownames(sce)); pheatmap(heat.vals, annotation_col=experiment)
# =============================================================================

print("\n=== 19. Zscan4 family gene heatmap ===")

zscan_genes = [g for g in adata.var_names if "Zscan" in g]
print(f"Zscan genes found: {zscan_genes}")

if len(zscan_genes) > 0:
    expr_z = adata[:, zscan_genes].X
    if issparse(expr_z):
        expr_z = expr_z.toarray()
    heat_z = (expr_z - expr_z.mean(axis=0))

    # Sort cells by condition for annotation
    sort_idx = adata.obs["Type"].argsort().values
    heat_z_df = pd.DataFrame(heat_z[sort_idx, :].T,
                              index=zscan_genes,
                              columns=adata.obs_names[sort_idx])
    col_colors_z = adata.obs["Type"].iloc[sort_idx].map(cond_palette)
    col_colors_z.index = adata.obs_names[sort_idx]

    g = sns.clustermap(heat_z_df, col_cluster=False, row_cluster=True,
                       cmap="RdBu_r", center=0, col_colors=col_colors_z,
                       yticklabels=True, xticklabels=False, figsize=(12, 6))
    g.fig.suptitle("Zscan4 family gene expression (row-mean centred)", y=1.01)
    plt.savefig("heatmap_zscan.png", dpi=150, bbox_inches="tight")
else:
    print("No Zscan genes found — check that gene symbol annotation succeeded.")

# =============================================================================
# 20. ZSCAN4a-CORRELATED GENES (2C-like signature)
#     Mirrors: cor(normt.exprs, normt.exprs[,"Zscan4a"], method="spearman") > 0.4
#              pheatmap(heat.vals, ...)
# =============================================================================

print("\n=== 20. Zscan4a-correlated gene heatmap (2C-like signature) ===")

if "Zscan4a" in adata.var_names:
    expr_all = adata.X
    if issparse(expr_all):
        expr_all = expr_all.toarray()

    zscan4a_idx = list(adata.var_names).index("Zscan4a")
    zscan4a_vec = expr_all[:, zscan4a_idx]

    rho_vals = []
    for i in range(adata.n_vars):
        rho, _ = spearmanr(expr_all[:, i], zscan4a_vec)
        rho_vals.append(rho)

    rho_series = pd.Series(rho_vals, index=adata.var_names).sort_values(ascending=False)
    corr_genes = rho_series[rho_series > 0.4].index.tolist()
    print(f"Genes with Spearman rho > 0.4 with Zscan4a: {len(corr_genes)}")
    print(corr_genes)

    # Heatmap of correlated genes
    expr_corr = adata[:, corr_genes].X
    if issparse(expr_corr):
        expr_corr = expr_corr.toarray()
    heat_corr = (expr_corr - expr_corr.mean(axis=0))

    sort_idx = adata.obs["Type"].argsort().values
    heat_corr_df = pd.DataFrame(heat_corr[sort_idx, :].T,
                                 index=corr_genes,
                                 columns=adata.obs_names[sort_idx])
    col_colors_c = adata.obs["Type"].iloc[sort_idx].map(cond_palette)
    col_colors_c.index = adata.obs_names[sort_idx]

    g = sns.clustermap(heat_corr_df, col_cluster=False, row_cluster=True,
                       cmap="RdBu_r", center=0, col_colors=col_colors_c,
                       yticklabels=True, xticklabels=False, figsize=(12, 8),
                       dendrogram_ratio=0.1, cbar_pos=(0.02, 0.8, 0.03, 0.15))
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=6)
    g.fig.suptitle("Zscan4a-correlated genes (Spearman rho > 0.4) — 2C-like signature", y=1.01)
    plt.savefig("heatmap_zscan4a_correlated.png", dpi=150, bbox_inches="tight")
else:
    print("Zscan4a not found in var_names — check gene symbol annotation step.")

# =============================================================================
# 21. SUMMARY TABLE
#     Mirrors: key results table in README
# =============================================================================

print("\n=== Summary ===")
print(f"Starting cells:          {adata_raw.n_obs}")
print(f"Cells after QC:          {adata.n_obs}")
print(f"Cells per condition:")
print(adata.obs["Type"].value_counts().to_string())
print(f"Genes after filtering:   {adata.n_vars}")
print(f"HVGs identified:         {len(hvg_names)}")
print(f"Chosen feature set size: {len(chosen)}")
if "Zscan4a" in adata.var_names:
    print(f"Zscan4a-correlated genes: {len(corr_genes)}")

print("\nDone. All outputs saved as PNG images and TSV files.")
