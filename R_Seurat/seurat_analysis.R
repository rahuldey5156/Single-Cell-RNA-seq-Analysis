# =============================================================================
# Single-Cell RNA-seq Analysis of Mouse ESCs — Seurat v5 (R)
# Reproduces: single_cell_es.Rmd (Bioconductor/scran/scater workflow)
# Dataset: E-MTAB-2600 | 869 cells, 3 conditions (serum, 2i, a2i)
# Input: counts_matrix.csv (genes x cells, Ensembl IDs as rownames)
# =============================================================================

options(Seurat.object.assay.version = "v5")

library(Seurat)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(gridExtra)
library(stringr)
library(org.Mm.eg.db)
library(AnnotationDbi)

SEED <- 100
set.seed(SEED)

# =============================================================================
# 1. LOAD DATA
# =============================================================================

cat("=== 1. Loading counts matrix ===\n")

counts_raw <- read.csv("counts_matrix.csv", row.names = 1, check.names = FALSE)
cat(sprintf("Raw matrix: %d genes x %d cells\n", nrow(counts_raw), ncol(counts_raw)))

counts_int <- round(counts_raw)

# =============================================================================
# 2. IDENTIFY ERCC SPIKES AND MITOCHONDRIAL GENES
# =============================================================================

cat("\n=== 2. Identifying ERCC spikes and mitochondrial genes ===\n")

is_ercc <- grepl("^ERCC", rownames(counts_int))
cat(sprintf("ERCC spikes: %d\n", sum(is_ercc)))

mito_ids <- c(
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
)
is_mito <- rownames(counts_int) %in% mito_ids
cat(sprintf("Mitochondrial genes: %d\n", sum(is_mito)))

# =============================================================================
# 3. CREATE SEURAT OBJECT
# =============================================================================

cat("\n=== 3. Creating Seurat object ===\n")

counts_main <- counts_int[!is_ercc, ]
counts_ercc <- counts_int[is_ercc, ]

seu <- CreateSeuratObject(
  counts       = counts_main,
  project      = "ES_scRNAseq",
  min.cells    = 0,
  min.features = 0
)

# Add ERCC as separate assay
seu[["ERCC"]] <- CreateAssayObject(counts = as.matrix(counts_ercc))

cat(sprintf("Seurat object: %d genes x %d cells\n", nrow(seu), ncol(seu)))

# =============================================================================
# 4. QC METRICS
# =============================================================================

cat("\n=== 4. Calculating QC metrics ===\n")

mito_genes_present <- rownames(seu)[rownames(seu) %in% mito_ids]
seu[["percent.mt"]] <- PercentageFeatureSet(seu, features = mito_genes_present)

total_counts        <- colSums(counts_int[, colnames(seu)])
ercc_counts         <- colSums(counts_ercc[, colnames(seu), drop = FALSE])
seu[["percent.ercc"]] <- (ercc_counts / total_counts) * 100

p1 <- ggplot(seu@meta.data, aes(x = nCount_RNA / 1e6)) +
  geom_histogram(bins = 20, fill = "grey80", colour = "white") +
  labs(x = "Library sizes (millions)", y = "Number of cells") +
  theme_classic()

p2 <- ggplot(seu@meta.data, aes(x = nFeature_RNA)) +
  geom_histogram(bins = 20, fill = "grey80", colour = "white") +
  labs(x = "Number of detected genes", y = "Number of cells") +
  theme_classic()

p3 <- ggplot(seu@meta.data, aes(x = percent.mt)) +
  geom_histogram(bins = 20, fill = "grey80", colour = "white") +
  labs(x = "Mitochondrial proportion (%)", y = "Number of cells") +
  theme_classic()

p4 <- ggplot(seu@meta.data, aes(x = percent.ercc)) +
  geom_histogram(bins = 20, fill = "grey80", colour = "white") +
  labs(x = "ERCC proportion (%)", y = "Number of cells") +
  theme_classic()

png("qc_before_filtering.png", width = 1400, height = 600, res = 150)
grid.arrange(p1, p2, p3, p4, ncol = 4)
dev.off()

# =============================================================================
# 5. CELL FILTERING — 3 MAD THRESHOLDS
# =============================================================================

cat("\n=== 5. Applying 3-MAD cell filters ===\n")

mad_filter <- function(x, nmads = 3, log = FALSE, direction = "both") {
  vals <- if (log) log(x + 1) else x
  med  <- median(vals, na.rm = TRUE)
  mad  <- median(abs(vals - med), na.rm = TRUE)
  lower <- med - nmads * mad
  upper <- med + nmads * mad
  if (log) lower <- exp(lower) - 1
  switch(direction,
    lower  = x < lower,
    higher = x > upper,
    both   = x < lower | x > upper
  )
}

libsize_drop <- mad_filter(seu$nCount_RNA,   log = TRUE,  direction = "lower")
feature_drop <- mad_filter(seu$nFeature_RNA, log = TRUE,  direction = "lower")
mito_drop    <- mad_filter(seu$percent.mt,   log = FALSE, direction = "higher")
spike_drop   <- mad_filter(seu$percent.ercc, log = FALSE, direction = "higher")

cat(sprintf("Removed by library size:  %d\n", sum(libsize_drop)))
cat(sprintf("Removed by feature count: %d\n", sum(feature_drop)))
cat(sprintf("Removed by mito:          %d\n", sum(mito_drop)))
cat(sprintf("Removed by ERCC spike:    %d\n", sum(spike_drop)))

keep_cells <- !(libsize_drop | feature_drop | mito_drop | spike_drop)
cat(sprintf("Cells remaining:          %d\n", sum(keep_cells)))

seu_raw <- seu
seu     <- seu[, keep_cells]
cat(sprintf("Post-QC Seurat: %d genes x %d cells\n", nrow(seu), ncol(seu)))

png("qc_violin_postfilter.png", width = 1400, height = 500, res = 150)
print(VlnPlot(seu, features = c("nCount_RNA","nFeature_RNA","percent.mt","percent.ercc"),
              ncol = 4, pt.size = 0.1) &
        theme(axis.text.x = element_blank()))
dev.off()

# =============================================================================
# 6. GENE ANNOTATION: Ensembl ID -> Gene Symbol
# =============================================================================

cat("\n=== 6. Annotating gene symbols ===\n")

non_ercc_ids <- rownames(seu)[!grepl("^ERCC", rownames(seu))]

anno <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys    = non_ercc_ids,
  keytype = "ENSEMBL",
  column  = "SYMBOL"
)

symbol_map <- setNames(anno$SYMBOL, anno$ENSEMBL)

old_names <- rownames(seu)
new_names <- old_names
non_ercc_idx <- !grepl("^ERCC", old_names)
new_names[non_ercc_idx] <- symbol_map[old_names[non_ercc_idx]]
new_names[is.na(new_names)] <- old_names[is.na(new_names)]
new_names <- make.names(new_names, unique = TRUE)

# Seurat v5 compatible gene renaming
counts_renamed <- LayerData(seu, assay = "RNA", layer = "counts")
rownames(counts_renamed) <- new_names

seu_renamed <- CreateSeuratObject(
  counts       = counts_renamed,
  project      = "ES_scRNAseq",
  min.cells    = 0,
  min.features = 0
)

# Transfer all metadata
seu_renamed@meta.data <- seu@meta.data

seu <- seu_renamed
cat(sprintf("Gene annotation done. %d genes.\n", nrow(seu)))

# =============================================================================
# 7. LOAD TARGETS FILE AND REMOVE BULK SAMPLES
# =============================================================================

cat("\n=== 7. Loading targets and removing bulk samples ===\n")

targets <- read.csv("E-MTAB-2600.targets.txt", sep = "\t")

targets_matched <- targets[targets$ERR %in% colnames(seu), ]
targets_matched <- targets_matched[match(colnames(seu), targets_matched$ERR), ]

seu$Type        <- targets_matched$Type
seu$ExtractName <- targets_matched$ExtractName

seu <- RenameCells(seu, new.names = targets_matched$ExtractName)

is_bulk <- grepl("bulk", colnames(seu))
cat(sprintf("Bulk samples found: %d\n", sum(is_bulk)))
seu <- seu[, !is_bulk]

cat(sprintf("Cells after bulk removal: %d\n", ncol(seu)))
print(table(seu$Type))

# =============================================================================
# 8. GENE FILTERING: expressed in >= 10 cells
# =============================================================================

cat("\n=== 8. Gene filtering (expressed in >= 10 cells) ===\n")

cat(sprintf("Genes before filtering: %d\n", nrow(seu)))

counts_mat   <- LayerData(seu, assay = "RNA", layer = "counts")
n_cells_expr <- rowSums(counts_mat > 0)
genes_keep   <- rownames(counts_mat)[n_cells_expr >= 10]
seu          <- seu[genes_keep, ]

cat(sprintf("Genes after filtering: %d\n", nrow(seu)))

ave_counts <- rowMeans(LayerData(seu, assay = "RNA", layer = "counts"))

png("gene_filtering.png", width = 1200, height = 500, res = 150)
par(mfrow = c(1, 2))
hist(log10(ave_counts + 1e-6), breaks = 100, col = "grey80",
     xlab = expression(Log[10] ~ "average count"), main = "")
abline(v = log10(1), col = "blue", lwd = 2, lty = 2)
smoothScatter(log10(ave_counts + 1e-6), n_cells_expr[genes_keep],
              xlab = expression(Log[10] ~ "average count"),
              ylab = "Number of expressing cells")
par(mfrow = c(1, 1))
dev.off()

# =============================================================================
# 9. NORMALISATION
# =============================================================================

cat("\n=== 9. Normalisation ===\n")

seu <- NormalizeData(seu, normalization.method = "LogNormalize",
                     scale.factor = 10000, verbose = FALSE)
cat("Log-normalisation applied\n")

seu$size_factor <- seu$nCount_RNA / median(seu$nCount_RNA)

png("size_factors.png", width = 800, height = 600, res = 150)
print(
  ggplot(seu@meta.data, aes(x = size_factor, y = nCount_RNA / 1e6)) +
    geom_point(size = 0.8, alpha = 0.5) +
    scale_x_log10() + scale_y_log10() +
    labs(x = "Size factor", y = "Library size (millions)",
         title = "Size factors vs library size") +
    theme_classic()
)
dev.off()

# =============================================================================
# 10. INITIAL PCA ON RAW LOG COUNTS (pre-normalisation)
# =============================================================================

cat("\n=== 10. Initial PCA on raw log counts ===\n")

seu_log_raw <- seu
raw_log <- log2(LayerData(seu_log_raw, assay = "RNA", layer = "counts") + 0.01)
seu_log_raw <- SetAssayData(seu_log_raw, assay = "RNA", layer = "data",
                             new.data = raw_log)

seu_log_raw <- FindVariableFeatures(seu_log_raw, nfeatures = 2000, verbose = FALSE)
seu_log_raw <- ScaleData(seu_log_raw, verbose = FALSE)
seu_log_raw <- RunPCA(seu_log_raw, npcs = 50, seed.use = SEED, verbose = FALSE)

png("pca_raw_log.png", width = 800, height = 700, res = 150)
print(DimPlot(seu_log_raw, reduction = "pca") +
        ggtitle("PCA — raw log counts (no normalisation)"))
dev.off()

# =============================================================================
# 11. CELL CYCLE SCORING
# =============================================================================

cat("\n=== 11. Cell cycle scoring ===\n")

s_genes_mouse   <- str_to_title(cc.genes$s.genes)
g2m_genes_mouse <- str_to_title(cc.genes$g2m.genes)

s_genes_use   <- s_genes_mouse[s_genes_mouse %in% rownames(seu)]
g2m_genes_use <- g2m_genes_mouse[g2m_genes_mouse %in% rownames(seu)]

seu <- CellCycleScoring(seu, s.features = s_genes_use,
                         g2m.features = g2m_genes_use,
                         set.ident = FALSE, seed = SEED)
seu$CellCycle <- seu$Phase

cat("Cell cycle distribution:\n")
print(table(seu$CellCycle))

png("cell_cycle_scores.png", width = 700, height = 600, res = 150)
print(
  ggplot(seu@meta.data, aes(x = S.Score, y = G2M.Score, colour = CellCycle)) +
    geom_point(size = 0.8, alpha = 0.7) +
    labs(x = "S score", y = "G2/M score", title = "Cell cycle scores") +
    theme_classic()
)
dev.off()

# =============================================================================
# 12. TOP 50 MOST HIGHLY EXPRESSED GENES
# =============================================================================

cat("\n=== 12. Top 50 most highly expressed genes ===\n")

norm_data  <- LayerData(seu, assay = "RNA", layer = "data")
mean_expr  <- rowMeans(norm_data)
top50      <- names(sort(mean_expr, decreasing = TRUE))[1:50]

top50_df <- data.frame(
  gene           = factor(top50, levels = rev(top50)),
  mean_logcounts = mean_expr[top50]
)

png("top50_expressed.png", width = 800, height = 1000, res = 150)
print(
  ggplot(top50_df, aes(x = gene, y = mean_logcounts)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(x = NULL, y = "Mean log-normalised expression",
         title = "Top 50 most highly expressed genes") +
    theme_classic(base_size = 8)
)
dev.off()

# =============================================================================
# 13. HIGHLY VARIABLE GENE DETECTION
# =============================================================================

cat("\n=== 13. Highly Variable Gene detection ===\n")

seu <- FindVariableFeatures(seu, selection.method = "vst",
                             nfeatures = 190, verbose = FALSE)
hvg_names <- VariableFeatures(seu)
cat(sprintf("HVGs identified: %d\n", length(hvg_names)))

png("hvg_plot.png", width = 900, height = 700, res = 150)
print(VariableFeaturePlot(seu) +
        ggtitle("Mean-variance trend — HVG selection") +
        theme_classic())
dev.off()

hvg_meta <- seu@assays$RNA@meta.data
if (!is.null(hvg_meta) && nrow(hvg_meta) > 0) {
  hvg_info <- hvg_meta[hvg_names, , drop = FALSE]
} else {
  hvg_info <- data.frame(gene = hvg_names)
}
write.table(hvg_info, file = "hsc_hvg.tsv", sep = "\t",
            quote = FALSE, col.names = NA)
cat("HVG table written to hsc_hvg.tsv\n")

# =============================================================================
# 14. CORRELATED HVG PAIRS (Spearman) — CHOSEN FEATURE SET
# =============================================================================

cat("\n=== 14. Correlated HVG pairs (Spearman) ===\n")

norm_data <- LayerData(seu, assay = "RNA", layer = "data")
hvg_expr  <- t(as.matrix(norm_data[hvg_names, ]))

gene_pairs <- combn(hvg_names, 2)

cat("Computing pairwise Spearman correlations...\n")
cor_vals <- apply(gene_pairs, 2, function(pair) {
  cor(hvg_expr[, pair[1]], hvg_expr[, pair[2]], method = "spearman")
})

n_cells <- nrow(hvg_expr)
t_stat  <- cor_vals * sqrt(n_cells - 2) / sqrt(1 - cor_vals^2)
p_vals  <- 2 * pt(-abs(t_stat), df = n_cells - 2)
fdr_vals <- p.adjust(p_vals, method = "BH")

cor_df <- data.frame(
  gene1  = gene_pairs[1, ],
  gene2  = gene_pairs[2, ],
  rho    = cor_vals,
  pvalue = p_vals,
  FDR    = fdr_vals,
  stringsAsFactors = FALSE
)
cor_df <- cor_df[order(cor_df$FDR), ]

write.table(cor_df, file = "hsc_cor.tsv", sep = "\t",
            quote = FALSE, row.names = FALSE)
cat("Correlation table written to hsc_cor.tsv\n")
print(head(cor_df))

sig_cor <- cor_df[cor_df$FDR <= 0.01, ]
chosen  <- unique(c(sig_cor$gene1, sig_cor$gene2))
chosen  <- chosen[!is.na(chosen) & chosen %in% rownames(seu)]
cat(sprintf("Chosen genes (FDR <= 0.01): %d\n", length(chosen)))

# =============================================================================
# 15. SCALE DATA ON CHOSEN FEATURES
# =============================================================================

cat("\n=== 15. Scaling data ===\n")
seu <- ScaleData(seu, features = chosen, verbose = FALSE)

# =============================================================================
# 16. PCA ON CHOSEN FEATURE SET
# =============================================================================

cat("\n=== 16. PCA on chosen feature set ===\n")

seu <- RunPCA(seu, features = chosen, npcs = 50,
              reduction.name = "pca_chosen",
              seed.use = SEED, verbose = FALSE)

pct_var <- seu@reductions$pca_chosen@stdev^2 /
           sum(seu@reductions$pca_chosen@stdev^2) * 100
cat(sprintf("Variance explained by PC1+PC2: %.1f%%\n", sum(pct_var[1:2])))

png("pca_cell_cycle.png", width = 800, height = 700, res = 150)
print(DimPlot(seu, reduction = "pca_chosen", group.by = "CellCycle",
              pt.size = 0.8) +
        ggtitle("PCA (chosen HVGs) — Cell cycle") + theme_classic())
dev.off()

png("pca_condition.png", width = 800, height = 700, res = 150)
print(DimPlot(seu, reduction = "pca_chosen", group.by = "Type", pt.size = 0.8,
              cols = c("2i" = "#377EB8", "a2i" = "#4DAF4A", "serum" = "#E41A1C")) +
        ggtitle("PCA (chosen HVGs) — Culture condition") + theme_classic())
dev.off()

# =============================================================================
# 17. HEATMAP OF 100 RANDOM CHOSEN HVGs
# =============================================================================

cat("\n=== 17. Heatmap — 100 random chosen HVGs ===\n")

set.seed(SEED)
sample_genes <- sample(chosen, size = min(100, length(chosen)))

norm_data  <- LayerData(seu, assay = "RNA", layer = "data")
expr_mat   <- as.matrix(norm_data[sample_genes, ])
heat_vals  <- expr_mat - rowMeans(expr_mat)

cell_order <- order(seu$Type)
ann_col <- data.frame(Condition = seu$Type[cell_order],
                       row.names = colnames(seu)[cell_order])
ann_colours <- list(Condition = c("2i"="#377EB8","a2i"="#4DAF4A","serum"="#E41A1C"))

png("heatmap_hvg100.png", width = 1200, height = 1400, res = 150)
pheatmap(heat_vals[, cell_order],
         cluster_cols = FALSE, cluster_rows = TRUE,
         annotation_col = ann_col, annotation_colors = ann_colours,
         show_colnames = FALSE, fontsize_row = 6,
         color = colorRampPalette(rev(brewer.pal(11,"RdBu")))(100),
         main = "100 random chosen HVGs (row-mean centred)")
dev.off()

# =============================================================================
# 18. tSNE PERPLEXITY SCAN (5, 10, 15, 20)
# =============================================================================

cat("\n=== 18. tSNE perplexity scan ===\n")

perplexities <- c(5, 10, 15, 20)
cond_cols    <- c("2i"="#377EB8","a2i"="#4DAF4A","serum"="#E41A1C")

plots_sum  <- list()
plots_type <- list()

for (perp in perplexities) {
  set.seed(SEED)
  rname <- paste0("tsne_p", perp)
  seu <- RunTSNE(seu, reduction = "pca_chosen",
                 dims = 1:min(30, ncol(seu@reductions$pca_chosen)),
                 perplexity = perp, reduction.name = rname,
                 seed.use = SEED, verbose = FALSE)

  plots_sum[[as.character(perp)]] <-
    FeaturePlot(seu, features = "nCount_RNA", reduction = rname, pt.size = 0.5) +
    ggtitle(paste("Perplexity", perp)) + theme_classic()

  plots_type[[as.character(perp)]] <-
    DimPlot(seu, group.by = "Type", reduction = rname,
            pt.size = 0.5, cols = cond_cols) +
    ggtitle(paste("Perplexity", perp)) + theme_classic()
}

png("tsne_perplexity_libsize.png", width = 1400, height = 1200, res = 150)
grid.arrange(grobs = plots_sum, ncol = 2,
             top = "tSNE perplexity scan — library size")
dev.off()

png("tsne_perplexity_condition.png", width = 1400, height = 1200, res = 150)
grid.arrange(grobs = plots_type, ncol = 2,
             top = "tSNE perplexity scan — culture condition")
dev.off()

# =============================================================================
# 19. FINAL tSNE AT PERPLEXITY 10
# =============================================================================

cat("\n=== 19. Final tSNE at perplexity=10 ===\n")

set.seed(SEED)
seu <- RunTSNE(seu, reduction = "pca_chosen",
               dims = 1:min(30, ncol(seu@reductions$pca_chosen)),
               perplexity = 10, reduction.name = "tsne",
               seed.use = SEED, verbose = FALSE)

png("tsne_condition.png", width = 800, height = 700, res = 150)
print(DimPlot(seu, reduction = "tsne", group.by = "Type",
              pt.size = 0.8, cols = cond_cols) +
        ggtitle("tSNE (perplexity=10) — Culture condition") + theme_classic())
dev.off()

png("tsne_cell_cycle.png", width = 800, height = 700, res = 150)
print(DimPlot(seu, reduction = "tsne", group.by = "CellCycle", pt.size = 0.8) +
        ggtitle("tSNE (perplexity=10) — Cell cycle") + theme_classic())
dev.off()

png("tsne_libsize.png", width = 800, height = 700, res = 150)
print(FeaturePlot(seu, features = "nCount_RNA", reduction = "tsne", pt.size = 0.8) +
        ggtitle("tSNE (perplexity=10) — Library size") + theme_classic())
dev.off()

# =============================================================================
# 20. ZSCAN4 FAMILY HEATMAP
# =============================================================================

cat("\n=== 20. Zscan4 family gene heatmap ===\n")

zscan_genes <- grep("Zscan", rownames(seu), value = TRUE)
cat(sprintf("Zscan genes found: %s\n", paste(zscan_genes, collapse = ", ")))

if (length(zscan_genes) > 0) {
  norm_data <- LayerData(seu, assay = "RNA", layer = "data")
  expr_z    <- as.matrix(norm_data[zscan_genes, ])
  heat_z    <- expr_z - rowMeans(expr_z)

  cell_order <- order(seu$Type)
  ann_col_z  <- data.frame(Condition = seu$Type[cell_order],
                            row.names = colnames(seu)[cell_order])

  png("heatmap_zscan.png", width = 1200, height = 700, res = 150)
  pheatmap(heat_z[, cell_order],
           cluster_cols = FALSE, cluster_rows = TRUE,
           annotation_col = ann_col_z, annotation_colors = ann_colours,
           show_colnames = FALSE, fontsize_row = 8,
           color = colorRampPalette(rev(brewer.pal(11,"RdBu")))(100),
           main = "Zscan4 family gene expression (row-mean centred)")
  dev.off()
}

# =============================================================================
# 21. ZSCAN4a-CORRELATED GENES (2C-like signature)
# =============================================================================

cat("\n=== 21. Zscan4a-correlated genes (2C-like signature) ===\n")

if ("Zscan4a" %in% rownames(seu)) {
  norm_data   <- as.matrix(LayerData(seu, assay = "RNA", layer = "data"))
  zscan4a_vec <- norm_data["Zscan4a", ]

  rho_vals <- apply(norm_data, 1, function(g) {
    cor(g, zscan4a_vec, method = "spearman")
  })

  rho_sorted  <- sort(rho_vals, decreasing = TRUE)
  corr_genes  <- names(rho_sorted[rho_sorted > 0.4])
  cat(sprintf("Genes with Spearman rho > 0.4 with Zscan4a: %d\n", length(corr_genes)))
  cat(paste(corr_genes, collapse = ", "), "\n")

  corr_present <- corr_genes[corr_genes %in% rownames(seu)]
  expr_corr    <- norm_data[corr_present, ]
  heat_corr    <- expr_corr - rowMeans(expr_corr)

  cell_order   <- order(seu$Type)
  ann_col_c    <- data.frame(Condition = seu$Type[cell_order],
                              row.names = colnames(seu)[cell_order])

  png("heatmap_zscan4a_correlated.png", width = 1200, height = 1000, res = 150)
  pheatmap(heat_corr[, cell_order],
           cluster_cols = FALSE, cluster_rows = TRUE,
           annotation_col = ann_col_c, annotation_colors = ann_colours,
           show_colnames = FALSE, fontsize_row = 6,
           color = colorRampPalette(rev(brewer.pal(11,"RdBu")))(100),
           main = "Zscan4a-correlated genes (rho > 0.4) — 2C-like signature")
  dev.off()
} else {
  cat("Zscan4a not found — check gene symbol annotation step.\n")
}

# =============================================================================
# 22. SUMMARY
# =============================================================================

cat("\n=== Summary ===\n")
cat(sprintf("Starting cells:            %d\n", ncol(seu_raw)))
cat(sprintf("Cells after QC:            %d\n", ncol(seu)))
cat("Cells per condition:\n")
print(table(seu$Type))
cat(sprintf("Genes after filtering:     %d\n", nrow(seu)))
cat(sprintf("HVGs identified:           %d\n", length(hvg_names)))
cat(sprintf("Chosen feature set size:   %d\n", length(chosen)))
if ("Zscan4a" %in% rownames(seu)) {
  cat(sprintf("Zscan4a-correlated genes:  %d\n", length(corr_genes)))
}
cat("\nDone. All plots saved as PNG files.\n")
