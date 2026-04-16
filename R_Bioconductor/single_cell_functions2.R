print("Just set global variable and functions here...")

  print("Set globals")
  #load colour scales
  COLOURS_CONT <- scale_colour_gradient2()
  COLOURS_GLOBAL <-c("Blue","Yellow","Red","Green","Orange","Plum2","Purple","Brown","Cyan","Green4","Red4","Coral","Burlywood4")
  COLOURS_DIS <-scale_fill_manual(values = COLOURS_GLOBAL,aesthetics=c("fill","colour"))
  COLOURS_TWO <-scale_fill_manual(values = COLOURS_GLOBAL[c(1,5)],aesthetics=c("fill","colour"))
  MY_SEED = 100
  FONTSIZE <- theme(axis.text=element_text(size=8), axis.title=element_text(size=12))


# -----------------------------------------------------------------------------
# run_tsne_perplexity_scan()
#
# Runs tSNE at multiple perplexity values and returns a list of ggplot objects.
# This replaces repetitive copy-pasted perplexity scan code blocks in the Rmd.
#
# Arguments:
#   sce_obj      - SingleCellExperiment object
#   feature_set  - character vector of genes to use for tSNE
#   colour_by    - metadata column name to colour points by (e.g. "sum", "Type")
#   perplexities - numeric vector of perplexity values to test (default: 5,10,15,20)
#   labels       - optional character vector of plot titles (one per perplexity value)
#
# Returns:
#   A list of ggplot objects, one per perplexity value
# -----------------------------------------------------------------------------
run_tsne_perplexity_scan <- function(sce_obj, feature_set, colour_by,
                                      perplexities = c(5, 10, 15, 20),
                                      labels = NULL) {
  plots <- lapply(seq_along(perplexities), function(i) {
    p <- perplexities[i]
    set.seed(MY_SEED)
    sce_obj <- runTSNE(sce_obj, feature_set = feature_set, perplexity = p)
    title <- if (!is.null(labels)) labels[i] else paste("Perplexity", p)
    plotTSNE(sce_obj, colour_by = colour_by) +
      FONTSIZE +
      ggtitle(title)
  })
  return(plots)
}
  
  
  

  
  
  
















