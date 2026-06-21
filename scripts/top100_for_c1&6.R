## ------------------------------------------------------------
## Print top N marker genes for clusters with few markers
## ------------------------------------------------------------

markersGS <- readRDS("brain.microglia.filter/markersGS_7.2.rds")
markerList <- getMarkers(markersGS)

target_clusters <- c("C1", "C6")
top_n <- 100
top_markers_list <- list()

for (cluster in target_clusters) {
  df <- markerList[[cluster]]
  
  if (is.null(df)) {
    cat("Cluster", cluster, "is NULL. Skipped.\n")
    next
  }
  
  # Automatically select sorting column
  if ("FDR" %in% colnames(df)) {
    sort_col <- "FDR"
  } else if ("Pval" %in% colnames(df)) {
    sort_col <- "Pval"
  } else {
    stop(paste("No FDR or Pval column found for cluster", cluster))
  }
  
  # Sort by significance
  df_sorted <- df[order(df[[sort_col]]), ]
  
  # Select top N genes (or all if fewer than top_n)
  top_genes <- head(df_sorted, min(nrow(df_sorted), top_n))
  top_markers_list[[cluster]] <- top_genes
  
  # Print top 10 genes for quick inspection
  cat("\nTop genes for", cluster, ":\n")
  print(head(top_genes[, c("name", sort_col)], 10))
}

## ------------------------------------------------------------
## View top genes (unsorted) for all clusters
## ------------------------------------------------------------

for (cluster in names(markerList)) {
  df <- markerList[[cluster]]
  
  if (is.null(df)) {
    cat("Cluster", cluster, "is NULL. Skipped.\n")
    next
  }
  
  # Extract first 10 genes without sorting
  top_genes <- head(df$name, 10)
  
  cat("\nTop 10 genes for", cluster, ":\n")
  print(top_genes)
}