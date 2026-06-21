# ===============================================
# GO Enrichment Analysis for Top 300 Marker Genes
# ===============================================

# Load required packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(tidyr)  # Added for additional data manipulation

# 1. Load and filter marker gene list
markerList <- readRDS("brain.microglia.filter/markerList.rds")
filtered_markerList <- Filter(function(x) !is.null(x) && nrow(x) > 0, markerList)

# Check if any clusters remain after filtering
if (length(filtered_markerList) == 0) {
  stop("No valid clusters found in markerList")
}

# 2. Combine all clusters into one data.frame
combined_df <- bind_rows(lapply(names(filtered_markerList), function(clu) {
  df <- filtered_markerList[[clu]]
  df$cluster <- clu
  df
})) %>% as.data.frame()

# Check if combined_df has data
if (nrow(combined_df) == 0) {
  stop("Combined data frame is empty")
}

# 3. Select top 300 marker genes by absolute Log2FC or p-value adjusted?
# Option A: By Log2FC (original approach)
top300 <- combined_df %>%
  arrange(desc(Log2FC)) %>%
  distinct(name, .keep_all = TRUE) %>%
  slice_head(n = min(300, nrow(.)))  # Handle cases with <300 genes

# Optional: Add a message about how many genes were selected
message(paste("Selected", nrow(top300), "unique genes for enrichment analysis"))

# 4. Convert gene symbols to Entrez IDs
gene.df <- tryCatch({
  bitr(
    top300$name,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Hs.eg.db
  )
}, error = function(e) {
  stop("Gene ID conversion failed: ", e$message)
})

# Filter to keep only successfully converted genes
valid_genes <- top300 %>% filter(name %in% gene.df$SYMBOL)
message(paste("Successfully converted", nrow(gene.df), "out of", nrow(top300), "genes"))

# Check if we have enough genes for enrichment
if (nrow(gene.df) < 10) {
  stop("Insufficient genes converted to Entrez IDs for enrichment analysis")
}

# 5. GO enrichment (Biological Process) with additional parameters
ego <- enrichGO(
  gene          = gene.df$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pvalueCutoff  = 0.05,      # Stricter cutoff
  qvalueCutoff  = 0.1,       # Stricter cutoff
  readable      = TRUE,
  minGSSize     = 10,        # Minimum gene set size
  maxGSSize     = 500        # Maximum gene set size
)

# Check if enrichment found any terms
if (is.null(ego) || nrow(ego) == 0) {
  warning("No enrichment terms found with current cutoffs. Consider relaxing pvalueCutoff")
  # Optionally try with original cutoffs
  ego <- enrichGO(
    gene          = gene.df$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pvalueCutoff  = 0.1,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
}

# 6. Process and filter GO result
ego_df <- as.data.frame(ego)

# Take top 20 terms by adjusted p-value
top_terms <- ego_df %>%
  arrange(p.adjust) %>%
  slice_head(n = min(20, nrow(.)))  # Handle cases with <20 terms

# Add message about number of terms found
message(paste("Found", nrow(ego_df), "enriched terms, showing top", nrow(top_terms)))

# 7. Enhanced bubble chart plot
p <- ggplot(top_terms, aes(x = GeneRatio, y = reorder(Description, GeneRatio))) +
  geom_point(aes(size = Count, color = p.adjust), alpha = 0.8) +
  scale_color_gradient(
    low = "red", 
    high = "blue", 
    trans = "log10",
    name = "Adjusted\np-value"
  ) +
  scale_size_continuous(
    name = "Gene Count",
    range = c(3, 10)
  ) +
  labs(
    title  = paste("GO Enrichment Analysis - Top", nrow(top300), "Marker Genes"),
    subtitle = paste("Biological Process |", nrow(gene.df), "genes analyzed"),
    x      = "Gene Ratio",
    y      = NULL  # Cleaner look
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray50"),
    legend.position = "right",
    legend.box = "vertical",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", linetype = "dashed")
  )

# Print the plot object
print(p)

# 8. Save plot with high resolution
ggsave(
  "top300_marker_gene_GO_bubble.pdf", 
  plot = p, 
  width = 10, 
  height = 7,
  dpi = 300,
  device = "pdf"
)

# Also save as PNG for quick viewing
ggsave(
  "top300_marker_gene_GO_bubble.png", 
  plot = p, 
  width = 10, 
  height = 7,
  dpi = 300,
  device = "png"
)

# 9. Export enrichment results as CSV
write.csv(ego_df, "GO_enrichment_results.csv", row.names = FALSE)

# 10. Optional: Create a dotplot using clusterProfiler's built-in function
p_dot <- dotplot(ego, showCategory = 20) + 
  ggtitle("GO Enrichment Analysis (Top 300 Marker Genes)")
ggsave("GO_enrichment_dotplot.pdf", plot = p_dot, width = 10, height = 8)

# Session info for reproducibility
sink("GO_analysis_session_info.txt")
sessionInfo()
sink()

message("\nAnalysis completed successfully!")
message(paste("Output files generated:", 
              "\n- top300_marker_gene_GO_bubble.pdf/png",
              "\n- GO_enrichment_results.csv",
              "\n- GO_enrichment_dotplot.pdf",
              "\n- GO_analysis_session_info.txt"))