# ===============================
# 1. Load Libraries
# ===============================
library(Matrix)
library(monocle3)

# ===============================
# 2. Load Input Data
# ===============================
expr_mat <- readRDS("expression_matrix.Rds")
gene_meta <- read.csv("gene_annotation.csv", header = TRUE)
cell_meta <- read.csv("cell_metadata.csv", header = TRUE)

# Ensure row/column names are properly set
rownames(gene_meta) <- gene_meta$gene_short_name
rownames(cell_meta) <- cell_meta$cell_id

# ===============================
# 3. Construct Monocle3 Object
# ===============================
cds <- new_cell_data_set(
  expr_mat,
  cell_metadata = cell_meta,
  gene_metadata = gene_meta
)

# ===============================
# 4. Preprocessing
# ===============================
cds <- preprocess_cds(cds, num_dim = 50)

# Optional: visualize variance explained
plot_pc_variance_explained(cds)

# ===============================
# 5. Dimensionality Reduction & Clustering
# ===============================
cds <- reduce_dimension(cds, reduction_method = "UMAP")
cds <- cluster_cells(cds)

# ===============================
# 6. Learn Trajectory Graph
# ===============================
cds <- learn_graph(cds)

# ===============================
# 7. Order Cells (Pseudotime)
# ===============================
# Option 1: Manual root selection (recommended)
# Replace "your_root_cluster" with known starting state
# cds <- order_cells(cds, root_cells = colnames(cds)[cds@clusters$UMAP$clusters == "your_root_cluster"])

# Option 2: Interactive selection
cds <- order_cells(cds)

# ===============================
# 8. Visualization
# ===============================
plot_cells(
  cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)