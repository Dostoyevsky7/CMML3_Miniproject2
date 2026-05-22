# =========================================================
# scDesign3 simulation script
# PBMC3K synthetic datasets
# Design:
#   2 separability settings:
#     1. high separability: B_cell, CD14_Monocyte, NK_cell
#     2. low separability: Naive_CD4_T, CD4_T_memory, CD8_T
#
#   3 composition settings:
#     1. balanced: 1000 / 1000 / 1000
#     2. moderate imbalance: 1800 / 900 / 300
#     3. extreme imbalance: 2400 / 540 / 60
#
#   Each condition is generated with multiple technical/simulation replicates.
# =========================================================


# =========================================================
# 0. Load packages
# =========================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(scDesign3)
  library(SingleCellExperiment)
  library(scran)
  library(ggplot2)
})


# =========================================================
# 1. Parameters
# =========================================================

# Input 10X directory
data_dir <- "E:/CMML/miniproject2/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19/"

# Output directory
out_dir <- "scdesign3_synthetic_replicates"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

rds_dir <- file.path(out_dir, "rds_files")
plot_dir <- file.path(out_dir, "quick_umap_plots")

dir.create(rds_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Number of synthetic replicates per condition
n_reps <- 5

# Number of HVGs used for scDesign3 fitting
n_hvg_scdesign3 <- 1000

# Total synthetic cells per dataset
ncell_target <- 3000

# Use one core for stability
n_cores <- 1


# =========================================================
# 2. Import PBMC3K data and create Seurat object
# =========================================================

pbmc.data <- Read10X(data.dir = data_dir)

seu <- CreateSeuratObject(
  counts = pbmc.data,
  project = "PBMC3K"
)


# =========================================================
# 3. Standard Seurat preprocessing
# =========================================================

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 2500 &
    percent.mt < 5
)

seu <- NormalizeData(seu, verbose = FALSE)

seu <- FindVariableFeatures(
  seu,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, verbose = FALSE)

seu <- FindNeighbors(seu, dims = 1:10, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.5, verbose = FALSE)

seu <- RunUMAP(seu, dims = 1:10, verbose = FALSE)


# Optional marker check plot
p_marker <- FeaturePlot(
  seu,
  features = c("CD3D", "MS4A1", "LYZ", "NKG7")
)

ggsave(
  filename = file.path(plot_dir, "PBMC_marker_FeaturePlot.png"),
  plot = p_marker,
  width = 10,
  height = 8,
  dpi = 300
)


# =========================================================
# 4. Find markers and manually annotate clusters
# =========================================================

markers <- FindAllMarkers(
  seu,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

top30_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 30)

write.csv(
  markers,
  file = file.path(out_dir, "PBMC3K_all_markers.csv"),
  row.names = FALSE
)

write.csv(
  top30_markers,
  file = file.path(out_dir, "PBMC3K_top30_markers_per_cluster.csv"),
  row.names = FALSE
)


# Manual annotation
# This assumes your clustering generates clusters 0-8.
# If your cluster numbers change, this label vector must be checked again.
new_labels <- c(
  "Naive_CD4_T",      # 0
  "CD14_Monocyte",    # 1
  "CD4_T_memory",     # 2
  "B_cell",           # 3
  "CD8_T",            # 4
  "FCGR3A_Monocyte",  # 5
  "NK_cell",          # 6
  "Dendritic_cell",   # 7
  "Platelet"          # 8
)

if (length(new_labels) != length(levels(seu))) {
  stop(
    "The number of manual labels does not match the number of Seurat clusters. ",
    "Please check levels(seu) and update new_labels."
  )
}

names(new_labels) <- levels(seu)

seu <- RenameIdents(seu, new_labels)
seu$celltype <- factor(Idents(seu))

print(table(seu$celltype))

p_celltype <- DimPlot(
  seu,
  group.by = "celltype",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Annotated PBMC3K reference")

ggsave(
  filename = file.path(plot_dir, "PBMC3K_annotated_celltypes.png"),
  plot = p_celltype,
  width = 8,
  height = 6,
  dpi = 300
)


# =========================================================
# 5. Define high- and low-separability cell-type groups
# =========================================================

# High separability: different major lineages
high_types <- c("B_cell", "CD14_Monocyte", "NK_cell")

# Low separability: T-cell lineage
low_types <- c("Naive_CD4_T", "CD4_T_memory", "CD8_T")

seu_high <- subset(seu, subset = celltype %in% high_types)
seu_low  <- subset(seu, subset = celltype %in% low_types)

seu_high$celltype <- droplevels(seu_high$celltype)
seu_low$celltype  <- droplevels(seu_low$celltype)

cat("\nHigh-separability reference cell numbers:\n")
print(table(seu_high$celltype))

cat("\nLow-separability reference cell numbers:\n")
print(table(seu_low$celltype))


# =========================================================
# 6. Prepare SingleCellExperiment objects for scDesign3
# =========================================================

prep_for_scdesign3 <- function(seu_obj, n_hvg = 1000) {
  
  DefaultAssay(seu_obj) <- "RNA"
  
  seu_obj <- NormalizeData(seu_obj, verbose = FALSE)
  
  seu_obj <- FindVariableFeatures(
    seu_obj,
    selection.method = "vst",
    nfeatures = n_hvg,
    verbose = FALSE
  )
  
  hvgs <- VariableFeatures(seu_obj)
  seu_obj <- subset(seu_obj, features = hvgs)
  
  sce <- as.SingleCellExperiment(seu_obj)
  
  colData(sce)$cell_type <- factor(seu_obj$celltype)
  
  return(sce)
}

sce_high <- prep_for_scdesign3(
  seu_obj = seu_high,
  n_hvg = n_hvg_scdesign3
)

sce_low <- prep_for_scdesign3(
  seu_obj = seu_low,
  n_hvg = n_hvg_scdesign3
)

cat("\nLevels in sce_high:\n")
print(levels(colData(sce_high)$cell_type))

cat("\nLevels in sce_low:\n")
print(levels(colData(sce_low)$cell_type))


# =========================================================
# 7. Define target composition
# =========================================================

balanced_counts <- c(1000, 1000, 1000)
moderate_counts <- c(1800, 900, 300)
extreme_counts  <- c(2400, 540, 60)

stopifnot(sum(balanced_counts) == ncell_target)
stopifnot(sum(moderate_counts) == ncell_target)
stopifnot(sum(extreme_counts) == ncell_target)


# =========================================================
# 8. Fit scDesign3 reference model
# =========================================================

fit_scdesign3_reference <- function(sce_ref, n_cores = 1) {
  
  data_obj <- construct_data(
    sce = sce_ref,
    assay_use = "counts",
    celltype = "cell_type",
    pseudotime = NULL,
    spatial = NULL,
    other_covariates = NULL,
    ncell = ncol(sce_ref),
    corr_by = "1"
  )
  
  marginal_obj <- fit_marginal(
    data = data_obj,
    predictor = "gene",
    mu_formula = "cell_type",
    sigma_formula = "1",
    family_use = "nb",
    n_cores = n_cores,
    usebam = FALSE
  )
  
  copula_obj <- fit_copula(
    sce = sce_ref,
    assay_use = "counts",
    marginal_list = marginal_obj,
    family_use = "nb",
    copula = "gaussian",
    n_cores = n_cores,
    input_data = data_obj$dat
  )
  
  return(
    list(
      data_obj = data_obj,
      marginal_obj = marginal_obj,
      copula_obj = copula_obj
    )
  )
}

cat("\nFitting scDesign3 model for high-separability reference...\n")
fit_high <- fit_scdesign3_reference(
  sce_ref = sce_high,
  n_cores = n_cores
)

cat("\nFitting scDesign3 model for low-separability reference...\n")
fit_low <- fit_scdesign3_reference(
  sce_ref = sce_low,
  n_cores = n_cores
)

saveRDS(
  fit_high,
  file = file.path(out_dir, "fit_high_reference_scDesign3.rds")
)

saveRDS(
  fit_low,
  file = file.path(out_dir, "fit_low_reference_scDesign3.rds")
)


# =========================================================
# 9. Functions for synthetic data generation
# =========================================================

make_target_covariate <- function(
    sce_ref,
    target_counts,
    ncell = 3000,
    seed = 123
) {
  
  stopifnot(sum(target_counts) == ncell)
  
  tmp <- construct_data(
    sce = sce_ref,
    assay_use = "counts",
    celltype = "cell_type",
    pseudotime = NULL,
    spatial = NULL,
    other_covariates = NULL,
    ncell = ncell,
    corr_by = "1"
  )
  
  new_cov <- as.data.frame(tmp$newCovariate)
  
  lvls <- levels(colData(sce_ref)$cell_type)
  
  if (length(target_counts) != length(lvls)) {
    stop("target_counts length does not match the number of cell_type levels.")
  }
  
  target_labels <- rep(lvls, times = target_counts)
  
  set.seed(seed)
  target_labels <- sample(target_labels, length(target_labels), replace = FALSE)
  
  new_cov$cell_type <- factor(target_labels, levels = lvls)
  
  return(new_cov)
}


simulate_one_dataset <- function(
    sce_ref,
    fit_obj,
    target_counts,
    seed = 123,
    n_cores = 1,
    ncell = 3000
) {
  
  target_cov <- make_target_covariate(
    sce_ref = sce_ref,
    target_counts = target_counts,
    ncell = ncell,
    seed = seed
  )
  
  para_obj <- extract_para(
    sce = sce_ref,
    marginal_list = fit_obj$marginal_obj,
    n_cores = n_cores,
    family_use = "nb",
    new_covariate = target_cov,
    data = fit_obj$data_obj$dat
  )
  
  set.seed(seed)
  
  new_count <- simu_new(
    sce = sce_ref,
    assay_use = "counts",
    mean_mat = para_obj$mean_mat,
    sigma_mat = para_obj$sigma_mat,
    zero_mat = para_obj$zero_mat,
    quantile_mat = NULL,
    copula_list = fit_obj$copula_obj$copula_list,
    n_cores = n_cores,
    family_use = "nb",
    input_data = fit_obj$data_obj$dat,
    new_covariate = target_cov,
    important_feature = fit_obj$copula_obj$important_feature,
    filtered_gene = fit_obj$data_obj$filtered_gene
  )
  
  return(
    list(
      counts = new_count,
      meta = target_cov
    )
  )
}


# =========================================================
# 10. Define all simulation conditions
# =========================================================

condition_list <- list(
  
  S1_high_balanced = list(
    separability = "high",
    composition = "balanced",
    sce_ref = sce_high,
    fit_obj = fit_high,
    target_counts = balanced_counts,
    seed_start = 1000
  ),
  
  S2_high_moderate = list(
    separability = "high",
    composition = "moderate",
    sce_ref = sce_high,
    fit_obj = fit_high,
    target_counts = moderate_counts,
    seed_start = 2000
  ),
  
  S3_high_extreme = list(
    separability = "high",
    composition = "extreme",
    sce_ref = sce_high,
    fit_obj = fit_high,
    target_counts = extreme_counts,
    seed_start = 3000
  ),
  
  S4_low_balanced = list(
    separability = "low",
    composition = "balanced",
    sce_ref = sce_low,
    fit_obj = fit_low,
    target_counts = balanced_counts,
    seed_start = 4000
  ),
  
  S5_low_moderate = list(
    separability = "low",
    composition = "moderate",
    sce_ref = sce_low,
    fit_obj = fit_low,
    target_counts = moderate_counts,
    seed_start = 5000
  ),
  
  S6_low_extreme = list(
    separability = "low",
    composition = "extreme",
    sce_ref = sce_low,
    fit_obj = fit_low,
    target_counts = extreme_counts,
    seed_start = 6000
  )
)


# =========================================================
# 11. Generate technical/simulation replicates
# =========================================================

all_sim_info <- list()

for (condition_name in names(condition_list)) {
  
  cfg <- condition_list[[condition_name]]
  
  cat("\n=================================================\n")
  cat("Generating condition:", condition_name, "\n")
  cat("Separability:", cfg$separability, "\n")
  cat("Composition:", cfg$composition, "\n")
  cat("Target counts:\n")
  print(cfg$target_counts)
  cat("=================================================\n")
  
  for (rep_i in seq_len(n_reps)) {
    
    seed_i <- cfg$seed_start + rep_i
    
    cat("\nReplicate", rep_i, "using seed", seed_i, "\n")
    
    sim_obj <- simulate_one_dataset(
      sce_ref = cfg$sce_ref,
      fit_obj = cfg$fit_obj,
      target_counts = cfg$target_counts,
      seed = seed_i,
      n_cores = n_cores,
      ncell = ncell_target
    )
    
    sim_obj$meta$condition <- condition_name
    sim_obj$meta$separability <- cfg$separability
    sim_obj$meta$composition <- cfg$composition
    sim_obj$meta$replicate <- rep_i
    sim_obj$meta$seed <- seed_i
    
    obj_name <- paste0(condition_name, "_rep", rep_i)
    
    saveRDS(
      sim_obj,
      file = file.path(rds_dir, paste0(obj_name, ".rds"))
    )
    
    celltype_table <- as.data.frame(table(sim_obj$meta$cell_type))
    colnames(celltype_table) <- c("cell_type", "n_cells")
    celltype_table$condition <- condition_name
    celltype_table$separability <- cfg$separability
    celltype_table$composition <- cfg$composition
    celltype_table$replicate <- rep_i
    celltype_table$seed <- seed_i
    
    all_sim_info[[obj_name]] <- celltype_table
    
    cat("Cell-type composition for", obj_name, ":\n")
    print(table(sim_obj$meta$cell_type))
  }
}

sim_info_df <- bind_rows(all_sim_info)

write.csv(
  sim_info_df,
  file = file.path(out_dir, "synthetic_dataset_celltype_composition.csv"),
  row.names = FALSE
)


# =========================================================
# 12. Convert one synthetic object to Seurat
# =========================================================

sim_to_seurat <- function(sim_obj, project_name) {
  
  seu_sim <- CreateSeuratObject(
    counts = sim_obj$counts,
    project = project_name
  )
  
  seu_sim$celltype <- factor(sim_obj$meta$cell_type)
  seu_sim$condition <- sim_obj$meta$condition
  seu_sim$separability <- sim_obj$meta$separability
  seu_sim$composition <- sim_obj$meta$composition
  seu_sim$replicate <- sim_obj$meta$replicate
  seu_sim$seed <- sim_obj$meta$seed
  
  return(seu_sim)
}


# =========================================================
# 13. Optional quick UMAP plots for all synthetic replicates
# =========================================================

quick_plot <- function(seu_obj, plot_title = NULL) {
  
  seu_obj <- NormalizeData(seu_obj, verbose = FALSE)
  
  seu_obj <- FindVariableFeatures(
    seu_obj,
    selection.method = "vst",
    nfeatures = 1000,
    verbose = FALSE
  )
  
  seu_obj <- ScaleData(seu_obj, verbose = FALSE)
  seu_obj <- RunPCA(seu_obj, verbose = FALSE)
  
  seu_obj <- RunUMAP(
    seu_obj,
    dims = 1:20,
    verbose = FALSE
  )
  
  p <- DimPlot(
    seu_obj,
    group.by = "celltype",
    label = TRUE,
    repel = TRUE
  ) +
    ggtitle(plot_title)
  
  return(p)
}


# Generate quick UMAP plots for every replicate
for (condition_name in names(condition_list)) {
  
  for (rep_i in seq_len(n_reps)) {
    
    obj_name <- paste0(condition_name, "_rep", rep_i)
    rds_file <- file.path(rds_dir, paste0(obj_name, ".rds"))
    
    sim_obj <- readRDS(rds_file)
    seu_sim <- sim_to_seurat(sim_obj, project_name = obj_name)
    
    p <- quick_plot(
      seu_obj = seu_sim,
      plot_title = obj_name
    )
    
    ggsave(
      filename = file.path(plot_dir, paste0(obj_name, "_UMAP_celltype.png")),
      plot = p,
      width = 7,
      height = 6,
      dpi = 300
    )
  }
}


# =========================================================
# 14. Save session information
# =========================================================

sink(file.path(out_dir, "sessionInfo.txt"))
sessionInfo()
sink()

cat("\nAll synthetic datasets generated successfully.\n")
cat("RDS files saved in:\n")
cat(rds_dir, "\n")
cat("Quick UMAP plots saved in:\n")
cat(plot_dir, "\n")





p_recall <- ggplot(recall_df, aes(x = true_celltype, y = recall, fill = method)) +
  geom_col(position = position_dodge(width = 0.8)) +
  facet_grid(separability ~ composition, scales = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Per-cell-type recall across scenarios",
    x = "True cell type",
    y = "Recall"
  )
