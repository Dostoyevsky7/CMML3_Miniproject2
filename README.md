# CMML3_Miniproject2

# scRNA-seq Clustering Benchmark Using Synthetic PBMC Data

This repository contains the R scripts used for the mini-project report:

**Benchmarking clustering methods on synthetic single-cell RNA-seq data with controlled cluster separability and cell-type imbalance**

## Project overview

Single-cell RNA sequencing (scRNA-seq) is widely used to identify cell types and cell states through clustering analysis. However, clustering performance can be influenced by the underlying data structure, including how well cell populations are separated and whether cell-type proportions are balanced.

This project uses **scDesign3** to generate synthetic PBMC-like scRNA-seq datasets with known ground-truth cell-type labels. The synthetic datasets were designed to test how two factors affect clustering reliability:

1. **Cluster separability**
   - High separability
   - Low separability

2. **Cluster size composition**
   - Balanced
   - Moderately imbalanced
   - Extremely imbalanced

Four clustering methods were benchmarked:

- Hierarchical clustering
- k-means clustering
- Leiden clustering
- Louvain clustering

The methods were evaluated using:

- Adjusted Rand Index (ARI)
- Normalized Mutual Information (NMI)
- Silhouette score
- UMAP visualization

## Repository structure

```text
.
├── scripts/
│   ├── 01_prepare_reference_data.R
│   ├── 02_generate_synthetic_data_scDesign3.R
│   ├── 03_clustering_benchmark.R
│   ├── 04_plot_benchmark_results.R
│   └── 05_plot_umap_examples.R
│
├── data/
│   └── README.md
│
├── results/
│   ├── clustering_metrics/
│   ├── confusion_matrices/
│   ├── umap_true_labels/
│   ├── umap_predicted_clusters/
│   └── figures/
│
└── README.md
