# Single-Cell RNA-seq Analysis of Mouse Spleen Using Seurat

## Overview

This project performs an **end-to-end single-cell RNA sequencing (scRNA-seq) analysis** using the **Seurat package in R**. The dataset used is a **10x Genomics mouse spleen scRNA-seq dataset**. The goal is to identify cellular heterogeneity, detect clusters of cells, and assign biological identities based on marker gene expression.

The workflow includes:

* Data loading from a `.h5` file
* Quality control and filtering
* Normalization and feature selection
* Dimensionality reduction (PCA and UMAP)
* Cell clustering
* Marker gene identification
* Differential gene expression analysis

This workflow follows a standard **QC → clustering → marker discovery → biological interpretation pipeline** for single-cell analysis. 

---

# Results and Interpretation

## 1. PCA and Elbow Plot Interpretation

The **Elbow plot** shows the amount of variance explained by each principal component.

Interpretation:

* The steep decrease occurs until around **PC 10–15**.
* After PC15, the variance explained stabilizes.
* Therefore **15 principal components were selected** for downstream analysis.

This indicates that the **major biological variation in the dataset is captured within the first 15 PCs**.

---

# 2. UMAP Visualization and Clustering

The **UMAP plot** visualizes cells in two dimensions based on gene expression similarity.

Interpretation:

* The analysis identified **~17 clusters of cells**.
* Each cluster represents **a group of transcriptionally similar cells**.
* Cells within a cluster likely represent **a specific immune cell population or subtype** in the spleen.
* The clusters are well separated, indicating **clear transcriptional differences between cell populations**.

Since the spleen contains multiple immune cells, this clustering likely represents populations such as:

* **T cells**
* **B cells**
* **Monocytes / macrophages**
* **Other immune cell subsets**

---

# 3. Marker Gene Heatmap

The **heatmap of top marker genes per cluster** shows genes that are highly expressed in specific clusters.

Interpretation:

* Each vertical block represents a **cluster**.
* Each row represents a **gene**.
* Yellow indicates **high expression**, while purple indicates **low expression**.

Observation:

* Distinct clusters show **unique gene expression signatures**.
* These marker genes allow identification of **cell type identities**.

This confirms that the clustering represents **biologically meaningful cell populations rather than random groupings**.

---

# 4. Dot Plot of Marker Genes

The dot plot shows expression of known immune marker genes:

| Gene  | Cell Type               |
| ----- | ----------------------- |
| Cd3d  | T cells                 |
| Ms4a1 | B cells                 |
| Lyz2  | Monocytes / macrophages |

Interpretation:

* **Cd3d** is strongly expressed in cluster 1 → indicates **T cell population**.
* **Ms4a1** shows high expression in cluster 2 → indicates **B cell population**.
* **Lyz2** is expressed in multiple clusters → indicates **myeloid lineage cells such as monocytes or macrophages**.

Dot size indicates **percentage of cells expressing the gene**, and color intensity indicates **average expression level**.

---

# 5. Feature Plot Interpretation

Feature plots visualize gene expression across the UMAP embedding.

Observations:

### Cd3d

* Highly expressed in a specific cluster.
* Confirms the presence of **T cells**.

### Ms4a1

* Expression localized in another cluster.
* Confirms **B cell population**.

These spatial expression patterns validate the **biological identity of clusters**.

---

# 6. Differential Expression Analysis

Differential expression analysis was performed between **cluster 0 and cluster 1**.

Purpose:

* Identify genes that distinguish these two cell populations.

Interpretation:

* Genes with high **log fold change** are more expressed in one cluster.
* These genes may represent **functional differences between immune cell types**.

The results were saved as:

```
DE_results.csv
```

These genes could represent **cell-type-specific markers or functional regulators**.

---

# Biological Insight

This analysis demonstrates that **mouse spleen tissue contains multiple immune cell populations**, which can be separated based on transcriptional profiles.

Key findings:

* Approximately **17 transcriptionally distinct clusters** were identified.
* Marker genes confirmed **major immune cell types** including:

  * **T cells (Cd3d)**
  * **B cells (Ms4a1)**
  * **Myeloid cells (Lyz2)**
* Differential gene expression further highlights **molecular differences between cell populations**.

This highlights the power of **single-cell RNA sequencing to resolve cellular heterogeneity within immune tissues**.

---

# Output Files

The analysis generated the following outputs:

```
cluster_markers.csv     # Marker genes for each cluster
DE_results.csv          # Differential expression results
UMAP plot               # Cell clustering visualization
Heatmap                 # Top marker genes
DotPlot                 # Marker gene expression
FeaturePlot             # Spatial gene expression
```

---
# Project Structure
scRNA-seq-Seurat-analysis
│
├── script
│   └── analysis.R
│
├── figures
│   ├── elbow_plot.png
│   ├── umap_clusters.png
│   ├── heatmap.png
│   ├── dotplot.png
│   └── featureplot.png
│
├── tables
│   ├── cluster_markers.csv
│   └── DE_results.csv
│
├── report
│   └── scRNAseq_report.pdf
│
└── README.md
# Tools Used

* **R**
* **Seurat**
* **tidyverse**
* **10x Genomics dataset**

---
