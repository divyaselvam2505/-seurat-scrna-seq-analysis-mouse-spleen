library(Seurat)
library(SeuratDisk)
install.packages("tidyverse")
library(tidyverse)
data<-Read10X_h5("C:\Users\jeeva\OneDrive\Documents\10k_mouse_spleen_scFFPE_singleplex_10k_mouse_spleen_scFFPE_singleplex_count_sample_filtered_feature_bc_matrix.h5")
library(Seurat)
library(SeuratDisk)

data <- Read10X_h5("C:/Users/jeeva/OneDrive/Documents/10k_mouse_spleen_scFFPE_singleplex_10k_mouse_spleen_scFFPE_singleplex_count_sample_filtered_feature_bc_matrix.h5")
str(data)
seurat.obj <- CreateSeuratObject(
  counts = data,
  project = "Mouse_Spleen",
  min.cells = 3,
  min.features = 200
)
seurat.obj
seurat.obj[["percent.mt"]] <- PercentageFeatureSet(seurat.obj, pattern = "^mt-")

VlnPlot(seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
        seurat.obj <- subset(seurat.obj,
                             subset = nFeature_RNA > 200 &
                               nFeature_RNA < 2500 &
                               percent.mt < 5)
seurat.obj <- NormalizeData(seurat.obj)

seurat.obj <- FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(seurat.obj)
seurat.obj <- ScaleData(seurat.obj, features = all.genes)

seurat.obj <- RunPCA(seurat.obj, features = VariableFeatures(object = seurat.obj))

ElbowPlot(seurat.obj)

seurat.obj <- FindNeighbors(seurat.obj, dims = 1:15)
seurat.obj <- FindClusters(seurat.obj, resolution = 0.5)


seurat.obj <- RunUMAP(seurat.obj, dims = 1:15)

DimPlot(seurat.obj, reduction = "umap")


markers <- FindAllMarkers(
  seurat.obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.csv(markers,"cluster_markers.csv")

top10 <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

DoHeatmap(seurat.obj, features = top10$gene)

DotPlot(seurat.obj, features = c("Cd3d","Ms4a1","Lyz2"))

FeaturePlot(seurat.obj,
            features = c("Cd3d","Ms4a1"))

de <- FindMarkers(
  seurat.obj,
  ident.1 = 0,
  ident.2 = 1
)

head(de)

write.csv(de,"DE_results.csv")
