library(Matrix)
library(Seurat)
library(qlcMatrix)
library(pheatmap)
library(ggplot2)
library(gridExtra)
library(stringr)
library(dplyr)
library(tidyr)
library(reshape2)
library(celda)
library(TENxPBMCData)
library(scater)
library(scCustomize)


# Set markers
markers <- list(Tcell_Markers = c("CD3E", "CD3D"),
                Bcell_Markers = c("CD79A", "CD79B", "MS4A1"),
                Monocyte_Markers = c("S100A8", "S100A9", "LYZ"),
                NKcell_Markers = "GNLY")

# Download the pbmc 4k dataset from the network
tmpDir = tempdir(check = TRUE)
download.file("https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz", 
              destfile = file.path(tmpDir, "tod.tar.gz"))
download.file("https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_filtered_gene_bc_matrices.tar.gz", 
              destfile = file.path(tmpDir, "toc.tar.gz"))
untar(file.path(tmpDir, "tod.tar.gz"), exdir = tmpDir)
untar(file.path(tmpDir, "toc.tar.gz"), exdir = tmpDir)

cellMatrixLocation = paste0(tmpDir,"/filtered_gene_bc_matrices/GRCh38")
fullMatrixLocation = paste0(tmpDir,"/raw_gene_bc_matrices/GRCh38")

# "cellMatrix" represents the matrix that has been filtered and only retains "cell" elements,
# while "fullMatrix" also includes empty droplets.
cellMatrix     = Read10X(cellMatrixLocation)
fullMatrix     = Read10X(fullMatrixLocation)

# preprocess
raw_seurat <- CreateSeuratObject(counts = cellMatrix, project = "PBMC4k_raw")
raw_seurat <- NormalizeData(raw_seurat)
raw_seurat <- FindVariableFeatures(raw_seurat)
raw_seurat <- ScaleData(raw_seurat)
raw_seurat <- RunPCA(raw_seurat, features = VariableFeatures(raw_seurat))
raw_seurat <- FindNeighbors(raw_seurat)
raw_seurat <- FindClusters(raw_seurat)
raw_seurat <- RunUMAP(raw_seurat, dims = 1:20)

# umap
DimPlot(
  object    = raw_seurat,
  reduction = "umap",
  label     = TRUE,
  label.size= 5,
  pt.size   = 0.5
) + ggtitle("UMAP: Cluster labels")

# Draw marker UMAP
features <- c("CD3D", "CD3E", "GNLY",
              "LYZ", "S100A8", "S100A9",
              "CD79A", "CD79B", "MS4A1")

FeaturePlot(
  object = raw_seurat,
  features = c("CD3D", "CD3E", "GNLY",
               "LYZ", "S100A8", "S100A9",
               "CD79A", "CD79B", "MS4A1"),
  reduction = "umap",
  cols = c("lightgrey", "red"),
  pt.size = 0.5,
  ncol = 3
)

# Annotate and draw the proportion of marker cells
groupClusters <- list(Tcells = c(0,1,2,4), Bcells = c(5,7), Monocytes = c(3,6,9,10), NKcells = 8)
p <- plotMarkerPercentageSeurat(raw_seurat,
                                markers = markers,
                                groupClusters = groupClusters,
                                assayName = c("RNA"))
print(p)

# Draw a violin plot for the expression of the marker gene
# p <- plotMarkerViolin(raw_seurat, markers, groupClusters)
# print(p)

#### FastCAR ####

# Used for estimating parameter selection
# ambProfile = describe.ambient.RNA.sequence(fullMatrix = fullMatrix, 
#                                            start = 10, 
#                                            stop = 500, 
#                                            by = 10, 
#                                            contaminationChanceCutoff = 0.05)
# 
# plot.ambient.profile(ambProfile)
# 
# 
# correctionEffectProfile = describe.correction.effect(fullMatrix, cellMatrix, 50, 250, 25, 0.05)
# 
# plot.correction.effect.chance(correctionEffectProfile)
# 
# 
# plot.correction.effect.removal(correctionEffectProfile)
# 
# emptyDropletCutoff = recommend.empty.cutoff(ambProfile)

# Remove Ambient RNA based on parameters
emptyDropletCutoff        = 250 
contaminationChanceCutoff = 0.0005

ambientProfile = determine.background.to.remove(fullMatrix, emptyDropletCutoff, contaminationChanceCutoff)
FastCARMatrix     = remove.background(cellMatrix, ambientProfile)

# preprocess
FastCAR_seurat = CreateSeuratObject(FastCARMatrix)
FastCAR_seurat <- NormalizeData(FastCAR_seurat)
FastCAR_seurat <- FindVariableFeatures(FastCAR_seurat)
FastCAR_seurat <- ScaleData(FastCAR_seurat)
FastCAR_seurat <- RunPCA(FastCAR_seurat, features = VariableFeatures(FastCAR_seurat))
FastCAR_seurat <- FindNeighbors(FastCAR_seurat)
FastCAR_seurat <- FindClusters(FastCAR_seurat)
FastCAR_seurat <- RunUMAP(FastCAR_seurat, dims = 1:20)

# umap
DimPlot(
  object    = FastCAR_seurat,
  reduction = "umap",
  label     = TRUE,
  label.size= 5,
  pt.size   = 0.5
) + ggtitle("UMAP: Cluster labels")

features <- c("CD3D", "CD3E", "GNLY",
              "LYZ", "S100A8", "S100A9",
              "CD79A", "CD79B", "MS4A1")

FeaturePlot(
  object = FastCAR_seurat,
  features = c("CD3D", "CD3E", "GNLY",
               "LYZ", "S100A8", "S100A9",
               "CD79A", "CD79B", "MS4A1"),
  reduction = "umap",
  cols = c("lightgrey", "red"),
  pt.size = 0.5,
  ncol = 3
)

### Annotations under different parameters
# # 50 cutoff
# groupClusters <- list(Tcells = c(1,2,3,4), Bcells = c(5,6), Monocytes = c(0,8,9,11), NKcells = 7)


# # 150 cutoff
# groupClusters <- list(Tcells = c(1,2,3,4), Bcells = c(5,6), Monocytes = c(0,8,9,11,12), NKcells = 7)

# 250 cutoff
groupClusters <- list(Tcells = c(0,1,3,5,12), Bcells = c(4,7), Monocytes = c(2,6,9,10), NKcells = 8)

p <- plotMarkerPercentageSeurat(FastCAR_seurat,
                                markers = markers,
                                groupClusters = groupClusters,
                                assayName = c("RNA","decontX"))
print(p)


# p <- plotMarkerViolin(FastCAR_seurat, markers, groupClusters)
# print(p)
#### FastCAR ####

#### DecontX ####
# Convert raw data into SingleCellExperiment objects
raw_counts <- GetAssayData(raw_seurat, assay = "RNA", slot = "counts")
sce <- SingleCellExperiment(assays = list(counts = raw_counts))

# Use DecontX to remove ambient RNA
sce <- decontX(sce)

# umap
umap <- reducedDim(sce, "decontX_UMAP")
plotDimReduceCluster(x = sce$decontX_clusters,
                     dim1 = umap[, 1], dim2 = umap[, 2])

plotDecontXContamination(sce)

# normalization
sce <- logNormCounts(sce)
plotDimReduceFeature(as.matrix(logcounts(sce)),
                     dim1 = umap[, 1],
                     dim2 = umap[, 2],
                     features = c("CD3D", "CD3E", "GNLY",
                                  "LYZ", "S100A8", "S100A9",
                                  "CD79A", "CD79B", "MS4A1"),
                     exactMatch = TRUE)

# Annotate and draw the proportion of marker cells
cellTypeMappings <- list(Tcells = 2, Bcells = 5, Monocytes = 1, NKcells = 6)
plotDecontXMarkerPercentage(sce,
                            markers = markers,
                            groupClusters = cellTypeMappings,
                            assayName = "counts")

plotDecontXMarkerPercentage(sce,
                            markers = markers,
                            groupClusters = cellTypeMappings,
                            assayName = "decontXcounts")

plotDecontXMarkerPercentage(sce,
                            markers = markers,
                            groupClusters = cellTypeMappings,
                            assayName = c("counts", "decontXcounts"))

plotDecontXMarkerExpression(sce,
                            markers = markers[["Monocyte_Markers"]],
                            groupClusters = cellTypeMappings,
                            ncol = 3)

sce <- logNormCounts(sce,
                     exprs_values = "decontXcounts",
                     name = "decontXlogcounts")
plotDecontXMarkerExpression(sce,
                            markers = markers[["Monocyte_Markers"]],
                            groupClusters = cellTypeMappings,
                            ncol = 3,
                            assayName = c("logcounts", "decontXlogcounts"))
#### DecontX ####

#### cellbender ####

# Convert the matrix into a type that Cellbender can recognize
# write10xCounts(
#   path = "pbmc4k.h5",
#   x = fullMatrix,
#   version = "3",
#   type = "HDF5"
# )
# Or you can download it in the terminal by using the following command
#curl -O https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices_h5.h5

# Run the command in the Linux terminal
# cellbender remove-background \
# --input pbmc4k.h5 \
# --output tiny_output.h5 \
# --expected-cells 4340 \
# --total-droplets-included 25000 \
# --cuda

# Subsequently, the PyTables package was used in the terminal to further modify the h5 file.
# ptrepack --complevel 5 tiny_output_filtered.h5:/matrix tiny_output_filtered_seurat.h5:/matrix

# Read the output .h5 file
cellbender_data <- Read10X_h5(filename = "tiny_output_filtered_seurat.h5", use.names = TRUE)
cellbenderseurat <- CreateSeuratObject(counts = cellbender_data)
cellbenderseurat <- NormalizeData(cellbenderseurat)
cellbenderseurat <- FindVariableFeatures(cellbenderseurat)
cellbenderseurat <- ScaleData(cellbenderseurat)
cellbenderseurat <- RunPCA(cellbenderseurat, features = VariableFeatures(object = cellbenderseurat))
cellbenderseurat <- FindNeighbors(cellbenderseurat)
cellbenderseurat <- FindClusters(cellbenderseurat)
cellbenderseurat <- RunUMAP(cellbenderseurat, dims = 1:20)

DimPlot(
  object    = cellbenderseurat,
  reduction = "umap",
  label     = TRUE,
  label.size= 5,
  pt.size   = 0.5
) + ggtitle("UMAP: Cluster labels")

features <- c("CD3D", "CD3E", "GNLY",
              "LYZ", "S100A8", "S100A9",
              "CD79A", "CD79B", "MS4A1")

FeaturePlot(
  object = cellbenderseurat,
  features = c("CD3D", "CD3E", "GNLY",
               "LYZ", "S100A8", "S100A9",
               "CD79A", "CD79B", "MS4A1"),
  reduction = "umap",
  cols = c("lightgrey", "red"),
  pt.size = 0.5,
  ncol = 3
)

groupClusters <- list(Tcells = c(0,2,3,4), Bcells = c(5,6,13), Monocytes = c(1,7,9,10,17), NKcells = 8)
p <- plotMarkerPercentageSeurat(cellbenderseurat,
                                markers = markers,
                                groupClusters = groupClusters,
                                assayName = c("RNA"))
# print(p)
# 
# p <- plotMarkerViolin(cellbenderseurat, markers, groupClusters)
# print(p)


#### Draw a summary chart based on the parameters shown in the above chart ####
# 1. Define four 4×4 matrices, with the row and column names all being "Cell Types"
cell_types <- c("Bcells","Monocytes","NKcells","Tcells")

raw_mat <- matrix(
  c(98,28,25,25,
    89,100,76,88,
    16,18,98,22,
    21,27,46,98),
  nrow = 4, byrow = TRUE,
  dimnames = list(cell_types, cell_types)
)
fastcar_mat <- matrix(
  c(100,28,27,25,
    0,96,0,0,
    0,0,96,7,
    12,14,42,93),
  nrow = 4, byrow = TRUE,
  dimnames = list(cell_types, cell_types)
)
cellbender_mat <- matrix(
  c(100,21,20,19,
    32,100,45,57,
    10,7,98,15,
    9,12,45,96),
  nrow = 4, byrow = TRUE,
  dimnames = list(cell_types, cell_types)
)
decontx_mat <- matrix(
  c(100,1,3,0,
    0,100,0,0,
    0,0,98,9,
    0,0,17,89),
  nrow = 4, byrow = TRUE,
  dimnames = list(cell_types, cell_types)
)

# 2. Melt each matrix into a growth table and add a "Method" column
df_raw     <- melt(raw_mat,     varnames = c("RowType","ColType"), value.name = "Percent") %>% mutate(Method = "Raw")
df_fastcar <- melt(fastcar_mat, varnames = c("RowType","ColType"), value.name = "Percent") %>% mutate(Method = "FastCAR")
df_cellb   <- melt(cellbender_mat, varnames = c("RowType","ColType"), value.name = "Percent") %>% mutate(Method = "Cellbender")
df_decontx <- melt(decontx_mat, varnames = c("RowType","ColType"), value.name = "Percent") %>% mutate(Method = "DecontX")

# 3. Merge the data of all methods
df_all <- bind_rows(df_raw, df_fastcar, df_cellb, df_decontx) %>%
  mutate(
    Method  = factor(Method,  levels = c("Raw","FastCAR","Cellbender","DecontX")),
    RowType = factor(RowType, levels = cell_types),
    ColType = factor(ColType, levels = cell_types)
  )

# 4. Draw a 4×4 grid, with each cell being a small bar chart (Raw, FastCAR, Cellbender, DecontX)
ggplot(df_all, aes(x = Method, y = Percent, fill = Method)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%d%%", Percent)), 
            vjust = -0.3, size = 2.5) +
  facet_grid(RowType ~ ColType) +
  scale_fill_manual(values = c(
    Raw       = "#1b9e77",
    FastCAR   = "#d95f02",
    Cellbender= "#7570b3",
    DecontX   = "#e7298a"
  )) +
  labs(x = NULL, y = "Percent", title = "Method Comparison Across Cell‐Type Pairs") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x       = element_text(angle = 45, hjust = 1, size = 6),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    strip.background  = element_rect(fill = "grey90", color = NA),
    strip.text        = element_text(size = 8),
    legend.position   = "none"
  )

