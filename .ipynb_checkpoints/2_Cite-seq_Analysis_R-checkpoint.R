# Load Seurat file
load("Swarbrick_metastatic_lymphnode.seurat.Rdata")
# Preprocessing: Remove the error cell (last one)
citeCells <- subset(citeCells, cells =colnames(citeCells)[1:3411] )

# Load Autoencoder file
ae_data = read.csv("autoencoder_50dim_test.csv")
row.names(ae_data) = row.names(citeCells@meta.data)
ae_data = as.matrix(ae_data)

# Import AE to Seurat object
citeCells[["ae"]] <- CreateDimReducObject(embeddings = ae_data, key = "AE_", assay = DefaultAssay(citeCells))

# Clustering
citeCells <- FindNeighbors(citeCells,reduction="ae",dims = 1:50)
citeCells <- FindClusters(citeCells, resolution = 0.2)

# Run Umap and Visualize
citeCells = RunUMAP(object = citeCells,reduction.key = "aeumap_",
                    reduction = "ae",dims = 1:50,n.components = 2)

DimPlot(citeCells, reduction = "umap",label = T, pt.size = 1) + NoLegend()

################################ DE Analysis #########################################
library(tidyverse)

# Find top 10 and visualize RNA markers
citeCells.rna.markers <- FindAllMarkers(citeCells, max.cells.per.ident = 100, only.pos = TRUE)
top10 <- citeCells.rna.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(citeCells,features = top10$gene)

# Find top 10 and visualize ADT markers
DefaultAssay(citeCells) <- "ADT"
citeCells.adt.markers <- FindAllMarkers(citeCells, max.cells.per.ident = 100, only.pos = TRUE)
top10_adt <- citeCells.adt.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(citeCells,features = top10_adt$gene)


################################ Annotation #########################################

library(SingleR)
counts <- GetAssayData(citeCells, slot = "data")

edit_genename <- function(x){
  x <- sub("GRCh38-","",x)
  x <- sub("mm10---","",x)
  return(x)
}

row.names(counts) = edit_genename(row.names(counts))
singler = CreateSinglerObject(counts, annot = NULL, "OZ", min.genes = 0,
                              technology = "10X", species = "Human", citation = "",
                              ref.list = list(), normalize.gene.length = F, variable.genes = "de",
                              fine.tune = F, do.signatures = F, clusters = citeCells@active.ident, do.main.types = T, 
                              reduce.file.size = T, numCores = SingleR.numCores)

human <- singler$singler[1][[1]]$SingleR.clusters.main$labels
singler = CreateSinglerObject(counts, annot = NULL, "OZ", min.genes = 0,
                              technology = "10X", species = "Mouse", citation = "",
                              ref.list = list(), normalize.gene.length = F, variable.genes = "de",
                              fine.tune = F, do.signatures = F, clusters = citeCells@active.ident, do.main.types = T, 
                              reduce.file.size = T, numCores = SingleR.numCores)

mouse <- singler$singler[1][[1]]$SingleR.clusters.main$labels

# Look at the heatmap, I manually annotate based on the RNA markers start with mm_ is 
# Mouse cell, GRCh38_ is Human cell
# Then we got these labels:

new.cluster.ids <- c("HG_T_cells_1","MM_Macrophages","HG_B_cell","HG_T_cells_2","MM_Monocytes","HG_Endothelial_cells",
                     "HG_Monocyte", "MM_Epithelial cells", "HG_NK_cell","HG_Pre-B_cell_CD34-")
names(new.cluster.ids) <- levels(citeCells)
citeCells <- RenameIdents(citeCells, new.cluster.ids)

# Visualize it again
DimPlot(citeCells, reduction = "umap",label = T, pt.size = 1) + NoLegend()