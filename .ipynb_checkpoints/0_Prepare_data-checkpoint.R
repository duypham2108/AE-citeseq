# Load Seurat file
load("Swarbrick_metastatic_lymphnode.seurat.Rdata")

# Prepare input data for Model

# Get count matrix of markers
# RNA matrix
rna_matrix = t(as.matrix(citeCells@assays$SCT@data))
temp = as.data.frame(rna_matrix)
temp = temp[,citeCells@assays$SCT@var.features]
rna_matrix = as.matrix(temp)
rna_matrix = cbind(rna_matrix,citeCells@active.ident)

# ADT matrix
citeCells <- NormalizeData(citeCells, assay = "ADT", normalization.method = "genesCLR")
proteomic_matrix = t(as.matrix(citeCells@assays$ADT@data))
proteomic_matrix = cbind(proteomic_matrix,citeCells@active.ident)

write.table(rna_matrix,"scRNAseq.txt",row.names = F, col.names = F,sep="\t")
write.table(proteomic_matrix,"scProteomics.txt",row.names = F, col.names = F,sep="\t")