# Load required packages

# Processing
library(Seurat)
library(dplyr)
library(cowplot)
library(patchwork)
library(ggplot2)
library(sctransform)
# Visualization
library(RColorBrewer)
library(rcartocolor)


########## Preprocessing 
# performed similarily for all samples

# "mock_24" is the example here, can be replaced by
# "moi0_01_24"
# "moi0_1_24"
# "moi01_24"
# "moi10_24"
# "moi0_01_6"
# "moi0_1_6"
# "moi01_6"
# "moi10_6"

### Filtering similar for all samples: nFeature_RNA >700 & percent.mt < 20 & 
#   nFeature_RNA < 7000


mock_24 <- Read10X("./Feature_Matrices/raw_feature_bc_matrix")
mock_24 <- CreateSeuratObject(mock_24, 
                              project = "chikv.mock.24", 
                              min.cells = 3, min.features = 200, 
                              names.delim="-")
# Show data
head(mock_24@meta.data,5)
# Save mock_24.rds
saveRDS(mock_24, file = "./RDS_single/mock_24.RDS")

# Determine percent MT DNA reads as well as nFeatures
# Add column to metadata describing percent mt (very high % MT DNA = dead cell)
mock_24[["percent.mt"]] <- PercentageFeatureSet(mock_24, pattern = "^MT-")

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Create single plots for quality control and save them
nfeature_rna_mock <- VlnPlot(mock_24, features = "nFeature_RNA") +
  geom_hline(yintercept=700, linetype="dashed", color = "red") +
  geom_hline(yintercept=7000, linetype="dashed", color = "red")+
  theme(legend.position = "none")

ncount_rna_mock <- VlnPlot(mock_24, features = "nCount_RNA") +
  theme(legend.position = "none")

percent_mt_mock <- VlnPlot(mock_24, features = "percent.mt") +
  geom_hline(yintercept=20, linetype="dashed", color = "red") +
  theme(legend.position = "none")

qc_mock_24 <- plot_grid(nfeature_rna_mock, ncount_rna_mock, percent_mt_mock,
                        ncol = 3)

qc_mock_24

ggsave("./Preprocessing/qc_mock_24.png", plot = qc_mock_24, width = 20, 
       height = 15, units = "cm")

# Filter cells based on criteria
mock_24sel <-subset(mock_24, 
                    subset= nFeature_RNA >700 & 
                      percent.mt < 20 & 
                      nFeature_RNA < 7000)

# Score cell cycle
mock_24sel <- CellCycleScoring(mock_24, s.features = s.genes, 
                               g2m.features = g2m.genes, set.ident = TRUE)

head(mock_24sel[[]])

# Filter out cells that have CHIKV genes falsely assigned
# This is due to "index hopping" on the Illumina 
# Only for Mock samples
mock_24sel[["percent.chikv"]] <- PercentageFeatureSet(mock_24sel, 
                                                      pattern = "^CHIKV-")

mock_24sel <-subset(mock_24sel, 
                    subset= percent.chikv <= 0)

#check the percentage of CHIKV reads (should be 0)
percent_chikv <- VlnPlot(mock_24sel, features = "percent.chikv") +
  theme(legend.position = "none")
percent_chikv

#m SCTtransform data and regress out cell cycle scores
mock_24sel <- SCTransform(mock_24sel, 
                          vars.to.regress = c("S.Score", "G2M.Score"),
                          verbose = FALSE)



# Run PCA and generate 30 pcs
mock_24sel <- RunPCA(mock_24sel, npcs=30, verbose = FALSE)

# Select how many PCs are needed
ElbowPlot(mock_24sel)

# Generate umap and cluster cells
mock_24clus <- RunUMAP(mock_24sel, reduction = "pca", dims = 1:20, 
                       verbose = FALSE)
mock_24clus <- FindNeighbors(mock_24clus, reduction = "pca", dims = 1:20, 
                             verbose = FALSE)
mock_24clus <- FindClusters(mock_24clus, resolution = 0.8, verbose = FALSE)

# Save RDS
DefaultAssay(mock_24clus) <- "SCT"
saveRDS(mock_24clus, file = "./RDS_single/mock_24clus.RDS")

