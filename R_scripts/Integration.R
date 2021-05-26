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


#### Load all required RDS objects. Here, I have used replicate 2 and 4 from
#### different runs, so I name them by replicate, timepoint, and MOI

# Adding metadata for infection status to the samples
# Infection status (moi*timepoint)
rep2_mock_6$infection <- 'mock_6'
rep2_moi01_6$infection <- 'moi01_6'
rep2_moi10_6$infection <- 'moi10_6'
rep4_mock_6$infection <- 'mock_6'
rep4_moi0_01_6$infection <- 'moi0_01_6'
rep4_moi0_1_6$infection <- 'moi0_1_6'

rep2_mock_24$infection <- 'mock_24'
rep2_moi01_24$infection <- 'moi01_24'
rep2_moi10_24$infection <- 'moi10_24'
rep4_mock_24$infection <- 'mock_24'
rep4_moi0_01_24$infection <- 'moi0_01_24'
rep4_moi0_1_24$infection <- 'moi0_1_24'


# timepoint
rep2_mock_6$timepoint <- '6hpi'
rep2_moi01_6$timepoint <- '6hpi'
rep2_moi10_6$timepoint <- '6hpi'
rep4_mock_6$timepoint <- '6hpi'
rep4_moi0_01_6$timepoint <- '6hpi'
rep4_moi0_1_6$timepoint <- '6hpi'

rep2_mock_24$timepoint <- '24hpi'
rep2_moi01_24$timepoint <- '24hpi'
rep2_moi10_24$timepoint <- '24hpi'
rep4_mock_24$timepoint <- '24hpi'
rep4_moi0_01_24$timepoint <- '24hpi'
rep4_moi0_1_24$timepoint <- '24hpi'

# MOI
rep2_mock_6$moi <- 'mock'
rep2_moi01_6$moi <- 'moi01'
rep2_moi10_6$moi <- 'moi10'
rep4_mock_6$moi <- 'mock'
rep4_moi0_01_6$moi <- 'moi0_01'
rep4_moi0_1_6$moi <- 'moi0_1'

rep2_mock_24$moi <- 'mock'
rep2_moi01_24$moi <- 'moi01'
rep2_moi10_24$moi <- 'moi10'
rep4_mock_24$moi <- 'mock'
rep4_moi0_01_24$moi <- 'moi0_01'
rep4_moi0_1_24$moi <- 'moi0_1'


# Change DefaultAssay to SCT
DefaultAssay(rep2_mock_24) <- "SCT"
DefaultAssay(rep2_moi01_24) <- "SCT"
DefaultAssay(rep2_moi10_24) <- "SCT"
DefaultAssay(rep2_mock_6) <- "SCT"
DefaultAssay(rep2_moi01_6) <- "SCT"
DefaultAssay(rep2_moi10_6) <- "SCT"
DefaultAssay(rep4_mock_24) <- "SCT"
DefaultAssay(rep4_moi0_1_24) <- "SCT"
DefaultAssay(rep4_moi0_01_24) <- "SCT"
DefaultAssay(rep4_mock_6) <- "SCT"
DefaultAssay(rep4_moi0_1_6) <- "SCT"
DefaultAssay(rep4_moi0_01_6) <- "SCT"

###### Integration ################################################################


# Set RAM capacity to maximum
options(future.globals.maxSize = 1000 * 1024^2)
memory.limit(9999999999)

# SCTransform pipeline integration
rep2_4_6hpi.features <- SelectIntegrationFeatures(object.list = c(rep2_mock_6,
                                                                  rep2_moi01_6,
                                                                  rep2_moi10_6,
                                                                  rep4_mock_6,
                                                                  rep4_moi0_01_6,
                                                                  rep4_moi0_1_6), 
                                                  nfeatures = 3000)

rep2_4_6hpi.list <- PrepSCTIntegration(object.list = c(rep2_mock_6,
                                                       rep2_moi01_6,
                                                       rep2_moi10_6,
                                                       rep4_mock_6,
                                                       rep4_moi0_01_6,
                                                       rep4_moi0_1_6), 
                                       anchor.features = rep2_4_6hpi.features, 
                                       verbose = TRUE)

rep2_4_6hpi.anchors <- FindIntegrationAnchors(object.list = rep2_4_6hpi.list, 
                                              normalization.method = "SCT", 
                                              anchor.features = rep2_4_6hpi.features, 
                                              verbose = TRUE)

integrated_6hpi <- IntegrateData(anchorset = rep2_4_6hpi.anchors, 
                                 normalization.method = "SCT", 
                                 verbose = TRUE)


# Saving the integrated file without clustering
saveRDS(integrated_6hpi, file = "./RDS_integrated/rep2_4_integrated_6hpi.RDS")

rm(list= c("rep2_4_6hpi.anchors", "rep2_4_6hpi.features", "rep2_4_6hpi.list",
           "rep2_mock_6", "rep2_moi01_6", "rep2_moi10_6", "rep4_mock_6", 
           "rep4_moi0_01_6", "rep4_moi0_1_6"))

# SCTransform pipeline integration, 24 hpi
rep2_4_24hpi.features <- SelectIntegrationFeatures(object.list = c(rep2_mock_24,
                                                                   rep2_moi01_24,
                                                                   rep2_moi10_24,
                                                                   rep4_mock_24,
                                                                   rep4_moi0_01_24,
                                                                   rep4_moi0_1_24), 
                                                   nfeatures = 3000)

rep2_4_24hpi.list <- PrepSCTIntegration(object.list = c(rep2_mock_24,
                                                        rep2_moi01_24,
                                                        rep2_moi10_24,
                                                        rep4_mock_24,
                                                        rep4_moi0_01_24,
                                                        rep4_moi0_1_24), 
                                        anchor.features = rep2_4_24hpi.features, 
                                        verbose = TRUE)

rep2_4_24hpi.anchors <- FindIntegrationAnchors(object.list = rep2_4_24hpi.list, 
                                               normalization.method = "SCT", 
                                               anchor.features = rep2_4_24hpi.features, 
                                               verbose = TRUE)

integrated_24hpi <- IntegrateData(anchorset = rep2_4_24hpi.anchors, 
                                  normalization.method = "SCT", 
                                  verbose = TRUE)


# Saving the integrated file without clustering
saveRDS(integrated_24hpi, file = "./RDS_integrated/rep2_4_integrated_24hpi.RDS")


###### Clustering and visualization #########################################################
#### REPEAT FOR 24 hpi ######

# Run the standard workflow for visualization and clustering
integrated_6hpi <- RunPCA(integrated_6hpi, npcs = 50, 
                          verbose = FALSE)

ElbowPlot(integrated_6hpi, 50)

# t-SNE and Clustering
integrated_6hpi_clus <- RunUMAP(integrated_6hpi, 
                                reduction = "pca", dims = 1:20)

integrated_6hpi_clus <- FindNeighbors(integrated_6hpi_clus, 
                                      reduction = "pca", dims = 1:20)

integrated_6hpi_clus <- FindClusters(integrated_6hpi_clus, 
                                     resolution = 0.5) # Resolution for 24 hpi: 0.8

# UMAP colored by infection status
Idents(integrated_6hpi_clus) <- "infection"

integrated_6hpi_umap_infection <-DimPlot(integrated_6hpi_clus,
                                         reduction = "umap", label = FALSE,
                                         pt.size = 1, label.size = 6, 
                                         split.by = "timepoint",
                                         cols = c("#DD8D29", "#46ACC8", "#972D15",
                                                  "#c94a53", "#be5168", "#a34974")
) + NoAxes() + theme(legend.text=element_text(size=16))

integrated_6hpi_umap_infection

# ggsave("./Preprocessing/integrated_all_infection.png", 
#        plot = integrated_all_umap_infection, 
#        width = 90, height = 15, units = "cm")


# UMAP colored by cluster
Idents(integrated_6hpi_clus) <- "integrated_snn_res.0.8"

integrated_6hpi_umap_cluster <-DimPlot(integrated_6hpi_clus,
                                       reduction = "umap", label = FALSE,
                                       pt.size = 1, label.size = 6, 
                                       cols = c("#51574a", "#447c69", "#74c493",
                                                "#8e8c6d", "#e4bf80", "#e9d78e",
                                                "#e2975d", "#f19670", "#e16552",
                                                "#c94a53", "#be5168", "#a34974",
                                                "#993767", "#65387d", "#4e2472",
                                                "#9163b6", "#e279a3", "#e0598b",
                                                "#7c9fb0", "#5698c4", "#9abf88")
) + NoAxes() + theme(legend.text=element_text(size=16))

integrated_6hpi_umap_cluster

# ggsave("./Preprocessing/integrated_all_cluster.png", 
#        plot = integrated_all_umap_cluster, 
#        width = 20, height = 15, units = "cm")

# Changing the Default Assay type to SCT for visualization
DefaultAssay(integrated_6hpi_clus) <- "SCT"

# Saving clustered file
saveRDS(integrated_6hpi_clus, file = "./RDS_integrated/integrated_6hpi_clus.RDS")




