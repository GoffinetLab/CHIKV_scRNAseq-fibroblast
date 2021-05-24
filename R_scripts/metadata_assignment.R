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


####################### 24 hpi ###########################

# Adding Metadata for CHIKV RNA level

# creating subsets for infection status and low/high CHIKV RNA combined

Idents(integrated_24hpi) <- "moi"

high_chikv_moi10 <- WhichCells(object = integrated_24hpi, 
                               expression = `CHIKV-sp` >= 4.5, idents = "moi10") 
high_chikv_moi01 <- WhichCells(object = integrated_24hpi, 
                               expression = `CHIKV-sp` >= 4, idents = "moi01") 
high_chikv_moi0_1 <- WhichCells(object = integrated_24hpi, 
                                expression = `CHIKV-sp` >= 3.5, idents = "moi0_1") 
high_chikv_moi0_01 <- WhichCells(object = integrated_24hpi, 
                                 expression = `CHIKV-sp` >= 2, idents = "moi0_01") 

high_chikv <- c(high_chikv_moi10, high_chikv_moi01, high_chikv_moi0_1,
                high_chikv_moi0_01)

low_chikv_moi10 <- WhichCells(object = integrated_24hpi, 
                              expression = `CHIKV-sp` < 4.5 & `CHIKV-sp` > 0,
                              idents = "moi10")
low_chikv_moi01 <- WhichCells(object = integrated_24hpi, 
                              expression = `CHIKV-sp` < 4 & `CHIKV-sp` > 0,
                              idents = "moi01")
low_chikv_moi0_1 <- WhichCells(object = integrated_24hpi, 
                               expression = `CHIKV-sp` < 3.5 & `CHIKV-sp` > 0,
                               idents = "moi0_1")
low_chikv_moi0_01 <- WhichCells(object = integrated_24hpi, 
                                expression = `CHIKV-sp` < 2 & `CHIKV-sp` > 0,
                                idents = "moi0_01")

low_chikv <- c(low_chikv_moi10, low_chikv_moi01, low_chikv_moi0_1,
               low_chikv_moi0_01)

no_chikv <- WhichCells(object = integrated_24hpi, 
                       expression = `CHIKV-sp` == 0, 
                       idents = c("moi01", "moi10", "moi0_1", "moi0_01"))

mock_24 <- WhichCells(object = integrated_24hpi, 
                      expression = `CHIKV-sp` == 0, idents = "mock")

all_chikv <- WhichCells(object = integrated_24hpi, 
                        expression = `CHIKV-sp` > 0)




# calculating the number of cells in each subset
length(high_chikv)
length(low_chikv)
length(no_chikv)
length(mock_24)
length(all_chikv)



# creating metadata for each cell with the information which MOI and CHIKV RNA low/high

integrated_24hpi[["chikrnalvl"]] <- ifelse(rownames(integrated_24hpi@meta.data) %in% 
                                             high_chikv, "high_chikv",
                                           ifelse(rownames(integrated_24hpi@meta.data) %in% low_chikv, "low_chikv",
                                                  ifelse(rownames(integrated_24hpi@meta.data) %in% no_chikv, "no_chikv",
                                                         ifelse(rownames(integrated_24hpi@meta.data) %in% all_chikv, "all_chikv",
                                                                ifelse(rownames(integrated_24hpi@meta.data) %in% mock_24, "mock",
                                                                       "undefined")))))


DefaultAssay(integrated_24hpi) <- "SCT"

# Check for CHIKV RNA amount in each subset
Idents(integrated_24hpi) <- "moi"

VlnPlot(integrated_24hpi, features = "CHIKV-sp", split.by = "chikrnalvl",
        cols = c("#f4a261", "#f4a261", "#2a9d8f", "#2a9d8f", "#2a9d8f",
                 "#2a9d8f", "#264653", "#264653", "#264653", "#264653"),
        pt.size = 0) + geom_hline(yintercept = 3.5)



# Add interferon module score
# Score cells depending on their expression of the Reactome IFN signalling genes

DefaultAssay(integrated_24hpi) <- "SCT"
Idents(integrated_24hpi) <- "chikrnalvl"

ifn_sig <- read.csv(file="Reactome_interferon_signalling.csv")
ifn_sig <- as.vector(ifn_sig[,1])
ifn_sig<- list(c(print(ifn_sig)))
ifn_sig

integrated_24hpi <- AddModuleScore(object = integrated_24hpi, 
                                   features = ifn_sig, ctrl = 100, 
                                   name = 'ifn_sig', assay="SCT")

head(integrated_24hpi[[]])

VlnPlot(integrated_24hpi, features = "ifn_sig1", split.by = "chikrnalvl",  
        cols = c("#f4a261", "#f4a261", "#2a9d8f", "#2a9d8f", "#2a9d8f",
                 "#2a9d8f", "#264653", "#264653", "#264653", "#264653"),
        pt.size = 0)

# Save RDS with metadata
saveRDS(integrated_24hpi, file = "integrated_24hpi.RDS")






######################## 6 hpi #######################


# Adding Metadata for CHIKV RNA level

# creating subsets for infection status and low/high CHIKV RNA combined
# Adding Metadata for CHIKV RNA level

# creating subsets for infection status and low/high CHIKV RNA combined

Idents(integrated_6hpi) <- "moi"

high_chikv_moi10 <- WhichCells(object = integrated_6hpi, 
                               expression = `CHIKV-sp` >= 4, idents = "moi10") 
high_chikv_moi01 <- WhichCells(object = integrated_6hpi, 
                               expression = `CHIKV-sp` >= 3, idents = "moi01") 
high_chikv_moi0_1 <- WhichCells(object = integrated_6hpi, 
                                expression = `CHIKV-sp` >= 2.5, idents = "moi0_1") 
high_chikv_moi0_01 <- WhichCells(object = integrated_6hpi, 
                                 expression = `CHIKV-sp` >= 1.5, idents = "moi0_01") 

high_chikv <- c(high_chikv_moi10, high_chikv_moi01, high_chikv_moi0_1,
                high_chikv_moi0_01)

low_chikv_moi10 <- WhichCells(object = integrated_6hpi, 
                              expression = `CHIKV-sp` < 4 & `CHIKV-sp` > 0,
                              idents = "moi10")
low_chikv_moi01 <- WhichCells(object = integrated_6hpi, 
                              expression = `CHIKV-sp` < 3 & `CHIKV-sp` > 0,
                              idents = "moi01")
low_chikv_moi0_1 <- WhichCells(object = integrated_6hpi, 
                               expression = `CHIKV-sp` < 2.5 & `CHIKV-sp` > 0,
                               idents = "moi0_1")
low_chikv_moi0_01 <- WhichCells(object = integrated_6hpi, 
                                expression = `CHIKV-sp` < 1.5 & `CHIKV-sp` > 0,
                                idents = "moi0_01")

low_chikv <- c(low_chikv_moi10, low_chikv_moi01, low_chikv_moi0_1,
               low_chikv_moi0_01)

no_chikv <- WhichCells(object = integrated_6hpi, 
                       expression = `CHIKV-sp` == 0, 
                       idents = c("moi01", "moi10", "moi0_1", "moi0_01"))

mock_6 <- WhichCells(object = integrated_6hpi, 
                     expression = `CHIKV-sp` == 0, idents = "mock")

all_chikv <- WhichCells(object = integrated_6hpi, 
                        expression = `CHIKV-sp` > 0)




# calculating the number of cells in each subset
length(high_chikv)
length(low_chikv)
length(no_chikv)
length(mock_6)
length(all_chikv)



# creating metadata for each cell with the information which MOI and CHIKV RNA low/high

integrated_6hpi[["chikrnalvl"]] <- ifelse(rownames(integrated_6hpi@meta.data) %in% 
                                            high_chikv, "high_chikv",
                                          ifelse(rownames(integrated_6hpi@meta.data) %in% low_chikv, "low_chikv",
                                                 ifelse(rownames(integrated_6hpi@meta.data) %in% no_chikv, "no_chikv",
                                                        ifelse(rownames(integrated_6hpi@meta.data) %in% all_chikv, "all_chikv",
                                                               ifelse(rownames(integrated_6hpi@meta.data) %in% mock_6, "mock",
                                                                      "undefined")))))


DefaultAssay(integrated_6hpi) <- "SCT"

# Check for CHIKV RNA amount in each subset
Idents(integrated_6hpi) <- "moi"

VlnPlot(integrated_6hpi, features = "CHIKV-sp", 
        split.by = "chikrnalvl",
        cols = c("#f4a261", "#f4a261", "#2a9d8f", "#2a9d8f", "#2a9d8f",
                 "#2a9d8f", "#264653", "#264653", "#264653", "#264653"),
        pt.size = 0) + geom_hline(yintercept = 4)


# Add interferon module score
# Score cells depending on their expression of the Reactome IFN signalling genes

DefaultAssay(integrated_6hpi) <- "SCT"
Idents(integrated_6hpi) <- "chikrnalvl"

ifn_sig <- read.csv(file="Reactome_interferon_signalling.csv")
ifn_sig <- as.vector(ifn_sig[,1])
ifn_sig<- list(c(print(ifn_sig)))
ifn_sig

integrated_6hpi <- AddModuleScore(object = integrated_6hpi, 
                                  features = ifn_sig, ctrl = 100, 
                                  name = 'ifn_sig', assay="SCT")

head(integrated_6hpi[[]])

VlnPlot(integrated_6hpi, features = "ifn_sig1", split.by = "chikrnalvl",  
        cols = c("#f4a261", "#f4a261", "#2a9d8f", "#2a9d8f", "#2a9d8f",
                 "#2a9d8f", "#264653", "#264653", "#264653", "#264653"),
        pt.size = 0)

# save RDS with metadata
saveRDS(integrated_6hpi, file = "integrated_6hpi.RDS")
