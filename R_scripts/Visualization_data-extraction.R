
######### All code has been used in a similar way to produce the 6 hpi figures


# Load required packages

# Processing
library(Seurat)
library(dplyr)
library(cowplot)
library(patchwork)
library(ggplot2)
library(tidyverse)
# Visualization
library(RColorBrewer)
library(rcartocolor)
library(svglite)

options(future.globals.maxSize = 10000 * 1024^2)
memory.limit(9999999999)

# load integrated Seurat object (integrated_24hpi or integrated_6hpi)


#Defining the color scheme

cols <- carto_pal(7, "Geyser") #rev(brewer.pal(9, 'GnBU')) reverses order
newcol <- colorRampPalette(cols) #extrapolate inbetween colors
ncols <- 100 #how many new colors?
cols2 <- newcol(ncols) #apply the function to get 100 colours

cols_single <- brewer.pal(9, "YlOrRd")[4:8] #rev(brewer.pal(9, 'GnBU')) reverses order
newcol_single <- colorRampPalette(cols_single) #extrapolate inbetween colors
ncols <- 100 #how many new colors?
cols2_single <- newcol_single(ncols) #apply the function to get 100 colours

######### UMAPS ###########

# UMAP colored by sample
# UMAP colored by MOI
integrated_24hpi$infection <- factor(x = integrated_24hpi$infection,
                               levels = c("mock_24", "moi0_01_24", "moi0_1_24",
                                          "moi01_24", "moi10_24"))

Idents(integrated_24hpi) <- "infection"

integrated_24hpi_umap_moi <-DimPlot(integrated_24hpi,
                                      reduction = "umap", label = FALSE,
                                      pt.size = 0.5, label.size = 6, 
                                      cols = c("#888888","#018181", "#B5C8A9",
                                               "#EDBB89", "#CA562D")
) + NoAxes() + theme(legend.text=element_text(size=16))

integrated_24hpi_umap_moi


# UMAP grid of CHIKV expression
DefaultAssay(integrated_24hpi) <- "SCT"

col.sc.chikv <- scale_color_gradientn(colours = cols2, limits=c(0.0001,10))


VlnPlot(integrated_24hpi, features = "CHIKV-sp")

umap_integrated_24hpi_chikv <- FeaturePlot(integrated_24hpi, features = "CHIKV-sp",
                                           split.by = "moi",
                                           label = FALSE, combine = FALSE, 
                                           pt.size = 1)

for(i in 1:length(umap_integrated_24hpi_chikv)) {
  umap_integrated_24hpi_chikv[[i]] <- umap_integrated_24hpi_chikv[[i]] + 
    NoAxes() + labs(title ="") + 
    theme(legend.position = "none") + 
    theme(panel.border =element_blank()) + 
    theme(plot.background = element_rect(fill = "transparent", color = NA)) +
    col.sc.chikv
} 


grid_splitbychikv_CHIKV_all <- cowplot::plot_grid(plotlist = 
                                                    umap_integrated_24hpi_chikv, 
                                                  ncol = 5)

grid_splitbychikv_CHIKV_all


# UMAPs of IFN module score
col.sc.ifn <- scale_color_gradientn(colours = cols2, limits=c(-0.1,0.65))

umap_integrated_24hpi_ifn <- FeaturePlot(integrated_24hpi, features = "ifn_sig1", 
                                         split.by = "moi",
                                         label = FALSE, combine = FALSE, 
                                         pt.size = 1)

for(i in 1:length(umap_integrated_24hpi_ifn)) {
  umap_integrated_24hpi_ifn[[i]] <- umap_integrated_24hpi_ifn[[i]] + NoAxes() + 
    NoAxes() + labs(title ="") + 
    theme(legend.position = "none") + 
    theme(panel.border =element_blank()) + 
    theme(plot.background = element_rect(fill = "transparent", color = NA)) +
    col.sc.ifn
} 

grid_splitbychikv_ifn_all <- cowplot::plot_grid(plotlist = 
                                                  umap_integrated_24hpi_ifn, 
                                                ncol = 5)

grid_splitbychikv_ifn_all



# UMAPs of Cofactors
cols_single <- brewer.pal(9, "YlOrRd")[4:8] 
newcol_single <- colorRampPalette(cols_single) #extrapolate inbetween colors
ncols <- 100 #how many new colors?
cols2_single <- newcol_single(ncols) #apply the function to get 100 colours
col.sc.cofac <- scale_color_gradientn(colours = cols2_single, 
                                      limits=c(0.0001, 6.2))


umap_integrated_24hpi_cofac <- FeaturePlot(integrated_24hpi, 
                                           features = c("MXRA8", "FHL1", "FCGR3A",
                                                        "FURIN", "VIM", "COL3A1"),
                                           split.by = "moi",
                                           label = FALSE, combine = FALSE, 
                                           pt.size = 0.5)

for(i in 1:length(umap_integrated_24hpi_cofac)) {
  umap_integrated_24hpi_cofac[[i]] <- umap_integrated_24hpi_cofac[[i]] + 
    NoAxes() + labs(title ="") + 
    theme(legend.position = "none") + 
    theme(panel.border =element_blank()) + 
    theme(plot.background = element_rect(fill = "transparent", color = NA)) +
    col.sc.cofac
} 



grid_integrated_24hpi_cofac <- cowplot::plot_grid(plotlist = umap_integrated_24hpi_cofac,
                                                  ncol = 6)

grid_integrated_24hpi_cofac

# UMAPs of ISGs
cols_single <- brewer.pal(9, "YlOrRd")[4:8] #rev(brewer.pal(9, 'GnBU')) reverses order
newcol_single <- colorRampPalette(cols_single) #extrapolate inbetween colors
ncols <- 100 #how many new colors?
cols2_single <- newcol_single(ncols) #apply the function to get 100 colours
col.sc.isg <- scale_color_gradientn(colours = cols2_single, 
                                      limits=c(0.0001, 6))


#VlnPlot(integrated_6hpi, features = "ISG15")
                                      

umap_integrated_6hpi_isg <- FeaturePlot(integrated_6hpi, 
                                         features = c("IFIT1", "MX1", "IFITM3",
                                                      "ISG15", "STAT1", "IFNB1"),
                                         split.by = "moi",
                                         label = FALSE, combine = FALSE, 
                                         pt.size = 0.5)

for(i in 1:length(umap_integrated_6hpi_isg)) {
  umap_integrated_6hpi_isg[[i]] <- umap_integrated_6hpi_isg[[i]] + 
    NoAxes() + labs(title ="") + 
    theme(legend.position = "none") + 
    theme(panel.border =element_blank()) + 
    theme(plot.background = element_rect(fill = "transparent", color = NA)) +
    col.sc.isg
} 



grid_integrated_6hpi_isg <- cowplot::plot_grid(plotlist = umap_integrated_6hpi_isg,
                                                  ncol = 6)

#grid_integrated_6hpi_cofac

ggsave("./UMAPs/integrated_6hpi_isg.png",
       plot = grid_integrated_6hpi_isg,
       width = 30, height = 30, units = "cm", dpi = 300)



########### VIOLIN PLOTS ##########

# reorder
integrated_24hpi$moi <- factor(x = integrated_24hpi$moi,
                               levels = c('mock', 'moi0_01', "moi0_1", 
                                          "moi01", "moi10"))
integrated_24hpi$chikrnalvl <- factor(x = integrated_24hpi$chikrnalvl, 
                                      levels = c('mock', "no_chikv", 
                                                 "low_chikv", "high_chikv"))

#VlnPlot for CHIKV expression level with cutoff indication
Idents(integrated_24hpi) <- "moi"

VlnPlot_chikv_24 <- VlnPlot(integrated_24hpi, features = "CHIKV-sp", 
        #split.by = "chikrnalvl",
        cols = c("#000000", "#000000", "#000000", "#000000", "#000000" ),
        pt.size = 0) + geom_hline(yintercept = c(2, 3.5, 4, 4.5)) #6 hpi: 1.5, 2.5, 3, 4

VlnPlot_chikv_24


#VlnPlot for IFN module score by MOI
Idents(integrated_24hpi) <- "moi"

Vln_ifnsig_24 <- VlnPlot(integrated_24hpi, features = "ifn_sig1", 
                         split.by = "chikrnalvl",
                         cols = c("#C8CE9C",
                                  "#5DBE82", "#298882", "#113C43",
                                  "#5DBE82",  "#298882", "#113C43",
                                  "#2A9D80", "#113C43"), 
                         pt.size = 0) + 
  theme(legend.position = "right") +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5,
               position = position_dodge(width = 0.9))

Vln_ifnsig_24




########### EXTRACT DATA FOR ANALYSIS ##########

# Fetch data for interferon module score analysis
chikv_ifnsig <- FetchData(integrated_24hpi, 
                          vars= c("timepoint", "chikrnalvl", "moi", "ifn_sig1"))

write.csv2(chikv_ifnsig, file = "chikv_ifnsig_all.csv")


## Extraction of data for volcano plots

# DEG testing
# no_chikv vs c(low_chikv, high_chikv)
Idents(integrated_24hpi) <- "chikrnalvl"
no_all_chikv <- FindMarkers(integrated_24hpi,
                            ident.1 = c("low_chikv", "high_chikv"),
                            ident.2 = "no_chikv")

head(no_all_chikv)
write.csv2(no_all_chikv, file = "Data output/no_vs_all_chikv24.csv")


# no_chikv vs low_chikv
no_low_chikv <- FindMarkers(integrated_24hpi,
                            ident.1 = "low_chikv",
                            ident.2 = "no_chikv")

head(no_low_chikv)
write.csv2(no_low_chikv, file = "Data output/no_vs_low_chikv24.csv")

# no_chikv vs high_chikv
no_high_chikv <- FindMarkers(integrated_24hpi,
                            ident.1 = "high_chikv",
                            ident.2 = "no_chikv")

head(no_high_chikv)
write.csv2(no_high_chikv, file = "Data output/no_vs_high_chikv24.csv")

# low_chikv vs high_chikv
low_high_chikv <- FindMarkers(integrated_24hpi,
                            ident.1 = "high_chikv",
                            ident.2 = "low_chikv")

head(low_high_chikv)
write.csv2(low_high_chikv, file = "Data output/low_vs_high_chikv24.csv")



## Extraction of data for heatmaps
# Fetch average expression data for heatmaps

# ISGs
genelist_bins_24 <- FetchData(integrated_24hpi, 
                          vars= c("chikv_bins", "CHIKV-sp", "ISG15", "IFI6",
                                  "MX1", "MX2", "IFIT1", "IFIT2", "IFIT3",
                                  "IFIT5", "OASL", "IFITM1", "IFITM2", "IFITM3",
                                  "CGAS","CAD", "IRF1", "IRF2", "IRF3",
                                  "IRF4", "IRF7", "IRF9", "STAT1", "STAT2", 
                                  "STAT6", "IFNB1", "IFNL1", "IFNL2", "IFNL3", 
                                  "TNF", "IL6", "IL1B"))

genelist_avg_bins24 <- genelist_bins_24 %>% 
  group_by(chikv_bins) %>%
  summarise(CHIKV = mean(`CHIKV-sp`), ISG15 = mean(ISG15),
            IFI6 = mean(IFI6), MX1 = mean(MX1),
            MX2 = mean(MX2), IFIT1 = mean(IFIT1),
            IFIT2 = mean(IFIT2), IFIT3 = mean(IFIT3),
            IFIT5 = mean(IFIT5), OASL = mean(OASL),
            IFITM1 = mean(IFITM1), IFITM2 = mean(IFITM2),
            IFITM3 = mean(IFITM3), cGAS = mean(CGAS),
            IRG1 = mean(CAD), IRF1 = mean(IRF1),
            IRF2 = mean(IRF2), IRF3 = mean(IRF3),
            IRF4 = mean(IRF4), IRF7 = mean(IRF7),
            IRF9 = mean(IRF9), STAT1 = mean(STAT1),
            STAT2 = mean(STAT2), STAT6 = mean(STAT6),
            IFNB1 = mean(IFNB1), IFNL1 = mean(IFNL1),
            IFNL2 = mean(IFNL2), IFNL3 = mean(IFNL3), 
            TNF = mean(TNF), IL6 = mean(IL6), 
            IL1B = mean(IL1B))


write.csv2(genelist_avg_bins24, file = "Data output/genelist_ISG_perbin24.csv")

# Marker genes
genelist_bins24_marker <- FetchData(integrated_24hpi, 
                                 vars= c("chikv_bins", "VIM", "COL3A1", 
                                         "ADGRE1", "CD14", "FCGR3A", "GAPDH", 
                                         "MAPK1", "MXRA8", "FHL1", "PHB", 
                                         "AXL", "DHX9", "SAMHD1", "FURIN"))

genelist_avg_marker_bins24 <- genelist_bins24_marker %>% 
  group_by(chikv_bins) %>%
  summarise(VIM = mean(VIM), COL3A1 = mean(COL3A1),
            EMR1 = mean(ADGRE1), CD14 = mean(CD14), 
            GAPDH = mean(GAPDH), 
            MAPK1 = mean(MAPK1), MXRA8 = mean(MXRA8), 
            FHL1 = mean(FHL1), PHB = mean(PHB),
            AXL = mean(AXL), DHX9 = mean(DHX9),
            SAMHD1 = mean(SAMHD1), 
            FURIN = mean(FURIN))


write.csv2(genelist_avg_marker_bins24, file = "genelist_marker_bins24.csv")

# RA genes
genelist_bins24_ra <- FetchData(integrated_24hpi, 
                                 vars= c("chikv_bins", "CXCL5", "CXCL8", 
                                         "ANPEP", "CCL5", "MMP3", "MMP9", 
                                         "MMP14", "ADAMTS5", "FGF2", "PDPN", 
                                         "NGF", "FAP"))

genelist_avg_ra_bins24 <- genelist_bins24_ra %>% 
  group_by(chikv_bins) %>%
  summarise(CXCL5 = mean(CXCL5), IL8 = mean(CXCL8),
            CD13 = mean(ANPEP), RANTES = mean(CCL5), 
            MMP3 = mean(MMP3), 
            MMP9 = mean(MMP9), MMP14 = mean(MMP14), 
            ADAMTS5 = mean(ADAMTS5), FGF2 = mean(FGF2),
            PDPN = mean(PDPN), NGF = mean(NGF),
            FAP = mean(FAP))


write.csv2(genelist_avg_ra_bins24, file = "genelist_RA_bins24.csv")


# DEG testing for bins

# bystander vs bin_low
Idents(integrated_24hpi) <- "chikv_bins"
bystander_binlow <- FindMarkers(integrated_24hpi,
                            ident.1 = "bin_low", 
                            ident.2 = "bystander")

write.csv2(bystander_binlow, file = "Data output/bystander_binlow24.csv")

# bystander vs bin_high
Idents(integrated_24hpi) <- "chikv_bins"
bystander_binhigh <- FindMarkers(integrated_24hpi,
                            ident.1 = "bin_high", 
                            ident.2 = "bystander")

write.csv2(bystander_binhigh, file = "Data output/bystander_binhigh24.csv")

# binlow vs bin_high
Idents(integrated_24hpi) <- "chikv_bins"
binlow_binhigh <- FindMarkers(integrated_24hpi,
                            ident.1 = "bin_high", 
                            ident.2 = "bin_low")

write.csv2(binlow_binhigh, file = "Data output/binlow_binhigh24.csv")


# Fetch expression for all genes from Interferon module score for 
# correlation analysis
genelist_all_ifnsig24 <- FetchData(integrated_24hpi, 
                                   vars= c("chikrnalvl", "CHIKV-sp", "AAAS","ABCE1",
                                           "ADAR","ARIH1","B2M",
                                           "BST2","CAMK2A","CAMK2B","CAMK2D",
                                           "CAMK2G","CD44","CIITA","DDX58","EGR1",
                                           "EIF2AK2","EIF4A1","EIF4A2","EIF4A3",
                                           "EIF4E","EIF4E2","EIF4E3","EIF4G1",
                                           "EIF4G2","EIF4G3","FCGR1A","FCGR1B",
                                           "FLNA","FLNB","GBP1","GBP2","GBP3",
                                           "GBP4","GBP5","GBP6","GBP7","HERC5",
                                           "HLA-A","HLA-B","HLA-C","HLA-DPA1",
                                           "HLA-DPB1","HLA-DQA1","HLA-DQA2",
                                           "HLA-DQB1","HLA-DQB2","HLA-DRA",
                                           "HLA-DRB1","HLA-DRB3","HLA-DRB4",
                                           "HLA-DRB5","HLA-E","HLA-F","HLA-G",
                                           "HLA-H","ICAM1","IFI27","IFI30","IFI35",
                                           "IFI6","IFIT1","IFIT2","IFIT3","IFITM1",
                                           "IFITM2","IFITM3","IFNA1","IFNA10",
                                           "IFNA13","IFNA14","IFNA16","IFNA17",
                                           "IFNA2","IFNA21","IFNA4","IFNA5",
                                           "IFNA6","IFNA7","IFNA8","IFNAR1",
                                           "IFNAR2","IFNB1","IFNG","IFNGR1",
                                           "IFNGR2","IP6K2","IRF1","IRF2","IRF3",
                                           "IRF4","IRF5","IRF6","IRF7","IRF8",
                                           "IRF9","ISG15","ISG20","JAK1","JAK2",
                                           "KPNA1","KPNA2","KPNA3","KPNA4","KPNA5",
                                           "KPNA7","KPNB1","MAPK3","MID1","MT2A",
                                           "MX1","MX2","NCAM1","NDC1","NEDD4",
                                           "NUP107","NUP133","NUP153","NUP155",
                                           "NUP160","NUP188","NUP205","NUP210",
                                           "NUP214","NUP35","NUP37","NUP42",
                                           "NUP43","NUP50","NUP54","NUP58","NUP62",
                                           "NUP85","NUP88","NUP93","NUP98","OAS1",
                                           "OAS2","OAS3","OASL","PDE12","PIAS1",
                                           "PIN1","PLCG1","PML","POM121","POM121C",
                                           "PPM1B","PRKCD","PSMB8","PTAFR","PTPN1",
                                           "PTPN11","PTPN2","PTPN6","RAE1","RANBP2",
                                           "RNASEL","RPS27A","RSAD2","SAMHD1",
                                           "SEC13","SEH1L","SOCS1","SOCS3","SP100",
                                           "STAT1","STAT2","SUMO1","TPR","TRIM10",
                                           "TRIM14","TRIM17","TRIM2","TRIM21",
                                           "TRIM22","TRIM25","TRIM26","TRIM29","TRIM3",
                                           "TRIM31","TRIM34","TRIM35","TRIM38","TRIM45",
                                           "TRIM46","TRIM48","TRIM5","TRIM6","TRIM62",
                                           "TRIM68","TRIM8","TYK2","UBA52","UBA7",
                                           "UBB","UBC","UBE2E1","UBE2L6","UBE2N",
                                           "USP18","USP41","VCAM1","XAF1"))

write.csv2(genelist_all_ifnsig24, file = "Data output/genelist_all_ifnsig24.csv")



# Calculate amount of CHIKV genes in each cluster
integrated_24hpi[["percent.chikv"]] <- PercentageFeatureSet(integrated_24hpi, 
                                                            pattern = "^CHIKV-")

df <- FetchData(integrated_24hpi, vars = c("percent.chikv", "moi"))

write.csv2(df, file = "Data output/percentchikv_percluster_24hpi.csv")
