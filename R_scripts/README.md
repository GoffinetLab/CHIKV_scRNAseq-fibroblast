These are the scripts used to preprocess, integrate, and analyze the files created by a 10X 3' GEM experiment with Chikungunya virus infected primary human synovial fibroblasts.

Preprocessing.R includes the standard workflow for CreateSeuratObject, quality control and filtering cells, cell cycle regression, PCA and clustering (RunUMAP, FindNeighbors, FindClusters) and saving the resulting RDS object. All samples were preprocessed similarily.

Integration.R assigns metadata information about MOI, timepoint, and and a combination of both ("infection") and uses the SCT workflow to integrate all samples from each timepoint together (resulting in two RDS objects, one for 6 hpi and one for 24 hpi). It performs the standard workflow mentioned above for PCA and clustering after integration. 

Metadata_assignment.R categorizes the cells depending on the amount of CHIKV in the cell, and groups them into bystander cells (no CHIKV, but from infected cultures), low CHIKV, and high CHIKV expression, next to the mock-infected cells. Additionally, it creates an interferon module score based on the REACTOME database pathway "Interferon Signaling" (R-HSA-913531) to score each cell for its expression of interferon signaling genes.

Visualization_data-extraction.R finally creates UMAPs and ViolinPlots to visualize the results. Additionally, it calculates differentially expressed genes between groups and extract data for further analysis and visualization in other programs such as GraphPad Prism.


The raw and processed data can be downloaded here: TO FOLLOW
The publication for this dataset is currently under review and will be posted as soon as it is published.
