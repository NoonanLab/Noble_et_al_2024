#Use WGCNA to call co-expression modules 

#Load Packages

#Seurat requires old version of Matrix
#remotes::install_version("Matrix", version = "1.6-1")
#library(Matrix)
#packageVersion("Matrix")

#MUST use Seurat V4; Seurat V5 doesn't work 
#remotes::install_version("Seurat", "4.3.0")
#remotes::install_version("SeuratObject", "4.1.3")

library(Seurat)
packageVersion("Seurat")

require(dplyr)
require(tidyverse)
require(ggplot2)
require(tidyr)

require(patchwork)
require(gridExtra)
require(cowplot)
require(scales)
require(reshape2)
require(readr)

library(impute)
library(WGCNA)

#Install hdWGCNA
#devtools::install_github('smorabit/hdWGCNA', ref='dev')
library(hdWGCNA)





#Set
theme_set(theme_cowplot())
set.seed(913)

#Load data
cts <- read_rds(file = 'Merged_Singles.rds')

#Normalize
cts <- NormalizeData(cts) %>% FindVariableFeatures() %>% ScaleData()

#Add PCA to object
cts <- RunPCA(cts)
cts <- FindNeighbors(cts, dims = 1:30)
cts <- FindClusters(cts, resolution = 0.6)
cts <- RunUMAP(cts, dims = 1:30)

#Check UMAP 
DimPlot(cts, label=TRUE, reduction = 'umap') +
  umap_theme()  + NoLegend()


#Step 1: Set Up Seurat Object for WGCNA 

#Multiple hdWGCNA experiments can be run on the same object and are stored in the Misc slot 
#This is a downstream step and cannot be subsetted
#We can use a number metrics to pick the gene list
#variable: use genes which are in VariableFeatures
#fraction: use gene expressed in a certain fraction of cells (can be used with group.by)
#custom: input a custom gene list
#In our case, we am using fraction -- we will use a higher percentage since we have high-depth data. 10% 

cts <- SetupForWGCNA(
  cts,
  gene_select = 'fraction',
  fraction = 0.1, 
  wgcna_name = 'ModCall1'
)

#Step 2: Construct Metacells 
#Uses the k-Nearest Neighbors algorithm to group together similar cells. 
#This reduces sparcity compared to the original matrix; WGCNA methods are sensitive to sparsity 

#We want to only construct metacells from the same biological origin (i.e. replicate) so we must group.by this. 
#We also want to construct metacells for each cell type separately. 

#The number of cells to be aggregated 'k' should be based on the size of the data. 
#In their example, they have 30k cells with 8188 in each sample and use k = 25. 
#These samples are much larger, we will use k = 50. 


#Construct metacells
cts <- MetacellsByGroups(
  seurat_obj = cts,
  group.by = c('seurat_clusters', 'rep'),
  k = 50, 
  max_shared = 10, 
  ident.group = 'seurat_clusters' #Sets the Idents of the Metacell Seurat Object
)


#Normalize the metacell expression matrix
cts <- NormalizeMetacells(cts)

#STEP 3: Co-Expression Network Analysis 

#In their example, they restrict their analysis to one group of cells. 
#In our case, the cells are all quite similar, they are all early RG-like cells whose main difference is cell cycle. 
#We would not have an objective criteria to select one cluster over another, so we will select them all 

cts <- SetDatExpr(
  cts,
  group_name = c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13'), 
  group.by = 'seurat_clusters', 
  assay = 'RNA', 
  slot = 'data'
)

#Select soft-power threshold 

#This is an important step. 
#Gene-gene correlation adjacency matrix to infer co-expression. 
#The correlations are raised to a power to reduce noise,
#this reduces weak connections and retains storng connections. 
#So we must determine the correct power. 

#We perform a parameter sweep. 
#The co-expression network should have a scale-free topology
#Thus we test for this. 

cts <- TestSoftPowers(cts, networkType = 'signed')
plot_list <- PlotSoftPowers(cts)

#assemble plots with patchwork 
wrap_plots(plot_list, ncol = 2)

#Return soft power of 20 

cts <- ConstructNetwork( 
  cts, soft_power = 20, 
  setDatExpr = F, 
  tom_name = 'all')

#STEP 4: Module Eigengenes and Connectivity 

#Compute harmonized module eigengenes 
#MEs are a metric to summarize the gene expression profile of an entire module. 
#MEs are computed by PCA on the subset of the gene expression matrix comprising each module. 
#The first PC of these are the MEs. 

#We can apply Harmony to batch correct MEs to yield harmonized module eigengenes (hMEs)

#Run ScaleData first or else harmony throws an error 
cts <- ScaleData(cts, features=VariableFeatures(cts))

cts <- ModuleEigengenes(
  cts,
  group.by.vars = 'rep'
)

#The ME matrices are stored as a matrix: Rows = Cells, Columns = Modules 
#GetMEs will retrieve harmonized MEs by default 

#Columns = Cell; Row = Value PC1 for each Module
hMEs <- GetMEs(cts)
write.table(hMEs, file = 'harmonized_ModuleEigengenes.csv', sep = ',', quote = F, col.names = T, row.names = T)

#Compute module connectivity
#We want to focus on hub genes
#So we want to know eigengene-based connectivity (kME) of each gene. 
#Computes kME values in the full single-cell dataset 
#Calculates pairwise correlations between genes and module eigenegenes 

cts <- ModuleConnectivity(cts)

#Now we can plot the genes in each module ranked by kME 

PlotKMEs(cts, ncol = 4)

#Get module assignment table 
modules <- GetModules(cts)
write.table(modules, file = 'modules_table.csv', sep = ',', quote = F, col.names = T, row.names = F)

#Extract hub genes with top N hubs 
hub_df <- GetHubGenes(cts, n_hubs = 15)
write.table(hub_df, file = 'ModuleHubs.csv', sep = ',', quote = F, col.names = T, row.names = F)

saveRDS(cts, file = 'Singles_hdWGCNA.rds')

#Compute hub gene signature scores 
#We often compute a score for the overall signature of a set of genes. 
#It's an alternative way of summarizing the expression of a module. 

cts <- ModuleExprScore(cts, n_genes = 25, method = 'Seurat')


##Visualization 

#Module Feature Plots 
plot_list <- ModuleFeaturePlot(cts, features = 'hMEs', order = T)
wrap_plots(plot_list, ncol=4)

#Module Correlations 
#Visiaulize the correlation between each module based on hMEs, MEs, or hub gene scores using corrplot 
ModuleCorrelogram(cts)


#Plot modules by Ident 
mods <- colnames(hMEs); mods <- mods[mods != 'grey']

#add hMEs to Seurat meta-data: 

cts@meta.data <- cbind(cts@meta.data, hMEs)

#We are interested in HAR function, thus: 

Idents(cts) <- 'har'

DotPlot(cts, features = mods, group.by = 'har') +
  coord_flip() +
  RotatedAxis() + 
  scale_color_gradient2(high = 'red', mid = 'grey95', low = 'blue')


PlotDendrogram(cts, main='hdWGCNA Dendrogram')

####Network Analysis #####
require(igraph)
require(uwot)

#Creates a folder of PDFs with module networks for each module 
ModuleNetworkPlot(cts)

#Create a network plot combinign all hub genes together 
# hubgene network

hub5 <- HubGeneNetworkPlot(
  cts,
  mods = 'all', 
  n_hubs = 10, n_other=0,
  edge_prop = 0.5#, 
  #return_graph = T
)

saveRDS(hub, file = 'HubNetwork_n8.rds')
saveRDS(hub5, file = 'HubNetwork_n5.rds')
saveRDS(hub6, file = 'HubNetwork_n6.rds')

#UMAP of module genes

#Add UMAP for modules to Seurat object
#Use n_hubs (number of hub genes), n_neighbors, neighbors parameter for UMAP
cts <- RunModuleUMAP(cts, n_hubs = 10, n_neighbors = 15, min_dist = 0.1)

#Get hub gene UMAP table
require(ggrepel)

dumap_df <- GetModuleUMAP(cts)

#Adjust input to include hub labels
dumap_df_input <-  dumap_df %>% 
  as.data.frame() %>%
  mutate(label = NA) %>%
  rowwise() %>%
  mutate(label = replace(label, hub == 'hub', gene))

#plot Module UMAP with labeled hubs 
ggplot(dumap_df_input, aes(x = UMAP1, y = UMAP2, label = label)) + 
  geom_point(
    color = dumap_df$color, #color = module color 
    size = dumap_df$kME*2 #size of each point = intramodular connectivity
  ) +
  geom_text_repel(size = 3, max.overlaps = 15) +
  umap_theme()

#Alternative plotting
g<- ModuleUMAPPlot(
  cts,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=4 ,# how many hub genes to plot per module?
  vertex.label.cex = 2, #font size
  keep_grey_edges=FALSE
)

###Gene Enrichment Analysis ####
#install.packages('enrichR')
require(enrichR)
#BiocManager::install("GeneOverlap")
require(GeneOverlap)

#Pull enrichr databases
dbs <- c('GO_Biological_Process_2021','GO_Cellular_Component_2021','GO_Molecular_Function_2021')

#Run Enricher
cts <- RunEnrichr(
  cts, 
  dbs = dbs, 
  max_genes = Inf
)

#Get results table
enrich_df <- GetEnrichrTable(cts)

write.table(enrich_df, file = 'module_enrichment_GO.tsv', sep = '\t', quote = F, row.names = F, col.names = T)


#Visualize 

#Go Term Plots
EnrichrBarPlot(cts, 
               outdir = 'enrichr_plots', 
               n_terms = 10, 
               plot_size = c(5,7),
               logscale = T)

#Table of gene counts per module 
mods.genes <- read_csv(file = 'modules_table.csv')

mod.genecount <- table(modules$module) %>% as.data.frame() 
