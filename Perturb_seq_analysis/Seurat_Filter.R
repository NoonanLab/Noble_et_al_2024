#Load packages
require(tidyverse)
require(ggplot2)
#Seurat Packages
library(Seurat)

#ScdblFinder equires old version of Matrix
remotes::install_version("Matrix", version = "1.6-1")
library(Matrix)
packageVersion("Matrix")

#Seurat equires new version of Matrix
remotes::install_version("Matrix", version = "1.6-4")
library(Matrix)
packageVersion("Matrix")


#Install scDblFinder
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("scDblFinder")

require(scDblFinder)

#Load files
cts1 <- readRDS(file = 'Batch1.rds')
cts2 <- readRDS(file = 'Batch2.rds')
cts3 <- readRDS(file = 'Batch3.rds')

#Convert
cts1 <- as.SingleCellExperiment(cts1)
cts2 <- as.SingleCellExperiment(cts2)
cts3 <- as.SingleCellExperiment(cts3)


cts1 <- scDblFinder(cts1, dbr = 0.05, verbose = T, samples = cts1$lane)
cts2 <- scDblFinder(cts2, dbr = 0.05, verbose = T, samples = cts2$lane)
cts3 <- scDblFinder(cts3, dbr = 0.05, verbose = T, samples = cts3$lane)

#Return to Seurat Object format
cts1 <- as.Seurat(cts1, data = NULL)
cts2 <- as.Seurat(cts2, data = NULL)
cts3 <- as.Seurat(cts3, data = NULL)

#Get Data
b1.dbl.calls <- FetchData(cts1, vars = c('scDblFinder.class', 'scDblFinder.score'))
b2.dbl.calls <- FetchData(cts2, vars = c('scDblFinder.class', 'scDblFinder.score'))
b3.dbl.calls <- FetchData(cts3, vars = c('scDblFinder.class', 'scDblFinder.score'))

#Write tables
write.table(b1.dbl.calls, file = 'batch1_dbl_calls.csv', sep = ',', quote = F, row.names = T)
write.table(b2.dbl.calls, file = 'batch2_dbl_calls.csv', sep = ',', quote = F, row.names = T)
write.table(b3.dbl.calls, file = 'batch3_dbl_calls.csv', sep = ',', quote = F, row.names = T)

#Remove Doublets
cts1 <- subset(cts1, scDblFinder.class == 'singlet')
cts2 <- subset(cts2, scDblFinder.class == 'singlet')
cts3 <- subset(cts3, scDblFinder.class == 'singlet')


#Filter by Feature counts, total counts, and mitochondrial reads
#List Mitochondrial genes 
DefaultAssay(cts1) <- 'RNA'
DefaultAssay(cts2) <- 'RNA'
DefaultAssay(cts3) <- 'RNA'

#Calculate percent mitochondrial reads per cell
cts1[['percent.mt']] <- PercentageFeatureSet(cts1, pattern = '^MT-') #This dataset has Mitochondrial genes beginning with MTND
cts2[['percent.mt']] <- PercentageFeatureSet(cts2, pattern = '^MT-') #This dataset has Mitochondrial genes beginning with MTND
cts3[['percent.mt']] <- PercentageFeatureSet(cts3, pattern = '^MT-') #This dataset has Mitochondrial genes beginning with MTND

#Calculate 2nd SD nFeature RNA in data 
feat.min.1 <- round(mean(cts1$nFeature_RNA) - 2 * sd(cts1$nFeature_RNA), digits = -2)
feat.min.2 <- round(mean(cts2$nFeature_RNA) - 2 * sd(cts2$nFeature_RNA), digits = -2)
feat.min.3 <- round(mean(cts3$nFeature_RNA) - 2 * sd(cts3$nFeature_RNA), digits = -2)

#Filter data by cell quality 
cts1 <-subset(cts1, subset = nFeature_RNA > feat.min.1 & percent.mt < 10)
cts2 <-subset(cts2, subset = nFeature_RNA > feat.min.2 & percent.mt < 10)
cts3 <-subset(cts3, subset = nFeature_RNA > feat.min.3 & percent.mt < 10)

#Save RDS
saveRDS(cts1, file = 'batch1_filtered.rds')
saveRDS(cts2, file = 'Batch2_filtered.rds')
saveRDS(cts3, file = 'Batch3_filtered.rds')

#Saved Batch 2 as Batch1 and Batch 1 as batch1; go fix

#Plots
Idents(cts) <- 'rep'

#Violin by batch
VlnPlot(cts, 
        pt.size = 0,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)

#Scatter by batch
plot1 <- FeatureScatter(cts, feature1 = "nCount_RNA", feature2 = "percent.mt", shuffle = T, raster = T)
plot2 <- FeatureScatter(cts, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", shuffle = T, raster = T)
plot3 <- FeatureScatter(cts, feature1 = 'nCount_RNA', feature2 = 'nCount_crispr', shuffle = T, raster = T)

CombinePlots(plots = list(plot1, plot2, plot3))
