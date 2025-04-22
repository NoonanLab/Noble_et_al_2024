#Read Vartrix output 
setwd("~/gibbs/Valid_p53")

library(Seurat)
library(Matrix)
library(stringr)
library(tidyverse)

#Matrix of alt allele (p53 mutant)
alt <- readMM('output') %>% as.matrix() %>% as.data.frame()
#Matrix of ref allele (wt)
ref <- readMM('output_ref') %>% as.matrix %>% as.data.frame()
#Get barcode sequences/names
barcodes <- read.table('/gpfs/gibbs/pi/noonan/man59/Sequencing/H3_Perturb_Rep2/CellRanger/Lane_1/outs/filtered_feature_bc_matrix/barcodes.tsv.gz')

#Construct alt and ref matrix 
colnames(alt) <- barcodes$V1 
colnames(ref) <- barcodes$V1 

alt.frame <- alt %>%
  t() %>% 
  as.data.frame() %>%
  rownames_to_column(var = 'barcode') %>%
  rename('alt_counts' = 'V1')

ref.frame <- ref %>%
  t() %>% 
  as.data.frame() %>%
  rownames_to_column(var = 'barcode') %>%
  rename('ref_counts' = 'V1')

comb.frame <- ref.frame %>%
  left_join(alt.frame, by = c('barcode' = 'barcode'))
  

