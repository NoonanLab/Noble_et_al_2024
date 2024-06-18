#Perturb-seq differential expression by target

library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(readr)

#Load seurat object 
cts <- read_rds(file = 'Singles_hdWGCNA.rds')

#Set Idents to HARs 
Idents(cts) <- 'har'

#Specify NTC Cells
ntc.list <- c('NTC1', 'NTC2', 'NTC3', 'NTC4', 'NTC5', 'NTC-10x')
ntc.cells <- subset(cts, subset = guide %in% ntc.list)
ntc.vec <- colnames(ntc.cells)

#List idents
ident.list <- data.frame(cts@meta.data$har) %>%
  table() %>%
  as.data.frame() %>% 
  filter(!cts.meta.data.har == 'Non-Targeting') %>%
  pull(cts.meta.data.har) %>%
  as.vector() 

#Run DE using FindMarkers against NTC cells
de.table <- lapply(ident.list, function(x){
  de <- FindMarkers(cts, ident.1 = x, ident.2 = ntc.vec, min.pct = 0.01, logfc.threshold = 0.05)
  de <- de %>% mutate(test = x,
                      p.adj = p.adjust(p_val, method = 'fdr')) %>%
    rownames_to_column(var = 'gene')
  return(de)
})

#Bind
de.df <- do.call(rbind, de.table) 


#Summarize results
de.summary <- de.df %>%
  filter(p.adj < 0.05) %>% 
  mutate(down = NA) %>%
  mutate(up = NA) %>% 
  rowwise() %>% 
  mutate(down = replace(down, avg_log2FC < 0, gene), 
         up = replace(up, avg_log2FC > 0, gene)) %>%
  ungroup() %>%
  group_by(test) %>% 
  arrange(p.adj) %>% 
  summarise(perturbed.down = (paste(na.omit(down), collapse = ',')),
            perturbed.up = paste(na.omit(up), collapse = ',')) %>%
  mutate(perturbed.down = na.omit(perturbed.down)) %>%
  rename('har' = 'test')

write.table(de.summary, file = 'Perturbation_DE_Summary_byHAR.txt', sep = '\t')  

#Summarize overal number of significant perturbations by directionality 
sig.de <- de.df %>%
  filter(p.adj < 0.05)

gene.summary.down <- table(((sig.de) %>% filter(avg_log2FC < 0) %>% select(gene, test))$gene) %>% as.data.frame()
gene.summary.up <- table(((sig.de) %>% filter(avg_log2FC > 0) %>% select(gene, test))$gene) %>% as.data.frame()


