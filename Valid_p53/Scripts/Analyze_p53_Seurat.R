#Single Cell Analysis of p53 mutant cells 
#GOALS
#1. Format the VarTrix Results as MetaData for the Seurat Object
#2. Visualize the mutant cells 
#3. Perform differential expression of p53 mutant versus wt cells

#Load Packages
require(Matrix)
require(tidyverse)
require(stringr)
require(ggplot2)
require(cowplot)
require(Seurat)
require(DESeq2)
require(ggrepel)
require(ggpmisc)


#Load Data
cts <- readRDS('/gpfs/gibbs/pi/noonan/man59/Sequencing/Perturb/Align/H3_Perturb_Full_Agg_UnNorm/outs/count/Singles_hdWGCNA.rds')

mut.df <- read_tsv('Vartrix_CBC_Results.txt') %>%
  dplyr::select(new.barcode, sample, batch, total_counts, alt_rate, wt_rate, status)

#Get CBCs from Seurat as DF, plust additional medadata 
cbc.frame <- FetchData(cts, vars = c('nCount_RNA', 'rep', 'lane', 'har', 'seurat_clusters')) %>%
  rownames_to_column(var = 'barcode') 

#Join calls to this 
cbc.frame <- cbc.frame %>%
  left_join(mut.df, by = c('barcode' = 'new.barcode')) #inner join keeps only rows from X that match

#Write results
write.table(cbc.frame, file = 'Results_Filtered_CBCs_p53.txt', sep = '\t', quote = F, row.names = F, col.names = T)

#Create a metadata frame 
mut.metadata <- cbc.frame %>%
  dplyr::select(barcode, status) %>%
  column_to_rownames(var = 'barcode')

#Add metadata to seurat object 
cts <- AddMetaData(cts, mut.metadata)

#Subset cells with calls 
mut.cells <- WhichCells(cts, expression = status == 'p53.mutant')
wt.cells <- WhichCells(cts, expression = status == 'wt')

#Visualize location of mutant cells
DimPlot(cts, cells.highlight = list(mut.cells, wt.cells), cols.highlight = c('#00BFC7','#F8766D'))
DimPlot(cts, cells.highlight = wt.cells)

#Perform differential expression between these cells 
de <- FindMarkers(cts, ident.1 = mut.cells, ident.2 = wt.cells, test.use = 'DESeq2') %>%
  mutate(FDR = p.adjust(p_val, method = 'BH'))

##Visualize the results 
#Volcano Plot of DE Genes
volcano.frame <- de %>%
  rownames_to_column(var = 'gene') %>% 
  mutate(sig = 'no') %>%
  mutate(sig = replace(sig, FDR < 0.05, 'yes')) %>%
  mutate(direction = 'down') %>%
  mutate(direction = replace(direction, avg_log2FC > 0, 'up')) %>%
  mutate(color = 'ns') %>%
  mutate(color = replace(color, direction == 'up' & sig == 'yes', 'up'),
        color = replace(color, direction == 'down' & sig == 'yes', 'down')) %>%
  mutate(label = NA) %>%
  rowwise() %>%
  mutate(label = replace(label, sig == 'yes', gene))

#Set Colors
volcano.colors <- c('#90e0ef', 'grey80', '#ffb703')
names(volcano.colors) <- c("down", "ns", "up")

#Plot 
ggplot(volcano.frame, aes(x = avg_log2FC, y = -log(FDR,10), color = color, label = label)) + 
  geom_point() + 
  scale_color_manual(values = volcano.colors) + 
  geom_hline(yintercept = -log(0.05, 10), color = 'black') + 
  geom_vline(xintercept = 0, linetype = 'dotted', color = 'grey30') + 
  geom_text_repel(color = 'black', size = 3) + 
  scale_y_continuous(expand = expansion(mult = c(0,2))) +
  scale_x_continuous(limits = symmetric_limits) + 
  theme_cowplot() + 
  panel_border(color = 'black') + 
  theme(aspect.ratio = 1:1)
  

#Make new barplots showing mutation rate in filtered cells 
rate.table <- cbc.frame %>%
  group_by(batch, status) %>%
  summarise(count = n(),
            depth = sum(total_counts)) %>%
  dplyr::select(batch, status, count) %>%
  pivot_wider(values_from = count, names_from = status) %>%
  mutate(p53.mutant = replace_na(p53.mutant, 0)) %>%
  mutate(mut.rate = (p53.mutant / (p53.mutant + wt))*100)

#Write results 
write.table(rate.table, file = 'Vartrix_Batch_Summary_FilteredCells.txt', sep = '\t', quote = F, row.names = F, col.names = T)

#BarPlot Summary 
barplot.df <- rate.table %>%
  dplyr::select(batch, p53.mutant, wt, mut.rate) %>%
  pivot_longer(cols = 2:3, names_to = 'status')

#Cell Counts
a <- ggplot(barplot.df, aes(x = batch, y = value, fill = status)) + 
  geom_bar(stat = 'identity', position = 'dodge', width = .7) + 
  geom_text(aes(label=value), vjust=0, position = position_dodge(width=.7)) +
  xlab("Batch") + 
  ylab("nCells Detected") + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) + 
  theme_cowplot() + 
  theme(aspect.ratio = 1:1)

#Mutation Rate
b <- ggplot(
  (barplot.df %>% dplyr::select(batch, mut.rate) %>% group_by(batch) %>% summarise(mean(mut.rate))), 
  aes(x = batch, y = `mean(mut.rate)`, group = batch)) + 
  geom_bar(stat = 'identity', width = .7) + 
  geom_text(aes(label= round(`mean(mut.rate)`, 2), vjust=-1)) +
  xlab("Batch") + 
  ylab("p53 Mutation Detection Rate") + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) + 
  theme_cowplot() + 
  theme(aspect.ratio = 1:1)

require(patchwork)
a + b



