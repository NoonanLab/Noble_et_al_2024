#To format a table of outputs so that barcodes can be neatly updated to their names in the final Seurat object
require(Matrix)
require(tidyverse)
require(stringr)
require(ggplot2)
require(cowplot)

df <- data.frame(file.name = list.files(pattern = 'vartrix_output*'))

sample.df <- df %>%
  rowwise() %>%
  #extract sample information
  mutate(rep = str_split_1(file.name, '_')[3],
         lane = str_split_1(file.name, '_')[4],
         allele = str_split_1(file.name, '_')[5]) %>%
  #usable sample tag
  mutate(sample = paste0(rep, '_', lane)) %>%
  #to get correct barcode suffixes, 8 has to be added to Batch2 and 16 has to be added to Batch3 
  mutate(lane.num = 0) %>%
  mutate(lane.num = replace(lane.num, rep == 'Rep2', 8),
         lane.num = replace(lane.num, rep == 'Rep3', 16)) %>%
  mutate(new.lane.num = as.numeric(str_split_1(lane, 'L')[2]) + lane.num)

#Analyze VarTrix

#Get list of samples 
samples <- unique(sample.df$sample)

#Initialize results
results <- list() 

for(sample in samples){
  ref.file <- paste0("vartrix_output_", sample, "_ref")
  alt.file <- paste0("vartrix_output_", sample, "_alt")
  rep <- str_split_1(sample, '_')[1]
  lane <- str_split_1(sample, '_L')[2]
  barcode.path <- paste0('/gpfs/gibbs/pi/noonan/man59/Sequencing/Perturb/Align/', rep, '_L', lane, '/outs/filtered_feature_bc_matrix/barcodes.tsv.gz')
  
  #Matrix of alt allele (p53)
  alt <- readMM(alt.file) %>% as.matrix() %>% as.data.frame()
  #Matrix of ref allele (wt)
  ref <- readMM(ref.file) %>% as.matrix %>% as.data.frame()
  #Get list of experiment barcodes
  barcodes <- read.table(barcode.path)
  
  #Construct alt and ref matrix 
  colnames(alt) <- barcodes$V1 
  colnames(ref) <- barcodes$V1 
  
  #Make into data_frame
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
    left_join(alt.frame, by = c('barcode' = 'barcode')) %>% 
    mutate(lane.num = 0) %>%
    mutate(lane.num = replace(lane.num, rep == 'Rep2', 8),
           lane.num = replace(lane.num, rep == 'Rep3', 16)) %>%
    mutate(new.lane.num = as.numeric(lane) + lane.num) %>%
    rowwise() %>%
    mutate(new.barcode = paste0(str_split_1(barcode, '-')[1],'-', new.lane.num)) %>%
    mutate(sample = sample)
  
  results[[length(results) + 1]] <- comb.frame
}

#Bind results into single frame
results.df <- do.call(rbind, results)


#Calculate proportions 
rate.df <- results.df %>%
  rowwise() %>%
  mutate(batch = str_split_1(sample, '_')[1]) %>%
  mutate(total_counts = ref_counts + alt_counts) %>% 
  mutate(alt_rate = alt_counts/total_counts,
         wt_rate = ref_counts/total_counts) %>%
  mutate(status = 'wt') %>%
  mutate(status = replace(status, total_counts == 0, 'no.call'),
         status = replace(status, total_counts > 0 & alt_rate > 0.2, 'p53.mutant'))

#Write results 
write.table(rate.df, file = 'Vartrix_CBC_Results.txt', sep = '\t', quote = F, row.names = F, col.names = T)


#Summarize 
rate.table <- rate.df %>%
  group_by(batch, status) %>%
  summarise(count = n(),
            depth = sum(total_counts)) %>%
  ungroup() %>% group_by(batch) %>%
  mutate(batch.depth = sum(depth)) %>%
  dplyr::select(batch, status, count, batch.depth) %>%
  pivot_wider(values_from = count, names_from = status) %>%
  mutate(p53.mutant = replace_na(p53.mutant, 0)) %>%
  mutate(mut.rate = (p53.mutant / (p53.mutant + wt))*100)

#Write results 
write.table(rate.table, file = 'Vartrix_Batch_Summary.txt', sep = '\t', quote = F, row.names = F, col.names = T)

#SUMMARY FIGURES

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

