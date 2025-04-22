#Analyze Batch 3, which is unaffected by the p53 mutation, against the others. 

#GOALS
#1. Compare Batch3 to others by Differential Expression
#2. Re-Run Linear Regression on Module Expression using only Batch 3

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
require(forcats)

#Load Data
cts <- readRDS('/gpfs/gibbs/pi/noonan/man59/Sequencing/Perturb/Align/H3_Perturb_Full_Agg_UnNorm/outs/count/Singles_hdWGCNA.rds')

#Subset cells by batch into different character vectors of CBCs

#Batch 3 NTC and Hits
b3.ntc.cells <- FetchData(cts, vars = c('rep', 'har')) %>%
  filter(rep == 'rep3', har == 'Non-Targeting') %>% 
  rownames_to_column(var = 'barcode') %>%
  pull(barcode)



#This isn't very helpful since we are just identifying generic batch effects. 
#Perhaps subset those with specific guides and then perform DE on Batch 1 and 3, for example. 

##Re-Run Linear Regression 

#Fetch HAR assignments 
har.table <- FetchData(cts, vars = c('har', 'rep', 'nFeature_RNA')) %>% 
  rownames_to_column(var = 'cell') %>%
  filter(rep == 'rep3')

#Load eigenvalues 
mod.tb <- read.table(file = '/gpfs/gibbs/pi/noonan/man59/Sequencing/Perturb/Align/H3_Perturb_Full_Agg_UnNorm/outs/count/WGCNA/harmonized_ModuleEigengenes.csv', sep = ',') %>% 
  rownames_to_column(var = 'cell') %>% 
  right_join(har.table, by = c('cell' = 'cell'))

#Get column names (modules) as list, remove non-module colnames
cols <- colnames(mod.tb)[-c(1,14,19:21)]
cols

#Format as list of formulas for linear regression
#Test module response as function of HAR and REP (batch)
#This includes no 'rep' since it's only Rep3
forms <- paste(cols, " ~nFeature_RNA + relevel(factor(har), ref = 'Non-Targeting')")
forms

#Run Linear Model
test.table <- lapply(forms, function(x){
  model <- lm(x, mod.tb)
  coef_df <- summary(model)$coefficient %>% #Extract p-Values 
    as.data.frame() %>% 
    mutate(mod = x) %>% 
    mutate(p.adj = p.adjust(`Pr(>|t|)`, method = 'fdr')) %>% #Add column specifying formula
    rownames_to_column(var = 'tests') 
  return(coef_df)
})

#Load counts and names 
mod.names <- read.csv('/gpfs/gibbs/pi/noonan/man59/Sequencing/Perturb/Align/H3_Perturb_Full_Agg_UnNorm/outs/count/WGCNA/module_labels.csv')

#Format the table to be readable
comb.df <- do.call(rbind, test.table) %>% 
  rowwise() %>% 
  mutate(module = str_split(mod, ' ')[[1]][1],
         har.name = str_replace(tests, 'relevel\\(.+\\)', '')) %>%
  ungroup() %>% 
  dplyr::select(module, har.name, Estimate, `Pr(>|t|)`, p.adj) %>%
  filter(!grepl('Intercept|nFeature', har.name)) %>%
  left_join(mod.names, by = c('module' = 'Var1')) %>%
  dplyr::select(mod.name, Term, har.name, Estimate, `Pr(>|t|)`, `p.adj`, new.color, mod.cat) %>%
  rename(`Pr(>|t|)` = 'p.value')

#Write results
write.table(comb.df, file = 'Module_LR_Batch3_Reanalysis.txt', sep = '\t', quote = F, col.names = T, row.names = F)

#Compare results to original by correlation 

#Load original LR and format the same way as above
orig.df <- read_tsv('/gpfs/gibbs/pi/noonan/man59/Sequencing/Perturb/Align/H3_Perturb_Full_Agg_UnNorm/outs/count/WGCNA/LM_Addv5_FullResults.tsv') %>%
  rowwise() %>% 
  mutate(module = str_split(mod, ' ')[[1]][1],
         har.name = str_replace(tests, 'relevel\\(.+\\)', '')) %>%
  ungroup() %>% 
  dplyr::select(module, har.name, Estimate, `Pr(>|t|)`, p.adj) %>%
  filter(!grepl('rep|Intercept|nFeature', har.name)) %>% #This one includes a Rep feature to account for
  left_join(mod.names, by = c('module' = 'Var1')) %>%
  dplyr::select(mod.name, Term, har.name, Estimate, `Pr(>|t|)`, `p.adj`, new.color, mod.cat) %>%
  rename(`Pr(>|t|)` = 'p.value')

#Merge 
joined.df <- orig.df %>%
  left_join(comb.df, by = c('mod.name' = 'mod.name', 'har.name' = 'har.name'))

#Compare 
ggplot(joined.df, aes(x = Estimate.x, y = Estimate.y)) + 
  geom_point(size = 0.5)

####Compare to Batch 1 
##Re-Run Linear Regression 

#Fetch HAR assignments 
har.table <- FetchData(cts, vars = c('har', 'rep', 'nFeature_RNA')) %>% 
  rownames_to_column(var = 'cell') %>%
  filter(rep %in% c('rep1', 'rep2'))

#Load eigenvalues 
mod.tb <- read.table(file = '/gpfs/gibbs/pi/noonan/man59/Sequencing/Perturb/Align/H3_Perturb_Full_Agg_UnNorm/outs/count/WGCNA/harmonized_ModuleEigengenes.csv', sep = ',') %>% 
  rownames_to_column(var = 'cell') %>% 
  right_join(har.table, by = c('cell' = 'cell'))

#Get column names (modules) as list, remove non-module colnames
cols <- colnames(mod.tb)[-c(1,14,19:21)]
cols

#Format as list of formulas for linear regression
#Test module response as function of HAR and REP (batch)
#This includes no 'rep' since it's only Rep3
forms <- paste(cols, " ~nFeature_RNA + relevel(factor(har), ref = 'Non-Targeting') + rep")
forms

#Run Linear Model
test.table <- lapply(forms, function(x){
  model <- lm(x, mod.tb)
  coef_df <- summary(model)$coefficient %>% #Extract p-Values 
    as.data.frame() %>% 
    mutate(mod = x) %>% 
    mutate(p.adj = p.adjust(`Pr(>|t|)`, method = 'fdr')) %>% #Add column specifying formula
    rownames_to_column(var = 'tests') 
  return(coef_df)
})

#Load counts and names 
mod.names <- read.csv('/gpfs/gibbs/pi/noonan/man59/Sequencing/Perturb/Align/H3_Perturb_Full_Agg_UnNorm/outs/count/WGCNA/module_labels.csv')

#Format the table to be readable
comb.df.b1 <- do.call(rbind, test.table) %>% 
  rowwise() %>% 
  mutate(module = str_split(mod, ' ')[[1]][1],
         har.name = str_replace(tests, 'relevel\\(.+\\)', '')) %>%
  ungroup() %>% 
  dplyr::select(module, har.name, Estimate, `Pr(>|t|)`, p.adj) %>%
  filter(!grepl('Intercept|nFeature', har.name)) %>%
  left_join(mod.names, by = c('module' = 'Var1')) %>%
  dplyr::select(mod.name, Term, har.name, Estimate, `Pr(>|t|)`, `p.adj`, new.color, mod.cat) %>%
  rename(`Pr(>|t|)` = 'p.value')

#Write results
write.table(comb.df.b1, file = 'Module_LR_Batch1_Reanalysis.txt', sep = '\t', quote = F, col.names = T, row.names = F)

#Compare Batch 1 and 3 
#Merge 
joined.df.batch <- comb.df.b1 %>%
  left_join(comb.df, by = c('mod.name' = 'mod.name', 'har.name' = 'har.name'))

#Compare 
library(ggpubr)
library(viridis)
library(ggpointdensity)

ggplot(joined.df.batch, aes(x = Estimate.x, y = Estimate.y)) + 
  geom_pointdensity(size = 0.5) + 
  scale_color_viridis() + 
  geom_smooth(method = "lm", 
              formula = y ~ x) + 
  stat_cor(label.x.npc = "left",
           label.y.npc = "top",
           size = 3) + 
  xlab('LR Effect Size (Original)') + 
  ylab('LR Effect Size (Batch3)') +
  scale_x_continuous(expand = expansion(mult = c(0,0))) + 
  theme_cowplot() + 
  panel_border(color = 'black') + 
  theme(aspect.ratio = 1:1)


