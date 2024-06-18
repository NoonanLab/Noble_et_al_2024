#hdWGCNA linear regression by perturbation 
require('Seurat')
require('tidyverse')

####Format Data####

#load Seurat object 
cts <- read_rds(file = 'Singles_hdWGCNA.rds')

#Fetch HAR assignments 
har.table <- FetchData(cts, vars = c('har', 'rep', 'nFeature_RNA')) %>% 
  rownames_to_column(var = 'cell')

#Load module eigenvalues 
mod.tb <- read.table(file = 'harmonized_ModuleEigengenes.csv', sep = ',') %>% 
  rownames_to_column(var = 'cell') %>% 
  left_join(har.table, by = c('cell' = 'cell'))


####Run LM####

#Get column names (modules) as list, remove non-module colnames
cols <- colnames(mod.tb)[-c(1,14,19:21)]
cols

#Paste formula to test module response as function of HAR and REP (batch)
forms <- paste(cols, " ~nFeature_RNA + relevel(factor(har), ref = 'Non-Targeting') + rep")
forms


#Run Linear Model
test.table <- lapply(forms, function(x){
  model <- lm(x, mod.tb)
  coef_df <- summary(model)$coefficient %>% #Extract p-Values 
    as.data.frame() %>% 
    mutate(mod = x) %>% 
    mutate(p.adj = p.adjust(`Pr(>|t|)`, method = 'fdr')) %>% 
    rownames_to_column(var = 'tests') 
  return(coef_df)
})

#Bind models
comb.df <- do.call(rbind, test.table) 
