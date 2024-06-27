#Linear Regression - by guide


#hdWGCNA linear regression by perturbation 
require('Seurat')
require('tidyverse')


####Format Data####

#Load Seurat object 
cts <- read_rds(file = 'Singles_hdWGCNA.rds')

#Fetch guide assignments 
har.table <- FetchData(cts, vars = c('guide', 'rep', 'nFeature_RNA')) %>% 
  rownames_to_column(var = 'cell') %>%
  mutate(guide2 = guide) %>%
  mutate(guide2 = replace(guide2, guide2 %in% c('NTC1', 'NTC2', 'NTC3', 'NTC4', 'NTC5', 'NTC-10x'), 'NTC'))

#Load eigenvalues 
mod.tb <- read.table(file = 'harmonized_ModuleEigengenes.csv', sep = ',') %>% 
  rownames_to_column(var = 'cell') %>% 
  left_join(har.table, by = c('cell' = 'cell'))

#Get column names (modules) as list, remove non-module colnames
cols <- colnames(mod.tb)[-c(1,14,19:22)]
cols

#Format as list of formulas for linear regression
#Test module response as function of HAR and REP (batch)
forms <- paste(cols, " ~nFeature_RNA + relevel(factor(guide2), ref = 'NTC') + rep")
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

#Bind
comb.df <- do.call(rbind, test.table) 



#Summarize results
mod.names <- read.csv(file = 'module_labels.csv')

sum.comb.df <- comb.df %>% 
  filter(!grepl('rep|Intercept|nFeature_RNA', tests)) %>%
  rowwise() %>%
  mutate(mod = str_split(mod, ' ')[[1]][1]) %>%
  mutate(mod = str_split(mod, ' ')[[1]][1]) %>%
  mutate(har.name = str_replace(tests, 'relevel\\(.+\\)', '')) %>% 
  mutate(guide.num = 'A') %>%
  mutate(guide.num = replace(guide.num, endsWith(har.name, '.1'), 'B')) %>%
  mutate(har = sub('[.]1$', '', har.name)) %>% 
  mutate(comparison = paste0(har,mod)) %>%
  left_join(mod.names, by = c('mod' = 'Var1'))








