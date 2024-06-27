#Perturb-seq Guide Assignment 
#Split by batch/replicate


#Make DF of guide counts 
rep1.guides <- as.data.frame(as.matrix(cts1[['crispr']]@counts))
rep2.guides <- as.data.frame(as.matrix(cts2[['crispr']]@counts))
rep3.guides <- as.data.frame(as.matrix(cts3[['crispr']]@counts))

#Replace counts with matrix to fill 0s 
rep1.matrix <- as.data.frame(t(as.matrix(cts1[['crispr']]@counts))) %>%
  replace(is.na(.), 0)
rep2.matrix <- as.data.frame(t(as.matrix(cts2[['crispr']]@counts))) %>%
  replace(is.na(.), 0)
rep3.matrix <- as.data.frame(t(as.matrix(cts3[['crispr']]@counts))) %>%
  replace(is.na(.), 0)

#Transpose to dataframe 
#Apply quantile 

#Rep1 
rep1.threshold <- rep1.matrix %>% 
  t() %>%
  as.data.frame()
rep1.threshold <- rep1.threshold %>% 
  mutate(threshold = apply(rep1.threshold, 1, quantile, probs = c(0.996)))
thresh1 <- as.data.frame(rep1.threshold$threshold)

#Rep2 
rep2.threshold <- rep2.matrix %>% 
  t() %>%
  as.data.frame()
rep2.threshold <- rep2.threshold %>% 
  mutate(threshold = apply(rep2.threshold, 1, quantile, probs = c(0.996)))
thresh2 <- as.data.frame(rep2.threshold$threshold)

#Rep3 
rep3.threshold <- rep3.matrix %>% 
  t() %>%
  as.data.frame()
rep3.threshold <- rep3.threshold %>% 
  mutate(threshold = apply(rep3.threshold, 1, quantile, probs = c(0.996)))
thresh3 <- as.data.frame(rep3.threshold$threshold)

#Make full matrix of thresholds 
#This is to subtract from counts 
threshold.m1 <- rep1.threshold$threshold %>%
  rep(times = ncol(rep1.guides)) %>%
  matrix(ncol = ncol(rep1.guides)) 

threshold.m2 <- rep2.threshold$threshold %>%
  rep(times = ncol(rep2.guides)) %>%
  matrix(ncol = ncol(rep2.guides)) 

threshold.m3 <- rep3.threshold$threshold %>%
  rep(times = ncol(rep3.guides)) %>%
  matrix(ncol = ncol(rep3.guides))

# subtract UMI threshold values from guides matrix
# adding and subtracting is to allow cells meeting the threshold to remain
rep1.guides.filter <- rep1.guides + 1 #add 1 to counts
rep1.guides.filter <- rep1.guides.filter - threshold.m1
rep1.guides.filter <- rep1.guides.filter - 1 #subtract
rep1.guides.filter[rep1.guides.filter < 0] <- 0 # set negative values to zero

rep2.guides.filter <- rep2.guides + 1 #add 1 to counts
rep2.guides.filter <- rep2.guides.filter - threshold.m2
rep2.guides.filter <- rep2.guides.filter - 1 #subtract
rep2.guides.filter[rep2.guides.filter < 0] <- 0 # set negative values to zero

rep3.guides.filter <- rep3.guides + 1 #add 1 to counts
rep3.guides.filter <- rep3.guides.filter - threshold.m3
rep3.guides.filter <- rep3.guides.filter - 1 #subtract
rep3.guides.filter[rep3.guides.filter < 0] <- 0 # set negative values to zero

# convert guide matrix to tibble
rep1.tb <- tibble(cell = colnames(rep1.guides),
                  x0 = apply(rep1.guides, 2, c, simplify = F),            # raw guide counts
                  x  = apply(rep1.guides.filter, 2, c, simplify = F)) %>% # filtered counts
  rowwise() %>%
  mutate(ug.raw = sum(x0 > 0),     # number of unique guides, before filtering
         ug.filtered = sum(x > 0)) # number of unique guides, after filtering

rep2.tb <- tibble(cell = colnames(rep2.guides),
                  x0 = apply(rep2.guides, 2, c, simplify = F),            # raw guide counts
                  x  = apply(rep2.guides.filter, 2, c, simplify = F)) %>% # filtered counts
  rowwise() %>%
  mutate(ug.raw = sum(x0 > 0),     # number of unique guides, before filtering
         ug.filtered = sum(x > 0)) # number of unique guides, after filtering

rep3.tb <- tibble(cell = colnames(rep3.guides),
                  x0 = apply(rep3.guides, 2, c, simplify = F),            # raw guide counts
                  x  = apply(rep3.guides.filter, 2, c, simplify = F)) %>% # filtered counts
  rowwise() %>%
  mutate(ug.raw = sum(x0 > 0),     # number of unique guides, before filtering
         ug.filtered = sum(x > 0)) # number of unique guides, after filtering

#Compute stats
rep1.stats <- rep1.tb %>%
  # z-score rows (cells)
  mutate(x.scale = list(scale(x)[, 1])) %>%
  # sort z-scores, highest to lowest
  mutate(x.scale.sort = list(sort(x.scale, decreasing = T, index.return = T))) %>%
  # stats
  mutate(total.counts = sum(x)) %>%
  # z-score and counts for top N guides (N = 5)
  mutate(topN.zscore = list(x.scale.sort$x[1:5]),
         topN.counts = list(x[x.scale.sort$ix[1:5]]),
         topN.guides = list(names(x.scale.sort$x[1:min(5, ug.filtered)]))) %>%
  # top guide stats
  mutate(top.guide     = topN.guides[1],    # top guide
         top.count.pct = topN.counts[1] / total.counts, # percent of total counts
         top.count     = topN.counts[1],    # counts for top guide
         count.2       = topN.counts[2],    # counts for 2nd highest guide
         count.diff  = topN.counts[1] - topN.counts[2], # count difference
         count.prop  = topN.counts[2]/topN.counts[1],   # proportion between top 2 guides
         zscore.max  = topN.zscore[1],  # zscore for top guide
         zscore.2    = topN.zscore[2],  # zscore for 2nd highest guide
         zscore.diff = topN.zscore[1] - topN.zscore[2]) %>% # zscore difference
  # collapse lists to string
  mutate(topN.zscore = str_c(round(topN.zscore, 1), collapse = ', '),
         topN.counts = str_c(topN.counts, collapse = ', '),
         topN.guides = str_c(topN.guides, collapse = ', ')) %>%
  # clean up
  relocate(cell, starts_with('ug'), total.counts, 
           top.guide, top.count.pct:zscore.diff) %>%
  dplyr::select(-starts_with('x'))

rep2.stats <- rep2.tb %>%
  # z-score rows (cells)
  mutate(x.scale = list(scale(x)[, 1])) %>%
  # sort z-scores, highest to lowest
  mutate(x.scale.sort = list(sort(x.scale, decreasing = T, index.return = T))) %>%
  # stats
  mutate(total.counts = sum(x)) %>%
  # z-score and counts for top N guides (N = 5)
  mutate(topN.zscore = list(x.scale.sort$x[1:5]),
         topN.counts = list(x[x.scale.sort$ix[1:5]]),
         topN.guides = list(names(x.scale.sort$x[1:min(5, ug.filtered)]))) %>%
  # top guide stats
  mutate(top.guide     = topN.guides[1],    # top guide
         top.count.pct = topN.counts[1] / total.counts, # percent of total counts
         top.count     = topN.counts[1],    # counts for top guide
         count.2       = topN.counts[2],    # counts for 2nd highest guide
         count.diff  = topN.counts[1] - topN.counts[2], # count difference
         count.prop  = topN.counts[2]/topN.counts[1],   # proportion between top 2 guides
         zscore.max  = topN.zscore[1],  # zscore for top guide
         zscore.2    = topN.zscore[2],  # zscore for 2nd highest guide
         zscore.diff = topN.zscore[1] - topN.zscore[2]) %>% # zscore difference
  # collapse lists to string
  mutate(topN.zscore = str_c(round(topN.zscore, 1), collapse = ', '),
         topN.counts = str_c(topN.counts, collapse = ', '),
         topN.guides = str_c(topN.guides, collapse = ', ')) %>%
  # clean up
  relocate(cell, starts_with('ug'), total.counts, 
           top.guide, top.count.pct:zscore.diff) %>%
  dplyr::select(-starts_with('x'))

#Rep 3 Stats
rep3.stats <- rep3.tb %>%
  # z-score rows (cells)
  mutate(x.scale = list(scale(x)[, 1])) %>%
  # sort z-scores, highest to lowest
  mutate(x.scale.sort = list(sort(x.scale, decreasing = T, index.return = T))) %>%
  # stats
  mutate(total.counts = sum(x)) %>%
  # z-score and counts for top N guides (N = 5)
  mutate(topN.zscore = list(x.scale.sort$x[1:5]),
         topN.counts = list(x[x.scale.sort$ix[1:5]]),
         topN.guides = list(names(x.scale.sort$x[1:min(5, ug.filtered)]))) %>%
  # top guide stats
  mutate(top.guide     = topN.guides[1],    # top guide
         top.count.pct = topN.counts[1] / total.counts, # percent of total counts
         top.count     = topN.counts[1],    # counts for top guide
         count.2       = topN.counts[2],    # counts for 2nd highest guide
         count.diff  = topN.counts[1] - topN.counts[2], # count difference
         count.prop  = topN.counts[2]/topN.counts[1],   # proportion between top 2 guides
         zscore.max  = topN.zscore[1],  # zscore for top guide
         zscore.2    = topN.zscore[2],  # zscore for 2nd highest guide
         zscore.diff = topN.zscore[1] - topN.zscore[2]) %>% # zscore difference
  # collapse lists to string
  mutate(topN.zscore = str_c(round(topN.zscore, 1), collapse = ', '),
         topN.counts = str_c(topN.counts, collapse = ', '),
         topN.guides = str_c(topN.guides, collapse = ', ')) %>%
  # clean up
  relocate(cell, starts_with('ug'), total.counts, 
           top.guide, top.count.pct:zscore.diff) %>%
  dplyr::select(-starts_with('x'))

#Score Cells

#Proportion Threshold 
count.prop.threshold <- 0.5

#Make calls 
rep1.calls <- rep1.stats %>% ungroup() %>%
  mutate(guide.call = case_when(ug.filtered == 0 ~ 'none_empty',
                                ug.filtered == 1 ~ 'single_true',
                                count.prop < count.prop.threshold ~ 'single_biased',
                                .default = 'multi')) %>%
  mutate(guide.call = factor(guide.call,
                             levels = c('single_true', 'single_biased',
                                        'multi', 'none_empty'))) %>%
  # clean up
  mutate(count.prop = case_match(count.prop, Inf ~ NA, .default = count.prop)) %>%
  arrange(guide.call, ug.filtered)

rep2.calls <- rep2.stats %>% ungroup() %>%
  mutate(guide.call = case_when(ug.filtered == 0 ~ 'none_empty',
                                ug.filtered == 1 ~ 'single_true',
                                count.prop < count.prop.threshold ~ 'single_biased',
                                .default = 'multi')) %>%
  mutate(guide.call = factor(guide.call,
                             levels = c('single_true', 'single_biased',
                                        'multi', 'none_empty'))) %>%
  # clean up
  mutate(count.prop = case_match(count.prop, Inf ~ NA, .default = count.prop)) %>%
  arrange(guide.call, ug.filtered)

rep3.calls <- rep3.stats %>% ungroup() %>%
  mutate(guide.call = case_when(ug.filtered == 0 ~ 'none_empty',
                                ug.filtered == 1 ~ 'single_true',
                                count.prop < count.prop.threshold ~ 'single_biased',
                                .default = 'multi')) %>%
  mutate(guide.call = factor(guide.call,
                             levels = c('single_true', 'single_biased',
                                        'multi', 'none_empty'))) %>%
  # clean up
  mutate(count.prop = case_match(count.prop, Inf ~ NA, .default = count.prop)) %>%
  arrange(guide.call, ug.filtered)

rep1.summary <- as.data.frame(table(rep1.calls$guide.call))
rep2.summary <- as.data.frame(table(rep2.calls$guide.call))
rep3.summary <- as.data.frame(table(rep3.calls$guide.call))


#Save calls
write.table(rep1.calls, file = 'GuideCalls_Rep1.tsv', sep = '\t', quote = F, col.names = T, row.names = F)
write.table(rep2.calls, file = 'GuideCalls_Rep2.tsv', sep = '\t', quote = F, col.names = T, row.names = F)
write.table(rep3.calls, file = 'GuideCalls_Rep3.tsv', sep = '\t', quote = F, col.names = T, row.names = F)

#Add Metadata
#Calls rep1 vector
call.num.1 <- rep1.calls %>%
  dplyr::select(cell, guide.call) %>%
  column_to_rownames(var = 'cell')
#Calls rep2 vector
call.num.2 <- rep2.calls %>%
  dplyr::select(cell, guide.call) %>%
  column_to_rownames(var = 'cell')
#Calls rep3 vector
call.num.3 <- rep3.calls %>%
  dplyr::select(cell, guide.call) %>%
  column_to_rownames(var = 'cell')

#Add merged calls to dataset
cts1 <- AddMetaData(cts1, call.num.1, col.name = 'num_guide')
cts2 <- AddMetaData(cts2, call.num.2, col.name = 'num_guide')
cts3 <- AddMetaData(cts3, call.num.3, col.name = 'num_guide')

#Subset for single-guide cells 
cts1 <- subset(cts1, subset = num_guide == 'single_true' | num_guide == 'single_biased')
cts2 <- subset(cts2, subset = num_guide == 'single_true' | num_guide == 'single_biased')
cts3 <- subset(cts3, subset = num_guide == 'single_true' | num_guide == 'single_biased')


ncol(cts1) + ncol(cts2) + ncol(cts3)

#Get guide calls
call.guide.1 <- rep1.calls %>%
  dplyr::select(cell, top.guide) %>%
  column_to_rownames(var = 'cell')
call.guide.2 <- rep2.calls %>%
  dplyr::select(cell, top.guide) %>%
  column_to_rownames(var = 'cell')
call.guide.3 <- rep3.calls %>%
  dplyr::select(cell, top.guide) %>%
  column_to_rownames(var = 'cell')

#Add guide metadata
cts1 <- AddMetaData(cts1, call.guide.1, col.name = 'guide')
cts2 <- AddMetaData(cts2, call.guide.2, col.name = 'guide')
cts3 <- AddMetaData(cts3, call.guide.3, col.name = 'guide')

#Format HAR Metadta 
df.hars.1 <- FetchData(cts1, vars = 'guide') %>% 
  mutate(temp = str_replace(guide, '2xHAR[.]1', 'hold')) %>%
  mutate(temp2 = str_replace(temp, "[.]1", '')) %>% #avoids regex usage of period 
  mutate(har = str_replace(temp2, 'hold', '2xHAR.1')) %>% 
  mutate(har = str_replace(har, 'NTC1', 'Non-Targeting')) %>% 
  mutate(har = str_replace(har, 'NTC2', 'Non-Targeting')) %>% 
  mutate(har = str_replace(har, 'NTC3', 'Non-Targeting')) %>% 
  mutate(har = str_replace(har, 'NTC4', 'Non-Targeting')) %>% 
  mutate(har = str_replace(har, 'NTC5', 'Non-Targeting')) %>% 
  mutate(har = str_replace(har, 'NTC-10x', 'Non-Targeting')) %>% 
  select(har)

df.hars.2 <- FetchData(cts2, vars = 'guide') %>% 
  mutate(temp = str_replace(guide, '2xHAR[.]1', 'hold')) %>%
  mutate(temp2 = str_replace(temp, "[.]1", '')) %>% #avoids regex usage of period 
  mutate(har = str_replace(temp2, 'hold', '2xHAR.1')) %>% 
  mutate(har = str_replace(har, 'NTC1', 'Non-Targeting')) %>% 
  mutate(har = str_replace(har, 'NTC2', 'Non-Targeting')) %>% 
  mutate(har = str_replace(har, 'NTC3', 'Non-Targeting')) %>% 
  mutate(har = str_replace(har, 'NTC4', 'Non-Targeting')) %>% 
  mutate(har = str_replace(har, 'NTC5', 'Non-Targeting')) %>% 
  mutate(har = str_replace(har, 'NTC-10x', 'Non-Targeting')) %>% 
  select(har)

df.hars.3 <- FetchData(cts3, vars = 'guide') %>% 
  mutate(temp = str_replace(guide, '2xHAR[.]1', 'hold')) %>%
  mutate(temp2 = str_replace(temp, "[.]1", '')) %>% #avoids regex usage of period 
  mutate(har = str_replace(temp2, 'hold', '2xHAR.1')) %>% 
  mutate(har = str_replace(har, 'NTC1', 'Non-Targeting')) %>% 
  mutate(har = str_replace(har, 'NTC2', 'Non-Targeting')) %>% 
  mutate(har = str_replace(har, 'NTC3', 'Non-Targeting')) %>% 
  mutate(har = str_replace(har, 'NTC4', 'Non-Targeting')) %>% 
  mutate(har = str_replace(har, 'NTC5', 'Non-Targeting')) %>% 
  mutate(har = str_replace(har, 'NTC-10x', 'Non-Targeting')) %>% 
  select(har)

cts1 <- AddMetaData(cts1, df.hars.1, col.name = 'har')  
cts2 <- AddMetaData(cts2, df.hars.2, col.name = 'har')  
cts3 <- AddMetaData(cts3, df.hars.3, col.name = 'har')  

saveRDS(cts1, 'batch1_singles.rds')
saveRDS(cts2, 'batch2_singles.rds')
saveRDS(cts3, 'batch3_singles.rds')

list <- c('cts1', 'cts2', 'cts3')


cts.full <- merge(cts1, cts2)
cts.full <- merge(cts.full, cts3)

saveRDS(cts.full, 'Merged_Singles.rds')
