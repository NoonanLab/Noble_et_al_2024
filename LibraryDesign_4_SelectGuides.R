#Select the 2 highest-scoring guides per region whose centers do not overlap by <10bp. 

require(tidyverse)

# Load Top HARs and Filtered Guides
hars.score <- read_tsv('Outs/hars_annotation_scores.txt')
simple.guides <- read_tsv('Outs/guides_merged.txt')

hars.top500 <- hars.score %>% slice_max(Total, n=500) %>% pull(HAR)

# Get guides for top 100 HARs
get.guides <- simple.guides %>% filter(HAR %in% hars.top500) %>%
  mutate(CENTER_Hg38 = (START_Hg38 + END_Hg38)/2) # calculate midpoint

# Get the top guide for each HAR
get.guides.top1 <- get.guides %>% group_by(HAR) %>% slice_max(MIT_Human, n=1)

# Compare other guides to top guides
get.guides.next2 <- get.guides %>% 
  left_join(get.guides.top1, by='HAR', suffix=c('', '.y')) %>%
  # compute distance between midpoints, keep if distance > 10bp
  mutate(CENTER_dist = abs(CENTER_Hg38 - CENTER_Hg38.y)) %>%
  filter(CENTER_dist > 10) %>%
  # get top 2 guides per HAR
  group_by(HAR) %>% slice_max(MIT_Human, n=1)

# Output guides
output.guides <- bind_rows(get.guides.top1, get.guides.next2) %>%
  arrange(HAR, -MIT_Human) %>%
  select(-ends_with('.y')) # drop columns from output
  
  
# Filter out singles 
filter.output <- output.guides %>% 
	add_count(HAR) %>%
	filter(n!=1) %>%
	select(-n)

# Collect missing HARs
missing.hars <- data.frame(setdiff(hars.top500, filter.output$HAR)) 

# Write table
write_tsv(filter.output, 'topguides_500hars.txt')
write_tsv(missing.hars, 'dropped_500hars.txt')
