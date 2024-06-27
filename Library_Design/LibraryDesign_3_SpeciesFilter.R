#Filter guides such that human and chimpanzee sequences are identical
#Guides shall not differ by MIT specificity of >5 and Jost specificity > 0.15 across species 

require(tidyverse)

# Load filtered guides (output from flashfilter_j65.py)
human.guides <- read_tsv('hars_hg38_300bp_filtered_guides.txt')
chimp.guides <- read_tsv('hars_panTro6_300bp_filtered_guides.txt')

# Load HARs
human.hars.bed <- read_tsv('hars_hg38_300_noprom.bed', 
                           col_names = c('chrom', 'start', 'end', 'HAR'))

# Define contig (combine columns)
human.hars <- human.hars.bed %>%
  transmute(contig.x = str_c(chrom, ':', start, '-', end), HAR = HAR)

# Intersect Human and Chimp guides ---------------------------------------------
intersect.guides <- human.guides %>% inner_join(chimp.guides, by='target') %>%
  # add HAR column
  left_join(human.hars, by='contig.x') %>%
  # split contig by ':'
  rowwise() %>%
  mutate(contig.x.split = str_split(contig.x, ':'),
         contig.y.split = str_split(contig.y, ':')) %>%
  # split again by '-'
  mutate(coord.x = str_split(contig.x.split[2], '-'),
         coord.y = str_split(contig.y.split[2], '-')) %>%
  # output columns
  transmute(SEQUENCE = target,
            CHROM_Hg38 = contig.x.split[1],
            START_Hg38 = as.numeric(coord.x[1]) + start.x,
            END_Hg38   = as.numeric(coord.x[1]) + stop.x,
            CHROM_panTro6 = contig.y.split[1],
            START_panTro6 = as.numeric(coord.y[1]) + start.y,
            END_panTro6   = as.numeric(coord.y[1]) + stop.y,
            HAR = HAR,
            MIT_Human = Hsu2013.x,
            MIT_Chimp = Hsu2013.y,
            Jost_maxOT_Human = JostCRISPRi_maxOT.x,
            Jost_maxOT_Chimp = JostCRISPRi_maxOT.y,
            Jost_guide_Human = JostCRISPRi_specificityscore.x,
            Jost_guide_Chimp = JostCRISPRi_specificityscore.y) %>%
  # rearrange columns
  relocate(-starts_with(c('MIT', 'Jost'))) %>%
  arrange(SEQUENCE)


# Filter guides ----------------------------------------------------------------
filter.guides <- intersect.guides %>% 
  # calculate difference between species scores
  mutate(MIT_diff = MIT_Human - MIT_Chimp,
         Jost_diff = Jost_guide_Human - Jost_guide_Chimp) %>%
  # species-efficacy-disparity
  filter(abs(MIT_diff) < 5 & abs(Jost_diff) < 0.15) %>% 
  filter(Jost_guide_Human >= 0.2 & Jost_guide_Chimp >= 0.2) %>%
  # specificity > off-target
  filter(Jost_guide_Human > Jost_maxOT_Human) %>%  
  filter(Jost_guide_Chimp > Jost_maxOT_Chimp)


# Output filtered guides -------------------------------------------------------
output.guides <- filter.guides %>% select(-ends_with('diff')) # drop columns from output
output.guides.simple <- output.guides %>% select(-starts_with('Jost')) # drop columns from output

write_tsv(output.guides, 'guides_300bp_j65.txt')
write_tsv(output.guides.simple, 'guides_300bp_full.txt')
