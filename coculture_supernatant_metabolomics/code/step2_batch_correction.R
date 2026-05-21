##
# Apply batch correction to normalized (median) data set (Trapecar Lab)
# Results from supernatant of co-cultured cells from different organs
#
# By Oscar Ospina
#
# Created: Sep 07, 2025
# Modified: May 21, 2026
#

library('tidyverse')

# Read normalized data
med_dat_dedup = readRDS('./data/median_normalized_only_id_intx_dup_sum.RDS')

# Read modified metadata
meta_dat = readr::read_delim('./data/Updated_Metadata_Table_MOD.csv', name_repair='minimal')

# Use ComBat for batch effect removal and do PCA
# Apply log2 x+1
med_dat_dedup_bc = sva::ComBat(dat=log2(med_dat_dedup + 1), batch=meta_dat[['experiment']])

# Save uncorrected and corrected data
saveRDS(med_dat_dedup_bc, './data/median_normalized_only_id_intx_dup_sum_ComBat.RDS')

# Perform PCA (pre and post bacth correction)
pca_med = prcomp(scale(t(med_dat_dedup)))
pca_med_bc = prcomp(scale(t(med_dat_dedup_bc)))

# Make data frame for PCA plot
pca_data_df = data.frame(pca_med[['x']][, c('PC1', 'PC2')]) %>% 
  rownames_to_column(var='mod_sample_name') %>%
  rename(PC1_MED=PC1, PC2_MED=PC2) %>%
  left_join(., data.frame(pca_med_bc[['x']][, c('PC1', 'PC2')]) %>% 
              rownames_to_column(var='mod_sample_name') %>%
              rename(PC1_MED_ComBat=PC1, PC2_MED_ComBat=PC2), by='mod_sample_name') %>%
  filter(!(mod_sample_name %in% c("untargeted_processing_blank01",
                                  "untargeted_processing_blank02"))) %>%
  left_join(., meta_dat, by='mod_sample_name') %>%
  select(-c('sample_name')) %>%
  mutate(group=case_when(group == "High + Metformin" ~ "High/Metformin",
                         group == "High + Semaglutide" ~ "High/Semaglutide",
                         TRUE ~ group)) %>%
  mutate(group=factor(group, levels=c("Low", "Mid", "High", 
                                      "High/Metformin", "High/Semaglutide")))

p5 = ggplot(pca_data_df, aes(x=PC1_MED, y=PC2_MED)) +
  geom_point(aes(color=group), size=2) +
  labs(color='', title='Median') +
  ggrepel::geom_text_repel(aes(label=replicate)) +
  scale_color_manual(values=col_pal) +
  guides(color=guide_legend(nrow=2)) +
  theme_bw() +
  theme(legend.position='bottom')

p6 = ggplot(pca_data_df, aes(x=PC1_MED, y=PC2_MED)) +
  geom_point(aes(color=experiment), size=2) +
  labs(color='', title='Median') +
  ggrepel::geom_text_repel(aes(label=replicate)) +
  scale_color_manual(values=col_pal) +
  guides(color=guide_legend(nrow=2)) +
  theme_bw() +
  theme(legend.position='bottom')

p7 = ggplot(pca_data_df, aes(x=PC1_MED_ComBat, y=PC2_MED_ComBat)) +
  geom_point(aes(color=group), size=2) +
  labs(color='', title='Median+ComBat') +
  ggrepel::geom_text_repel(aes(label=replicate)) +
  scale_color_manual(values=col_pal) +
  guides(color=guide_legend(nrow=2)) +
  theme_bw() +
  theme(legend.position='bottom')

p8 = ggplot(pca_data_df, aes(x=PC1_MED_ComBat, y=PC2_MED_ComBat)) +
  geom_point(aes(color=experiment), size=2) +
  labs(color='', title='Median+ComBat') +
  ggrepel::geom_text_repel(aes(label=replicate)) +
  scale_color_manual(values=col_pal) +
  guides(color=guide_legend(nrow=2)) +
  theme_bw() +
  theme(legend.position='bottom')

print(ggpubr::ggarrange(p5, p6, p7, p8, ncol=4, nrow=2))

