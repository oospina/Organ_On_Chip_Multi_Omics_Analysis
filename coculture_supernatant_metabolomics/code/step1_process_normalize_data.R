##
# Process metabolomics data set (Trapecar Lab)
# Results from supernatant of co-cultured cells from different organs
#
# By Oscar Ospina
#
# Created: Sep 07, 2025
# Modified: May 21, 2026
#

library('tidyverse')

# Specify file paths
pos_fp = './data/Export_POS_Trapecar_Untargeted_Metabolomics.csv'
neg_fp = './data/Export_NEG_Trapecar_Untargeted_Metabolomics.csv'
meta_fp = './data/Updated_Metadata_Table.csv'

# Read sample meta data
# Remove '_POS' token from sample names for easier match with both POS/NEG
# Also remove "ProcessingBlank*" samples
meta_dat = read_delim(meta_fp, delim=',') %>%
  mutate(mod_sample_name=gsub('_POS', '', sample_name), .after=1) %>%
  filter(!grepl('ProcessingBlank', mod_sample_name))

# Read intensity data 
# Clean sample names to make it easier to match with sample metadata 
pos_dat = data.table::fread(pos_fp, data.table=FALSE)
colnames(pos_dat) = gsub("_POS|\\.mzXML\\ Peak\\ height", '', colnames(pos_dat))
neg_dat = data.table::fread(neg_fp, data.table=FALSE)
colnames(neg_dat) = gsub("_NEG|\\.mzXML\\ Peak\\ height", '', colnames(neg_dat))

rm(pos_fp, neg_fp, meta_fp) # Clean env

# Find columns in intx data that are not sample names (or missing metadata)
no_intx_cols_pos = grep(paste0(meta_dat[[2]], collapse='|'), colnames(pos_dat), value=TRUE, invert=TRUE)
no_intx_cols_neg = grep(paste0(meta_dat[[2]], collapse='|'), colnames(neg_dat), value=TRUE, invert=TRUE)

# Keep metabolite annotation data separated from intensity data
pos_annot_df = pos_dat[, no_intx_cols_pos] %>% 
  janitor::clean_names() %>% mutate(row_id=paste0('row_', row_id))
neg_annot_df = neg_dat[, no_intx_cols_neg] %>% janitor::clean_names() %>% 
  janitor::clean_names() %>% mutate(row_id=paste0('row_', row_id))

# Remove annotation columns from intensity data
pos_dat = pos_dat %>% select(-all_of(no_intx_cols_pos))
neg_dat = neg_dat %>% select(-all_of(no_intx_cols_neg))

rm(no_intx_cols_pos, no_intx_cols_neg) # Clean env

# Match order of samples in intensity data with metadata
pos_dat = pos_dat[, meta_dat[['mod_sample_name']] ]
neg_dat = neg_dat[, meta_dat[['mod_sample_name']] ]

if(identical(colnames(pos_dat), meta_dat[['mod_sample_name']]) &
   identical(colnames(neg_dat), meta_dat[['mod_sample_name']])){
  # Clean sample names
  pos_dat = pos_dat %>% janitor::clean_names()
  neg_dat = neg_dat %>% janitor::clean_names()
  # Add clean sample names to metadata
  meta_dat[['mod_sample_name']] = colnames(pos_dat)
} else{
  stop('SAMPLES ARE NOT IN THE SAME ORDER!!!')
}

# Add rownames to intensity data
rownames(pos_dat) = pos_annot_df[['row_id']]
rownames(neg_dat) = neg_annot_df[['row_id']]

# Perform median normalization
pos_dat_med = sweep(pos_dat, 2, apply(pos_dat, 2, median), FUN = "/") * median(apply(pos_dat, 2, median))
neg_dat_med = sweep(neg_dat, 2, apply(neg_dat, 2, median), FUN = "/") * median(apply(neg_dat, 2, median))

# Filter out peaks detected in few (<= 50%) samples
pos_dat_med = pos_dat_med[rowSums(pos_dat_med == 0) < round(ncol(pos_dat_med) * 0.5, 0), ]
neg_dat_med = neg_dat_med[rowSums(neg_dat_med == 0) < round(ncol(neg_dat_med) * 0.5, 0), ]

# Put all normalized data together
# First, modify row names to identify channel
rownames(pos_dat_med) = paste0(rownames(pos_dat_med), '_POS')
rownames(neg_dat_med) = paste0(rownames(neg_dat_med), '_NEG')

pos_annot_df = pos_annot_df %>% mutate(row_id_mod=paste0(row_id, '_POS'), .after=1)
neg_annot_df = neg_annot_df %>% mutate(row_id_mod=paste0(row_id, '_NEG'), .after=1)

# Then, combine data
med_dat = rbind(pos_dat_med, neg_dat_med)
annot_df = rbind(pos_annot_df, neg_annot_df)

rm(pos_dat, neg_dat) # Clean env

# Save metabolite annotation data
saveRDS(annot_df, './data/metabolite_annotation_data.RDS')

rm(pos_dat_med, neg_dat_med, pos_annot_df, neg_annot_df) # Clean env

# Filter out peaks without annotation
annot_df_ident = annot_df %>% filter(row_identity_all_i_ds != '')
med_dat_ident = med_dat[annot_df_ident[['row_id_mod']], ]

rm(med_dat, annot_df) # Clean env

# Average duplicate annotated peaks
unique_mets = unique(annot_df_ident[['row_identity_all_i_ds']])

## THEN WITH MEDIAN NORMALIZATION
med_ls = lapply(unique_mets, function(i){
  # Get all rows with the same annotation
  row_ids = annot_df_ident[["row_id_mod"]][annot_df_ident[['row_identity_all_i_ds']] == i]
  rows = med_dat_ident[row_ids, , drop=FALSE]
  # Get the mean intensity for each sample
  mean_intx = as.data.frame(t(as.matrix(colMeans(rows))))
  rownames(mean_intx) = paste0(row_ids, collapse='_')
  
  return(mean_intx)
})
med_dat_dedup = do.call(rbind, med_ls)

rm(med_dat_ident, med_ls, unique_mets) # Clean env

# Save normalized data
saveRDS(med_dat_dedup, './data/median_normalized_only_id_intx_dup_sum.RDS')

# Extract variables from sample names and further clean metadata
meta_dat = meta_dat %>%
  mutate(clean_sample_name=gsub('untargeted_', '', mod_sample_name) %>%
           gsub('int_', '', .), .after='mod_sample_name') %>%
  select(-c('raw_type', 'mode'))

# Save modified metadata
write.csv(meta_dat, './data/Updated_Metadata_Table_MOD.csv', row.names=FALSE, quote=TRUE)

