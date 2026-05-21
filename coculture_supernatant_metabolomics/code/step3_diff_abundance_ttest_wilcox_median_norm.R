##
# Differential abundance tests using T-tests and Wilcoxon tests
# USING MEDIAN NORMALIZED DATA
#
# By Oscar Ospina
#
# Created: Sep 07, 2025
# Modified: May 21, 2026
#

library('tidyverse')
library('ComplexHeatmap')

# Read metadata
sample_meta = read_delim('./data/Updated_Metadata_Table_MOD.csv') %>%
  filter(group != 'Media Control')

# Read batch corrected intensity data
intx_mtx = readRDS('./data/median_normalized_only_id_intx_dup_sum_ComBat.RDS')

# Read metabolite annotation data
annot_df = readRDS('./data/metabolite_annotation_data.RDS')

# Crate table with all comparisons to be made
compars = data.frame()
excl = c()
for(i in unique(sample_meta[['group']])){
  excl = c(excl, i)
  compars = rbind(compars, 
                  expand.grid(group1=i, group2=unique(sample_meta[['group']][!(sample_meta[['group']] %in% excl)])))
}

rm(excl) # Clean env

# Perform DA tests (all combinations)
res_ls = lapply(1:nrow(compars), function(j){
  # Select samples for each group
  grp_1 = compars[j, 'group1']
  grp_2 = compars[j, 'group2']
  grp_1 = sample_meta[['mod_sample_name']][sample_meta[['group']] == grp_1]
  grp_2 = sample_meta[['mod_sample_name']][sample_meta[['group']] == grp_2]
  
  # Test each peak
  df_tmp = data.frame()
  for(i in rownames(intx_mtx)){
    intx1 = as.vector(intx_mtx[i, grp_1])
    intx2 = as.vector(intx_mtx[i, grp_2])
    
    ttest = t.test(intx1, intx2)
    wilcox = wilcox.test(intx1, intx2)
    
    metab_id = unlist(ifelse(grepl('row_[_0-9POSNEG]+_row', i), 
                             str_split(i, pattern='_(?=row)'),
                             i))
    metab_id = paste0(unique(annot_df[['row_identity_all_i_ds']][annot_df[['row_id_mod']] %in% metab_id]), collapse='; ')
    
    df_tmp = bind_rows(df_tmp,
                       data.frame(
                         mol_id=i,
                         mol_annot=metab_id,
                         log2FC=mean(intx1, na.rm=TRUE) - mean(intx2, na.rm=TRUE),
                         t_stat=as.vector(ttest[['statistic']]),
                         t_pval=as.vector(ttest[['p.value']]),
                         wilcox_stat=as.vector(wilcox[['statistic']]),
                         wilcox_pval=as.vector(wilcox[['p.value']]),
                         compar=paste0(compars[[1]][j], '_vs_', compars[[2]][j])
                       )
    )
  }
  df_tmp = df_tmp %>%
    mutate(t_pval_adj=p.adjust(t_pval, method='BH'), .after='t_pval') %>%
    mutate(wilcox_pval_adj=p.adjust(wilcox_pval, method='BH'), .after='wilcox_pval')
  
  return(df_tmp)
})

# Put all results in single sample
res_df = do.call(bind_rows, res_ls) %>%
  arrange(compar, wilcox_pval_adj, abs(log2FC))

# Save results to file
write.csv(res_df, './results/diff_abundance_ttest_wilcox_median_norm.csv', quote=TRUE, row.names=FALSE)

# Save to Excel by comparison
res_df_split = split(res_df, res_df$compar)
names(res_df_split) = gsub('[ \\+]+', '', names(res_df_split)) %>% gsub('ormin|utide', '', .)
openxlsx::write.xlsx(res_df_split, './results/diff_abundance_ttest_wilcox_by_compar_median_norm.xlsx')

