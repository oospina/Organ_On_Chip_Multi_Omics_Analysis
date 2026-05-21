##
# Gene set enrichment analysis
#
# By Oscar Ospina
#
# Created: Jul 24, 2025
# Modified: May 21, 2026
#

library('tidyverse')
library('fgsea')

# Get hallmark pathways
fps = list.files('./data/', pattern='\\.gmt', full.names=TRUE)
pws_ls = lapply(c('kegg', 'go.bp', 'reactome'), function(i){
  pws_raw = readLines(grep(i, fps, value=TRUE))
  pws = lapply(pws_raw, function(i){
    pw_tmp = unlist(strsplit(i, split='\\t'))
    pw_name_tmp = pw_tmp[1]
    pw_genes_tmp = pw_tmp[-c(1:2)]
    return(list(pw_name=pw_name_tmp,
                pw_genes=pw_genes_tmp))
  })
  rm(pws_raw)
  pws_names = c()
  for(i in 1:length(pws)){
    pws_names = append(pws_names, pws[[i]][['pw_name']])
    pws[[i]] = pws[[i]][['pw_genes']]
  }
  names(pws) = pws_names
  
  return(pws)
})
names(pws_ls) = c('kegg', 'go.bp', 'reactome')

rm(fps) # Clean env

# Save pathways in R object
saveRDS(pws_ls, './data/gene_set_list_kegg_gobp_reactome.RDS')

# Get size of pathway
pws_sz = lapply(names(pws_ls), function(i){
  ls_tmp = lapply(pws_ls[[i]], function(j){
    return(length(j))
  })
  return(ls_tmp)
})
names(pws_sz) = c('kegg', 'go.bp', 'reactome')

# Save pathway sizes in R object
saveRDS(pws_sz, './data/gene_set_sizes_list_kegg_gobp_reactome.RDS')

# Read DE results
deg_ls = readRDS('./results/diff_expression_tests_all_compars.RDS')
compars = names(deg_ls)

# Put all DGE tests in single data frame and add rank statistic
## Reverse logFC for controls
## Also calculate weighted ranks
deg_ls = lapply(compars, function(i){
  df_tmp = deg_ls[[i]] %>%
    dplyr::select('gene_symbol', 'log2FoldChange', 'pvalue') %>%
    filter(!is.na(pvalue)) %>%
    add_column(compar=i) %>%
    separate(col=compar, into=c('tissue_type', 'compar1'), sep='_', extra='merge', remove=FALSE) %>%
    separate(col=compar1, into=c('group1', 'group2'), sep='_vs_', remove=TRUE) %>%
    mutate(rank_stat_grp1_w= -log10(pvalue) * -log2FoldChange) %>% # Revert logFC since "group1" is the baseline/ctrl
    mutate(rank_stat_grp2_w= -log10(pvalue) * log2FoldChange) %>%
    mutate(rank_stat_grp1= -1 * log2FoldChange) %>% # Revert logFC since "group1" is the baseline/ctrl
    mutate(rank_stat_grp2= log2FoldChange)
  
  return(df_tmp)
})
names(deg_ls) = compars

# Create named vectors with ranks (only logFC)
ranks_sets = lapply(compars, function(i){
  df_tmp = deg_ls[[i]]
  ls_tmp_1 = sort(setNames(df_tmp[['rank_stat_grp1']], nm=df_tmp[['gene_symbol']]), decreasing=TRUE)
  ls_tmp_2 = sort(setNames(df_tmp[['rank_stat_grp2']], nm=df_tmp[['gene_symbol']]), decreasing=TRUE)
  ls_tmp_3 = sort(setNames(df_tmp[['rank_stat_grp1_w']], nm=df_tmp[['gene_symbol']]), decreasing=TRUE)
  ls_tmp_4 = sort(setNames(df_tmp[['rank_stat_grp2_w']], nm=df_tmp[['gene_symbol']]), decreasing=TRUE)
  ls_tmp = list(ls_tmp_1, ls_tmp_2, ls_tmp_3, ls_tmp_4)
  names(ls_tmp) = c(paste0(unique(df_tmp[['tissue_type']]), '_', unique(df_tmp[['group1']])),
                    paste0(unique(df_tmp[['tissue_type']]), '_', unique(df_tmp[['group2']])),
                    paste0(unique(df_tmp[['tissue_type']]), '_', unique(df_tmp[['group1']]), '_w'),
                    paste0(unique(df_tmp[['tissue_type']]), '_', unique(df_tmp[['group2']]), '_w'))
  return(ls_tmp)
})
names(ranks_sets) = compars

# Run FGSEA
fgsea_ls = lapply(compars, function(i){
  ls_tmp1 = lapply(names(pws_ls), function(p){
    df_tmp = fgsea(pathways=pws_ls[[p]], stats=ranks_sets[[i]][[1]], minSize=3, scoreType='pos', nproc=4)
    return(df_tmp)
  })
  names(ls_tmp1) = names(pws_ls)
  
  ls_tmp2 = lapply(names(pws_ls), function(p){
    df_tmp = fgsea(pathways=pws_ls[[p]], stats=ranks_sets[[i]][[2]], minSize=3, scoreType='pos', nproc=4)
    return(df_tmp)
  })
  names(ls_tmp2) = names(pws_ls)
  
  ls_tmp3 = lapply(names(pws_ls), function(p){
    df_tmp = fgsea(pathways=pws_ls[[p]], stats=ranks_sets[[i]][[3]], minSize=3, scoreType='pos', nproc=4)
    return(df_tmp)
  })
  names(ls_tmp3) = names(pws_ls)
  
  ls_tmp4 = lapply(names(pws_ls), function(p){
    df_tmp = fgsea(pathways=pws_ls[[p]], stats=ranks_sets[[i]][[4]], minSize=3, scoreType='pos', nproc=4)
    return(df_tmp)
  })
  names(ls_tmp4) = names(pws_ls)
  
  ls_tmp = list(ls_tmp1, ls_tmp2, ls_tmp3, ls_tmp4)
  names(ls_tmp) = names(ranks_sets[[i]])
  
  return(ls_tmp)
})
names(fgsea_ls) = compars

rm(pws_ls, deg_ls, compars, ranks_sets) # Clean env

# Save scores
saveRDS(fgsea_ls, './results/fgsea_enrichment_scores_kegg_cc_bp.RDS')

