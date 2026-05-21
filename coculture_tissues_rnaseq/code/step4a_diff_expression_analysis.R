##
# Differential gene expression analysis
#
# By Oscar Ospina
#
# Created: Jun 26, 2025
# Modified: May 21, 2026
#

library('tidyverse')
library('DESeq2')
library('gridExtra')
library('ComplexHeatmap')

# Read DESeq object
deseq_obj = readRDS('./data/deseq_object_sva_included.RDS')

# Add group column to DESeq object
colData(deseq_obj) = cbind(colData(deseq_obj),
                           group=factor(paste0(colData(deseq_obj)[['culture']], '_',
                                               colData(deseq_obj)[['diet_tx']])))

# Set order of groups for comparisons
grp_order = c("isolation_Normal", "isolation_Fasting", "isolation_Western",
              "interaction_Normal", "interaction_Fasting", "interaction_Western",
              "interaction_Western_Semaglutide", "interaction_Western_Metformin")

# Set comparisons to be made
compars = as.data.frame(colData(deseq_obj)) %>%
  dplyr::select('tissue_type', 'group') %>%
  distinct() %>%
  mutate(group=factor(group, levels=c(grp_order)))
compars = expand_grid(as.vector(unique(compars[['tissue_type']])),
                      as.data.frame(t(combn(as.vector(levels(compars[['group']])), 2))))
colnames(compars) = c('tissue', 'grp1', 'grp2')
compars = compars %>%
  dplyr::filter(!(str_detect('Semaglutide$|Metformin$', grp1) & !str_detect('Semaglutide$|Metformin$', grp2))) %>%
  dplyr::filter(str_detect(grp1, '^isolation') & str_detect(grp2, '^isolation') |
                  str_detect(grp1, '^interaction') & str_detect(grp2, '^interaction') |
                  str_match(grp1, 'Fasting$|Normal$|Western$') ==
                  str_match(grp2, 'Fasting$|Normal$|Western$|Semaglutide$|Metformin$'))

# Add columns to comparison DF to store number of samples
compars[['grp1_n']] = NA
compars[['grp2_n']] = NA

ls_names = c() # List names
# Perform DE analysis
deg_ls = lapply(1:nrow(compars), function(i){
  # Identify samples in relevant comparison
  ttype = as.vector(unlist(compars[i, 1]))
  grps = as.vector(unlist(compars[i, 2:3]))
  samples_grp1 = rownames(colData(deseq_obj)[colData(deseq_obj)[['group']] == grps[1] & colData(deseq_obj)[['tissue_type']] == ttype, ])
  samples_grp2 = rownames(colData(deseq_obj)[colData(deseq_obj)[['group']] == grps[2] & colData(deseq_obj)[['tissue_type']] == ttype, ])
  
  # Sample number of samples in comparison
  compars[i, 4] <<- length(samples_grp1)
  compars[i, 5] <<- length(samples_grp2)
  
  # Save order of comparison
  grp_lvls = grp_order[grp_order %in% grps]
  ls_names <<- append(ls_names, paste0(ttype, '_', grp_lvls[1], '_vs_', grp_lvls[2]))
  
  res_df = NULL
  # MAke sure groups have enough samples to run tests
  if(length(samples_grp1) > 1 & length(samples_grp2) > 1){
    # Subset DESeq object
    deseq_obj_tmp = deseq_obj[, c(samples_grp1, samples_grp2)]
    
    # Make sure factors in DESeq object are correct
    deseq_obj_tmp$group = factor(deseq_obj_tmp$group, levels=grp_lvls)
    
    # Keep genes expressed at least in half of the samples
    genes_keep = rowSums(assay(deseq_obj_tmp) == 0) <= round(ncol(assay(deseq_obj_tmp))*0.5, 0)
    deseq_obj_tmp = deseq_obj_tmp[genes_keep, ]
    
    # Run DE analysis
    design(deseq_obj_tmp) = ~experiment + group
    deseq_obj_tmp = DESeq(deseq_obj_tmp)
    
    # Get p-values
    res_df = as.data.frame(results(deseq_obj_tmp, contrast=c('group', grps[2], grps[1])))
    
    res_df = res_df %>%
      rownames_to_column('gene_symbol') %>%
      arrange(padj, desc(log2FoldChange))
  }
  
  return(res_df)
})
names(deg_ls) = ls_names %>%
  gsub("eraction|lation", "", .) %>%
  gsub("Western_Semaglutide", "Semag", .) %>%
  gsub("Western_Metformin", "Metfo", .) %>%
  gsub("ing|al|ern", "", .) %>%
  gsub("Large_Intestine", "LargInt", .) %>%
  gsub("Skelet_Muscle", "SkMuscl", .)
deg_ls = deg_ls[!unlist(lapply(deg_ls, is.null))]

# Save results
saveRDS(deg_ls, './results/diff_expression_tests_all_compars.RDS')

