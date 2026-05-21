##
# Differential gene expression analysis among tissue types
#
# By Oscar Ospina
#
# Created: Jun 29, 2025
# Modified: May 21, 2026
#

library('tidyverse')
library('DESeq2')
library('gridExtra')
library('ComplexHeatmap')

# Read DESeq object
deseq_obj = readRDS('./data/deseq_object_sva_included.RDS')

# Extract tissue types
ttypes = as.vector(unique(deseq_obj@colData[['tissue_type']]))

# # Subset object to "normal" samples
samples_sub = as.data.frame(deseq_obj@colData) %>% filter(culture == 'isolation' & diet_tx == 'Normal')
deseq_obj_sub = deseq_obj[, rownames(samples_sub)]

# Perform DE analysis
deg_ls = lapply(ttypes, function(i){
  # Make new variable to compare tissue against others
  deseq_obj_tmp = deseq_obj_sub
  deseq_obj_tmp@colData$tissue_grp = factor(ifelse(deseq_obj_tmp@colData$tissue_type == i, i, 'other'))
  
  # Keep genes expressed at least in half of the samples
  genes_keep = rowSums(assay(deseq_obj_tmp) == 0) <= round(ncol(assay(deseq_obj_tmp))*0.5, 0)
  deseq_obj_tmp = deseq_obj_tmp[genes_keep, ]
  
  # Set new model
  design(deseq_obj_tmp) = ~experiment + tissue_grp
  
  # Run DE analysis
  deseq_obj_tmp = DESeq(deseq_obj_tmp)
  
  # Get p-values
  res_df = as.data.frame(results(deseq_obj_tmp, contrast=c('tissue_grp', i, 'other')))
  
  res_df = res_df %>%
    rownames_to_column('gene_symbol') %>%
    arrange(padj, desc(log2FoldChange))
  
  return(res_df)
})
names(deg_ls) = ttypes

# Save results
saveRDS(deg_ls, './results/diff_expression_tests_among_tissues.RDS')

