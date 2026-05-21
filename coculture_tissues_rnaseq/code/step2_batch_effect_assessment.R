##
# Assess batch effects on data
#
# By Oscar Ospina
#
# Created: Jun 26, 2025
# Modified: May 21, 2026
#

library('tidyverse')
library('DESeq2')
library('gridExtra')

# Read DESeq object
deseq_obj = readRDS('./data/deseq_object_sva_included.RDS')

# Make sample sets to look at
sample_sets = list(all=rownames(colData(deseq_obj)),
                   adipose=rownames(as.data.frame(colData(deseq_obj)) %>% filter(tissue_type == 'Adipose')),
                   brain=rownames(as.data.frame(colData(deseq_obj)) %>% filter(tissue_type == 'Brain')),
                   intestine=rownames(as.data.frame(colData(deseq_obj)) %>% filter(tissue_type == 'Large_Intestine')), 
                   pancreas=rownames(as.data.frame(colData(deseq_obj)) %>% filter(tissue_type == 'Pancreas')),
                   liver=rownames(as.data.frame(colData(deseq_obj)) %>% filter(tissue_type == 'Liver')),
                   muscle=rownames(as.data.frame(colData(deseq_obj)) %>% filter(tissue_type == 'Skeletal_Muscle')))

# Make PCA plots to assess effect of batch term inclusion
pca_ls = lapply(names(sample_sets), function(i){
  # Extract samples
  deseq_obj_tmp = deseq_obj[, sample_sets[[i]]]
  
  # Identify high variance genes from log2(VST) counts
  high_var_genes = apply(log2(assay(vst(deseq_obj_tmp, blind=FALSE))), 1, sd) %>%
    sort(decreasing=TRUE) %>%
    head(5000)
  
  # Calculate PCs to assess correction effect
  ## Without correction
  counts_noncorr = assay(vst(deseq_obj_tmp, blind=TRUE))
  pca_res = prcomp(t(counts_noncorr[names(high_var_genes), ]), center=TRUE, scale=TRUE)
  expl_var = round(as.vector(summary(pca_res)[['importance']][2, 1:2])*100, 2)
  ## With correction
  counts_corr = vst(deseq_obj_tmp, blind=FALSE)
  counts_corr = limma::removeBatchEffect(assay(counts_corr), batch=counts_corr$experiment)
  pca_res_corr = prcomp(t(counts_corr[names(high_var_genes), ]), center=TRUE, scale=TRUE)
  expl_var_corr = round(as.vector(summary(pca_res_corr)[['importance']][2, 1:2])*100, 2)
  
  # Plot PCs before and after correction
  pca_df = as.data.frame(pca_res[['x']][, 1:2]) %>%
    rownames_to_column('sample_name') %>%
    mutate(sva='Uncorrected') %>%
    bind_rows(., as.data.frame(pca_res_corr[['x']][, 1:2]) %>%
                rownames_to_column('sample_name') %>%
                mutate(sva='Corrected')) %>%
    mutate(sva=factor(sva, levels=c('Uncorrected', 'Corrected'))) %>%
    left_join(., as.data.frame(colData(deseq_obj_tmp)) %>% 
                rownames_to_column('sample_name'), by='sample_name') %>% 
    pivot_longer(cols=c("experiment", "culture", "tissue_origin", 
                        "number_of_samples_in_pool", "diet_tx"),
                 values_to='value', names_to='variable') %>%
    mutate(variable=paste0(sva, '_', variable)) %>%
    split(., f=.$variable)
  
  pca_p = list()
  pca_p[[1]] = ggplot(pca_df$Uncorrected_diet_tx) +
    geom_point(aes(x=PC1, y=PC2, color=value, shape=tissue_type), size=2) +
    guides(color=guide_legend(override.aes=list(size=3))) +
    labs(x=paste0('PC1 (', expl_var[1], ')'), y=paste0('PC2 (', expl_var[2], ')')) +
    theme_bw() +
    facet_grid(~variable)
  
  pca_p[[2]] = pca_p[[1]] %+% pca_df$Corrected_diet_tx + labs(x=paste0('PC1 (', expl_var_corr[1], ')'), y=paste0('PC2 (', expl_var_corr[2], ')'))
  pca_p[[3]] = pca_p[[1]] %+% pca_df$Uncorrected_culture + labs(x=paste0('PC1 (', expl_var[1], ')'), y=paste0('PC2 (', expl_var[2], ')'))
  pca_p[[4]] = pca_p[[1]] %+% pca_df$Corrected_culture + labs(x=paste0('PC1 (', expl_var_corr[1], ')'), y=paste0('PC2 (', expl_var_corr[2], ')'))
  pca_p[[5]] = pca_p[[1]] %+% pca_df$Uncorrected_experiment + labs(x=paste0('PC1 (', expl_var[1], ')'), y=paste0('PC2 (', expl_var[2], ')'))
  pca_p[[6]] = pca_p[[1]] %+% pca_df$Corrected_experiment + labs(x=paste0('PC1 (', expl_var_corr[1], ')'), y=paste0('PC2 (', expl_var_corr[2], ')'))
  pca_p[[7]] = pca_p[[1]] %+% pca_df$Uncorrected_number_of_samples_in_pool + labs(x=paste0('PC1 (', expl_var[1], ')'), y=paste0('PC2 (', expl_var[2], ')'))
  pca_p[[8]] = pca_p[[1]] %+% pca_df$Corrected_number_of_samples_in_pool + labs(x=paste0('PC1 (', expl_var_corr[1], ')'), y=paste0('PC2 (', expl_var_corr[2], ')'))
  pca_p[[9]] = pca_p[[1]] %+% pca_df$Uncorrected_tissue_origin + labs(x=paste0('PC1 (', expl_var[1], ')'), y=paste0('PC2 (', expl_var[2], ')'))
  pca_p[[10]] = pca_p[[1]] %+% pca_df$Corrected_tissue_origin + labs(x=paste0('PC1 (', expl_var_corr[1], ')'), y=paste0('PC2 (', expl_var_corr[2], ')'))
  
  return(pca_p)
})

