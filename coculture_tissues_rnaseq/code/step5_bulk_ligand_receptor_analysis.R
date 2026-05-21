##
# Ligand-receptor inference from bulk RNA seq
#
# By Oscar Ospina
#
# Created: Dec 22, 2025
# Modified: May 21, 2025
#

library('tidyverse')
library('DESeq2')
library('BulkSignalR')
library('ComplexHeatmap')

# Set parallelization
cl = parallel::makeCluster(8)
parallel::clusterSetRNGStream(cl, iseed=12345)
doParallel::registerDoParallel(cl)

# Read DESeq object
deseq_obj = readRDS('./data/deseq_object_sva_included.RDS')

# Subset object to "normal" samples
samples_keep = as.data.frame(deseq_obj@colData) %>% filter(diet_tx == 'Normal')
deseq_obj_sub = deseq_obj[, rownames(samples_keep)]
rm(deseq_obj, samples_keep) # Clean env

# Extract tissue types
cult = as.vector(unique(deseq_obj_sub$culture))

datamod_ls = lapply(cult, function(i){
  # Subset DESeq to relevant tissue
  samples_sub = rownames(colData(deseq_obj_sub))[colData(deseq_obj_sub)$culture == i]
  
  # Prepare expression matrix
  deseq_tmp = deseq_obj_sub[, samples_sub]
  expr_mtx = vst(deseq_tmp, blind=TRUE)
  expr_mtx = limma::removeBatchEffect(assay(expr_mtx), batch=expr_mtx$experiment)
  expr_mtx = expr_mtx[apply(expr_mtx, 1, sd) != 0, ] # learnParameters crashes if no variation
  
  # Create BulkSignalR model
  data_model = BSRDataModel(counts=expr_mtx,
                            species="hsapiens",
                            normalize=FALSE,
                            log.transformed=TRUE,
                            method='vst')
  
  return(data_model)
  
})
names(datamod_ls) = cult

# Make LR inferences
bsr_ls = lapply(names(datamod_ls), function(i){
  # Update parameters
  bsr_tmp = learnParameters(datamod_ls[[i]],
                            verbose=TRUE,
                            min.positive=5,
                            quick=FALSE,
                            plot.folder="./results/",
                            filename=paste0("bulksignalr_parameters_", i))
  
  # Infer L-R interactions
  bsr_inf = BSRInference(bsr_tmp,
                         min.cor=0.3,
                         min.positive=5,
                         fdr.proc='BH',
                         reference="REACTOME-GOBP")
  
  return(bsr_inf)
})
names(bsr_ls) = cult

saveRDS(bsr_ls, './results/bulksignalr_inferences.RDS')

# Extract and save LR interactions
lr_ls = lapply(names(bsr_ls), function(i){
  lr_df = LRinter(bsr_ls[[i]])
  lr_df = lr_df %>% arrange(qval)
  
  return(lr_df)
})
names(lr_ls) = names(bsr_ls)

# AGGREGATE pathways
aggr_lr_ls = lapply(names(bsr_ls), function(i){
  bsr_inf_aggr = reduceToPathway(bsr_ls[[i]])
  #bsr_inf_aggr = reduceToBestPathway(bsr_inf_aggr)
  return(bsr_inf_aggr)
})
names(aggr_lr_ls) = names(bsr_ls)

# Save AGGREGATED pathways to file
lapply(names(bsr_ls), function(i){
  lr_aggr_df = LRinter(aggr_lr_ls[[i]]) 
  lr_aggr_df = lr_aggr_df %>% arrange(qval)
  
  # Save AGGREGATED ligand receptor interactions to file
  openxlsx::write.xlsx(lr_aggr_df, paste0("./results/enriched_LR_reactome_gobp_aggregated_", i, '.xlsx'))
})

# Get L-R scores per sample
scores_ls = lapply(names(aggr_lr_ls), function(i){
  #  bsr_inf_aggr_sign = reduceToPathway(bsr_ls[[i]])
  bsr_inf_aggr_sign = BSRSignature(aggr_lr_ls[[i]], qval.thres=0.01)
  
  scores_lr = scoreLRGeneSignatures(datamod_ls[[i]], 
                                    bsr_inf_aggr_sign,
                                    name.by.pathway=FALSE, 
                                    rownames.LRP=TRUE,
                                    abs.z.score=FALSE)
  
  return(scores_lr)
})
names(scores_ls) = names(aggr_lr_ls)

# Make heatmaps
## Make color palette
col_pal = c(exp_2="orange", exp_3="#C71585")

# Make heatmaps of L-R interaction scores by tissue
## Extract matrices 
mtx_ls = lapply(as.vector(unique(colData(deseq_obj_sub)$tissue_type)), function(i){
  samples_sub = rownames(colData(deseq_obj_sub))[colData(deseq_obj_sub)$tissue_type == i]
  lr_mtx = as.data.frame(scores_ls[['isolation']][, grep('ISO_', samples_sub, value=TRUE)]) %>%
    rownames_to_column('gene_sign') %>%
    full_join(., as.data.frame(scores_ls[['interaction']][, grep('INT_', samples_sub, value=TRUE)]) %>%
                rownames_to_column('gene_sign'), by='gene_sign') %>%
    column_to_rownames('gene_sign') %>%
    as.matrix() %>%
    t() %>% scale %>% t()
  
  lr_mtx = as.data.frame(lr_mtx) %>% 
    rownames_to_column('lr_signature')
  
  return(lr_mtx)
}) 
names(mtx_ls) = as.vector(unique(colData(deseq_obj_sub)$tissue_type))

# Save scores by tissue
openxlsx::write.xlsx(mtx_ls, file=paste0('./results/tissue_level_LR_scores_reactome_gobp_aggregated.xlsx'))

# Stop parallelization env
parallel::stopCluster(cl)

