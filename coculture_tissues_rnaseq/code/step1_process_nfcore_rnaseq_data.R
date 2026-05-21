##
# Load STAR-Salmon-generated counts and preprocess for DGE analysis
#
# By Oscar Ospina
#
# Created: Jun 26, 2025
# Modified: May 21, 2026
#

library('tidyverse')
library('DESeq2')

# Read transcript to gene name table (from nfcore/rnaseq/salmon)
fp = './nfcore_rnaseq_run/nfcore_rnaseq_output/star_salmon/tx2gene.tsv'
tx2gene = readr::read_delim(fp, delim='\t')
# Collapse ambiguous transcript IDs to single gene symbol
tx2gene = tx2gene %>% mutate(gene_id=str_replace(gene_id, "_[0-9]{1,2}$", ''))
tx2gene = tx2gene[, c(1:2)]

# Pass salmon quant files to tximport
fp = list.files('./nfcore_rnaseq_run/nfcore_rnaseq_output/star_salmon/', recursive=TRUE, pattern='quant\\.sf', full.names=TRUE)
names(fp) = stringr::str_extract(fp, '[INTISO]{3}_[A-Z]_[0-9A-Z]{1,2}_[EX23]{3}')
txi = tximport::tximport(fp, type='salmon', tx2gene=tx2gene)
rm(fp, tx2gene) # Clean env

# Save tximport object
saveRDS(txi, './data/imported_salmon_counts_tximport.RDS')

# Read and process sample meta data
fp = './data/RNA_SampleList.xlsx'
sample_meta = readxl::read_excel(fp, .name_repair='minimal') %>% 
  janitor::clean_names() %>%
  mutate(sample_name=str_replace(sample_name, "EXP", "EX")) %>%
  mutate(tissue_type=factor(str_replace(tissue_type, " ", "_"))) %>%
  mutate(diet_tx=factor(case_when(drug != 'None' ~ paste0(diet, '_', drug), TRUE ~ diet))) %>%
  mutate(culture=ifelse(str_detect(sample_name, 'ISO'), 'isolation', 'interaction')) %>%
  mutate(experiment=factor(paste0('exp_', experiment))) %>% 
  mutate(tissue_origin=factor(tissue_origin)) %>%
  mutate(number_of_samples_in_pool=factor(number_of_samples_in_pool)) %>%
  column_to_rownames('sample_name')
rm(fp) # Clean env

# Make sure sample names order match in sample meta and counts
sample_meta = sample_meta[colnames(txi$counts), ]

# Create initial DESeq object (no batch correction)
deseq_obj = DESeqDataSetFromTximport(txi=txi,
                                     colData=sample_meta,
                                     design=~diet_tx)

# Keep genes expressed at least in half of the samples
genes_keep = rowSums(assay(deseq_obj) == 0) <= round(ncol(assay(deseq_obj))*0.5, 0)
deseq_obj = deseq_obj[genes_keep, ]
rm(genes_keep, txi) # Clean env

# Calculate library size factors
deseq_obj = estimateSizeFactors(deseq_obj)

# Estimate SVs
sva_res = sva::svaseq(dat=assay(vst(deseq_obj, blind=TRUE)),
                      mod=model.matrix(~diet_tx, data=sample_meta),
                      mod0=model.matrix(~1, data=sample_meta))

# Extract SVs to add to DESeq object
sva_mtx = sva_res[['sv']]
colnames(sva_mtx) = paste0('SV', 1:ncol(sva_mtx))
colData(deseq_obj) = S4Vectors::DataFrame(cbind(colData(deseq_obj), sva_mtx))

# Plot SVs
sv_df = as.data.frame(colData(deseq_obj)) %>% 
  pivot_longer(cols=c("experiment", "culture", "tissue_type", 
                      "tissue_origin", "number_of_samples_in_pool", "diet_tx"),
               values_to='value', names_to='variable') %>%
  split(., f=.$variable)

sv_p1 = ggplot(sv_df$experiment) +
  geom_point(aes(x=SV1, y=SV2, color=value), size=2) +
  guides(color=guide_legend(override.aes=list(size=3))) +
  facet_grid(~variable)

sv_p2 = sv_p1 %+% sv_df$tissue_type
sv_p3 = sv_p1 %+% sv_df$diet_tx
sv_p4 = sv_p1 %+% sv_df$culture
sv_p5 = sv_p1 %+% sv_df$tissue_origin
sv_p6 = sv_p1 %+% sv_df$number_of_samples_in_pool

gridExtra::grid.arrange(sv_p1, sv_p2, sv_p3, sv_p4, sv_p5, sv_p6)

rm(sv_p1, sv_p2, sv_p3, sv_p4, sv_p5, sv_p6, sva_res, sva_mtx, sv_df) # Clean env

# Save DESeq object
saveRDS(deseq_obj, './data/deseq_object_sva_included.RDS')

