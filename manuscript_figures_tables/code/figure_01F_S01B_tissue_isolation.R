##
# Code to generate Figures 1F and S1B of the manuscript
#
# By Oscar Ospina
#
# Created: Aug 18, 2025
# Modified: May 20, 2026
#

library('DESeq2')
library('tidyverse')
library('ComplexHeatmap')

# Load color palette
col_pal = readRDS('./data/color_palettes.RDS')

# Read DESeq object
deseq_obj = readRDS('../coculture_tissues_rnaseq/data/deseq_object_sva_included.RDS')

# Read DEGs among "mid" diet tissues in isolation
fp = '../coculture_tissues_rnaseq/results/diff_expression_tests_among_tissues.RDS'
iso_degs = readRDS(fp)
rm(fp) # Clean env

## Make annotation table
ann_df = as.data.frame(colData(deseq_obj)) %>% 
  select(c('tissue_type', 'culture')) %>%
  mutate(tissue_type=case_when(tissue_type == 'Skeletal_Muscle' ~ 'Muscle',
                               tissue_type == 'Large_Intestine' ~ 'Colon',
                               TRUE ~ tissue_type)) %>%
  arrange(tissue_type, culture)

# Make sure correct model is specified
design(deseq_obj) = ~experiment + tissue_type

# Subset DEseq object to relevant samples (Mid; ISO)
samples_iso = grep('^ISO_[A-Z]_C[0-9]+', colnames(deseq_obj), value = TRUE)
deseq_iso_obj = deseq_obj[, samples_iso]

rm(samples_iso) # Clean env

# DE genes to highlight
## Selected with help of GPT 4o as well as co-authors selection from list of DEGs
genes_tissue_high = c(
  # Adipose
  "PLIN1", "APOC1", "CFD", "FZD4", "LIPE",
  # Brain
  "GFAP", "NEFL", "SYP", "TTYH1", "PTPRZ1",
  # Large Intestine
  "CDX1", "MUC12", "MUC6", "FOXL1", "SPDEF",
  # Liver
  "TTR", "APOE", "TM6SF2", "PNPLA3", 'ALB', 'TTR', 
  'APOA1', 'APOA2', 'HNF4A', 'CYP34A', 'CYP2E1', 'ASGR1', 'SERPINA1', 'KRT18',
  # Pancreas
  "CPA1", "CEL", "PCSK1", "CHGA", "PCSK1N",
  # Skeletal Muscle
  "MYH3", "MYH7", "TTN", "ACTN2", "ENO3"
)

# Identify top 100 DEGs within each tissue (by FDR and logFC)
degs_iso_df = lapply(names(iso_degs), function(i){
  df_tmp = iso_degs[[i]] %>% add_column(tissue=i)
}) %>% bind_rows() %>%
  filter(padj < 0.05 & log2FoldChange > 1) %>%
  filter(!str_detect(gene_symbol, '^LOC[0-9]+')) %>%
  group_by(tissue) %>%
  arrange(padj, desc(log2FoldChange)) %>%
  ungroup()

# Subset metadata
ann_iso_df = ann_df %>% filter(rownames(.) %in% colnames(deseq_iso_obj))

# Get normalized counts and scale
hm_iso_mtx = vst(deseq_iso_obj, blind=FALSE)
hm_iso_mtx = limma::removeBatchEffect(assay(hm_iso_mtx), batch=deseq_iso_obj$experiment)
hm_iso_mtx = hm_iso_mtx[unique(degs_iso_df[['gene_symbol']]), rownames(ann_iso_df)]
hm_iso_mtx_scl = scale(t(hm_iso_mtx)) # Scale by row (gene)


################################################################################
# FIGURE_S01_B

# Calculate PCs
pca_iso_obj = prcomp(hm_iso_mtx_scl)
var_iso_exp = as.vector(summary(pca_iso_obj)[['importance']][2, c(1,2)])
pca_iso_df = as.data.frame(pca_iso_obj[['x']][, c(1:2)]) %>%
  rownames_to_column('samplename') %>%
  left_join(., ann_iso_df %>% rownames_to_column('samplename'), by='samplename')

pca_iso_p = ggplot(pca_iso_df, aes(x=PC1, y=PC2)) +
  geom_point(aes(color=tissue_type), shape=16, size=5) +
  labs(x=paste0('PC1 (', round(var_iso_exp[1]*100, 1), '%)'),
       y=paste0('PC2 (', round(var_iso_exp[2]*100, 1), '%)')) +
  ggtitle("PCA of DEGs (log(FC >1), FDR < 0.05) among tissues\n(Control diet; only ISO)") +
  scale_color_manual(values=col_pal) +
  guides(color=guide_legend(title='Tissue type')) +
  theme(legend.key=element_blank(), legend.text=element_text(size=15), legend.title=element_text(size=15),
        axis.text=element_text(size=15), axis.title=element_text(size=15),
        title=element_text(size=12),
        axis.line.x=element_line(linewidth=0), axis.line.y=element_line(linewidth=0),
        panel.background=element_rect(fill=NULL, color='black'))
# graphics.off()
# pdf('./results/FIGURE_S01_B.pdf', height=5, width=5)
# print(pca_iso_p)
# dev.off()
rm(pca_iso_obj, var_iso_exp, pca_iso_df, pca_iso_p, hm_iso_mtx, ann_df, degs_iso_df, deseq_iso_obj, iso_degs) # Clean env


################################################################################
# FIGURE_01_F

# Create modified legend
lgd_iso_list = list(
  Legend(labels=unique(ann_iso_df[['tissue_type']]), 
         legend_gp=gpar(fill=col_pal[as.vector(unique(ann_iso_df[['tissue_type']]))]))
)

# Combine legends
combined_iso_legend = packLegend(lgd_iso_list[[1]], direction="vertical")

# Transpose matrix and match order of columns with meta data
hm_iso_mtx_scl = t(hm_iso_mtx_scl) # Scale by row (gene)
hm_iso_mtx_scl = hm_iso_mtx_scl[, match(rownames(ann_iso_df), colnames(hm_iso_mtx_scl)), drop=FALSE] 

# Generate heatmap annotation object
hm_iso_ann = HeatmapAnnotation(df=ann_iso_df[, !grepl('culture', colnames(ann_iso_df)), drop=FALSE], 
                                  show_legend=FALSE,
                                  col=list(tissue_type=col_pal[as.vector(unique(ann_iso_df[['tissue_type']]))])
)

# Make annotation to highlight gene names
row_iso_ann = rowAnnotation(foo=anno_mark(at=as.vector(na.omit(match(genes_tissue_high, rownames(hm_iso_mtx_scl)))),  
                                          labels=genes_tissue_high[genes_tissue_high %in% rownames(hm_iso_mtx_scl)],
                                          labels_gp=gpar(fontface="italic")))

# Generate plot
hm_iso_p = Heatmap(hm_iso_mtx_scl, 
                   name='Scaled gene\nexpression',
                   column_title='Top DEGs among tissues\ncontrol diet; Only ISO\n(batch corr.; log(FC) >1, FDR < 0.05)',
                   cluster_rows=TRUE, 
                   cluster_columns=TRUE, 
                   show_row_names=FALSE, 
                   show_column_names=FALSE, 
                   show_row_dend=FALSE,
                   col=circlize::colorRamp2(c(-4.5, -1.5, 0, 1.5, 4.5), colors=c('darkblue', 'blue', 'white', 'red', 'darkred')),
                   bottom_annotation=hm_iso_ann,
                   right_annotation=row_iso_ann)

hm_iso_p = draw(hm_iso_p, annotation_legend_list=NULL, heatmap_legend_list=combined_iso_legend,
                heatmap_legend_side="right", annotation_legend_side="right", padding=unit(c(2, 2, 2, 20), "mm"))

# graphics.off()
# pdf('./results/FIGURE_01_F_HEATMAP.pdf', width=5)
# print(hm_iso_p)
# dev.off()

