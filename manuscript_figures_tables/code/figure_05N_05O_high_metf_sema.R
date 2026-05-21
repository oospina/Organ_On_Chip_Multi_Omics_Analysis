##
# Code to generate Figures 5N, 5O
#
# By Oscar Ospina
#
# Created: Sep 09, 2025
# Modified: May 21, 2026
#

# Load libraries
library('tidyverse')
library('DESeq2')
library('ComplexHeatmap')

# Load DA analysis results
fp = '../coculture_supernatant_metabolomics/results/diff_abundance_ttest_wilcox_median_norm.csv'
res_df = read.csv(fp)
rm(fp) # Clean env

# Read batch corrected intensity data
fp = '../coculture_supernatant_metabolomics/data/median_normalized_only_id_intx_dup_sum_ComBat.RDS'
intx_mtx = readRDS(fp)
rm(fp) # Clean env

# Select metabolites to plot
mets_plot = res_df %>%
  filter(compar %in% c("High_vs_High + Metformin", "High_vs_High + Semaglutide")) %>%
  filter(abs(log2FC) > 0 & t_pval < 0.05) %>%
  select(c('mol_id', 'mol_annot')) %>%
  distinct()

# Read metadata and subset to relevant samples
fp = '../coculture_supernatant_metabolomics/data/Updated_Metadata_Table_MOD.csv'
df_annot = read_delim(fp) %>%
  filter(group %in% c("High", "High + Metformin", "High + Semaglutide"))  %>%
  select(c('mod_sample_name', 'clean_sample_name', 'group')) %>%
  column_to_rownames('clean_sample_name')
rm(fp) # Clean env

# Subset matrix to DA metabolites
intx_mtx_sub = intx_mtx[ mets_plot[['mol_id']], df_annot[['mod_sample_name']] ]
intx_mtx_df = as.data.frame(t(intx_mtx_sub)) %>% 
  rownames_to_column('mod_sample_name') %>%
  left_join(df_annot, ., by='mod_sample_name')

################################################################################
# FIGURE_05_O

# Create boxplots
bp_ls = lapply(unique(mets_plot[['mol_id']]), function(i){
  # Get annotation of molecule
  metab_name_tmp = mets_plot[['mol_annot']][mets_plot[['mol_id']] == i]
  
  # Make plot
  bp_p = ggplot(intx_mtx_df, aes(x=group, y=get(i), color=group)) +
    geom_boxplot(aes(fill=group), alpha=0.5) +
    geom_point() +
    scale_color_manual(values=col_pal) +
    scale_fill_manual(values=col_pal) +
    labs(x='', y='Median normalized intensity', title=metab_name_tmp) +
    theme_classic() +
    theme(panel.background=element_rect(fill=NA, color='black'),
          axis.line=element_line(color=NA),
          axis.text.x=element_text(angle=30, vjust=1, hjust=1),
          legend.background=element_blank(),
          legend.key=element_blank())
  
  return(bp_p)
})

# graphics.off()
# pdf('./results/FIGURE_05_O.pdf', height=10, width=10)
# print(ggpubr::ggarrange(plotlist=bp_ls, ncol=3, nrow=3))
# dev.off()

rm(intx_mtx_df, bp_ls) # Clean env


################################################################################
# FIGURE_05_N

# Metabolites to highlight
metab_show = c(
  "D-Pantothenic acid",
  "Biotin",
  "L-Isoleucine",
  "Leucine",
  "Adrenic acid",
  "Methyl acetoacetate;3-Methyl-2-oxobutanoic acid",
  "Succinic acid",
  "L-Lysine",
  "L-Sorbose;Myo-inositol;Mannose;Allose;D-Psicose;chiro-Inositol",
  "Pyridoxine",
  "L-Phenylalanine",
  "D-Glucuronic acid;D-(+)-Galacturonic acid",
  "Taurine",
  "Histidine",
  "Folic acid"
)

# Prepare intensity matrix
hm_mtx = t(scale(t(intx_mtx_sub[, match(df_annot[['mod_sample_name']], colnames(intx_mtx_sub))])))
hm_mtx = hm_mtx[match(mets_plot[['mol_id']], rownames(hm_mtx)), ]
colnames(hm_mtx) = rownames(df_annot)
rownames(hm_mtx) = mets_plot[['mol_annot']]

# Create heatmap annotation object
hm_annot = HeatmapAnnotation(df=df_annot %>% select(-c('mod_sample_name')),
                             col=list(group=col_pal[as.vector(unique(df_annot[['group']]))]),
                             show_legend=FALSE)

# Put legends together for vertical layout
lgd_list = list(
  Legend(labels=unique(df_annot[['group']]), 
         legend_gp=gpar(fill=col_pal[as.vector(unique(df_annot[['group']]))]))
)
comb_leg = packLegend(lgd_list[[1]], direction="vertical")

# Make annotation to highlight gene names
row_ann = rowAnnotation(foo=anno_mark(at=as.vector(na.omit(match(metab_show, rownames(hm_mtx)))),  
                                      labels=metab_show[metab_show %in% rownames(hm_mtx)],
                                      labels_gp=gpar(fontsize=8)))

# Generate heatmap
hm_long = Heatmap(hm_mtx,
                  clustering_distance_columns='canberra', 
                  clustering_method_columns='complete',
                  name='Scaled peak\nintensity',
                  cluster_rows=TRUE,
                  cluster_columns=TRUE,
                  show_row_names=FALSE,
                  show_column_names=FALSE,
                  show_row_dend=FALSE,
                  column_title='Differentially abundant metabolites\nHigh vs Semaglutide/Metformin\n(T-test; p-val < 0.05 and abs(log2FC) > 0)',
                  column_title_gp=gpar(fontsize=8),
                  bottom_annotation=hm_annot,
                  right_annotation=row_ann)

# graphics.off()
# pdf('./results/FIGURE_05_N.pdf', width=8, height=3)
# draw(hm_long, annotation_legend_list=NULL, heatmap_legend_list=comb_leg, 
#      heatmap_legend_side="right", annotation_legend_side="right", padding=unit(c(10, 2, 2, 10), "mm"))
# dev.off()

