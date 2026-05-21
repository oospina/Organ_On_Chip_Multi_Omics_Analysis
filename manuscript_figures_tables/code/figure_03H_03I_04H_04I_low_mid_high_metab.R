##
# Code to generate Figures 3H, 3I, 4H, 4I
#
# By Oscar Ospina
#
# Created: Sep 07, 2025
# Modified: May 21, 2026
#

# Load libraries
library('tidyverse')
library('DESeq2')
library('ComplexHeatmap')

# Load color palette
col_pal = readRDS('../manuscript_figures_tables/data/color_palettes.RDS')

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
  filter(compar %in% c("Mid_vs_Low", "Mid_vs_High", "Low_vs_High")) %>%
  filter(t_pval_adj < 0.05) %>%
  select(c('mol_id', 'mol_annot')) %>%
  distinct()

# Read metadata and subset to relevant samples
fp = '../coculture_supernatant_metabolomics/data/Updated_Metadata_Table_MOD.csv'
df_annot = read_delim(fp) %>%
  filter(group %in% c("Mid", "High", "Low"))  %>%
  mutate(group=factor(group, c('Low', 'Mid', 'High'))) %>%
  select(c('mod_sample_name', 'clean_sample_name', 'group')) %>%
  column_to_rownames('clean_sample_name')
rm(fp) # Clean env

# Subset matrix to DA metabolites
intx_mtx_sub = intx_mtx[ mets_plot[['mol_id']], df_annot[['mod_sample_name']] ]
intx_mtx_df = as.data.frame(t(intx_mtx_sub)) %>% 
  rownames_to_column('mod_sample_name') %>%
  left_join(df_annot, ., by='mod_sample_name')

################################################################################
# FIGURES_03_I AND FIGURES_04_I

# Create boxplots
## Low vs Mid
bp_ls1 = lapply(unique(mets_plot[['mol_id']]), function(i){
  # Get annotation of molecule
  metab_name_tmp = mets_plot[['mol_annot']][mets_plot[['mol_id']] == i]
  
  # Make plot
  intx_tmp = intx_mtx_df[intx_mtx_df[[1]] %in% c(df_annot %>% filter(group %in% c('Low', 'Mid')) %>% pull(mod_sample_name)), ]
  bp_p = ggplot(intx_tmp, aes(x=group, y=get(i), color=group)) +
    geom_boxplot(aes(fill=group), alpha=0.5) +
    geom_point() +
    scale_color_manual(values=col_pal) +
    scale_fill_manual(values=col_pal) +
    labs(x='', y='Median normalized intensity', title=metab_name_tmp) +
    theme_classic() +
    theme(panel.background=element_rect(fill=NA, color='black'),
          axis.line=element_line(color=NA),
          legend.background=element_blank(),
          legend.key=element_blank())
  
  return(bp_p)
})

## High vs Mid
bp_ls2 = lapply(unique(mets_plot[['mol_id']]), function(i){
  # Get annotation of molecule
  metab_name_tmp = mets_plot[['mol_annot']][mets_plot[['mol_id']] == i]
  
  # Make plot
  intx_tmp = intx_mtx_df[intx_mtx_df[[1]] %in% c(df_annot %>% filter(group %in% c('High', 'Mid')) %>% pull(mod_sample_name)), ]
  bp_p = ggplot(intx_tmp, aes(x=group, y=get(i), color=group)) +
    geom_boxplot(aes(fill=group), alpha=0.5) +
    geom_point() +
    scale_color_manual(values=col_pal) +
    scale_fill_manual(values=col_pal) +
    labs(x='', y='Median normalized intensity', title=metab_name_tmp) +
    theme_classic() +
    theme(panel.background=element_rect(fill=NA, color='black'),
          axis.line=element_line(color=NA),
          legend.background=element_blank(),
          legend.key=element_blank())
  
  return(bp_p)
})

# graphics.off()
# pdf('./results/FIGURE_03_I.pdf', height=10, width=8)
# print(ggpubr::ggarrange(plotlist=bp_ls1, ncol=3, nrow=3))
# dev.off()

# graphics.off()
# pdf('./results/FIGURE_04_I.pdf', height=10, width=8)
# print(ggpubr::ggarrange(plotlist=bp_ls2, ncol=3, nrow=3))
# dev.off()

rm(intx_mtx_df, bp_ls1, bp_ls2) # Clean env


################################################################################
# FIGURES_03_H

# Metabolites to highlight (assisted by GPT 4.1)
metab_show = c(
  "Biotin",
  "L-Alanine;D-Alanine",
  "L-Lysine",
  "Propionate",
  "Taurocholic acid (TCA)",
  "Gamma-linolenic acid",
  "Docosahexaenoic acid",
  "D-Pantothenic acid",
  "Xanthine",
  "Stearidonic acid",
  "D-(+)-Trehalose;D-Lactose;Maltose;D-(+)-Cellobiose;Isomaltose",
  "Nicotinamide",
  "Cortisol",
  "4-Pyridoxate",
  "Ursodeoxycholic acid;Chenodeoxycholic acid (CDCA);Deoxycholate;Isoursodeoxycholic acid"
)

# Subset to only Fasting and Control samples
hm_mtx = intx_mtx_sub[,  grepl('_low_|_mid_', colnames(intx_mtx_sub))]
df_annot1 = df_annot %>% filter(group %in% c('Low', 'Mid'))

# Prepare intensity matrix
hm_mtx = t(scale(t(hm_mtx[, match(df_annot1[['mod_sample_name']], colnames(hm_mtx))])))
hm_mtx = hm_mtx[match(mets_plot[['mol_id']], rownames(hm_mtx)), ]
colnames(hm_mtx) = rownames(df_annot1)
rownames(hm_mtx) = mets_plot[['mol_annot']]

# Create heatmap annotation object
hm_annot = HeatmapAnnotation(df=df_annot1 %>% select(-c('mod_sample_name')),
                             col=list(group=col_pal[c('Low', 'Mid')]),
                             show_legend=FALSE)

# Put legends together for vertical layout
lgd_list = list(
  Legend(labels=names(col_pal[c('Low', 'Mid')]), 
         legend_gp=gpar(fill=col_pal[c('Low', 'Mid')]))
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
                  column_title='Differentially abundant metabolites\nMid vs Low/Mid diet\n(T-test; FDR < 0.05)',
                  column_title_gp=gpar(fontsize=8),
                  bottom_annotation=hm_annot,
                  right_annotation=row_ann)

# graphics.off()
# pdf('./results/FIGURE_03_H.pdf', width=8, height=5)
# draw(hm_long, annotation_legend_list=NULL, heatmap_legend_list=comb_leg, 
#      heatmap_legend_side="right", annotation_legend_side="right", padding=unit(c(10, 2, 2, 10), "mm"))
# dev.off()

rm(list=grep('intx_mtx_sub|df_annot|col_pal|res_df|intx_mtx|metab_show|mets_plot', ls(), value=TRUE, invert=TRUE)) # Clean env


################################################################################
# FIGURES_04_H

# Subset to only Fasting and Control samples
hm_mtx = intx_mtx_sub[,  grepl('_high_|_mid_', colnames(intx_mtx_sub))]
df_annot2 = df_annot %>% filter(group %in% c('High', 'Mid'))

# Prepare intensity matrix
hm_mtx = t(scale(t(hm_mtx[, match(df_annot2[['mod_sample_name']], colnames(hm_mtx))])))
hm_mtx = hm_mtx[match(mets_plot[['mol_id']], rownames(hm_mtx)), ]
colnames(hm_mtx) = rownames(df_annot2)
rownames(hm_mtx) = mets_plot[['mol_annot']]

# Create heatmap annotation object
hm_annot = HeatmapAnnotation(df=df_annot2 %>% select(-c('mod_sample_name')),
                             col=list(group=col_pal[c('High', 'Mid')]),
                             show_legend=FALSE)

# Put legends together for vertical layout
lgd_list = list(
  Legend(labels=names(col_pal[c('High', 'Mid')]), 
         legend_gp=gpar(fill=col_pal[c('High', 'Mid')]))
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
                  column_title='Differentially abundant metabolites\nMid vs Mid/High diet\n(T-test; FDR < 0.05)',
                  column_title_gp=gpar(fontsize=8),
                  bottom_annotation=hm_annot,
                  right_annotation=row_ann)

# graphics.off()
# pdf('./results/FIGURE_04_H.pdf', width=8, height=5)
# draw(hm_long, annotation_legend_list=NULL, heatmap_legend_list=comb_leg, 
#      heatmap_legend_side="right", annotation_legend_side="right", padding=unit(c(10, 2, 2, 10), "mm"))
# dev.off()

