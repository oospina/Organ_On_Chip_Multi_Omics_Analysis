##
# Code to generate Figure S2 of the manuscript
#
# By Oscar Ospina
#
# Created: Aug 29, 2025
# Modified: May 21, 2025
#

library('tidyverse')

# Read FGSEA results
fp = '../coculture_tissues_rnaseq/results/fgsea_enrichment_scores_kegg_cc_bp.RDS'
fgsea_all_ls = readRDS(fp)

# Read gene set sizes
fp = '../coculture_tissues_rnaseq/data/gene_set_sizes_list_kegg_gobp_reactome.RDS'
pws_sz = readRDS(fp)

# Read gene sets to display selected by Dr. Trapecar
fp = './data/Selected_GSEA_Iso_vs_Int_all_tissues.xlsx'
tissues = readxl::excel_sheets(fp)

# Create key for spreadsheet names
key = c(Adipocytes="Adipose",  Muscle="SkMuscl", `Large intestine`= "LargInt", 
        Brain="Brain", Pancreas="Pancreas", Liver="Liver")


############### Normal samples; isolation vs interaction ############### 
# Select relevant scores
fgsea_ls = fgsea_all_ls[grep('_iso_Norm_vs_int_Norm', names(fgsea_all_ls), value=TRUE)]

# Extract relevant gene set names by tissue
select_pw = lapply(tissues, function(i){
  # Read gene set names for each tissue
  df_tmp = readxl::read_excel(fp, sheet=i)
  
  # Match line separating each selection (ISOLATION) and extract subsequent lines until empty line (NA)
  iso_idx = which(str_detect(df_tmp[[1]], 'IN ISOLATION \\(control\\)')) + 1
  iso_tmp = ''
  iso_pws = c()
  while(!is.na(iso_tmp)){
    iso_tmp = df_tmp[[1]][iso_idx]
    if(!is.na(iso_tmp)){
      iso_pws = append(iso_pws, iso_tmp) 
    }
    iso_idx = iso_idx + 1
  }
  rm(iso_idx, iso_tmp) # Clean env
  
  # Match line separating each selection (INTERACTION) and extract subsequent lines until empty line (NA)
  int_idx = which(str_detect(df_tmp[[1]], 'IN INTERACTION \\(control\\)')) + 1
  int_tmp = ''
  int_pws = c()
  while(!is.na(int_tmp)){
    int_tmp = df_tmp[[1]][int_idx]
    if(!is.na(int_tmp)){
      int_pws = append(int_pws, int_tmp) 
    }
    int_idx = int_idx + 1
  }
  rm(int_idx, int_tmp) # Clean env
  
  return(c(iso_pws, int_pws))
})
names(select_pw) = tissues

# Create dataframe with NES data
bub_df_ls = lapply(names(select_pw), function(i){
  key_tmp = as.vector(key[i])
  gsea_tmp = fgsea_ls[[grep(paste0('^', key_tmp), names(fgsea_ls), value=TRUE)]]
  gsea_tmp = as.data.frame(gsea_tmp[[paste0(key_tmp, '_iso_Norm_w')]][['go.bp']]) %>%
    dplyr::select('pathway', 'padj_iso'='padj', 'NES_iso'='NES', 'size_iso'='size') %>%
    left_join(., as.data.frame(gsea_tmp[[paste0(key_tmp, '_int_Norm_w')]][['go.bp']]) %>%
                dplyr::select('pathway', 'padj_int'='padj', 'NES_int'='NES', 'size_int'='size'), by='pathway') %>%
    filter(pathway %in% select_pw[[i]]) %>% # Subset to pathways selected by Dr. Trapecar
    left_join(., as.data.frame(t(data.frame(pws_sz[['go.bp']]))) %>%
                rownames_to_column('pathway') %>%
                dplyr::rename(pw_size=2), by='pathway') %>%
    mutate(size_iso=(size_iso/pw_size) * 100) %>%
    mutate(size_int=(size_int/pw_size) * 100) %>%
    dplyr::select(-c('pw_size')) %>%
    mutate(pathway=str_replace(pathway, '^GOBP_', '')) %>%
    mutate(pathway=case_when(is.na(padj_iso) | is.na(padj_int) ~ pathway,
                             padj_iso >= 0.05 & padj_int >= 0.05 ~ pathway,
                             padj_iso < 0.05 & padj_int < 0.05 ~ paste0(pathway, '_1.2'),
                             padj_iso < 0.05 & padj_int >= 0.05 ~ paste0(pathway, '_1'),
                             padj_iso >= 0.05 & padj_int < 0.05 ~ paste0(pathway, '_2'))) %>%
    select(-c('padj_iso', 'padj_int')) %>%
    pivot_longer(cols=c('NES_iso', 'NES_int', 'size_iso', 'size_int', ), names_to=c(".value", "group"), names_sep = "_") %>%
    mutate(group=ifelse(group == 'iso', 'Isolation', 'Interaction')) %>%
    mutate(group=factor(group, levels=c('Isolation', 'Interaction')))
  
  return(gsea_tmp)
})
names(bub_df_ls) = names(select_pw)

# Generate plots
bub_p = lapply(names(bub_df_ls), function(i){
  # Create plot
  bp = ggplot(bub_df_ls[[i]]) +
    geom_point(aes(x=NES, y=reorder(pathway, dplyr::desc(NES)), color=NES, size=size)) +
    ylab('') + xlab("Normalized enrichment score (NES)") +
    ggtitle(i) +
    khroma::scale_color_YlOrBr() +
    theme(panel.background=element_rect(color='grey10', fill=NULL)) +
    facet_wrap(~group)
})
names(bub_p) = as.vector(key[names(bub_df_ls)])

# Print plots to file
lapply(names(bub_p), function(i){
  graphics.off()
  pdf(paste0('./results/', i, '_NormalIsoVsInt_selected_fgsea_scores_go.bp.pdf'), width=12)
  print(bub_p[[i]])
  dev.off()
})

