##
# Code to generate Figures S1E,F of the manuscript
#
# By Oscar Ospina
#
# Created: Aug 21, 2025
# Modified: May 20, 2026
#

library('DESeq2')
library('tidyverse')
library('ComplexHeatmap')

# Load color palette
col_pal = readRDS('./data/color_palettes.RDS')

# Read gene sets to display selected by Dr. McGilvrey
fp = './data/Selected_GOBP_Liver_Islets.xlsx'
selected_pw = readxl::read_excel(fp, sheet=1) %>%
  janitor::clean_names() %>%
  add_column(tissue=readxl::excel_sheets(fp)[1]) %>%
  bind_rows(., readxl::read_excel(fp, sheet=2) %>%
              janitor::clean_names() %>%
              add_column(tissue=readxl::excel_sheets(fp)[2])) %>%
  mutate(tissue=dplyr::recode(tissue, "Islets"="Pancreas"))
rm(fp) # Clean env

# Read FGSEA results
fp = '../coculture_tissues_rnaseq/results/fgsea_enrichment_scores_kegg_cc_bp.RDS'
fgsea_ls = readRDS(fp)

# Select relevant scores
fgsea_ls = fgsea_ls[grep('_iso_Fast_vs_iso_West', names(fgsea_ls), value=TRUE)]

# Get pathways
## GMT files obtained from https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
fps = list.files('../coculture_tissues_rnaseq/data/', pattern='\\.gmt', full.names=TRUE)
pws_ls = lapply(c('go.bp'), function(i){
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
names(pws_ls) = c('go.bp')
rm(fps) # Clean env

# Get size of pathway
pws_sz = lapply(names(pws_ls), function(i){
  ls_tmp = lapply(pws_ls[[i]], function(j){
    return(length(j))
  })
  return(ls_tmp)
})
names(pws_sz) = c('go.bp')
rm(pws_ls) # Clean env

# Create dataframe with NES data
gsea_df = lapply(names(fgsea_ls), function(i){
  key_tmp = str_extract(i, '^[A-Za-z]+')
  gsea_tmp = fgsea_ls[[i]][paste0(key_tmp, c('_iso_Fast_w', '_iso_West_w'))]
  gsea_tmp = lapply(names(gsea_tmp), function(j){
    tmp = gsea_tmp[[j]][['go.bp']] %>% add_column(tissue=key_tmp, 
                                                  tx=gsub(paste0(key_tmp, '_'), '', j)) %>%
      mutate(diet=case_when(str_detect(tx, 'Fast') ~ 'Low', TRUE ~ 'High'))
    return(tmp)
  }) %>% bind_rows() %>%
    as.data.frame() %>%
    dplyr::select('pathway', 'padj', 'NES', 'size', 'tissue', 'diet') %>%
    filter(pathway %in% (selected_pw %>% filter(tissue == key_tmp) %>% pull(pathway))) %>% # Subset to selected pathways
    left_join(., as.data.frame(t(data.frame(pws_sz[['go.bp']]))) %>%
                rownames_to_column('pathway') %>%
                dplyr::rename(pw_size=2), by='pathway') %>%
    mutate(size_frac=(size/pw_size) * 100) %>%
    dplyr::select(-c('pw_size')) %>%
    mutate(pathway=str_replace(pathway, '^GOBP_', '')) %>%
    mutate(pathway=case_when(is.na(padj) ~ pathway,
                             padj >= 0.05 ~ pathway,
                             padj < 0.05 ~ paste0(pathway, '_*'))) %>%
    select(-c('padj'))
  
  return(gsea_tmp)
}) %>% bind_rows() %>% 
  mutate(tissue=factor(tissue, levels=c('Pancreas', 'Liver')),
         diet=factor(diet, levels=c('Low', 'High'))) %>%
  filter(NES > 1) # Positive enrichment

# Generate bubble plot
bp = ggplot(gsea_df) +
  geom_point(aes(x=NES, y=reorder(pathway, dplyr::desc(NES)), color=NES, size=size_frac)) +
  ylab('') + xlab("Normalized enrichment score (NES)") +
  #ggtitle(i) +
  scale_color_gradient2(low="#FDBB47", mid="#B74202", high="#662506", midpoint=3) +
  theme(panel.background=element_rect(color='grey10', fill=NULL)) +
  facet_wrap(~tissue*diet, dir='v', scales='free_y')
# graphics.off()
# pdf('./results/FIGURE_XX_ISO_LOW_VS_HIGH_X_1.pdf', width=10, height=10)
# print(bp)
# dev.off()

