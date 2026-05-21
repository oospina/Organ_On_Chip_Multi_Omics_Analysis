##
# Color palettes
#

# Experiment/batch
exp = c(Exp_2="orange", Exp_3="#C71585")

# Diets
diet = c(Mid="#752db0", Low="#b24d9e", High="#131e94",
         Metformin="#69ab5f", Semaglutide="#52a3c5")

# Tissues
tissue = c(Adipose="#000000", Brain="#E69F00", Colon="#56B4E9", 
           Pancreas="#009E73", Liver="#F0E442", Muscle="#0072B2")

# Culture type
culture = c(Interaction="#00CED1", Isolation="#DC143C")

# Put all palettes together
col_pal =  c(exp, diet, tissue, culture)

# Save palette
saveRDS(col_pal, '../data/color_palettes.RDS')

