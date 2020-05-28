#libraries & custom functions
suppressWarnings(library(data.table))
suppressWarnings(library(plyr))
suppressWarnings(library(dplyr))
suppressWarnings(library(stringr))
suppressWarnings(library(ggplot2))
suppressWarnings(library(reshape2))
suppressWarnings(library(forcats))
suppressWarnings(library(RColorBrewer))
library(cowplot)
library(purrr)
library(grid)
library(gridExtra)
source('/ifs/res/taylorlab/chavans/scripts-orphan/multiplot.R')

specify_decimal = function(x, k) format(round(x, k), nsmall=k)
"%ni%" = Negate("%in%")
curr_date = format(Sys.Date(),"%d-%b-%y")

ex_clin = fread('~/tempo-cohort-level/WES_metadata_040620.txt') %>% arrange(desc(TMB)); 
dim(ex_clin); head(ex_clin)
im_clin = fread('~/tempo-cohort-level/data_clinical_sample_011920.txt',skip = 4) %>% 
  select(SAMPLE_ID, CANCER_TYPE) %>% 
  distinct(.); 
dim(im_clin); head(im_clin)

theme_set(theme_classic(base_size = 12))

SigWES = ex_clin %>% filter_at(vars(starts_with("Signature")),any_vars(. >= 0.2))
mol_desc = ex_clin %>% 
           mutate(SigWESPresent = ifelse(DMP %in% SigWES$DMP,TRUE,FALSE)) %>% 
           select(DMP, 
                  Cancer_Type_Aggregate, 
                  MSIscore, 
                  SigWESPresent, 
                  MutationsPerSample, 
                  TMB,
                  TMBWES_NonIMgenes,
                  NMBWES,
                  NMBWES_NonIMGenes,
                  Purity_Reviewed,
                  WGD_status, 
                  FGA)
table(mol_desc$SigWESPresent)

mol_desc$DMP = factor(mol_desc$DMP, levels=mol_desc$DMP[order(mol_desc$Cancer_Type_Aggregate)])
colourCount = length(unique(mol_desc$Cancer_Type_Aggregate))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

cancer = ggplot(mol_desc, aes(x=reorder(DMP, mol_desc$DMP), y = "Cancer_Type", fill=Cancer_Type_Aggregate)) + 
  geom_bar(show.legend = F, stat = 'identity', width =1, na.rm = F) +
  scale_fill_manual(values = getPalette(colourCount)) + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank(), 
        axis.ticks.y = element_blank(), axis.title.y=element_blank(), 
        panel.grid = element_blank(), panel.background = element_blank(),
        legend.position="bottom", legend.direction="horizontal", legend.title = element_blank())
#cancer

msiexome = ggplot(mol_desc, aes(x=reorder(DMP, mol_desc$DMP), y = "MSIscore", fill=MSIscore)) + 
  geom_bar(show.legend = F, stat = 'identity', width = 1, na.rm = F) +
  scale_fill_gradient(low="cornsilk", high="darkgrey", na.value = 'grey80') + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank(), 
        axis.ticks.y = element_blank(), axis.title.y=element_blank(),
        panel.grid = element_blank(), panel.background = element_blank(),
        legend.position="bottom", legend.direction="horizontal")

sigexome = ggplot(mol_desc, aes(x=reorder(DMP, mol_desc$DMP), y = "SignatureWES", fill=SigWESPresent)) + 
  geom_bar(show.legend = F, stat = 'identity', width = 1, na.rm = F) +
  scale_fill_manual(values = c('cornsilk','darkgrey'), na.value = 'grey80') + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank(), 
        axis.ticks.y = element_blank(), axis.title.y=element_blank(),
        panel.grid = element_blank(), panel.background = element_blank(),
        legend.position="bottom", legend.direction="horizontal")

# ###
# coverage = ggplot(mol_desc, aes(x=reorder(DMP, mol_desc$DMP), y = CoverageWES, fill=CoverageWES)) + 
#         geom_bar(show.legend = F, stat = 'identity', width =1, na.rm = F) +
#         #scale_fill_gradient(low="white", high="forestgreen", na.value = 'grey80') + 
#         theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank(), 
#         axis.ticks.y = element_blank(), axis.title.y = element_text(size = 7),
#         panel.grid = element_blank(), panel.background = element_blank())
# 
# baitset = ggplot(mol_desc, aes(x=reorder(DMP, mol_desc$DMP), y = "BaitSet", fill=BaitSet)) + 
#         geom_bar(show.legend = F, stat = 'identity', width =1, na.rm = F) +
#         scale_fill_manual(values = c('cornsilk','grey','khaki','seashell'),na.value = 'grey80') + 
#         theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank(), 
#         axis.ticks.y = element_blank(), axis.title.y=element_blank(),
#         panel.grid = element_blank(), panel.background = element_blank(),
#         legend.position="bottom", legend.direction="horizontal")
# ###

purity = ggplot(mol_desc, aes(x=reorder(DMP, mol_desc$DMP), y = "Purity", fill=Purity_Reviewed)) + 
          geom_bar(show.legend = F, stat = 'identity', width =1, na.rm = F) +
          scale_fill_gradient(low="darkgrey", high="skyblue", na.value = 'grey80') + 
          theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank(), 
          axis.ticks.y = element_blank(), axis.title.y=element_blank(), 
          panel.grid = element_blank(), panel.background = element_blank(),
          legend.position="bottom", legend.direction="horizontal")

wgd = ggplot(mol_desc, aes(x=reorder(DMP, mol_desc$DMP), y = "WGD", fill=WGD_status)) + 
        geom_bar(show.legend = F, stat = 'identity', width =1, na.rm = F) +
        scale_fill_manual(values = c('cornsilk','darkgrey'), na.value = 'grey80') + 
        theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank(), 
        axis.ticks.y = element_blank(), axis.title.y=element_blank(),
        panel.grid = element_blank(), panel.background = element_blank(),
        legend.position="bottom", legend.direction="horizontal")

fga = ggplot(mol_desc, aes(x=reorder(DMP, mol_desc$DMP), y = "FGA", fill=FGA)) + 
        geom_bar(show.legend = F, stat = 'identity', width =1, na.rm = F) +
        scale_fill_gradient(low="darkgrey", high="skyblue", na.value = 'grey80') + 
        theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank(), 
        axis.ticks.y = element_blank(), axis.title.y=element_blank(),
        panel.grid = element_blank(), panel.background = element_blank(),
        legend.position="bottom", legend.direction="horizontal")

###
wesmutcnt = ggplot(mol_desc, aes(x=reorder(DMP, mol_desc$DMP), y = MutationsPerSample, fill=MutationsPerSample)) + 
  geom_bar(show.legend = F, stat = 'identity', width =1, na.rm = F) + scale_y_log10() + ylab('Mutation Count') +
  #scale_fill_gradient(low="white", high="skyblue", na.value = 'grey80') + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank(), 
        axis.ticks.y = element_blank(), axis.title.y = element_text(size = 9),
        panel.grid = element_blank(), panel.background = element_blank())

westmb = ggplot(mol_desc, aes(x=reorder(DMP, mol_desc$DMP), y = TMB, fill=TMB)) + 
         geom_bar(show.legend = F, stat = 'identity', width =1, na.rm = F) + scale_y_log10() + ylab('TMB') +
        #scale_fill_gradient(low="skyblue", high="skyblue", na.value = 'grey80') + 
        theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank(), 
              axis.ticks.y = element_blank(), axis.title.y = element_text(size = 9),
              panel.grid = element_blank(), panel.background = element_blank())

westmb_nim = ggplot(mol_desc, aes(x=reorder(DMP, mol_desc$DMP), y = TMBWES_NonIMgenes, fill=TMBWES_NonIMgenes)) + 
        geom_bar(show.legend = F, stat = 'identity', width =1, na.rm = F) + scale_y_log10() + ylab('TMB_NIM') +
        #scale_fill_gradient(low="white", high="dodgerblue", na.value = 'grey80') + 
        theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank(), 
              axis.ticks.y = element_blank(), axis.title.y = element_text(size = 9),
              panel.grid = element_blank(), panel.background = element_blank())

wesnmb = ggplot(mol_desc, aes(x=reorder(DMP, mol_desc$DMP), y = NMBWES, fill=NMBWES)) + 
        geom_bar(show.legend = F, stat = 'identity', width =1, na.rm = F) + scale_y_log10() + ylab('NMB') +
        #scale_fill_gradient(low="white", high="skyblue", na.value = 'grey80') + 
        theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank(), 
          axis.ticks.y = element_blank(), axis.title.y = element_text(size = 9),
          panel.grid = element_blank(), panel.background = element_blank())

wesnmb_nim = ggplot(mol_desc, aes(x=reorder(DMP, mol_desc$DMP), y = NMBWES_NonIMGenes, fill=NMBWES_NonIMGenes)) + 
  geom_bar(show.legend = F, stat = 'identity', width =1, na.rm = F) + scale_y_log10() + ylab('NMB_NIM') +
  #scale_fill_gradient(low="white", high="skyblue", na.value = 'grey80') + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank(), 
        axis.ticks.y = element_blank(), axis.title.y = element_text(size = 9),
        panel.grid = element_blank(), panel.background = element_blank())
###

plot = plot_grid(msiexome, sigexome, purity, wgd, fga, wesmutcnt,westmb, wesnmb,westmb_nim, wesnmb_nim, cancer, align = "hv", ncol = 1, rel_heights = c(0.15,0.15,0.15,0.15,0.15,0.25,0.25,0.25,0.25,0.25,0.15))
plot

ggsave('~/tempo-cohort-level/Figure1B_MolecularDescriptors.pdf', plot = plot, height=8,width =11)
       