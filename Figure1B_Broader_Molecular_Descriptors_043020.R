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
       

###############
### Version2 plot - As of 080220
###############
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
  mutate(SigWESPresent = ifelse(DMP %in% SigWES$DMP,TRUE,FALSE),
         PercStrongBindingNeo = 100*NeoAntigenCountPerSampleWES_SB/MutationsPerSample) %>% 
  select(DMP, 
         Cancer_Type_Aggregate, 
         MSIscore, 
         SigWESPresent, 
         MutationsPerSample, 
         TMB,
         PercStrongBindingNeo,
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

cancer = ggplot(mol_desc, aes(x=reorder(DMP, mol_desc$DMP), y = "Cancer", fill=Cancer_Type_Aggregate)) + ylab('Cancer') +
  geom_bar(show.legend = T, stat = 'identity', width =1, na.rm = F) +
  scale_fill_manual(values = getPalette(colourCount)) + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank(), axis.ticks.y = element_blank(), axis.title.y=element_blank(), panel.grid = element_blank(), panel.background = element_blank(),legend.position="bottom", legend.direction="horizontal", legend.title = element_blank())
#cancer

msiexome = ggplot(mol_desc, aes(x=reorder(DMP, mol_desc$DMP), y = "MSI", fill=MSIscore)) + 
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

purity = ggplot(mol_desc, aes(x=reorder(DMP, mol_desc$DMP), y = Purity_Reviewed, fill=Purity_Reviewed)) + 
  geom_bar(show.legend = F, stat = 'identity', width =1, na.rm = F) + ylab('Purity') + ylim(0,1) +
  #scale_fill_gradient(low="skyblue", high="skyblue", na.value = 'grey80') + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank(), 
        axis.ticks.y = element_blank(), axis.title.y = element_text(size = 9),
        panel.grid = element_blank(), panel.background = element_blank())

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

westmb = ggplot(mol_desc, aes(x=reorder(DMP, mol_desc$DMP), y = TMB+1, fill=TMB+1)) + 
  geom_bar(show.legend = F, stat = 'identity', width =1, na.rm = F) + 
  scale_y_log10() + 
  ylab('TMB+1') + 
  #ylim(0,400) +
  #scale_fill_gradient(low="skyblue", high="skyblue", na.value = 'grey80') + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_text(size = 9),
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
wespercstrongbinbingneo = ggplot(mol_desc, aes(x=reorder(DMP, mol_desc$DMP), y = PercStrongBindingNeo, fill=PercStrongBindingNeo)) +
  geom_bar(show.legend = F, stat = 'identity', width =1, na.rm = F) + ylab('% Strong Binding \n Neoantigens') + ylim(0,100) +
  #scale_fill_gradient(low="skyblue", high="skyblue", na.value = 'grey80') + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x=element_blank(), 
        axis.ticks.y = element_blank(), axis.title.y = element_text(size = 9),
        panel.grid = element_blank(), panel.background = element_blank())
# 
# plot = plot_grid(msiexome, sigexome, purity, wgd, fga, wesmutcnt,westmb, wesnmb,westmb_nim, wesnmb_nim, cancer, align = "hv", ncol = 1, rel_heights = c(0.15,0.15,0.15,0.15,0.15,0.25,0.25,0.25,0.25,0.25,0.15))
# plot
###

plot = plot_grid(westmb, msiexome, wgd, purity, wespercstrongbinbingneo, cancer, align = "hv", ncol = 1, rel_heights = c(0.15,0.10,0.10,0.15,0.20,0.30))
plot
ggsave('~/tempo-cohort-level/Figure1B_MolecularDescriptors_v2.pdf', plot = plot, height=8,width =11)

###################################################################################
##Supplemental Figure2 Violing  plot by Cancer types
###################################################################################
ex_clin = fread('~/tempo-cohort-level/WES_metadata_040620.txt') %>% arrange(desc(TMB)); 
dim(ex_clin); head(ex_clin)

im_clin1 = fread('~/tempo-cohort-level/IM_metadata_040620.txt') %>% 
  select(DMP, TMBIMPACT,  MSIIMPACT = MSIscore) 
dim(im_clin1); head(im_clin1)

im_clin2 = fread('~/tempo-cohort-level/data_clinical_sample_011920.txt',skip = 4) %>% 
  select(SAMPLE_ID, PurityIMPACT = TUMOR_PURITY) %>% filter(SAMPLE_ID %in% im_clin1$DMP)
  distinct(.); 
dim(im_clin2); head(im_clin2)

im_clin = full_join(im_clin1, im_clin2, by = c(DMP = 'SAMPLE_ID'))  

theme_set(theme_classic(base_size = 12))

SigWES = ex_clin %>% filter_at(vars(starts_with("Signature")),any_vars(. >= 0.2))
mol_desc = ex_clin %>% 
  mutate(SigWESPresent = ifelse(DMP %in% SigWES$DMP,TRUE,FALSE),
         PercStrongBindingNeo = 100*NeoAntigenCountPerSampleWES_SB/MutationsPerSample) %>% 
  select(DMP, 
         Cancer_Type_Aggregate, 
         MSIscore, 
         SigWESPresent, 
         MutationsPerSample, 
         TMB,
         Ploidy,
         PercStrongBindingNeo,
         TMBWES_NonIMgenes,
         NMBWES,
         NMBWES_NonIMGenes,
         Purity_Reviewed,
         WGD_status, 
         FGA)
table(mol_desc$SigWESPresent)

mol_desc1 = left_join(mol_desc, im_clin, by = 'DMP')

mol_desc$DMP = factor(mol_desc$DMP, levels=mol_desc$DMP[order(mol_desc$Cancer_Type_Aggregate)])
colourCount = length(unique(mol_desc$Cancer_Type_Aggregate))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

#TMB
tmb_w = ggplot(mol_desc, aes(x=Cancer_Type_Aggregate, y=TMB+1, fill=Cancer_Type_Aggregate)) + 
  scale_y_log10() + scale_fill_manual(values = getPalette(colourCount)) + 
  geom_violin(show.legend = F) + geom_jitter(size =0.1, show.legend = F) +
  theme(axis.text.x = element_text(angle = 90), axis.title.x = element_blank())

#MSI
msi_w = ggplot(mol_desc, aes(x=Cancer_Type_Aggregate, y=MSIscore, fill=Cancer_Type_Aggregate)) + 
  scale_fill_manual(values = getPalette(colourCount)) + 
  geom_violin(show.legend = F) + geom_jitter(size =0.1, show.legend = F) +
  theme(axis.text.x = element_text(angle = 90), axis.title.x = element_blank())

#Purity
purity_w = ggplot(mol_desc, aes(x=Cancer_Type_Aggregate, y=Purity_Reviewed, fill=Cancer_Type_Aggregate)) + 
  scale_fill_manual(values = getPalette(colourCount)) + scale_y_continuous(breaks=seq(0,1,0.2)) +
  geom_violin(show.legend = F) + geom_jitter(size =0.1, show.legend = F) +
  theme(axis.text.x = element_text(angle = 90), axis.title.x = element_blank())

#Ploidy
ploidy_w = ggplot(mol_desc, aes(x=Cancer_Type_Aggregate, y=Ploidy, fill=Cancer_Type_Aggregate)) + 
  scale_fill_manual(values = getPalette(colourCount)) + scale_y_continuous(breaks=seq(0,8,2)) +
  geom_violin(show.legend = F) + geom_jitter(size =0.1, show.legend = F) +
  theme(axis.text.x = element_text(angle = 90), axis.title.x = element_blank())

plot = plot_grid(tmb_w, msi_w, purity_w, ploidy_w, align = "hv",ncol=1,rel_heights = c(0.25,0.25,0.25,0.25))
#plot
ggsave('~/tempo-cohort-level/SuppFigure2_MolecularDescriptorsByCancerTypes_Voilins_WES.pdf', plot = plot, height=11,width =8)

  # fill=name allow to automatically dedicate a color for each group

#TMB
tmb_i = ggplot(mol_desc1, aes(x=Cancer_Type_Aggregate, y=TMBIMPACT+1, fill=Cancer_Type_Aggregate)) + 
  scale_y_log10() + scale_fill_manual(values = getPalette(colourCount)) + 
  geom_violin(show.legend = F) + geom_jitter(size =0.1, show.legend = F) +
  theme(axis.text.x = element_text(angle = 90), axis.title.x = element_blank())

#MSI
msi_i = ggplot(mol_desc1, aes(x=Cancer_Type_Aggregate, y=MSIIMPACT, fill=Cancer_Type_Aggregate)) + 
  scale_fill_manual(values = getPalette(colourCount)) + 
  geom_violin(show.legend = F) + geom_jitter(size =0.1, show.legend = F) +
  theme(axis.text.x = element_text(angle = 90), axis.title.x = element_blank())

#Purity
mol_desc0 = filter(mol_desc1, !is.na(PurityIMPACT), PurityIMPACT != "N/A") %>% mutate(.,PurityIMPACT = as.numeric(PurityIMPACT)/100)
table(is.na(mol_desc0$PurityIMPACT))
purity_i = ggplot(mol_desc0, aes(x=Cancer_Type_Aggregate, y=PurityIMPACT, fill=Cancer_Type_Aggregate)) + 
  scale_fill_manual(values = getPalette(colourCount)) + scale_y_continuous(breaks=seq(0,1,0.2)) +
  geom_violin(show.legend = F) + geom_jitter(size =0.1, show.legend = F) +
  theme(axis.text.x = element_text(angle = 90), axis.title.x = element_blank())

plot2 = plot_grid(tmb_i, msi_i, purity_i, align = "hv",ncol=1,rel_heights = c(0.25,0.25,0.25))
plot2
ggsave('~/tempo-cohort-level/SuppFigure2_MolecularDescriptorsByCancerTypes_Voilins_IM.pdf', plot = plot2, height=11,width =8)
