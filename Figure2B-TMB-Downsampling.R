library(RColorBrewer)
library(Rmisc)
library(purrr)
library(grid)
library(gridExtra)
library(gsheet)
library(ggpubr)
library(cowplot)
library(gridExtra)
source('/ifs/res/taylorlab/chavans/scripts-orphan/multiplot.R')

specify_decimal = function(x, k) format(round(x, k), nsmall=k)
"%ni%" = Negate("%in%")
curr_date = format(Sys.Date(),"%d-%b-%y")

ex_clin = fread('~/tempo-cohort-level/WES_metadata_040620.txt')
dim(ex_clin); head(ex_clin)

ex_clin2 = fread('~/tempo-cohort-level/WES_allTMBs_051020.txt') %>% 
  select(Tumor_Sample_Barcode, TMBWES_IMgenes) %>%
  mutate(TMBWES_IMgenes = ifelse(is.na(TMBWES_IMgenes),0,TMBWES_IMgenes))

ex_clin3 = inner_join(ex_clin, ex_clin2, by = c(Tumor_Sample_Barcode = "Tumor_Sample_Barcode"))

im_clin = fread('~/tempo-cohort-level/IM_metadata_040620.txt') %>% 
  arrange(desc(TMBIMPACT)) %>%
  select(CMO, TMBIMPACT, TMBIMPACT_Downsampled); 
dim(im_clin); head(im_clin)

im_ex_clin = inner_join(ex_clin3, im_clin, by = c(Tumor_Sample_Barcode = "CMO"))


im_ex_clin0 = im_ex_clin %>% 
  select(DMP, TMBWES = TMB, TMBIMPACT, TMBWES_IMgenes, TMBWES_NonIMgenes, TMBIMPACT_Downsampled, Purity_Reviewed) %>% 
  filter(TMBIMPACT<10, Purity_Reviewed >=0.5) %>% 
  mutate(adj.depth = TMBIMPACT_Downsampled/TMBIMPACT, 
         adj.genecontent = TMBWES_NonIMgenes/TMBWES_IMgenes, 
         adj.TMBIMPACT = TMBIMPACT*adj.depth*adj.genecontent) %>% 
  select(-c("adj.depth","adj.genecontent", "Purity_Reviewed"))
dim(im_ex_clin0) #679   7

my_comparisons <- list( c("TMBIMPACT","TMBIMPACT_Downsampled"), 
                        c("TMBWES_IMgenes","TMBWES_NonIMgenes"), 
                        c("TMBWES","TMBWES_IMgenes"), 
                        c("TMBIMPACT","TMBWES"), 
                        c("TMBWES","TMBWES_NonIMgenes"), 
                        c("TMBIMPACT","adj.TMBIMPACT"), 
                        c("TMBWES","adj.TMBIMPACT"), 
                        c("TMBIMPACT","TMBWES_IMgenes") )

pdatm = melt(im_ex_clin0, id = c("DMP"))
pdatm = pdatm %>% mutate(group = variable, TMB = value) %>% select(-c(variable,value))
head(pdatm)  
mm_tmbs_ = pdatm
mm_tmbs_$group = factor(mm_tmbs_$group, levels = c("TMBIMPACT","TMBIMPACT_Downsampled","TMBWES_IMgenes","TMBWES_NonIMgenes","adj.TMBIMPACT","TMBWES"))
mm_tmbs_ = mm_tmbs_ %>% mutate(seq = ifelse(group %like% "WES", "WES", "IM"))
head(mm_tmbs_)  

p14 = ggplot(mm_tmbs_, aes(x=group, y=TMB, fill = seq)) + theme_classic(base_size = 16) +
  geom_boxplot() + ylab("TMB in Purity >50% \n(TMBIMPACT <10 Samples Only)") + xlab("") + 
  scale_fill_manual(values = c('gray35','gray65')) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  
p14 + stat_compare_means(aes(label = ..p.signif..), method = "wilcox.test", comparisons = my_comparisons)