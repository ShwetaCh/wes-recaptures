#libraries & custom functions
suppressWarnings(library(data.table))
suppressWarnings(library(plyr))
suppressWarnings(library(dplyr))
suppressWarnings(library(stringr))
suppressWarnings(library(ggplot2))
suppressWarnings(library(reshape2))
suppressWarnings(library(forcats))
suppressWarnings(library(RColorBrewer))
library(purrr)
library(grid)
library(gridExtra)
source('/ifs/res/taylorlab/chavans/scripts-orphan/multiplot.R')

specify_decimal = function(x, k) format(round(x, k), nsmall=k)
"%ni%" = Negate("%in%")
curr_date = format(Sys.Date(),"%d-%b-%y")

ex_clin = fread('~/tempo-cohort-level/WES_metadata_040620.txt'); 
dim(ex_clin); head(ex_clin)

im_clin = fread('~/tempo-cohort-level/data_clinical_sample_011920.txt',skip = 4) %>% 
  select(SAMPLE_ID, CANCER_TYPE) %>% 
  distinct(.); 
dim(im_clin); head(im_clin)

#Exome
ct_count_ex = ex_clin %>% 
  group_by(Cancer_Type_Aggregate) %>% 
  summarise(n = n()) %>% 
  arrange(desc(-n))
order_needed = filter(ct_count_ex, Cancer_Type_Aggregate != "Other") 
order_needed_ = c("Other",order_needed$Cancer_Type_Aggregate)
ct_count_ex$Cancer_Type_Aggregate = factor(ct_count_ex$Cancer_Type_Aggregate, levels = order_needed_)

#ct_count_ex$Cancer_Type <- factor(ct_count_ex$Cancer_Type, levels = rev(c("Breast Carcinoma","Small Cell Lung Cancer","Prostate Cancer","Non-Small Cell Lung Cancer",
                                                             # "Germ Cell Tumor","Bladder Cancer","Endometrial Cancer","Renal Cell Carcinoma",
                                                             # "Pancreatic Cancer","Ovarian Cancer","Glioma","Head and Neck Carcinoma",
                                                             # "Colorectal Cancer","Cancer of Unknown Primary","Melanoma","Soft Tissue Sarcoma",
                                                             # "Biliary Cancer","Other")))

ex_ct_barplot = ggplot(ct_count_ex, aes(x=Cancer_Type_Aggregate, y=n),levels=ct_count_ex$Cancer_Type_Aggregate) + 
  geom_bar(stat="identity", fill='steelblue') + coord_flip() + 
  xlab('Tumor types in Recaptured-WES cohort') + ylab('Sample count') + 
  theme_classic(base_size = 12) + theme(legend.position = 'none')

#IMPACT
ct_counv_im = count(im_clin, CANCER_TYPE) %>% arrange(desc(n)) %>% mutate(CANCER_TYPE = ifelse(n<25, 'Other', CANCER_TYPE))
ct_count_im = ct_counv_im %>% group_by(CANCER_TYPE) %>% summarise(n = sum(n)) %>% arrange(desc(-n))
order_needed = filter(ct_count_im, CANCER_TYPE != "Other") 
order_needed_ = c("Other",order_needed$CANCER_TYPE)
ct_count_im$CANCER_TYPE = factor(ct_count_im$CANCER_TYPE, levels = order_needed_)

im_ct_barplot = ggplot(ct_count_im, aes(x=CANCER_TYPE, y=n)) + 
  geom_bar(stat="identity", fill='darkgray') + coord_flip() + 
  xlab('Tumor types in MSK-IMPACT-Clinical-Series') + ylab('Sample count') + 
  theme_classic(base_size = 12) + theme(legend.position = 'none')

#grid.arrange(ex_ct_barplot,im_ct_barplot, nrow = 1, newpage = F)
#ggsave('~/tempo-cohort-level/Figure1A_CancerTypeDistributions.pdf', plot = grid.arrange(ex_ct_barplot,im_ct_barplot, nrow = 1, newpage = F))

ct_count_ex$Cancer_Type_Aggregate = str_replace(ct_count_ex$Cancer_Type_Aggregate,"Carcinoma","Cancer")
ct_count_ex$Cancer_Type_Aggregate = str_replace(ct_count_ex$Cancer_Type_Aggregate,"Biliary Cancer","Hepatobiliary Cancer")
ct_count_im$CANCER_TYPE = str_replace(ct_count_im$CANCER_TYPE,"Carcinoma","Cancer")

comb = left_join(ct_count_im,ct_count_ex, by=c(CANCER_TYPE='Cancer_Type_Aggregate')) %>%
  mutate(n.y = ifelse(is.na(n.y),0,n.y),
         CANCER_TYPE1 = ifelse(n.y==0,'Other',CANCER_TYPE)) %>%
  select(IMcount = n.x, WEScount = n.y, CANCER_TYPE = CANCER_TYPE1) %>%
  group_by(CANCER_TYPE) %>%
  summarise(IMcount = sum(IMcount), WEScount = sum(WEScount)) %>% 
  arrange(desc(-IMcount))

comb.m = melt(comb, id = 'CANCER_TYPE')  
unique(comb.m$CANCER_TYPE)

order_needed = filter(ct_count_ex, Cancer_Type_Aggregate != "Other") 
order_needed_ = c("Other",order_needed$Cancer_Type_Aggregate)
comb.m$CANCER_TYPE <- factor(comb.m$CANCER_TYPE, levels = rev(order_needed_))

comb_barplot = ggplot(comb.m, aes(x=factor(CANCER_TYPE), y=value, fill=variable)) + 
  theme_classic(base_size = 12) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  scale_fill_manual(values=c('steelblue','darkgrey')) +
  scale_y_log10()+
  #coord_flip() + 
  #xlab('Tumor types in MSK-IMPACT-Clinical-Series') + 
  ylab('Sample count') + xlab('') +
  theme(legend.position = 'bottom', legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))
comb_barplot

ggsave('~/tempo-cohort-level/Figure1A_CancerTypeDistributions.pdf', plot = comb_barplot, width = 11, height = 8)
