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

##################

im_clin = fread('~/tempo-cohort-level/IM_metadata_040620.txt')
dim(im_clin)  #1636

url <- 'docs.google.com/spreadsheets/d/1ZQCZ-02b8VNNDL05w-FdUiTpIVJ0WAW0EwxCkgX3oM0'

# Scatter plot for clonality effect on TMB
ex_roslin <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE) %>%
  select(DMP = DMPID, TMBWES = TMBExome, TMBIMPACT, TMBWES_IMgenes = TMBExomeIMgenes, TMBWES_NonIMgenes = TMBExome_NonIMgenes, Purity_Reviewed = PurityExome, Cancer_Type, TMBIMPACTClonal) %>% 
  filter(DMP %in% im_clin$DMP) %>% 
  filter(Cancer_Type %ni% c("Melanoma","Non-Small Cell Lung Cancer","Small Cell Lung Cancer","Bladder Cancer","Colorectal Cancer")) %>%
  filter(TMBIMPACT >=9, TMBIMPACT <=11) %>% 
  select(DMP, TMBIMPACTClonal, TMBIMPACT, Purity_Reviewed,Cancer_Type) %>%
  arrange(as.numeric(TMBIMPACT))
dim(ex_roslin) #34
ex_roslin

im_ex_clin_roslin_ = ex_roslin %>% 
  mutate(TMBIMPACT = ifelse(is.na(TMBIMPACT),0,TMBIMPACT),
         TMBIMPACTClonal = ifelse(is.na(TMBIMPACTClonal),0,TMBIMPACTClonal)) %>%
  mutate(grouptmb = ifelse(TMBIMPACT>=10,'High-TMB','Low-TMB'))

im_ex_clin_roslin_$grouptmb = as.factor(im_ex_clin_roslin_$grouptmb)

overall = ggplot(im_ex_clin_roslin_, aes(x = TMBIMPACT, y = TMBIMPACTClonal,  color = grouptmb)) + 
  geom_point(pch = 21, size = 2, fill = 'lightgray') + 
  #scale_fill_manual(values = c('lightgreen','lightgray')) +
  scale_color_manual(values = c('red','blue')) + 
  #coord_fixed() + 
  xlab('TMBIMPACT') + ylab('TMBIMPACT-Clonal') + 
  xlim(9,11) +
  geom_vline(xintercept = 10, linetype = 'dashed', color = 'red') +
  geom_hline(yintercept = 10, linetype = 'dashed', color = 'red') +
  #labs(title = paste0("R^2 = ",specify_decimal(r2,3),"\n")) +
  ggtitle('TMBIMPACT in non-ICB cancer types
  limited to samples with TMB ~10') + 
  theme_classic(base_size = 14) +
  theme(legend.title = element_blank(), legend.position = 'right', plot.margin = unit(c(1,1,1,1), 'lines'))
overall

# Scatter plot for MSI effect on TMB
ex_roslin0 <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE) 
head(ex_roslin0); names(ex_roslin0); dim(ex_roslin0);

ex_roslin = ex_roslin0 %>%
  select(DMP = DMPID, TMBWES = TMBExome, TMBIMPACT, TMBWES_IMgenes = TMBExomeIMgenes, TMBWES_NonIMgenes = TMBExome_NonIMgenes, Purity_Reviewed = PurityExome, Cancer_Type, MSIIMPACT_Class, MSIIMPACT, TMBIMPACTClonal) %>% 
  filter(DMP %in% im_clin$DMP) %>% 
  #filter(Cancer_Type %ni% c("Melanoma","Non-Small Cell Lung Cancer","Small Cell Lung Cancer","Bladder Cancer","Colorectal Cancer")) %>%
  filter(TMBIMPACT >=5, TMBIMPACT <=15) %>%
  filter(MSIIMPACT_Class %in% c("Instable","Stable","Indeterminate")) %>%
  select(DMP, MSIIMPACT, MSIIMPACT_Class, TMBIMPACT, TMBIMPACTClonal, Purity_Reviewed, Cancer_Type) %>%
  arrange(desc(TMBIMPACT))
dim(ex_roslin) #54 #170
write.table(ex_roslin,'~/tempo-cohort-level/TMBAround10data.txt',quote=F, row.names=F, append=F, sep ="\t")

im_ex_clin_roslin_ = ex_roslin %>% 
  mutate(TMBIMPACT = ifelse(is.na(TMBIMPACT),0,TMBIMPACT),
         TMBIMPACTClonal = ifelse(is.na(TMBIMPACTClonal),0,TMBIMPACTClonal)) %>%
  mutate(grouptmb = ifelse(TMBIMPACT>=10,'High-TMB','Low-TMB')) 

im_ex_clin_roslin_$grouptmb = as.factor(im_ex_clin_roslin_$grouptmb)
im_ex_clin_roslin_$MSIIMPACT_Class = as.factor(im_ex_clin_roslin_$MSIIMPACT_Class)


overall = ggplot(im_ex_clin_roslin_, aes(x = TMBIMPACT, y = MSIIMPACT,  color = MSIIMPACT_Class)) + 
  geom_point(pch = 21, size = 2, fill = 'lightgray') + 
  #scale_fill_manual(values = c('lightgreen','lightgray')) +
  scale_color_manual(values = c('blue','red','green')) + 
  #coord_fixed() + 
  xlab('TMBIMPACT') + ylab('MSIIMPACT') + 
  xlim(5,15) +
  geom_vline(xintercept = 10, linetype = 'dashed', color = 'red') +
  geom_hline(yintercept = 10, linetype = 'dashed', color = 'red') +
  #labs(title = paste0("R^2 = ",specify_decimal(r2,3),"\n")) +
  ggtitle('MSIscores for ~10 TMBIMPACT samples') + 
  theme_classic(base_size = 14) +
  theme(legend.title = element_blank(), legend.position = 'right', plot.margin = unit(c(1,1,1,1), 'lines'))
overall

# Scatter plot for Neo effect on TMB
