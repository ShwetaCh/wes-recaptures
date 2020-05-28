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
library(gsheet)
source('/ifs/res/taylorlab/chavans/scripts-orphan/multiplot.R')

specify_decimal = function(x, k) format(round(x, k), nsmall=k)
"%ni%" = Negate("%in%")
curr_date = format(Sys.Date(),"%d-%b-%y")

surv = fread('~/tempo-cohort-level/top_surv_data_MSI_level1_2_3_050520.txt'); head(surv); dim(surv);
ex_clin0 = fread('~/tempo-cohort-level/WES_metadata_040620.txt') %>% 
  select(DMP, Tumor_Sample_Barcode, Cancer_Type_Aggregate); head(ex_clin0); dim(ex_clin0); 
im_clin0 = fread('~/tempo-cohort-level/IM_metadata_040620.txt') %>% 
  select(DMP, CMO, MSIscore, MSIscore_Class) %>% 
  distinct(.); dim(im_clin0); head(im_clin0)
##########################################################################################

surv_master = inner_join(surv, ex_clin0, by=c(DMPID = 'DMP')); dim(surv_master)
surv_master = inner_join(surv_master, im_clin0, by=c(DMPID = 'DMP')); dim(surv_master) #225

surv_master$Death_Code <- as.factor(surv_master$Death_Code)
surv_master$Progression_Code <- as.factor(surv_master$Progression_Code)
surv_master$Cancer_Type <- as.factor(surv_master$Cancer_Type_Aggregate)
##########################################################################################
surv_ = distinct(surv_master) %>% 
  filter(MSIGotTherapy %in% c('Indeterminate','High'), MSIscore_Class %in% c('Indeterminate','Instable')) %>% 
  distinct(.) 
dim(surv_) #56

surv_ = surv_ %>%
  mutate(DMPID1 = factor(DMPID, 
                         levels = surv_[order(surv_$MSIscore_Class,surv_$TotalTime),1]))
head(surv_)
surv_.m = melt(surv_ %>% select(DMPID1, TreatmentStopped, Death, TreatmentContinued, Progression),id.var=c("DMPID1")) #, WES_miMSS, DMP_miMSS
head(surv_.m)

p1 = ggplot(surv_, aes(x=DMPID1, y=TotalTime)) + 
  geom_bar(stat = "identity", aes(fill = as.factor(Cancer_Type_Aggregate), color = MSIscore_Class), width=0.7) + coord_flip() +   
  geom_point(data = surv_.m, aes(DMPID1, value, shape=variable )) + #, colour=variable
  scale_color_manual(values = c('darkgray','brown')) + 
  scale_fill_manual(values = getPalette(colourCount)) + 
  scale_shape_manual(values = c(10,7,9,8,6,25)) +
  scale_y_continuous(limits=c(-3,90), breaks = seq(0, 90, len = 16)) +
  xlab('Patients with MSI Instable / Indeterminate & recieved immunotherapy') + 
  ylab("Months from treatment start date") + 
  theme_classic(base_size=9) + 
  theme(legend.position = "right", legend.direction = "vertical", legend.title = element_blank())
p1

#--------------------------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------------------------

surv_2 = distinct(surv_master) %>% 
  filter(!is.na(LevelGotTherapy)) %>% 
  distinct(.)
dim(surv_2)
surv_2 = surv_2 %>%
  mutate(DMPID1 = factor(DMPID, levels = unique(surv_2[order(surv_2$LevelGotTherapy, surv_2$TotalTime),1])))
#,
head(surv_2)
dim(surv_2)

surv_2.m = melt(surv_2 %>% select(DMPID1, TreatmentStopped, Death, TreatmentContinued,Progression),id.var=c("DMPID1")) #, WES_miMSS, DMP_miMSS
head(surv_2.m)

p2 = ggplot(surv_2, aes(x=DMPID1, TotalTime, y=TotalTime)) + 
     geom_bar(stat = "identity", aes(fill = as.factor(Cancer_Type_Aggregate)), width=0.7) + #aes(fill = as.factor(Cancer_Type_Aggregate), color = Group)
  coord_flip() +   
     geom_point(data = surv_2.m, aes(DMPID1, value, shape=variable)) + #colour=variable, 
  scale_fill_manual(values = getPalette(colourCount)) +
  #scale_color_manual(values = c('darkgray','brown')) + 
       scale_shape_manual(values = c(10,7,9,8)) +
       scale_y_continuous(limits=c(-3,90), breaks = seq(0, 90, len = 16)) +
          xlab('Patients with OncoKB alterations & recieved targeted-therapy') + 
          ylab("Months from treatment start date") + 
            theme_classic(base_size=9) + 
            theme(legend.position = "right", legend.direction = "vertical", legend.title = element_blank())
p2

#######################
#plot_grid( p1,p2, align = "hv", ncol = 1, rel_heights = c(0.5,0.5))
ggsave('~/tempo-cohort-level/Figure1C_SwimmersMSIsensor.pdf', plot = p1 )
ggsave('~/tempo-cohort-level/Figure1C_SwimmersOncoKBLevels.pdf', plot = p2 )



