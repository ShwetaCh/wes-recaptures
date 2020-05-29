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

########### For Roslin data TMB boxplot to show adj.TMB

im_dnsmpl_tmb = fread('~/tempo-cohort-level/fixed.downsampledTMB_ForRoslin.txt') %>% 
  select(DMP, TMBIMPACT_Downsampled = TMBIMPACT.downsampled.fixed) %>%
  filter(DMP %in% im_ex_clin$DMP)
dim(im_dnsmpl_tmb) #1413

url <- 'docs.google.com/spreadsheets/d/1ZQCZ-02b8VNNDL05w-FdUiTpIVJ0WAW0EwxCkgX3oM0'
ex_roslin <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE) %>%
  select(DMP = DMPID, TMBWES = TMBExome, TMBIMPACT, TMBWES_IMgenes = TMBExomeIMgenes, TMBWES_NonIMgenes = TMBExome_NonIMgenes, Purity_Reviewed = PurityExome) %>% 
  filter(DMP %in% im_ex_clin$DMP)
dim(ex_roslin) #1636

im_ex_clin_roslin = inner_join(im_dnsmpl_tmb, ex_roslin, by = 'DMP')
dim(im_ex_clin_roslin) #1413 -- figure out why rest are missing?!

im_ex_clin0_roslin = im_ex_clin_roslin %>% 
  select(DMP, TMBWES, TMBIMPACT, TMBWES_IMgenes, TMBWES_NonIMgenes, TMBIMPACT_Downsampled, Purity_Reviewed) %>% 
  filter(TMBIMPACT<10, Purity_Reviewed >=0.5) %>% 
  mutate(adj.depth = TMBIMPACT_Downsampled/TMBIMPACT, 
         adj.genecontent = TMBWES_NonIMgenes/TMBWES_IMgenes, 
         adj.TMBIMPACT = TMBIMPACT*adj.depth*adj.genecontent) %>% 
  select(-c("adj.depth","adj.genecontent", "Purity_Reviewed"))
dim(im_ex_clin0_roslin) #601   7

pdatm = melt(im_ex_clin0_roslin, id = c("DMP"))
pdatm = pdatm %>% mutate(group = variable, TMB = value) %>% select(-c(variable,value))
head(pdatm)  
mm_tmbs_ = pdatm
mm_tmbs_$group = factor(mm_tmbs_$group, levels = c("TMBIMPACT","TMBIMPACT_Downsampled","TMBWES_IMgenes","TMBWES_NonIMgenes","adj.TMBIMPACT","TMBWES"))
mm_tmbs_ = mm_tmbs_ %>% mutate(seq = ifelse(group %like% "WES", "WES", "IM"))
head(mm_tmbs_)  
p15 = ggplot(mm_tmbs_, aes(x=group, y=TMB, fill = seq)) + theme_classic(base_size = 16) +
  geom_boxplot() + ylab("TMB in Purity >50% \n(TMBIMPACT <10 Samples Only)") + xlab("") + 
  scale_fill_nejm() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        plot.margin = unit(c(1,1,1,1), 'lines'))  
p15 + stat_compare_means(aes(label = ..p.signif..), method = "wilcox.test", comparisons = my_comparisons)

pdf("~/tempo-cohort-level/Figure2B_TMB_Roslin.pdf", paper = "a4r")
p15 + stat_compare_means(aes(label = ..p.signif..), method = "wilcox.test", comparisons = my_comparisons)
dev.off()


######### For Roslin showing TMB IMPACT vs. WES

im_clin = fread('~/tempo-cohort-level/IM_metadata_040620.txt') %>% select(DMP, TMBIMPACT)
im_ex_clin_roslin1 = inner_join(im_clin, ex_roslin, by = c('DMP')) %>% 
  select(DMP, TMBIMPACT=TMBIMPACT.x, TMBWES)

dim(im_ex_clin_roslin1) #1636
# Equation
label = lm_eqn(im_ex_clin_roslin,y=im_ex_clin_roslin$TMBIMPACT,x=im_ex_clin_roslin$TMBWES)
label
# R²
r2 = summary(lm(TMBWES~TMBIMPACT,data=im_ex_clin_roslin))$adj.r.squared
r2
#Spearman Correlation
corr = cor.test(im_ex_clin_roslin$TMBWES, im_ex_clin_roslin$TMBIMPACT,method = "spearman") 
#print(str(corr))
rho = as.numeric(corr$estimate)
rho
### scatter plot LINEAR scale Original

# Plot
im_ex_clin_roslin_ = im_ex_clin_roslin1 %>% 
  mutate(TMBIMPACT = ifelse(is.na(TMBIMPACT),0,TMBIMPACT),
         TMBIMPACT = TMBIMPACT + 1, TMBWES = TMBWES +1 ) %>%
  mutate(grouptmb = ifelse(TMBIMPACT>=20 | TMBWES>=20,'High-TMB','Low-TMB'))

im_ex_clin_roslin_$grouptmb = as.factor(im_ex_clin_roslin_$grouptmb)

overall = ggplot(im_ex_clin_roslin_, aes(x = TMBIMPACT, y = TMBWES)) + 
  geom_point(pch = 21, size = 2, fill = 'lightgray') + 
  #scale_fill_manual(values = c('lightgreen','lightgray')) +
  #scale_color_manual(values = c('red','blue')) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  geom_smooth(method=lm, se=FALSE, lty=1) + 
  scale_y_continuous(trans = 'log10',  limits = c(1,500), breaks = c(0,1,2,5,10,20,50,100,200,500)) + 
  scale_x_continuous(trans='log10', limits = c(1,500), breaks = c(0,1,2,5,10,20,50,100,200,500)) + 
  coord_fixed() + xlab('TMBIMPACT + 1') + ylab('TMBWES + 1') + 
  labs(title = paste0("R^2 = ",specify_decimal(r2,3),"\n")) +
  theme_classic(base_size = 14) +
  theme(legend.title = element_blank(), legend.position = 'right', plot.margin = unit(c(1,1,1,1), 'lines'))
#overall

r2_ = summary(lm(TMBWES~TMBIMPACT,data=filter(im_ex_clin_roslin_, grouptmb == 'Low-TMB')))$adj.r.squared
r2_
overall_2 = ggplot(im_ex_clin_roslin_, aes(x = TMBIMPACT, y = TMBWES, fill = grouptmb)) + 
  geom_point(pch = 21, size = 2) + 
  scale_fill_manual(values = c('lightgreen','lightgray')) +
  scale_color_manual(values = c('red','blue')) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  geom_smooth(method=lm, se=FALSE, lty=1, aes(group=grouptmb, color=grouptmb)) + 
  scale_y_continuous(trans = 'log10',  limits = c(1,500), breaks = c(0,1,2,5,10,20,50,100,200,500)) + 
  scale_x_continuous(trans='log10', limits = c(1,500), breaks = c(0,1,2,5,10,20,50,100,200,500)) + 
  coord_fixed() + xlab('TMBIMPACT + 1') + ylab('TMBExome + 1') + 
  labs(title = paste0("R^2 = ",specify_decimal(r2_,3),"\n")) +
  theme_classic(base_size = 14) +
  theme(legend.title = element_blank(),legend.position = 'bottom', plot.margin = unit(c(1,1,1,1), 'lines'))
#overall_2
grid.arrange(overall, overall_2, newpage = FALSE, ncol = 2)

# tmb_plot = ggplot(im_ex_clin_roslin, 
#                   aes(x = TMBIMPACT, y = TMBWES)) + 
#   #facet_wrap(~Algorithm) +
#   geom_point(alpha = 0.5, size = 3, shape = 16, color = "#08519c") +
#   geom_abline(col = "blue") +
#   scale_x_continuous(limits = c(0,150)) +
#   scale_y_continuous(limits = c(0,150)) +
#   labs(x = "IMPACT TMB", y = "Exome TMB") +
#   geom_point(alpha = 0.3, size = 3, shape = 16) + 
#   geom_abline(col = "black",slope=1, intercept =0) + 
#   geom_smooth(method=lm, se=FALSE, lty=2) +
#   coord_fixed() +
#   theme_classic(base_size = 16) +
#   theme(axis.text.x = element_text(colour = "blue")) +
#   theme(axis.text.y = element_text(colour = "blue")) +
#   theme(plot.title = element_text(colour = "blue", size = 14)) +
#   labs(title = paste0("Non-Synonymous Coding Mutations-based-TMB","\n",
#                       "		R^2 = ",specify_decimal(r2,3),"\n",
#                       "		Spearman Correlation = ",specify_decimal(rho,3),"\n",
#                       "		IMPACT_TMB = 1.2*Exome_TMB + 1.9"))

pdf("~/tempo-cohort-level/Figure2A_TMB_Compare_Roslin.pdf", paper = "a4r")
grid.arrange(overall, overall_2, newpage = FALSE, ncol = 2)
dev.off()


#########################################################################################################################################################################

url <- 'docs.google.com/spreadsheets/d/1ZQCZ-02b8VNNDL05w-FdUiTpIVJ0WAW0EwxCkgX3oM0'
ex_roslin_purity_depth <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE) 
ex_roslin_purity_depth0 = ex_roslin_purity_depth %>%
  select(DMP = DMPID, TMBWES = TMBExome, purity_bin, depth_bin, TMBWESClonal = TMBExomeClonal) %>% 
  filter(DMP %in% im_ex_clin$DMP)
dim(ex_roslin_purity_depth0) #1636


# > filter(ex_roslin_purity_depth0, is.na(purity_bin))
# DMP TMBWES purity_bin depth_bin TMBWESClonal
# 1 P-0020573-T01-IM6      0       <NA> (100,200]            0
# 2 P-0020769-T01-IM6      0       <NA> (100,200]            0
# 3 P-0021159-T01-IM6      0       <NA> (200,300]            0
# 4 P-0021256-T01-IM6      0       <NA> (100,200]            0
# 5 P-0028349-T01-IM6      0       <NA> (100,200]            0

##Panel B : Effect of purity, and clonality on TMB
#mm_se1 = ddply(mm_, c("purity_bin"), summarise, N = length(TMBExome),mean = mean(TMBExome),sd   = sd(TMBExome),se   = sd / sqrt(N))
#mm_se = summarySE(mm_, measurevar="TMBExome", groupvars=c("purity_bin"))
mm_se_low = filter(ex_roslin_purity_depth0, !is.na(purity_bin), TMBWES<10) %>% summarySE(., measurevar="TMBWES", groupvars=c("purity_bin"))
mm_se_high = filter(ex_roslin_purity_depth0, !is.na(purity_bin), TMBWES>=10) %>% summarySE(., measurevar="TMBWES", groupvars=c("purity_bin"))

p3 = ggplot(mm_se_low, aes(x=purity_bin, y=TMBWES)) + theme_classic(base_size = 16) +
  geom_point(size=5, shape=21, fill="white") + ylab("TMB WES\n(<10)") +
  scale_y_continuous(limits=c(0,5),breaks = c(1,2,3,4,5)) +
  scale_x_discrete(labels = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)) +
  geom_errorbar(aes(ymin=TMBWES-ci, ymax=TMBWES+ci), colour="black", width=.1)

p4 = ggplot(mm_se_high, aes(x=purity_bin, y=TMBWES)) + theme_classic(base_size = 16) +
  geom_point(size=5, shape=21, fill="white") + ylab("TMB WES\n(>=10)") +
  scale_y_continuous(limits=c(10,100),breaks = c(10,20,40,60,80,100)) +
  scale_x_discrete(labels = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)) +
  geom_errorbar(aes(ymin=TMBWES-ci, ymax=TMBWES+ci), colour="black", width=.1)

grid.arrange(p3,p4, ncol = 2, widths = c(0.5,0.5), newpage = F)

## Clonal TMB
mm_se_low_cl = filter(ex_roslin_purity_depth0, !is.na(purity_bin), TMBWESClonal<10) %>% summarySE(., measurevar="TMBWESClonal", groupvars=c("purity_bin"))
mm_se_high_cl = filter(ex_roslin_purity_depth0, !is.na(purity_bin), TMBWESClonal>=10) %>% summarySE(., measurevar="TMBWESClonal", groupvars=c("purity_bin"))

p5 = ggplot(mm_se_low_cl, aes(x=purity_bin, y=TMBWESClonal)) + theme_classic(base_size = 16) +
  geom_point(size=5, shape=21, fill="white") + ylab("Clonal TMB WES\n(<10)") +
  scale_y_continuous(limits=c(0,5),breaks = c(1,2,3,4,5)) +
  scale_x_discrete(labels = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)) +
  geom_errorbar(aes(ymin=TMBWESClonal-ci, ymax=TMBWESClonal+ci), colour="black", width=.1)

p6 = ggplot(mm_se_high_cl, aes(x=purity_bin, y=TMBWESClonal)) + theme_classic(base_size = 16) +
  geom_point(size=5, shape=21, fill="white") + ylab("Clonal TMB WES\n(>=10)") +
  scale_y_continuous(limits=c(0,125),breaks = c(20,40,60,80,100)) +
  scale_x_discrete(labels = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)) +
  geom_errorbar(aes(ymin=TMBWESClonal-ci, ymax=TMBWESClonal+ci), colour="black", width=.1)

grid.arrange(p5,p6, ncol = 2, widths = c(0.5,0.5), newpage = F)
#########################################################################################################################################################################
#Panel C : Effect of depth on TMB

#Hist with a line
p<-ggplot(ex_roslin_purity_depth0, aes(x=TMBWES)) + theme_classic(base_size = 16) +
  geom_histogram(color="black", fill="white", bins = 100) + scale_x_continuous(breaks = c(10,50,100,150,200,250,300,350))
p7 = p + geom_vline(aes(xintercept=10),color="red", lty=2);
p7

mm_se_low = filter(ex_roslin_purity_depth0, TMBWES<10, !is.na(depth_bin)) %>% summarySE(., measurevar="TMBWES", groupvars=c("depth_bin"))
mm_se_high = filter(ex_roslin_purity_depth0, TMBWES>=10, !is.na(depth_bin)) %>% summarySE(., measurevar="TMBWES", groupvars=c("depth_bin"))

#TMB by Sequencing depth (update for Exome depth)

p9 = ggplot(mm_se_low, aes(x=depth_bin, y=TMBWES)) + theme_classic(base_size = 16) +
  geom_bar(stat="identity") + ylab("TMB WES\n(<10)") +
  #scale_y_continuous(limits=c(0,5),breaks = c(1,2,3,4,5)) +
  #scale_x_discrete(labels = c(50,100,150,200,250,300,350)) +
  geom_errorbar(aes(ymin=TMBWES-ci, ymax=TMBWES+ci), colour="black", width=.1)

p10 = ggplot(mm_se_high, aes(x=depth_bin, y=TMBWES)) + theme_classic(base_size = 16) +
  geom_bar(stat="identity") + ylab("TMB WES\n(>=10)") +
  #scale_y_continuous(limits=c(0,120),breaks = c(10,20,40,60,80,100,120)) +
  #scale_x_discrete(labels = c(50,100,150,200,250,300,350)) +
  geom_errorbar(aes(ymin=TMBWES-ci, ymax=TMBWES+ci), colour="black", width=.1)

pdf('')
grid.arrange(p9,p10, ncol = 2, widths = c(0.5,0.5), newpage = F)

pdf("~/tempo-cohort-level/Figure2CDE_Purity_Depth_TMB_Roslin.pdf", paper = "a4")
grid.arrange(p3, p4, p5, p6, p9,p10, ncol = 2, nrow =3, widths = c(0.5,0.5), newpage = F)
dev.off()
