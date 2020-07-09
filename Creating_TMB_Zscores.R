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

url <- 'docs.google.com/spreadsheets/d/1ZQCZ-02b8VNNDL05w-FdUiTpIVJ0WAW0EwxCkgX3oM0'
ex_roslin <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE) %>%
  select(DMP = DMPID, TMBWES = TMBExome, TMBIMPACT, TMBWES_IMgenes = TMBExomeIMgenes, TMBWES_NonIMgenes = TMBExome_NonIMgenes, Purity_Reviewed = PurityExome) %>% 
  filter(DMP %in% im_ex_clin$DMP)
dim(ex_roslin) #1636

##################

######### For Roslin showing TMB IMPACT vs. WES

im_clin = fread('~/tempo-cohort-level/IM_metadata_040620.txt') %>% select(DMP, TMBIMPACT)
im_ex_clin_roslin1 = inner_join(im_clin, ex_roslin, by = c('DMP')) %>% 
  select(DMP, TMBIMPACT=TMBIMPACT.x, TMBWES)

dim(im_ex_clin_roslin1) #1636
# Equation
label = lm_eqn(im_ex_clin_roslin1,y=im_ex_clin_roslin1$TMBIMPACT,x=im_ex_clin_roslin1$TMBWES)
label
# R²
r2 = summary(lm(TMBWES~TMBIMPACT,data=im_ex_clin_roslin1))$adj.r.squared
r2
#Spearman Correlation
corr = cor.test(im_ex_clin_roslin1$TMBWES, im_ex_clin_roslin1$TMBIMPACT,method = "spearman") 
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
overall


### Tranformation + zscores
library(rcompanion)
T_tuk = transformTukey(im_ex_clin_roslin_$TMBIMPACT, plotit=FALSE)
plotNormalHistogram(T_tuk)
T_tuk_wes = transformTukey(im_ex_clin_roslin_$TMBWES, plotit=FALSE)
plotNormalHistogram(T_tuk_wes)

T_tuk_z = scale(T_tuk, center = T, scale = F)
plotNormalHistogram(T_tuk_z)
T_tuk_z_wes = scale(T_tuk_wes, center = T, scale = F)
plotNormalHistogram(T_tuk_z_wes)

##Plotting them both
#####
plot(density(T_tuk_z), col = 'red',xlim=c(-0.5,0.5),ylim=c(0,6)) 
lines(density(T_tuk_z_wes),col='blue')
#####

##Plotting them both again, different ways
x <- data.frame(IMPACT=im_ex_clin_roslin_$TMBIMPACT,WES=im_ex_clin_roslin_$TMBWES)
data<- melt(x)
ggplot(data,aes(x=value, fill=variable)) + geom_density(alpha=0.25)
ggplot(data,aes(x=value, fill=variable)) + geom_histogram(alpha=0.25)

x <- data.frame(IMPACT=T_tuk_z,WES=T_tuk_z_wes)
data<- melt(x)
ggplot(data,aes(x=value, fill=variable)) + geom_density(alpha=0.25)
ggplot(data,aes(x=value, fill=variable)) + geom_histogram(alpha=0.25)
ggplot(data,aes(x=variable, y=value, fill=variable)) + geom_boxplot()
