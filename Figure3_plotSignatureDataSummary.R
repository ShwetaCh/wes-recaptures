#written by Noah Friedman
#a modified version of https://github.com/taylor-lab/mutationSignatureWaterfallPlots/blob/master/plot_signature_waterfall_plot.R


library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(data.table); setDTthreads(6)
library(ggpubr)



#PLOTS THE COMBINED WATERFALL
plot_signature_waterfall <- function(df, title, 
                                           orderingValParam, primarySigValParam, primarySigNameParam,
                                           secondPredominantSigValParam, secondPredominantSigNameParam,
                                           thirdPredominantSigValParam, thirdPredominantSigNameParam, 
                                           mutBurdenParam,
                                           plotLevels){
  
  barColorPalette = c(
    "#00DFFF", "#FF0000", "#FF1493",
    "#FF9200", '#ffb347', '#ffd394', 
    "#FFF600", "#267574","#ADFF2F", "#00dd5d",
    "#2A52BE","#551A8B","pink", "#9acd32",
    "#D3D3D3","gray"
  )
  names(barColorPalette)= plotLevels
  
  #column renaming
  df <- my.rename(df, orderingValParam, "orderingVal")
  df <- my.rename(df, primarySigValParam, "primarySigVal")
  df <- my.rename(df, primarySigNameParam, "primarySigName")
  df <- my.rename(df, secondPredominantSigValParam, "secondPredominantSigVal")
  df <- my.rename(df, secondPredominantSigNameParam, "secondPredominantSigName")
  df <- my.rename(df, thirdPredominantSigValParam, "thirdPredominantSigVal")
  df <- my.rename(df, thirdPredominantSigNameParam, "thirdPredominantSigName")
  df <- my.rename(df, mutBurdenParam, "mutBurden")
  
  #SET TITLE
  nSufficientMuts <- dim(df[df$mutBurden == "False",])[1]
  finalTitle <- paste(title,  ': ', nSufficientMuts, '/', dim(df)[1], 'sufficient mutations')
  
  plt <- ggplot(df)+
    
    #The first bar of the predominant signature in the positive direction
    geom_bar(aes(x = reorder(Tumor_Sample_Barcode, orderingVal), y=primarySigVal, 
                 fill = factor(primarySigName, levels=plotLevels)), stat="identity")+
    
    #THIRD BAR (plot it first so it is covered by the second bar)
    geom_bar(aes(x = reorder(Tumor_Sample_Barcode, orderingVal), y=-secondPredominantSigVal - thirdPredominantSigVal, 
                 fill = factor(thirdPredominantSigName, levels=plotLevels)), stat = "identity")+ #color the lower signature column by which signature it is                                                             
    
    #The second bar of the second predominant signature in the negative direction
    geom_bar(aes(x = reorder(Tumor_Sample_Barcode, orderingVal), y=-secondPredominantSigVal, 
                 fill = factor(secondPredominantSigName, levels=plotLevels)), stat = "identity")+ #color the lower signature column by which signature it is                                                                                
    
    scale_fill_manual(values=barColorPalette)+
    theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks= element_blank(), axis.line = element_blank())+
    ylab('Signature Fraction')+
    ylim(-1,1)+
    guides(fill=guide_legend(title="Signature"))+
    ggtitle(finalTitle)+
    theme(legend.position = "none")
  
  mutBurdenRug <- plot_tmb_color_tiles(df, 'TMBExome', orderingValParam)
  alignedPlot <- plot_grid(plt, mutBurdenRug, nrow=2, rel_heights = c(1,.3))
  return(alignedPlot)
}

#RENAME COLUMN FUNCTIONS
my.rename <- function(df, old.name, new.name){ #R SUCKS!!!!!! heres my renaming function cause god forbid there would be an easy or intiuitive way to do this with R
  names(df)[names(df) == old.name] <- new.name 
  return(df)
}

plot_tmb_color_tiles <- function(df, tmbColParam, orderingValParam){
  df <- my.rename(df, tmbColParam, "tmbCol")
  df <- my.rename(df, orderingValParam, "orderingVal")
  ggplot(df, aes(reorder(Tumor_Sample_Barcode, orderingVal), y=0))+
    geom_tile(aes(fill=log(tmbCol)))+
    scale_fill_viridis_c(option='magma', direction = -1, limits=c(-2,8))+
    theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())+
    theme(legend.position = "none")
}

#TODO FIX SORTING ERROR THAT OCCURS WITH NON AGING DOMINANT CANCERs

df <- read.table('/Users/friedman/Desktop/exomeProjectHelp/plottingData.tsv', sep = '\t', header=TRUE)
#plot_signature_waterfall(df[df$Cancer_Type == 'Pancreatic Cancer',], 'alexandria')


plottingLevelsExome <- c('exomeSig_Aging', 'exomeSig_AID/APOBEC', 'exomeSig_BRCA1/2', 
                    'exomeSig_Smoking', 'exomeSig_Tobacco_chewing', 'exomeSig_Aflatoxin',
                    'exomeSig_UV', 'exomeSig_MMR/MSI', 'exomeSig_POLE', 'exomeSig_Signature_14', 
                    'exomeSig_TMZ/Alkylating', 'exomeSig_Signature_17', 'exomeSig_Signature_5', 'exomeSig_Signature_23',  'other')

plottingLevelsImpact <- c('impactSig_Aging', 'impactSig_AID/APOBEC', 'impactSig_BRCA1/2', 
                          'impactSig_Smoking', 'impactSig_Tobacco_chewing', 'impactSig_Aflatoxin',
                          'impactSig_UV', 'impactSig_MMR/MSI', 'impactSig_POLE', 'impactSig_Signature_14', 
                          'impactSig_TMZ/Alkylating', 'impactSig_Signature_17', 'impactSig_Signature_5', 'impactSig_Signature_23', 'other')

#MAKE MASSIVE ALIGNED PLOT
alignedPlot <- plot_grid(
  
  #BILIARY CANCER
  plot_signature_waterfall(df[df$Cancer_Type == 'Biliary Cancer',], 'Biliary Cancer -- IMPACT', 
                           'orderingValImpact', 'signatureOfInterestMagnitudeImpact', 'signatureOfInterestNameImpact',
                           'secondPredominantSigMagnitudeImpact', 'secondPredominantSigNameImpact',
                           'thirdPredominantSigMagnitudeImpact', 'thirdPredominantSigNameImpact',
                           'insufficientImpactMutBurden', plottingLevelsImpact
  ),
  plot_signature_waterfall(df[df$Cancer_Type == 'Biliary Cancer',], 'Biliary Cancer -- EXOME', 
                           'orderingValExome', 'signatureOfInterestMagnitudeExome', 'signatureOfInterestNameExome',
                           'secondPredominantSigMagnitudeExome', 'secondPredominantSigNameExome',
                           'thirdPredominantSigMagnitudeExome', 'thirdPredominantSigNameExome',
                           'insufficientExomeMutBurden', plottingLevelsExome
  ),
  
  #BLADDER CANCER
  plot_signature_waterfall(df[df$Cancer_Type == 'Bladder Cancer',], 'Bladder Cancer -- IMPACT', 
                           'orderingValImpact', 'signatureOfInterestMagnitudeImpact', 'signatureOfInterestNameImpact',
                           'secondPredominantSigMagnitudeImpact', 'secondPredominantSigNameImpact',
                           'thirdPredominantSigMagnitudeImpact', 'thirdPredominantSigNameImpact',
                           'insufficientImpactMutBurden', plottingLevelsImpact
  ),
  plot_signature_waterfall(df[df$Cancer_Type == 'Bladder Cancer',], 'Bladder Cancer -- EXOME', 
                           'orderingValExome', 'signatureOfInterestMagnitudeExome', 'signatureOfInterestNameExome',
                           'secondPredominantSigMagnitudeExome', 'secondPredominantSigNameExome',
                           'thirdPredominantSigMagnitudeExome', 'thirdPredominantSigNameExome',
                           'insufficientExomeMutBurden', plottingLevelsExome
  ),
  
  #BREAST CARCINOMA
  plot_signature_waterfall(df[df$Cancer_Type == 'Breast Carcinoma',], 'Breast Carcinoma -- IMPACT', 
                           'orderingValImpact', 'signatureOfInterestMagnitudeImpact', 'signatureOfInterestNameImpact',
                           'secondPredominantSigMagnitudeImpact', 'secondPredominantSigNameImpact',
                           'thirdPredominantSigMagnitudeImpact', 'thirdPredominantSigNameImpact',
                           'insufficientImpactMutBurden', plottingLevelsImpact
  ),
  plot_signature_waterfall(df[df$Cancer_Type == 'Breast Carcinoma',], 'Breast Carcinoma -- EXOME', 
                           'orderingValExome', 'signatureOfInterestMagnitudeExome', 'signatureOfInterestNameExome',
                           'secondPredominantSigMagnitudeExome', 'secondPredominantSigNameExome',
                           'thirdPredominantSigMagnitudeExome', 'thirdPredominantSigNameExome',
                           'insufficientExomeMutBurden', plottingLevelsExome
  ),
  
  #COLORECTAL CANCER
  plot_signature_waterfall(df[df$Cancer_Type == 'Colorectal Cancer',], 'Colorectal Cancer -- IMPACT', 
                           'orderingValImpact', 'signatureOfInterestMagnitudeImpact', 'signatureOfInterestNameImpact',
                           'secondPredominantSigMagnitudeImpact', 'secondPredominantSigNameImpact',
                           'thirdPredominantSigMagnitudeImpact', 'thirdPredominantSigNameImpact', 
                           'insufficientImpactMutBurden', plottingLevelsImpact
  ),
  plot_signature_waterfall(df[df$Cancer_Type == 'Colorectal Cancer',], 'Colorectal Cancer -- EXOME', 
                           'orderingValExome', 'signatureOfInterestMagnitudeExome', 'signatureOfInterestNameExome',
                           'secondPredominantSigMagnitudeExome', 'secondPredominantSigNameExome',
                           'thirdPredominantSigMagnitudeExome', 'thirdPredominantSigNameExome',
                           'insufficientExomeMutBurden', plottingLevelsExome
  ),
  
  #ENDOMETRIAL CANCER
  plot_signature_waterfall(df[df$Cancer_Type == 'Endometrial Cancer',], 'Endometrial Cancer -- IMPACT', 
                           'orderingValImpact', 'signatureOfInterestMagnitudeImpact', 'signatureOfInterestNameImpact',
                           'secondPredominantSigMagnitudeImpact', 'secondPredominantSigNameImpact',
                           'thirdPredominantSigMagnitudeImpact', 'thirdPredominantSigNameImpact', 
                           'insufficientImpactMutBurden', plottingLevelsImpact
  ),
  plot_signature_waterfall(df[df$Cancer_Type == 'Endometrial Cancer',], 'Endometrial Cancer -- EXOME', 
                           'orderingValExome', 'signatureOfInterestMagnitudeExome', 'signatureOfInterestNameExome',
                           'secondPredominantSigMagnitudeExome', 'secondPredominantSigNameExome',
                           'thirdPredominantSigMagnitudeExome', 'thirdPredominantSigNameExome',
                           'insufficientExomeMutBurden', plottingLevelsExome
  ),
  
  
  #Germ Cell Tumor
  plot_signature_waterfall(df[df$Cancer_Type == 'Germ Cell Tumor',], 'Germ Cell Tumor -- IMPACT', 
                           'orderingValImpact', 'signatureOfInterestMagnitudeImpact', 'signatureOfInterestNameImpact',
                           'secondPredominantSigMagnitudeImpact', 'secondPredominantSigNameImpact',
                           'thirdPredominantSigMagnitudeImpact', 'thirdPredominantSigNameImpact',
                           'insufficientImpactMutBurden', plottingLevelsImpact
  ),
  plot_signature_waterfall(df[df$Cancer_Type == 'Germ Cell Tumor',], 'Germ Cell Tumor -- EXOME', 
                           'orderingValExome', 'signatureOfInterestMagnitudeExome', 'signatureOfInterestNameExome',
                           'secondPredominantSigMagnitudeExome', 'secondPredominantSigNameExome',
                           'thirdPredominantSigMagnitudeExome', 'thirdPredominantSigNameExome',
                           'insufficientExomeMutBurden', plottingLevelsExome
  ),
  
  #Glioma
  plot_signature_waterfall(df[df$Cancer_Type == 'Glioma',], 'Glioma -- IMPACT', 
                           'orderingValImpact', 'signatureOfInterestMagnitudeImpact', 'signatureOfInterestNameImpact',
                           'secondPredominantSigMagnitudeImpact', 'secondPredominantSigNameImpact',
                           'thirdPredominantSigMagnitudeImpact', 'thirdPredominantSigNameImpact',
                           'insufficientImpactMutBurden', plottingLevelsImpact
  ),
  plot_signature_waterfall(df[df$Cancer_Type == 'Glioma',], 'Glioma -- EXOME', 
                           'orderingValExome', 'signatureOfInterestMagnitudeExome', 'signatureOfInterestNameExome',
                           'secondPredominantSigMagnitudeExome', 'secondPredominantSigNameExome',
                           'thirdPredominantSigMagnitudeExome', 'thirdPredominantSigNameExome',
                           'insufficientExomeMutBurden', plottingLevelsExome
  ),
  
  #Head and Neck Carcinoma
  plot_signature_waterfall(df[df$Cancer_Type == 'Head and Neck Carcinoma',], 'Head and Neck Carcinoma -- IMPACT', 
                           'orderingValImpact', 'signatureOfInterestMagnitudeImpact', 'signatureOfInterestNameImpact',
                           'secondPredominantSigMagnitudeImpact', 'secondPredominantSigNameImpact',
                           'thirdPredominantSigMagnitudeImpact', 'thirdPredominantSigNameImpact',
                           'insufficientImpactMutBurden', plottingLevelsImpact
  ),
  plot_signature_waterfall(df[df$Cancer_Type == 'Head and Neck Carcinoma',], 'Head and Neck Carcinoma -- EXOME', 
                           'orderingValExome', 'signatureOfInterestMagnitudeExome', 'signatureOfInterestNameExome',
                           'secondPredominantSigMagnitudeExome', 'secondPredominantSigNameExome',
                           'thirdPredominantSigMagnitudeExome', 'thirdPredominantSigNameExome',
                           'insufficientExomeMutBurden', plottingLevelsExome
  ),
  
  #MELANOMA
  plot_signature_waterfall(df[df$Cancer_Type == 'Melanoma',], 'Melanoma -- IMPACT', 
                           'orderingValImpact', 'signatureOfInterestMagnitudeImpact', 'signatureOfInterestNameImpact',
                           'secondPredominantSigMagnitudeImpact', 'secondPredominantSigNameImpact',
                           'thirdPredominantSigMagnitudeImpact', 'thirdPredominantSigNameImpact',
                           'insufficientImpactMutBurden', plottingLevelsImpact
  ),
  plot_signature_waterfall(df[df$Cancer_Type == 'Melanoma',], 'Melanoma -- EXOME', 
                           'orderingValExome', 'signatureOfInterestMagnitudeExome', 'signatureOfInterestNameExome',
                           'secondPredominantSigMagnitudeExome', 'secondPredominantSigNameExome',
                           'thirdPredominantSigMagnitudeExome', 'thirdPredominantSigNameExome',
                           'insufficientExomeMutBurden', plottingLevelsExome
  ),
  
  
  #Non-Small Cell Lung Cancer
  plot_signature_waterfall(df[df$Cancer_Type == 'Non-Small Cell Lung Cancer',], 'Non-Small Cell Lung Cancer -- IMPACT', 
                           'orderingValImpact', 'signatureOfInterestMagnitudeImpact', 'signatureOfInterestNameImpact',
                           'secondPredominantSigMagnitudeImpact', 'secondPredominantSigNameImpact',
                           'thirdPredominantSigMagnitudeImpact', 'thirdPredominantSigNameImpact',
                           'insufficientImpactMutBurden', plottingLevelsImpact
  ),
  plot_signature_waterfall(df[df$Cancer_Type == 'Non-Small Cell Lung Cancer',], 'Non-Small Cell Lung Cancer -- EXOME', 
                           'orderingValExome', 'signatureOfInterestMagnitudeExome', 'signatureOfInterestNameExome',
                           'secondPredominantSigMagnitudeExome', 'secondPredominantSigNameExome',
                           'thirdPredominantSigMagnitudeExome', 'thirdPredominantSigNameExome',
                           'insufficientExomeMutBurden', plottingLevelsExome
  ),
  
  #Small Cell Lung Cancer
  plot_signature_waterfall(df[df$Cancer_Type == 'Small Cell Lung Cancer',], 'Small Cell Lung Cancer -- IMPACT', 
                           'orderingValImpact', 'signatureOfInterestMagnitudeImpact', 'signatureOfInterestNameImpact',
                           'secondPredominantSigMagnitudeImpact', 'secondPredominantSigNameImpact',
                           'thirdPredominantSigMagnitudeImpact', 'thirdPredominantSigNameImpact',
                           'insufficientImpactMutBurden', plottingLevelsImpact
  ),
  plot_signature_waterfall(df[df$Cancer_Type == 'Small Cell Lung Cancer',], 'Small Cell Lung Cancer -- EXOME', 
                           'orderingValExome', 'signatureOfInterestMagnitudeExome', 'signatureOfInterestNameExome',
                           'secondPredominantSigMagnitudeExome', 'secondPredominantSigNameExome',
                           'thirdPredominantSigMagnitudeExome', 'thirdPredominantSigNameExome',
                           'insufficientExomeMutBurden', plottingLevelsExome
  ),
  
  #OVARIAN CANCER
  plot_signature_waterfall(df[df$Cancer_Type == 'Ovarian Cancer',], 'Ovarian Cancer -- IMPACT', 
                           'orderingValImpact', 'signatureOfInterestMagnitudeImpact', 'signatureOfInterestNameImpact',
                           'secondPredominantSigMagnitudeImpact', 'secondPredominantSigNameImpact',
                           'thirdPredominantSigMagnitudeImpact', 'thirdPredominantSigNameImpact',
                           'insufficientImpactMutBurden', plottingLevelsImpact
  ),
  plot_signature_waterfall(df[df$Cancer_Type == 'Ovarian Cancer',], 'Ovarian Cancer -- EXOME', 
                           'orderingValExome', 'signatureOfInterestMagnitudeExome', 'signatureOfInterestNameExome',
                           'secondPredominantSigMagnitudeExome', 'secondPredominantSigNameExome',
                           'thirdPredominantSigMagnitudeExome', 'thirdPredominantSigNameExome',
                           'insufficientExomeMutBurden', plottingLevelsExome
  ),
  
  #PANCREATIC
  plot_signature_waterfall(df[df$Cancer_Type == 'Pancreatic Cancer',], 'Pancreatic Cancer -- IMPACT', 
                           'orderingValImpact', 'signatureOfInterestMagnitudeImpact', 'signatureOfInterestNameImpact',
                           'secondPredominantSigMagnitudeImpact', 'secondPredominantSigNameImpact',
                           'thirdPredominantSigMagnitudeImpact', 'thirdPredominantSigNameImpact',
                           'insufficientImpactMutBurden', plottingLevelsImpact
  ),
  plot_signature_waterfall(df[df$Cancer_Type == 'Pancreatic Cancer',], 'Pancreatic Cancer -- EXOME', 
                           'orderingValExome', 'signatureOfInterestMagnitudeExome', 'signatureOfInterestNameExome',
                           'secondPredominantSigMagnitudeExome', 'secondPredominantSigNameExome',
                           'thirdPredominantSigMagnitudeExome', 'thirdPredominantSigNameExome',
                           'insufficientExomeMutBurden',plottingLevelsExome
  ),
  
  #PROSTATE
  plot_signature_waterfall(df[df$Cancer_Type == 'Prostate Cancer',], 'Prostate Cancer -- IMPACT', 
                           'orderingValImpact', 'signatureOfInterestMagnitudeImpact', 'signatureOfInterestNameImpact',
                           'secondPredominantSigMagnitudeImpact', 'secondPredominantSigNameImpact',
                           'thirdPredominantSigMagnitudeImpact', 'thirdPredominantSigNameImpact', 
                           'insufficientImpactMutBurden', plottingLevelsImpact
  ),
  plot_signature_waterfall(df[df$Cancer_Type == 'Prostate Cancer',], 'Prostate Cancer -- EXOME', 
                           'orderingValExome', 'signatureOfInterestMagnitudeExome', 'signatureOfInterestNameExome',
                           'secondPredominantSigMagnitudeExome', 'secondPredominantSigNameExome',
                           'thirdPredominantSigMagnitudeExome', 'thirdPredominantSigNameExome',
                           'insufficientExomeMutBurden',plottingLevelsExome
  ),
  
  
  #RENAL CELL
  plot_signature_waterfall(df[df$Cancer_Type == 'Renal Cell Carcinoma',], 'Renal Cell Carcinoma -- IMPACT', 
                           'orderingValImpact', 'signatureOfInterestMagnitudeImpact', 'signatureOfInterestNameImpact',
                           'secondPredominantSigMagnitudeImpact', 'secondPredominantSigNameImpact',
                           'thirdPredominantSigMagnitudeImpact', 'thirdPredominantSigNameImpact', 
                           'insufficientImpactMutBurden', plottingLevelsImpact
  ),
  plot_signature_waterfall(df[df$Cancer_Type == 'Renal Cell Carcinoma',], 'Renal Cell Carcinoma -- EXOME', 
                           'orderingValExome', 'signatureOfInterestMagnitudeExome', 'signatureOfInterestNameExome',
                           'secondPredominantSigMagnitudeExome', 'secondPredominantSigNameExome',
                           'thirdPredominantSigMagnitudeExome', 'thirdPredominantSigNameExome',
                           'insufficientExomeMutBurden',plottingLevelsExome
  ),
  
  #SOFT TISSUE SARCOMA
  plot_signature_waterfall(df[df$Cancer_Type == 'Soft Tissue Sarcoma',], 'Soft Tissue Sarcoma -- IMPACT', 
                           'orderingValImpact', 'signatureOfInterestMagnitudeImpact', 'signatureOfInterestNameImpact',
                           'secondPredominantSigMagnitudeImpact', 'secondPredominantSigNameImpact',
                           'thirdPredominantSigMagnitudeImpact', 'thirdPredominantSigNameImpact', 
                           'insufficientImpactMutBurden', plottingLevelsImpact
  ),
  plot_signature_waterfall(df[df$Cancer_Type == 'Soft Tissue Sarcoma',], 'Soft Tissue Sarcoma -- EXOME', 
                           'orderingValExome', 'signatureOfInterestMagnitudeExome', 'signatureOfInterestNameExome',
                           'secondPredominantSigMagnitudeExome', 'secondPredominantSigNameExome',
                           'thirdPredominantSigMagnitudeExome', 'thirdPredominantSigNameExome',
                           'insufficientExomeMutBurden', plottingLevelsExome
  ),
  
  #UNKNOWN
  plot_signature_waterfall(df[df$Cancer_Type == 'Cancer of Unknown Primary',], 'Cancer of Unknown Primary -- IMPACT', 
                           'orderingValImpact', 'signatureOfInterestMagnitudeImpact', 'signatureOfInterestNameImpact',
                           'secondPredominantSigMagnitudeImpact', 'secondPredominantSigNameImpact',
                           'thirdPredominantSigMagnitudeImpact', 'thirdPredominantSigNameImpact',
                           'insufficientImpactMutBurden', plottingLevelsImpact
  ),
  plot_signature_waterfall(df[df$Cancer_Type == 'Cancer of Unknown Primary',], 'Cancer of Unknown Primary -- EXOME', 
                           'orderingValExome', 'signatureOfInterestMagnitudeExome', 'signatureOfInterestNameExome',
                           'secondPredominantSigMagnitudeExome', 'secondPredominantSigNameExome',
                           'thirdPredominantSigMagnitudeExome', 'thirdPredominantSigNameExome',
                           'insufficientExomeMutBurden', plottingLevelsExome
  ),
  
  #Other
  plot_signature_waterfall(df[df$Cancer_Type == 'Other',], 'Other -- IMPACT', 
                           'orderingValImpact', 'signatureOfInterestMagnitudeImpact', 'signatureOfInterestNameImpact',
                           'secondPredominantSigMagnitudeImpact', 'secondPredominantSigNameImpact',
                           'thirdPredominantSigMagnitudeImpact', 'thirdPredominantSigNameImpact', 
                           'insufficientImpactMutBurden',plottingLevelsImpact
  ),
  plot_signature_waterfall(df[df$Cancer_Type == 'Other',], 'Other -- EXOME', 
                           'orderingValExome', 'signatureOfInterestMagnitudeExome', 'signatureOfInterestNameExome',
                           'secondPredominantSigMagnitudeExome', 'secondPredominantSigNameExome',
                           'thirdPredominantSigMagnitudeExome', 'thirdPredominantSigNameExome',
                           'insufficientExomeMutBurden',plottingLevelsExome
  ),
  nrow = 18, ncol=2
)


#########DUMMY CODE TO CREATE LEGENDS


#CREATE SIGNATURES LEGEND
barColorPalette = c(
  "#00DFFF", "#FF0000", "#FF1493",
  "#FF9200", '#ffb347', '#ffd394', 
  "#FFF600", "#267574","#ADFF2F", "#00dd5d",
  "#2A52BE","#551A8B","pink", "#9acd32",
  "#D3D3D3","gray"
)
plottingLevelsNeutral  <- c('Aging', 'AID/APOBEC', 'BRCA1/2', 
                            'Smoking', 'Tobacco_chewing', 'Aflatoxin',
                            'UV', 'MMR/MSI', 'POLE', 'Signature_14', 
                            'TMZ/Alkylating', 'Signature_17', 'Signature_5', 'Signature_23', 'other')

names(barColorPalette)= plottingLevelsNeutral
signaturesLegend <-
  get_legend(ggplot(df, aes(x=Tumor_Sample_Barcode, y=1, fill=factor(secondPredominantSigNameExome, levels=plottingLevelsNeutral)))+
  geom_bar(stat='identity')+
  scale_fill_manual(values=barColorPalette, drop=FALSE)+
  theme(legend.title = element_text(size = 30))+
  guides(fill=guide_legend(title="Signature")))

########LEGEND FOR TMB
tmbLegend <- get_legend(ggplot(df, aes(x=Tumor_Sample_Barcode, y=0))+
  geom_tile(aes(fill=log(TMBExome)))+
  theme(legend.title = element_text(size = 30))+
  scale_fill_viridis_c(option='magma', direction = -1, limits=c(-2,8)))


##PASTE ALL PLOTS TOGETHER
legendsGrid <- plot_grid(ggplot(), signaturesLegend, tmbLegend, ggplot(), nrow = 4)

alignedPlotWithLegend <- plot_grid(alignedPlot,legendsGrid, ncol = 2, rel_widths = c(1,.3))

alignedPlotWithLegendAndTitle <- plot_grid(ggplot()+ggtitle('SIGNATURES IN IMPACT AND MATCHED EXOME')+theme(plot.title = element_text(size = 40, face = "bold"))
                                           , alignedPlotWithLegend, nrow=2, rel_heights = c(.1,1))


#TODO FIX MARGINS AND TEXT SIZE
ggsave('~/Desktop/plot.pdf', plot=alignedPlotWithLegendAndTitle,  width = 20, height = 40, units = c("in"), limitsize = FALSE)


#
#
#
#
#######OTHER AUXILIARY PLOTS

df <- read.table('/Users/friedman/Desktop/exomeProjectHelp/plottingData.tsv', sep = '\t', header=TRUE)

df$exomeSig_BRCA1.2
df$TMBExome
df$exomeSig_BRCA1.2
#TMB in BRCA related cases 
ggplot(df[df$DominantSignatureExome == 'exomeSig_BRCA1/2',], aes(x = isHRDCancer, y = exomeConfidence_BRCA1.2))+
  geom_boxplot(fatten = NULL)+ #TODO FIX THIS TO PLOT THE MIDDLE LINE AS THE MEAN!!!
  stat_summary(fun.y = median, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = 0.75, size = 1, linetype = "solid")+
  labs(caption = 'format_sigs_for_exome_proj.ipynb  plotSignatureDataSummary.R')+
  stat_compare_means(label.y=1.2)+
  ggtitle('Confidence differences in BRCA dominant cases')

######################################

#MSI/AGING SIGNATURES

ggplot(df[df$DominantSignatureExome == 'exomeSig_MMR/MSI' | df$DominantSignatureExome == 'exomeSig_Aging',],
       aes(x=DominantSignatureExome, y=MSIExome))+
       geom_boxplot(fatten = NULL)+ #TODO FIX THIS TO PLOT THE MIDDLE LINE AS THE MEAN!!!
       stat_summary(fun.y = median, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = 0.75, size = 1, linetype = "solid")+
       stat_compare_means()+
       labs(caption = 'format_sigs_for_exome_proj.ipynb  plotSignatureDataSummary.R')

ggplot(df[df$DominantSignatureImpact == 'impactSig_MMR/MSI' | df$DominantSignatureImpact == 'impactSig_Aging',],
       aes(x=MSIExomeClass, y=TMBIMPACT))+
  geom_boxplot(fatten = NULL)+ #TODO FIX THIS TO PLOT THE MIDDLE LINE AS THE MEAN!!!
  stat_summary(fun.y = median, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = 0.75, size = 1, linetype = "solid")+
  labs(caption = 'format_sigs_for_exome_proj.ipynb  plotSignatureDataSummary.R')



##
##
##
##
##
df <- read.table('/Users/friedman/Desktop/exomeProjectHelp/rocPlotting.tsv', sep = '\t', header=TRUE)
p <- ggplot(df, aes(x = 1 - Specificity, y=Sensitivity, group=signature))+
  geom_line(aes(color=signature))+
  scale_color_manual(values = c("#00DFFF", "red", "#ff4d00", "#FF1493", "#267574", "orange", "yellow"))

ggsave('~/Desktop/plot.pdf', plot=p,  width = 5, height = 3, units = c("in"), limitsize = FALSE)



#####
########
#############
###################
########################
###############################
#########################
####################
###############
#########
######
#
#BRCA plots SECTION

library(ggridges)
df <- read.table('/Users/friedman/Desktop/exomeProjectHelp/brcaSigAnlysis.tsv', sep = '\t', header=TRUE)
df$brcaPresentImpact = ifelse(df$impactConfidence_BRCA1.2 >.9, "Present", "Absent")
df$brcaDominantImpact = ifelse(df$DominantSignatureImpact == 'impactSig_BRCA1/2', "Present", "Absent")
df$sigExistsImpact = ifelse(df$impactSig_BRCA1.2 > .1, "Present", "Absent")

emptyTheme <- theme(axis.line = element_blank(),
                    #axis.text = element_blank(),
                    axis.ticks.y = element_blank(),
                    #axis.title.y = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank())


p1 <- ggplot(df[df$exomeConfidence_BRCA1.2 > .9,])+
  stat_bin(aes(x=TMBExome, y=..count../sum(..count..), fill=brcaPresentImpact), bins=15)+
  ylab('Fraction of all cases')+
  ggtitle('Method 1: Present if:\nIMPACT Sig3 confidence > 90%')+
  theme(legend.position = "none")+
  theme(plot.title=element_text(size=10))+
  emptyTheme
#method1: Impact Confidence Sig3 > 90%

p2 <- ggplot(df[df$exomeConfidence_BRCA1.2 > .9,])+
  stat_bin(aes(x=TMBExome, y=..count../sum(..count..), fill=brcaDominantImpact), bins=15)+
  ylab('Fraction of all cases')+
  ggtitle('Method 2: Present if:\nIMPACT Sig3 is dominant signature')+
  theme(legend.position = "none")+
  ylab('.')+
  xlab('.')+
  theme(plot.title=element_text(size=10))+
  emptyTheme
#Impact Sig3 is Dominant

p3 <- ggplot(df[df$exomeConfidence_BRCA1.2 > .9,])+
  stat_bin(aes(x=TMBExome, y=..count../sum(..count..), fill=sigExistsImpact), bins=15)+
  ylab('Fraction of all cases')+
  ggtitle('Method 3: Present if:\nIMPACT Sig3 fraction >10%')+
  ylab('.')+
  xlab('.')+
  labs(fill='Signature detected\nby IMPACT')+
  theme(plot.title=element_text(size=10))+
  emptyTheme
#Impact Sig3 > 10%

emptyTheme <- theme(axis.line = element_blank(),
                    axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    axis.title = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank())

legend <- get_legend(p3)
#fullLegend <- plot_grid(legend, textGuide, nrow=2, rel_heights = c(.75,1))
p3 <- p3 + theme(legend.position = "none")
alignedPlot <- plot_grid(p1, p2, p3, legend, ncol = 4, rel_widths = c(1,1,1,.5))
alignedPlotWithTitle <- plot_grid(ggplot() + ggtitle('Inferability of BRCA Signature from IMPACT Sequencing')+
                                    theme(plot.title = element_text(hjust=0.5, face='bold')),
                                  alignedPlot, rel_heights = c(.1, 1), nrow=2)
ggsave('~/Desktop/plot.pdf', plot=alignedPlotWithTitle,  width = 10.5, height = 3, units = c("in"), limitsize = FALSE)

hrdCancers <- c('Breast Carcinoma', 'Pancreatic Cancer', 'Ovarian Cancer', 'Prostate Cancer')

#LABELS a Column with HRD cancer/not hrd cancer
myfxn <- function(var){
  if(var %in% hrdCancers){
    return(hrdCancers[[match(var, hrdCancers)]])
  }
  else{
    return('_Other')
  }
}
df$cancerName <- sapply(df$Cancer_Type, myfxn)
p1 <- ggplot(df[df$DominantSignatureExome == 'exomeSig_BRCA1/2',], aes(x=TMBExome, y=cancerName, fill=isHRDCancer))+
  geom_density_ridges()+
  theme(legend.position = "none")+
  ggtitle('TMB in BRCA Sig Dominant Cases')

p2 <- ggplot(df[df$DominantSignatureExome == 'exomeSig_BRCA1/2',], aes(x=exomeSig_BRCA1.2, y=cancerName, fill=isHRDCancer))+
  geom_density_ridges()+
  theme(legend.position = "none")+
  ggtitle('BRCA Sig Decomposition')


p1Impact <- ggplot(df[df$DominantSignatureImpact == 'impactSig_BRCA1/2',], aes(x=TMBIMPACT, y=cancerName, fill=isHRDCancer))+
  geom_density_ridges()+
  theme(legend.position = "none")+
  ggtitle('TMB in BRCA Sig Dominant Cases')

p2Impact <- ggplot(df[df$DominantSignatureImpact == 'impactSig_BRCA1/2',], aes(x=impactSig_BRCA1.2, y=cancerName, fill=isHRDCancer))+
  geom_density_ridges()+
  theme(legend.position = "none")+
  ggtitle('BRCA Sig Decomposition')

#p3 <- ggplot(df[df$DominantSignatu == 'exomeSig_BRCA1/2',], aes(x=TMBExome, y=cancerName, fill=isHRDCancer))

alignedP <- plot_grid(p1Impact, p2Impact, p1, p2, nrow=2, ncol=2)
alignedPWithCaption <- plot_grid(alignedP, ggplot() + labs(caption='plotSignatureDataSummary.R'), nrow=2, rel_heights = c(1,.1))

ggsave('~/Desktop/plot.pdf', plot=alignedPWithCaption,  width = 9, height = 8, units = c("in"), limitsize = FALSE)


#OTHER WAY TO DO THIS
ggplot()+
  geom_density(data=df[df$Cancer_Type == 'Ovarian Cancer' & df$DominantSignatureExome == 'exomeSig_BRCA1/2',], aes(x=exomeSig_BRCA1.2*TMBExome, group='Cancer_Type', colour='Cancer_Type'))+
  geom_density(data=df[df$Cancer_Type == 'Ovarian Cancer' & df$DominantSignatureExome == 'exomeSig_BRCA1/2',], aes(x=exomeSig_BRCA1.2*TMBExome), colour='red')+
  geom_density(data=df[df$Cancer_Type == 'Breast Carcinoma' & df$DominantSignatureExome == 'exomeSig_BRCA1/2',], aes(x=exomeSig_BRCA1.2*TMBExome), colour='pink')+
  geom_density(data=df[df$Cancer_Type == 'Prostate Cancer' & df$DominantSignatureExome == 'exomeSig_BRCA1/2',], aes(x=exomeSig_BRCA1.2*TMBExome), colour='green')+
  geom_density(data=df[df$Cancer_Type == 'Pancreatic Cancer' & df$DominantSignatureExome == 'exomeSig_BRCA1/2',], aes(x=exomeSig_BRCA1.2*TMBExome), colour='purple')+
  geom_density(data=df[df$isHRDCancer == 'False' & df$DominantSignatureExome == 'exomeSig_BRCA1/2',], aes(x=exomeSig_BRCA1.2*TMBExome), colour='black')+
  ggtitle('Distribution of TMB by ')

#
#
#
#
#
##############
#

plot_sig_concordance <- function(df, title, lineColor){
  p1 <- ggplot(df)+
    stat_summary_bin(aes(x = TMBExome, y=dominantSignaturesAgree), bins=10, colour=lineColor, size=1.25)+
    ggtitle(title)+
    ylab('Fraction Concordant')+
    scale_x_log10(limits=c(0.1,500))+
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.line.x = element_blank())+
    ylim(0,1)+
    theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank())
  p2 <- ggplot(df)+
    stat_bin(aes(x = TMBExome, y = ..count..), bins=10)+
    scale_x_log10(limits=c(0.1,500))+
    ylab('N Cases')+
    scale_y_log10(limits= c(1,300))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  alignedPlot <- plot_grid(p1, p2, nrow=2, rel_heights = c(1,.5))
  return(alignedPlot)
}

df <- read.table('/Users/friedman/Desktop/exomeProjectHelp/exomeImpactConcordance.tsv', sep = '\t', header=TRUE)
allSigs <- ggplot(df)+
  stat_summary_bin(aes(x = TMBExome, y=dominantSignaturesAgree), bins=10, size=1.25)+
  scale_x_log10(limits=c(0.1,500))+
  ggtitle('Concordance between IMPACT and Exome dominant signature (ALL SIGS)')+
  ylab('Fraction of Cases that Agree')+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

alignedPlot <- plot_grid(
  plot_sig_concordance(df[df$DominantSignatureExome == 'exomeSig_Aging',], 'Aging', "#00DFFF"),
  plot_sig_concordance(df[df$DominantSignatureExome == 'exomeSig_AID/APOBEC',], 'APOBEC', "#FF0000"),
  plot_sig_concordance(df[df$DominantSignatureExome == 'exomeSig_BRCA1/2',],'BRCA1/2', "#FF1493"),
  plot_sig_concordance(df[df$DominantSignatureExome == 'exomeSig_Smoking',], 'Smoking', "#FF9200"),
  plot_sig_concordance(df[df$DominantSignatureExome == 'exomeSig_Signature_5',], 'Signature 5', "pink"),
  plot_sig_concordance(df[df$DominantSignatureExome == 'exomeSig_MMR/MSI',], 'MMR/MSI', "#267574"),
  plot_sig_concordance(df[df$DominantSignatureExome == 'exomeSig_UV',], 'UV', "#FFF600"),
  plot_sig_concordance(df[df$DominantSignatureExome == 'exomeSig_TMZ/Alkylating',], 'TMZ/Alkylating', "#2A52BE"),
  nrow = 4, ncol=2
  )

finalPlot <- plot_grid(allSigs, alignedPlot, ggplot() + labs(caption ='plotSignatureDataSummary.R   inferability_of_signatures_impact_vs_exome.ipynb'),
                       nrow=3, rel_heights = c(.3,1,.05))

ggsave('~/Desktop/plot.pdf', plot=finalPlot,  width = 10, height = 20, units = c("in"), limitsize = FALSE)

############MMR/aging signature

df <- read.table('/Users/friedman/Desktop/exomeProjectHelp/exomeImpactConcordance.tsv', sep = '\t', header=TRUE)

p1 <- ggplot(df, aes(x = MSIExome))+
  geom_smooth(aes(y = exomeSig_Aging*TMBExome, colour='Aging Sig'))+
  geom_smooth(aes(y = exomeSig_MMR.MSI*TMBExome, colour='MMR Sig'))+
  geom_smooth(aes(y = exomeNonMMRorSig1Magnitude*TMBExome, colour='Non MMR/Aging Sigs'))+
  ggtitle('Exome')

p2 <- ggplot(df, aes(x = MSIIMPACT))+
  geom_smooth(aes(y = impactSig_Aging*TMBIMPACT, colour='Aging Sig'))+
  geom_smooth(aes(y = impactSig_MMR.MSI*TMBIMPACT, colour='MMR Sig'))+
  geom_smooth(aes(y = impactNonMMRorSig1Magnitude*TMBIMPACT, colour='Non MMR/Aging Sigs'))+
  ggtitle('IMPACT')

alignedPlots <- plot_grid(p1, p2, ncol=2)
fullPlot <- plot_grid(ggplot()+ ggtitle('Aging/MMR Sig Co-occurence NMUT ATTRIBUTED TO PROCESS'),alignedPlots, ggplot()+ labs(caption='plotSignatureDataSummary.R, inferability_of_signature.ipynb'),
                      nrow=3, rel_heights = c(.1,1,.1))

ggsave('~/Desktop/plot.pdf', plot=fullPlot,  width = 10, height = 5, units = c("in"), limitsize = FALSE)








#TODO MOVE
#TODO MOVE
#TODO MOVE
#TODO MOVE
#TODO MOVE
#TODO MOVE
#TODO MOVE
#TODO MOVE
#TODO MOVE
#TODO MOVE
#TODO MOVE
#TODO MOVE
#TODO MOVE
#TODO MOVE
#TODO MOVE
#TODO MOVE
#TODO MOVE
#TODO MOVE
#TODO MOVE
#TODO MOVE
#TODO MOVE

###############################################
#TODO move this elsewhere!!
plot_obs_expected <- function(df, title){
  plt <- ggplot(df, aes(x = Nmut))+
    stat_summary_bin(aes(y=nHotspotExpected), colour='red', alpha=0.5)+
    stat_summary_bin(aes(y=nHotspotObserved), colour='blue', alpha=0.5)+
    scale_x_log10()+
    ggtitle(title)+
    ylim(0,15)
  return(plt)
}

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/obsVsExpectedHotspots.tsv', sep = '\t', header=TRUE)
alignedPlot <- plot_grid(
          plot_obs_expected(df[df$cancer_type == 'Endometrial Cancer',], 'Endometrial'),
          plot_obs_expected(df[df$cancer_type == 'Colorectal Cancer',], 'Colorectal'),
          plot_obs_expected(df[df$cancer_type == 'Glioma',], 'Glioma'), ncol=3)

fullPlot <- plot_grid(
  ggplot()+ ggtitle('Observed Vs Expected Hotspot Mutation Burden '),
  alignedPlot,
  ggplot()+labs(caption='R--currently plotSignatureDataSummary, oncogenic_mut_prob_simulation_by_gene.ipynb'),
  nrow=3,
  rel_heights = c(.1,1,.1)
)

ggsave('~/Desktop/plot.pdf', plot=fullPlot,  width = 15, height = 5, units = c("in"), limitsize = FALSE)

######################
#######
######
####














#TODO NOAH move this shit to another area





df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/obsVsExpectedHotspots.tsv', sep = '\t', header=TRUE)

ggplot(df, aes(x = Nmut))+
  stat_summary_bin(aes(y=nExpectedOncogenic), colour='red', alpha=0.5)+
  stat_summary_bin(aes(y=nOncogenicObserved), colour='blue', alpha=0.5)+
  scale_x_log10()+
  ggtitle('')+
  ylim(0,55)


######
###
##

barColorPalette = c(
  "#0b6623", "#67b826",
  "#A4DBE8", "#3255A4",
  "#ffb347", '#CC5500',
  'red'
)
plottingLevelsNeutral  <- c('colorectalHyper', 'colorectalNormal', 
                            'endometrialNormal', 'endometrialHyper',
                            'gliomaNormal', 'gliomaHyper',
                            'theoreticalSusceptibility')

names(barColorPalette)= plottingLevelsNeutral

df <- read.table('/Users/friedman/Desktop/WORK/hotspotFractionOfMutationBurden.tsv', sep = '\t', header=TRUE)
p1 <- ggplot(df, aes(group = class, x=hotspotFrac, colour=class))+
  scale_x_log10()+
  geom_density(size=2)+
  scale_color_manual(values = barColorPalette)

p2 <- ggplot(df, aes(group = simplifiedClass, x=hotspotFrac, colour=simplifiedClass))+
  scale_x_log10()+
  geom_density(size=2)

alignedPlot <- plot_grid(p2, p1, ncol=2)
finalPlot <- plot_grid(ggplot() + ggtitle('Hypothetical and observed hotspot mutation burdens'),
                       alignedPlot,
                       ggplot() + labs(caption='plotSignatureDataSummary.R--temp, oncogenic_mut_prob_simulation_by_gene.ipynb'),
                       nrow=3, rel_heights = c(.1,1,.1)
                       )

ggsave('~/Desktop/plot.pdf', plot=finalPlot,  width = 12, height = 5, units = c("in"), limitsize = FALSE)











####################
##################
###############
############
########
######
####
##
#
#PLOT a simplistic version of dominant signatures for a figure

emptyTheme <- theme(axis.line = element_blank(),
                    axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    axis.title = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank())

plot_two_bars_and_title_for_cancer_type <- function(df, ctype, colorPal, hideLegend = TRUE){
  nCases <- dim(df)[1]
  title <- paste(ctype, '  n = ', nCases)
  
  exomeBar <- ggplot(df)+
    geom_bar(aes(x=1, y=.5, fill=dominantSigExomeAdj), stat='identity')+
    scale_fill_manual(values=colorPal, drop=FALSE)+
    emptyTheme+
    coord_flip()+
    ggtitle(' ')
  if(hideLegend == TRUE){
    exomeBar <- exomeBar + theme(legend.position = 'None')
  }
  impactBar <- ggplot(df)+
    geom_bar(aes(x=1, y=1, fill=dominantSigImpactAdj), stat='identity')+
    scale_fill_manual(values=colorPal, drop=FALSE)+
    emptyTheme+
    coord_flip()+
    ggtitle(title)
  if(hideLegend == TRUE){
    impactBar <- impactBar + theme(legend.position = 'None')
  }
  
  barPlots <- plot_grid(impactBar, exomeBar, ncol=2)
  fullPlot <- plot_grid(barPlots)
  
  return(fullPlot)
}

make_full_impact_vs_exome_plot <- function(df){
  
  #SET UP COLOR INFO
  barColorPalette = c(
    "#00DFFF", "#FF0000", "#FF1493",
    "#FF9200", '#FF6347', '#ffd394', 
    "#FFF600", "#267574","#ADFF2F", "#00dd5d",
    "#2A52BE","#551A8B","pink", "#9acd32",
    "#808080", 	"#F5F5F5"
  )
  plottingLevelsNeutral  <- c('Aging', 'AID/APOBEC', 'BRCA1/2', 
                              'Smoking', 'Tobacco_chewing', 'Aflatoxin',
                              'UV', 'MMR/MSI', 'POLE', 'Signature_14', 
                              'TMZ/Alkylating', 'Signature_17', 'Signature_5', 'Signature_23', 'other', '_NotEnoughMuts')
  
  names(barColorPalette)= plottingLevelsNeutral
  signaturesLegend <-
  get_legend(ggplot(df, aes(x=Tumor_Sample_Barcode, y=1, fill=factor(dominantSigExomeAdj, levels=plottingLevelsNeutral)))+
                  geom_bar(stat='identity')+
                  scale_fill_manual(values=barColorPalette, drop=FALSE)+
                  theme(legend.title = element_text(size = 30))+
                  guides(fill=guide_legend(title="Signature")))
  
  fullPlot <- plot_grid(
    plot_two_bars_and_title_for_cancer_type(df[df$Cancer_Type == "Biliary Cancer",], "Biliary Cancer", barColorPalette),
    plot_two_bars_and_title_for_cancer_type(df[df$Cancer_Type == "Bladder Cancer",], "Bladder Cancer", barColorPalette),
    plot_two_bars_and_title_for_cancer_type(df[df$Cancer_Type == "Breast Carcinoma",], "Breast Carcinoma", barColorPalette),
    plot_two_bars_and_title_for_cancer_type(df[df$Cancer_Type == "Cancer of Unknown Primary",], "Cancer of Unknown Primary", barColorPalette),
    plot_two_bars_and_title_for_cancer_type(df[df$Cancer_Type == "Colorectal Cancer",], "Colorectal Cancer", barColorPalette),
    plot_two_bars_and_title_for_cancer_type(df[df$Cancer_Type == "Endometrial Cancer",], "Endometrial Cancer", barColorPalette),
    plot_two_bars_and_title_for_cancer_type(df[df$Cancer_Type == "Germ Cell Tumor",], "Germ Cell Tumor", barColorPalette),
    plot_two_bars_and_title_for_cancer_type(df[df$Cancer_Type == "Glioma",], "Glioma", barColorPalette),
    plot_two_bars_and_title_for_cancer_type(df[df$Cancer_Type == "Head and Neck Carcinoma",], "Head and Neck Carcinoma", barColorPalette),
    plot_two_bars_and_title_for_cancer_type(df[df$Cancer_Type == "Melanoma",], "Melanoma", barColorPalette),
    plot_two_bars_and_title_for_cancer_type(df[df$Cancer_Type == "Non-Small Cell Lung Cancer",], "Non-Small Cell Lung Cancer", barColorPalette),
    plot_two_bars_and_title_for_cancer_type(df[df$Cancer_Type == "Other",], "Other", barColorPalette),
    plot_two_bars_and_title_for_cancer_type(df[df$Cancer_Type == "Ovarian Cancer",], "Ovarian Cancer", barColorPalette),
    plot_two_bars_and_title_for_cancer_type(df[df$Cancer_Type == "Pancreatic Cancer",], "Pancreatic Cancer", barColorPalette),
    plot_two_bars_and_title_for_cancer_type(df[df$Cancer_Type == "Prostate Cancer",], "Prostate Cancer", barColorPalette),
    plot_two_bars_and_title_for_cancer_type(df[df$Cancer_Type == "Renal Cell Carcinoma",], "Renal Cell Carcinoma", barColorPalette),
    plot_two_bars_and_title_for_cancer_type(df[df$Cancer_Type == 'Small Cell Lung Cancer',], 'Small Cell Lung Cancer', barColorPalette),
    plot_two_bars_and_title_for_cancer_type(df[df$Cancer_Type == 'Soft Tissue Sarcoma',], 'Soft Tissue Sarcoma', barColorPalette), ncol=1, nrow=19)
  fullPlotWLeg <- plot_grid(fullPlot, signaturesLegend, ncol=2, rel_widths = c(1,.2))
  
  alignedPlotWithTitle <- plot_grid(ggplot()+ggtitle('Impact Signatures                    vs                    Exome Signatures')+theme(plot.title = element_text(size = 30, face = "bold")),
                                    fullPlotWLeg, ggplot()+labs(caption='plotSignatureDataSummary.R  inferrability_of_signatures_impact_vs_exome.ipynb'),
                                    nrow=3, rel_heights = c(.05,1,.025))
  return(alignedPlotWithTitle)
}

print(levels(df$Cancer_Type))

df <- read.table('/Users/friedman/Desktop/exomeProjectHelp/sigsDataForFigurePlot.tsv', sep = '\t', header=TRUE)

plt <- make_full_impact_vs_exome_plot(df)
ggsave('~/Desktop/plot.pdf', plot=plt,  width = 15, height = 20, units = c("in"), limitsize = FALSE)






barColorPalette = c(
  "#00DFFF", "#FF0000", "#FF1493",
  "#FF9200", '#FF6347', '#ffd394', 
  "#FFF600", "#267574","#ADFF2F", "#00dd5d",
  "#2A52BE","#551A8B","pink", "#9acd32",
  "#808080", 	"#F5F5F5"
)
plottingLevelsNeutral  <- c('Aging', 'AID/APOBEC', 'BRCA1/2', 
                            'Smoking', 'Tobacco_chewing', 'Aflatoxin',
                            'UV', 'MMR/MSI', 'POLE', 'Signature_14', 
                            'TMZ/Alkylating', 'Signature_17', 'Signature_5', 'Signature_23', 'other', '_NotEnoughMuts')

names(barColorPalette)= plottingLevelsNeutral

plot_two_bars_and_title_for_cancer_type(df[df$Cancer_Type == 'Endometrial Cancer',], 'Endometrial Cancer', barColorPalette)



##############
######################
##########################
################################
####################################
##########################################
emptyTheme <- theme(axis.line = element_blank(),
                    #axis.text.x = element_blank(),
                    axis.ticks = element_blank(),
                    axis.title.x = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank())

df <- read.table('/Users/friedman/Desktop/exomeProjectHelp/sigBinsInfo.tsv', sep = '\t', header=TRUE)

p <- ggplot(df, aes(x=nbin.96, y=reorder(Signature, nbin.96), colour=fracStronglyFavored))+
  geom_point()+
  scale_colour_gradient(high='black', low='gray', name = "Fraction signature \nprobability at strongly\nfavored motifs")+
  emptyTheme+
  ylab('Signature')+
  xlab('N bins > 1/96 signature chance')+
  ggtitle('Signature Specificity of Stratton 30 Signatures')

ggsave('~/Desktop/plot.pdf', plot=p,  width = 5, height = 5, units = c("in"), limitsize = FALSE)


