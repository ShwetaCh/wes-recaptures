########################################################################################################################
### Load data
########################################################################################################################
  library(ggsci)
  library(ggplot2)
  library(ggpubr)
  library(cowplot)
  library(tidyr)
  library(data.table)
  library(plyr)
  library(dplyr)
  library(stringr)
  library(reshape2)
  library(purrr)
  library(forcats)
  library(grid)
  library(gridExtra)
'%nin%' = Negate('%in%')
#PREP RUN THIS FIRST /ifs/res/taylorlab/chavans/roslin_2.4_deliveries/all_comb_mafs/get_fillout_impact_called_exome_missed-parity-032619.R
########################################################################################################################

#CLINICAL
clin_data = fread('~/tempo-cohort-level/WES_metadata_040620.txt') %>% 
  select(CMO = Tumor_Sample_Barcode, 
         SAMPLE_ID = DMP, 
         CANCER_TYPE = Cancer_Type_Aggregate, 
         ONCOTREE_CODE = OncoTreeCode)
dim(clin_data)

#STATS
stats = fread('~/tempo-cohort-level/WES_IM_parity_stats.txt'); dim(stats)
#CMO_ID	Exome_Called	DMP_ID	IMPACT_Called	Exome_Recalled	Exome_Missed	Perc_Recalled	MissingVariants

#CLIN + STATS
big_wes = inner_join(stats, clin_data, by=c("DMP_ID"="SAMPLE_ID")) %>%
  filter(!is.na(Perc_Recalled)); dim(big_wes) #1473

#IMPACT MAF

big_wes_maf = fread('~/tempo-cohort-level/variants_mafanno.mut_purity.ccf.IMPACT.maf') %>% 
  select(-c(mutation_effect,oncogenic, LEVEL_1, LEVEL_2A, LEVEL_2B, LEVEL_3A, LEVEL_3B, LEVEL_4, LEVEL_R1, Highest_level))

DMPmapping = fread('~/tempo-cohort-level/CMO_DMP_mapping1636.txt')
oncokb_impact = fread('~/tempo-cohort-level/data_mutations_extended_somatic.oncokb.txt') %>% 
  filter(Tumor_Sample_Barcode %in% DMPmapping$DMP) %>%
  select(Tumor_Sample_Barcode,Highest_level)

length(intersect(oncokb_impact$Tumor_Sample_Barcode, DMPmapping$DMP)) #1473
length(setdiff(big_wes_maf$Tumor_Sample_Barcode1, oncokb_impact$Tumor_Sample_Barcode)) #0

big_wes_maf = big_wes_maf %>%  mutate(var_tag = str_c(Chromosome, ':', Start_Position, ':', Hugo_Symbol),
         Highest_level = mapvalues(Tumor_Sample_Barcode1, oncokb_impact$Tumor_Sample_Barcode, oncokb_impact$Highest_level)) %>%
  filter(Mutation_Status!='GERMLINE') 

########################################################################################################################
### Compare calls in impact/wes
### Then looking manually into mising variants
missingstuff = fread('~/tempo-cohort-level/RecallRate_WES_t_alt_WES_t_depth.txt')
write.table(
  filter(impact_mutations_big_wes,type %in% c('snv','indel'),wes_all_detected==T,wes_all_called==F) %>% 
  mutate(var_tag1 = str_c(Tumor_Sample_Barcode,':',Chromosome, ':', Start_Position, ':', Hugo_Symbol),
         WES_t_alt = mapvalues(var_tag1, missingstuff$var_tag1, missingstuff$WES_t_alt),
         WES_t_depth = mapvalues(var_tag1, missingstuff$var_tag1, missingstuff$WES_t_depth)) %>%  
  select(Tumor_Sample_Barcode,Hugo_Symbol,HGVSp_Short,t_alt_count,t_depth,n_alt_count,Variant_Classification,`is-a-hotspot`,oncokb_level, vaf_bin, purity, WES_t_alt, WES_t_depth,var_tag1), 
  '~/tempo-cohort-level/RecallRate_NotCalled_by_VAF_var_tag.txt', sep = "\t", row.names = FALSE, quote = FALSE)

write.table(
  filter(impact_mutations_big_wes,type %in% c('snv','indel'),wes_all_detected==F,wes_all_called==F) %>% 
  mutate(var_tag1 = str_c(Tumor_Sample_Barcode,':',Chromosome, ':', Start_Position, ':', Hugo_Symbol),
         WES_t_alt = mapvalues(var_tag1, missingstuff$var_tag1, missingstuff$WES_t_alt),
         WES_t_depth = mapvalues(var_tag1, missingstuff$var_tag1, missingstuff$WES_t_depth)) %>%  
  select(Tumor_Sample_Barcode,Hugo_Symbol,HGVSp_Short,t_alt_count,t_depth,n_alt_count,Variant_Classification,`is-a-hotspot`,oncokb_level, vaf_bin, purity, WES_t_alt, WES_t_depth,var_tag1), 
  '~/tempo-cohort-level/RecallRate_NotDetected_by_VAF.txt', sep = "\t", row.names = FALSE, quote = FALSE)

#!is.na(oncokb_level)

write.table(
  filter(impact_mutations_big_wes,type %in% c('snv','indel'),wes_all_detected==T,wes_all_called==F,!is.na(oncokb_level)) %>% 
    mutate(var_tag1 = str_c(Tumor_Sample_Barcode,':',Chromosome, ':', Start_Position, ':', Hugo_Symbol),
           WES_t_alt = mapvalues(var_tag1, missingstuff$var_tag1, missingstuff$WES_t_alt),
           WES_t_depth = mapvalues(var_tag1, missingstuff$var_tag1, missingstuff$WES_t_depth)) %>%  
    select(Tumor_Sample_Barcode,Hugo_Symbol,HGVSp_Short,t_alt_count,t_depth,n_alt_count,Variant_Classification,`is-a-hotspot`,oncokb_level, vaf_bin, purity, WES_t_alt, WES_t_depth,var_tag1), 
  '~/tempo-cohort-level/RecallRate_NotCalled_by_OncoKB_var_tag.txt', sep = "\t", row.names = FALSE, quote = FALSE)

write.table(
  filter(impact_mutations_big_wes,type %in% c('snv','indel'),wes_all_detected==F,wes_all_called==F,!is.na(oncokb_level)) %>% 
    mutate(var_tag1 = str_c(Tumor_Sample_Barcode,':',Chromosome, ':', Start_Position, ':', Hugo_Symbol),
           WES_t_alt = mapvalues(var_tag1, missingstuff$var_tag1, missingstuff$WES_t_alt),
           WES_t_depth = mapvalues(var_tag1, missingstuff$var_tag1, missingstuff$WES_t_depth)) %>%  
    select(Tumor_Sample_Barcode,Hugo_Symbol,HGVSp_Short,t_alt_count,t_depth,n_alt_count,Variant_Classification,`is-a-hotspot`,oncokb_level, vaf_bin, purity, WES_t_alt, WES_t_depth,var_tag1),
  '~/tempo-cohort-level/RecallRate_NotDetected_by_OncoKB.txt', sep = "\t", row.names = FALSE, quote = FALSE)
#filter(impact_mutations_big_wes,type=='snv',wes_all_detected==T,wes_all_called==F,vaf_bin=="5%") %>% select(Tumor_Sample_Barcode,Hugo_Symbol,HGVSp_Short,t_alt_count,t_depth,n_alt_count,Variant_Classification,`is-a-hotspot`,oncokb_level)
#filter(impact_mutations_big_wes,type=='snv',wes_all_detected==T,wes_all_called==F,oncokb_level=="4") %>% select(Tumor_Sample_Barcode,cmo_id,Hugo_Symbol,HGVSp_Short,t_alt_count,t_depth,n_alt_count,Variant_Classification,`is-a-hotspot`,oncokb_level) 
#samplelevel_maf = fread('/ifs/res/taylorlab/chavans/roslin_2.4_deliveries/all_samplelevel_mafs_Proj_07951_RSTUWZ/comb.RSTUWZ.sample.level.v2.maf')

########################################################################################################################
# Put together fillouts
  # unfil_variants_fillouts -- signed out
  # unfil_variants_fillout_all_impact -- all impact
fillouts = dir('/juno/work/ccs/chavans/wes-paper/all_fillouts_imcalls', pattern = '_fo.maf$', full.names = T)
fillouts = ldply(fillouts, fread) %>%
    mutate(var_tag = str_c(Chromosome, ':', Start_Position, ':', Hugo_Symbol))

fillouts_big_wes =dir('/juno/work/ccs/chavans/wes-paper/all_fillouts_imcalls', pattern = 'per_sample_output_fo.maf$', full.names = T)
head(fillouts_big_wes)

missingstuff = NA
impact_mutations_big_wes = map_dfr(big_wes$CMO_ID,function(id) {
    sample_info = filter(big_wes, CMO_ID==id)
    dmp_id = sample_info$DMP_ID
    sample_impact = filter(big_wes_maf, Tumor_Sample_Barcode1 == dmp_id)
    print(id)
    
    if (sample_info$Perc_Recalled==100) {
        sample_impact$wes_all_called = TRUE
        sample_impact$wes_all_detected = TRUE
        sample_impact$wes_not_detectable = FALSE
        sample_impact$cmo_id = id
        sample_impact
        #missingstuff = NA
    } else {
        fo = fread(grep(id, fillouts_big_wes, value = TRUE)) %>%
        filter(Tumor_Sample_Barcode == id) %>%
        mutate(var_tag = str_c(Chromosome, ':', Start_Position, ':', Hugo_Symbol),
               t_total_count = t_alt_count+t_ref_count)
        fo_detected = filter(fo, t_alt_count>0)
        fo_not_detectable = filter(fo, t_total_count<20 & t_alt_count==0)
        mutate(sample_impact,
               cmo_id = id,
               wes_all_called = var_tag %nin% fo$var_tag,
               wes_all_detected = var_tag %in% fo_detected$var_tag | wes_all_called == T,
               wes_not_detectable = var_tag %in% fo_not_detectable$var_tag)
        #missingstuff = paste(fo$Tumor_Sample_Barcode, fo$Hugo_Symbol,	fo$HGVSp_Short, fo$var_tag, fo$t_alt_count, fo$t_total_count, sep="|")
        #write.table(missingstuff, '~/tempo-cohort-level/RecallRate_WES_t_alt_WES_t_depth.txt', sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)      
        }
  }) %>% 
  mutate(type = plyr::mapvalues(Variant_Type, c('INS','DEL','SNP'), c('indel', 'indel', 'snv')),
              vaf = as.numeric(t_alt_count)/(as.numeric(t_ref_count)+as.numeric(t_alt_count))
         #missingstuff = missingstuff
         ) %>%
  filter(type %in% c('indel','snv'))

dim(impact_mutations_big_wes)#13496   188
# impact_mutations_big_wes$vaf_bin = cut(impact_mutations_big_wes$vaf, breaks = c(seq(0,.75,.05),1),
#                                labels = str_c(c(seq(5,75,5),100), '%'),

impact_mutations_big_wes$vaf_bin = cut(impact_mutations_big_wes$vaf, breaks = c(seq(0,1,.05)),
                               #labels = c(str_c(c(seq(5,70,5)), '%'),'>75%'),
                               labels = c(str_c(c(seq(5,100,5)), '%')),
                               ordered_result = T)
#impact_mutations_big_wes = replace_na(impact_mutations_big_wes, list(vaf_bin = '>75%'))

### add oncokb annotation
# impact_mutations = left_join(impact_mutations, select(impact_maf, Tumor_Sample_Barcode, var_tag, Highest_level)) %>%
impact_mutations_big_wes = mutate(impact_mutations_big_wes, oncokb_level = str_extract(Highest_level, '[0-9]{1}[A-Z]{0,1}'))
unique(impact_mutations_big_wes$Highest_level)
(unique(impact_mutations_big_wes$oncokb_level))

#impact_mutations_big_wes_saved = impact_mutations_big_wes

#head(impact_mutations_big_wes)
###################

#filter(.,cmo_id %nin% sample_mapping_redacted$CMO_Sample_ID_fixed) %>%
p1 = filter(impact_mutations_big_wes, type=='snv') %>%
      group_by(vaf_bin) %>%
    dplyr::summarise(total = n(),
                `Called mutations` = sum(wes_all_called==T)/n(),
              #`Detected mutations` = sum(wes_all_detected==T)/n(),
              `Detected mutations at coverage 20` = sum(wes_all_detected==T)/sum(wes_not_detectable==F)) %>%
    select(-total) %>%
    melt(id.vars = 'vaf_bin') %>%
    #summarySE(., measurevar="value", groupvars=c("vaf_bin")) %>%
    ggplot(., aes(x = as.factor(vaf_bin), y = value, group = variable, color = variable)) +
    geom_vline(xintercept = c(1,2), color = 'darkred', linetype = 'dashed') +
    annotate('text', x = c(1,2)+.4, y = .05, color = 'darkred', label = c('5%', '10%'), size = 4) +
    geom_hline(yintercept = 0.95, linetype = 'dashed', color = 'black') +
    # annotate('text', x = 20, y = 0.9520198-.02, color = 'black', label = 'Median') +
    geom_line(size = 1.5) +
    #geom_errorbar(aes(ymin = value - ci, ymax = value + ci)) +
    scale_color_jama() +
    scale_y_continuous(expand = c(0,0), limits = c(0,1.04)) +
    scale_x_discrete(breaks = c('10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')) +
    theme_classic() +
    theme(#aspect.ratio = 1,
      legend.position = c(1,0), legend.justification = c(1,0),legend.background=element_blank()) +
    labs(x = 'Variant allele fraction\n(bins of 10%)', y = 'Recall') +
    guides(color = guide_legend(title = '', keywidth = unit(2,'lines'))) 

p2 = filter(impact_mutations_big_wes, type=='snv') %>%
     dplyr::count(vaf_bin) %>%
    ggplot(., aes(x = vaf_bin, y = n)) +
    geom_col(fill = '#08519c') +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = '', y = 'Count') +
    theme_classic() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    labs(title = 'SNVs only')


p3 = filter(impact_mutations_big_wes,  type %in% c('snv','indel')) %>%
     group_by(., vaf_bin) %>%
    dplyr::summarise(total = n(),
              `Called mutations` = sum(wes_all_called==T)/n(),
              #`Detected mutations` = sum(wes_all_detected==T)/n(),
              `Detected mutations at coverage 20` = sum(wes_all_detected==T)/sum(wes_not_detectable==F)) %>%
  select(-total) %>%
    melt(id.vars = 'vaf_bin') %>%
    ggplot(., aes(x = as.factor(vaf_bin), y = value, group = variable, color = variable)) +
    geom_vline(xintercept = c(1,2), color = 'darkred', linetype = 'dashed') +
    annotate('text', x = c(1,2)+.4, y = .05, color = 'darkred', label = c('5%', '10%'), size = 4) +
    geom_hline(yintercept = 0.948045, linetype = 'dashed', color = 'black') +
    # annotate('label', x = 20, y = 0.948045-.02, color = 'black', label = 'Median') +
    geom_line(size = 1.5) +
    scale_color_jama() +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.04)) +
  scale_x_discrete(breaks = c('10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')) +
    theme_classic() +
    theme(legend.position = c(1,0), legend.justification = c(1,0), legend.background=element_blank()) +
    labs(x = 'Variant allele fraction\n(bins of 10%)', y = 'Recall') +
    guides(color = guide_legend(title = '', keywidth = unit(2,'lines')))

p4 = filter(impact_mutations_big_wes, type %in% c('snv','indel')) %>%
    dplyr::count(., vaf_bin, type) %>% 
    mutate(type = fct_recode(type, Indel = 'indel', SNV = 'snv')) %>%
    ggplot(., aes(x = vaf_bin, y = n, fill = type)) +
    geom_col() +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values = c('SNV'='#08519c', 'Indel'='#6baed6'), '') +
    labs(x = '', y = 'Count') +
    theme_classic() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = c(1,1), legend.justification = c(1,1), legend.background=element_blank()) +
    labs(title = 'SNVs + Indels')
#pdf('~/Desktop/variant-qc-vaf-bins_07951.pdf', w = 10, h = 6)
#par(mfrow=c(2,2))

grid.arrange(p2,p4,p1,p3,ncol = 2, heights = c(2,2), newpage = F)
ggsave('~/tempo-cohort-level/Figure1D_RecallRates_VAF.pdf', plot = grid.arrange(p2,p4,p1,p3,ncol = 2, heights = c(2,2), newpage = F))

#dev.off()


#all_samples = bind_rows(master, rename(big_wes, Cancer_Type = CANCER_TYPE))
all_samples = dplyr::rename(big_wes, Cancer_Type = CANCER_TYPE) 

master_ct_conv = dplyr::count(all_samples, Cancer_Type) %>%
    mutate(ct = ifelse(n<25, 'Other', Cancer_Type))

master_ct_count = group_by(master_ct_conv, ct) %>%
    dplyr::summarise(n = sum(n)) %>%
    ungroup() %>%
    mutate(ct = fct_reorder(ct, -n),
           ct = fct_relevel(ct, 'Other', after = Inf))
arrange(master_ct_conv,desc(n))

p6 = ggplot(master_ct_count, aes(x = ct, y = n)) +
  geom_col() +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = '', y = 'Sample count') +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
p5=p6
p7 = left_join(impact_mutations_big_wes, select(all_samples, cmo_id = CMO_ID, Cancer_Type)) %>%
     left_join(., select(master_ct_conv, Cancer_Type, ct, n)) %>%
    mutate(ct = fct_reorder(ct, -n),
         ct = fct_relevel(ct, 'Other', after = Inf)) %>%
    filter(type=='snv') %>%
    group_by(ct) %>%
    dplyr::summarise(total = n(),
              `Called mutations` = sum(wes_all_called==T),
              `False negative` = sum(wes_all_called==F)) %>%
    select(-total) %>%
    melt(id.vars = 'ct') %>%
    mutate(variable = fct_relevel(variable, 'False negative')) %>%
    ggplot(., aes(x = ct, y = value, fill = variable)) +
    geom_col(position = 'fill') +
    geom_hline(yintercept = 0.95, color = 'black', linetype = 'dashed') +
    scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = rev(c('#4393c3','#d6604d'))) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
          plot.margin = unit(c(0.75,.25,.25,.25), 'lines'), legend.position = 'top') +
  labs(x = '', y = 'Recall') +  labs(title = 'Only SNVs') +
    guides(fill = guide_legend(title = '', reverse = T))

p8 = left_join(impact_mutations_big_wes, select(all_samples, cmo_id = CMO_ID, Cancer_Type)) %>%
   left_join(., select(master_ct_conv, Cancer_Type, ct, n)) %>%
  mutate(ct = fct_reorder(ct, -n),
         ct = fct_relevel(ct, 'Other', after = Inf)) %>%
  filter(type %in% c('snv','indel')) %>%
  group_by(ct) %>%
  dplyr::summarise(total = n(),
            `Called mutations` = sum(wes_all_called==T),
            `False negative` = sum(wes_all_called==F)) %>%
  select(-total) %>%
  melt(id.vars = 'ct') %>%
  mutate(variable = fct_relevel(variable, 'False negative')) %>%
  ggplot(., aes(x = ct, y = value, fill = variable)) +
  geom_col(position = 'fill') +
  geom_hline(yintercept = 0.95, color = 'black', linetype = 'dashed') +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = rev(c('#4393c3','#d6604d'))) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
        plot.margin = unit(c(0.75,.25,.25,.25), 'lines'), legend.position = 'top') +
  labs(x = '', y = 'Recall') +  labs(title = 'SNVs and indels') +
  guides(fill = guide_legend(title = '', reverse = T))

grid.arrange(p5,p6,p7,p8,ncol = 2, heights = c(2,3), newpage = F)
ggsave('~/tempo-cohort-level/Figure1D_RecallRates_CancerTypes.pdf', plot = grid.arrange(p5,p6,p7,p8,ncol = 2, heights = c(2,3), newpage = F))


 ### By actionability
p7 = filter(impact_mutations_big_wes, oncokb_level!='', type=='snv') %>%
  dplyr::count(oncokb_level) %>%
  mutate(oncokb_level = fct_relevel(oncokb_level, rev(c('1','2','3A','3B','4')))) %>%
  ggplot(., aes(x = oncokb_level, y = n)) +
  geom_col() +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = '', y = 'Count', title = 'SNVs Only') +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

p8 = filter(impact_mutations_big_wes, oncokb_level!='', type=='snv') %>%
  group_by(oncokb_level) %>%
  dplyr::summarise(total = n(),
            `Called mutations` = sum(wes_all_called==T)/n(),
            #`Detected mutations` = sum(wes_all_detected==T)/n(),
            `Detected mutations cov.20` = sum(wes_all_detected==T)/sum(wes_not_detectable==F)) %>%
  select(-total) %>%
  ungroup() %>%
  mutate(oncokb_level = fct_relevel(oncokb_level, rev(c('1','2','3A','3B','4')))) %>%
  melt(id.vars = c('oncokb_level')) %>%
  ggplot(., aes(x = oncokb_level, y = value, fill = variable)) +
  geom_col(position = 'dodge') +
  scale_y_continuous(expand = c(0,0)) +
  geom_hline(yintercept = 1, color = 'black', linetype = 'dashed', size =1.5) +
  scale_fill_jama() +
  theme_classic() +
  theme(plot.margin = unit(c(0.75,.25,.25,.25), 'lines'), legend.position = 'top') +
  labs(x = '', y = 'Recall') +
  guides(fill = guide_legend(title = ''))

p9 = filter(impact_mutations_big_wes, oncokb_level!='', type %in% c('snv','indel')) %>%
  dplyr::count(oncokb_level) %>%
  mutate(oncokb_level = fct_relevel(oncokb_level, rev(c('1','2','3A','3B','4')))) %>%
  ggplot(., aes(x = oncokb_level, y = n)) +
  geom_col() +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = '', y = 'Count', title = 'SNVs and Indels') +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

p10 = filter(impact_mutations_big_wes, oncokb_level!='', type %in% c('snv','indel')) %>%
  group_by(oncokb_level) %>%
  dplyr::summarise(total = n(),
            `Called mutations` = sum(wes_all_called==T)/n(),
            #`Detected mutations` = sum(wes_all_detected==T)/n(),
            `Detected mutations cov.20` = sum(wes_all_detected==T)/sum(wes_not_detectable==F)) %>%
  select(-total) %>%
  ungroup() %>%
  mutate(oncokb_level = fct_relevel(oncokb_level, rev(c('1','2','3A','3B','4')))) %>%
  melt(id.vars = c('oncokb_level')) %>%
  ggplot(., aes(x = oncokb_level, y = value, fill = variable)) +
  geom_col(position = 'dodge') +
  scale_y_continuous(expand = c(0,0)) +
  geom_hline(yintercept = 1, color = 'black', linetype = 'dashed', size = 1.5) +
  scale_fill_jama() +
  theme_classic() +
  theme(plot.margin = unit(c(0.75,.25,.25,.25), 'lines'), legend.position = 'top') +
  labs(x = '', y = 'Recall') +
  guides(fill = guide_legend(title = ''))
#pdf('~/Desktop/variant-qc-oncokb.pdf', w = 8, h = 7)
grid.arrange(p7, p9, p8, p10, nrow = 2, newpage = F)
ggsave('~/tempo-cohort-level/Figure1D_RecallRates_OncoKBLevels.pdf', plot = grid.arrange(p7, p9, p8, p10, ncol = 2, heights = c(2,2), newpage = F))

#dev.off()

#oncokb_not_detected=filter(impact_mutations_big_wes,wes_all_called==FALSE,wes_all_detected==FALSE, oncokb_level=="4") %>% select(cmo_id,vaf_bin,Tumor_Sample_Barcode, HGVSp_Short, var_tag, Hugo_Symbol, t_alt_count, t_depth)

#filter(impact_mutations_big_wes,vaf_bin=="10%",wes_all_called==FALSE,wes_all_detected==TRUE) %>% select(cmo_id,Tumor_Sample_Barcode, HGVSp_Short, var_tag, Hugo_Symbol, t_alt_count, t_depth)

### by purity
impact_mutations_big_wes1 = impact_mutations_big_wes_saved
facets_purity = fread('~/tempo-cohort-level/WES_metadata_040620.txt') %>% select(Tumor_Sample_Barcode,Purity_Reviewed) %>%
  filter(Tumor_Sample_Barcode %in% DMPmapping$CMO)
dim(facets_purity)
dim(impact_mutations_big_wes1)
impact_mutations_big_wes1 = impact_mutations_big_wes1 %>% 
  mutate(.,Purity = plyr::mapvalues(cmo_id, 
                                    facets_purity$Tumor_Sample_Barcode, 
                                    as.numeric(facets_purity$Purity_Reviewed)))

impact_mutations_big_wes1$fp_bin = cut(as.numeric(impact_mutations_big_wes1$Purity), breaks = c(seq(0,1,.05)),
                                       #labels = c(str_c(c(seq(5,70,5)), '%'),'>75%'),
                                       labels = c(str_c(c(seq(5,100,5)),'%')),
                                       ordered_result = T)
dim(impact_mutations_big_wes1)
range(facets_purity$Purity_Reviewed)
#[1] 0.11 0.97
unique(impact_mutations_big_wes1$fp_bin)
#[1] 80%  40%  65%  55%  75%  95%  45%  60%  35%  25%  15%  30%  85% 
#[14] 70%  50%  100% 90%  20% 
#20 Levels: 5% < 10% < 15% < 20% < 25% < 30% < 35% < 40% < ... < 100%

#impact_mutations_big_wes = impact_mutations_big_wes %>% mutate(.,fp_bin_perc = plyr::mapvalues(cmo_id, facets_purity$Sample_ID,as.character(facets_purity$fp_bin_perc)))

#impact_mutations_big_wes_p = impact_mutations_big_wes_p %>% mutate(.,Purity = ifelse(is.na(Purity),1.1,Purity))#, !is.na(Purity)

df =  filter(impact_mutations_big_wes1, type %in% c('snv','indel')) %>%
  group_by(vaf_bin) %>%
  #arrange(fp_bin) %>%
  dplyr::summarise(total = n(),
                   `Called mutations` = sum(wes_all_called==T)/n(),
                   #`Detected mutations` = sum(wes_all_detected==T)/n(),
                   `Detected mutations at coverage 20` = sum(wes_all_detected==T)/sum(wes_not_detectable==F))
dim(df)
df

df1 =  filter(impact_mutations_big_wes1, type %in% c('snv','indel')) %>%
  group_by(fp_bin) %>%
  #arrange(fp_bin) %>%
  dplyr::summarise(total = n(),
            `Called mutations` = sum(wes_all_called==T)/n(),
            #`Detected mutations` = sum(wes_all_detected==T)/n(),
            `Detected mutations at coverage 20` = sum(wes_all_detected==T)/sum(wes_not_detectable==F))
dim(df1)
df1

df2 =  filter(impact_mutations_big_wes, type %in% c('snv','indel')) %>%
  group_by(oncokb_level) %>%
  #arrange(fp_bin) %>%
  dplyr::summarise(total = n(),
                   `Called mutations` = sum(wes_all_called==T)/n(),
                   #`Detected mutations` = sum(wes_all_detected==T)/n(),
                   `Detected mutations at coverage 20` = sum(wes_all_detected==T)/sum(wes_not_detectable==F))
dim(df2)
df2

df3 =  left_join(impact_mutations_big_wes, select(all_samples, cmo_id = CMO_ID, Cancer_Type)) %>%
  left_join(., select(master_ct_conv, Cancer_Type, ct, n)) %>%
  mutate(ct = fct_reorder(ct, -n),
         ct = fct_relevel(ct, 'Other', after = Inf)) %>%
  filter(type %in% c('snv','indel')) %>%
  group_by(ct) %>%
  dplyr::summarise(total = n(),
                   `Called mutations` = sum(wes_all_called==T),
                   `False negative` = sum(wes_all_called==F))
dim(df3)
df3

#impact_mutations_big_wes1 = filter(impact_mutations_big_wes1, grepl("%",fp_bin)) %>% filter(Tumor_Sample_Barcode %in% big_wes$DMP_ID)
dim(impact_mutations_big_wes1)

fpp = filter(impact_mutations_big_wes1, type %in% c('snv')) %>%
  group_by(fp_bin) %>%
  dplyr::summarise(total = n(),
            `Called mutations` = sum(wes_all_called==T)/n(),
            #`Detected mutations` = sum(wes_all_detected==T)/n(),
            `Detected mutations at coverage 20` = sum(wes_all_detected==T)/sum(wes_not_detectable==F)) %>%
  select(-total) %>%
  melt(id.vars = 'fp_bin') %>%
  ggplot(., aes(x = as.factor(fp_bin), y = value, group = variable, color = variable)) +
  #geom_vline(xintercept = c(1,2), color = 'darkred', linetype = 'dashed') +
  #annotate('text', x = c(1,2)+.4, y = .05, color = 'darkred', label = c('5%', '10%'), size = 2) +
  #geom_hline(yintercept = 0.9520198, linetype = 'dashed', color = 'black') +
  #annotate('text', x = 20, y = 0.9520198-.02, color = 'black', label = 'Median') +
  geom_line(size = 1.5) +
  scale_color_jama() +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.04)) +
  scale_x_discrete(breaks = c('10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')) +
        #labels = c('10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')) +
  theme_classic() +
  theme(#aspect.ratio = 1,
    legend.position = c(1,0), legend.justification = c(1,0), legend.background = element_blank()) +
  labs(x = 'Purity estimate\n(bins of 10%)', y = 'Recall') +
  guides(color = guide_legend(title = '', keywidth = unit(2,'lines'))) +
  labs(title = 'SNVs only')
  
fpp
fpp2 = filter(impact_mutations_big_wes1, type %in% c('snv','indel')) %>%
  group_by(fp_bin) %>%
  dplyr::summarise(total = n(),
                   `Called mutations` = sum(wes_all_called==T)/n(),
                   #`Detected mutations` = sum(wes_all_detected==T)/n(),
                   `Detected mutations at coverage 20` = sum(wes_all_detected==T)/sum(wes_not_detectable==F)) %>%
  select(-total) %>%
  melt(id.vars = 'fp_bin') %>%
  ggplot(., aes(x = as.factor(fp_bin), y = value, group = variable, color = variable)) +
  #geom_vline(xintercept = c(1,2), color = 'darkred', linetype = 'dashed') +
  #annotate('text', x = c(1,2)+.4, y = .05, color = 'darkred', label = c('5%', '10%'), size = 2) +
  #geom_hline(yintercept = 0.9520198, linetype = 'dashed', color = 'black') +
  #annotate('text', x = 20, y = 0.9520198-.02, color = 'black', label = 'Median') +
  geom_line(size = 1.5) +
  scale_color_jama() +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.04)) +
  scale_x_discrete(breaks = c('10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')) +
  #scale_x_discrete(labels = c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,'NA')) +
  theme_classic() +
  theme(#aspect.ratio = 1,
    legend.position = c(1,0), legend.justification = c(1,0), legend.background = element_blank()) +
  labs(x = 'Purity estimate\n(bins of 10%)', y = 'Recall') +
  guides(color = guide_legend(title = '', keywidth = unit(2,'lines')))+
  labs(title = 'SNVs and Indels')
fpp
fpp2

fpc = filter(impact_mutations_big_wes1, type %in% c('snv','indel')) %>%
  dplyr::count(., fp_bin, type) %>% 
  mutate(type = fct_recode(type, Indel = 'indel', SNV = 'snv')) %>%
  ggplot(., aes(x = fp_bin, y = n, fill = type)) +
  geom_col() +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(breaks = c('10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')) +
  scale_fill_manual(values = c('SNV'='#08519c', 'Indel'='#6baed6'), '') +
  labs(x = '', y = 'Mutations') +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = c(0.96,1.1), legend.justification = c(1,1)) 
#+labs(title = 'SNVs and Indels')

sc = filter(impact_mutations_big_wes1, type %in% c('snv')) %>% 
  select(Tumor_Sample_Barcode, fp_bin) %>% distinct(.) %>% 
  group_by(fp_bin) %>% 
  dplyr::summarise(total_samples = n()) %>%
  ggplot(., aes(x = fp_bin, y = total_samples)) +
  geom_col(fill = '#08519c') +
  labs(x = '', y = 'Samples') +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_x_discrete(breaks = c('10%','20%','30%','40%','50%','60%','70%','80%','90%','100%'))
#fpc
grid.arrange(sc, fpc, fpp, fpp2, nrow = 2, newpage = F)
ggsave('~/tempo-cohort-level/Figure1D_RecallRates_Purity.pdf', plot = grid.arrange(sc, fpc, fpp, fpp2, nrow = 2, newpage = F))

###################################
# Final figure with errors bars, new version as of 08/01/2020
###################################
df = fread('/Users/chavans/tempo-cohort-level/Recall_graph_w_errorbars_080120.txt')
head(df)

dff1 = df %>% 
  group_by(Bin) %>% 
  filter(Category == "VAF") %>% 
  mutate(se_Perc_called_mutations = Perc_called_mutations/sqrt(N), 
                                        se_low_c = Perc_called_mutations - se_Perc_called_mutations,  
                                        se_high_c = Perc_called_mutations + se_Perc_called_mutations,
                                        se_Perc_detected_mutations_20cov = Perc_detected_mutations_20cov/sqrt(N), 
                                        se_low_d = Perc_detected_mutations_20cov - se_Perc_detected_mutations_20cov,  
                                        se_high_d = Perc_detected_mutations_20cov + se_Perc_detected_mutations_20cov) %>%
  select(-c(Category, N,starts_with('se_l'),starts_with('se_h')))

dff1.m = melt(dff1, id.vars = c('Bin','se_Perc_detected_mutations_20cov','se_Perc_called_mutations'))

dff2.m = mutate(dff1.m, se = ifelse(variable == "Perc_called_mutations",se_Perc_called_mutations,
                                    ifelse(variable == "Perc_detected_mutations_20cov",se_Perc_detected_mutations_20cov,"")))
dff2.m = mutate(dff2.m, se = as.numeric(se), value = as.numeric(value))

ggplot(dff2.m, aes(x=Bin*100, y=value, group=variable, color=variable)) + 
  geom_line() +
  geom_pointrange(aes(ymin=value-(se), ymax=value+(se))) +
  # annotate('text', x = 20, y = 0.9520198-.02, color = 'black', label = 'Median') +
  geom_line(size = 1) +
  #geom_errorbar(,aes(ymin = value - ci, ymax = value + ci)) +
  scale_color_jama() +
  scale_y_continuous(expand = c(0,0), n.breaks=10, limits = c(0,1.15)) +
  scale_x_discrete(limits = c(10,20,30,40,50,60,70,80,90,100)) +
  theme_classic() +
  theme(#aspect.ratio = 1,
    legend.position = c(1,0), legend.justification = c(1,0),legend.background=element_blank()) +
  labs(x = 'Variant allele fraction\n(bins of 10%)', y = 'Recall') +
  guides(color = guide_legend(title = '', keywidth = unit(2,'lines'))) +
  geom_vline(xintercept = c(5,10), color = 'darkred', linetype = 'dashed') +
  annotate('text', x = c(7,12), y = .05, color = 'darkred', label = c('5%', '10%'), size = 4) +
  geom_hline(yintercept = 0.95, linetype = 'dashed', color = 'black')
  
###########################

dff1 = df %>% 
  group_by(Bin) %>% 
  filter(Category == "PURITY") %>% 
  mutate(se_Perc_called_mutations = Perc_called_mutations/sqrt(N), 
         se_low_c = Perc_called_mutations - se_Perc_called_mutations,  
         se_high_c = Perc_called_mutations + se_Perc_called_mutations,
         se_Perc_detected_mutations_20cov = Perc_detected_mutations_20cov/sqrt(N), 
         se_low_d = Perc_detected_mutations_20cov - se_Perc_detected_mutations_20cov,  
         se_high_d = Perc_detected_mutations_20cov + se_Perc_detected_mutations_20cov) %>%
  select(-c(Category, N,starts_with('se_l'),starts_with('se_h')))

dff1.m = melt(dff1, id.vars = c('Bin','se_Perc_detected_mutations_20cov','se_Perc_called_mutations'))

dff2.m = mutate(dff1.m, se = ifelse(variable == "Perc_called_mutations",se_Perc_called_mutations,
                                    ifelse(variable == "Perc_detected_mutations_20cov",se_Perc_detected_mutations_20cov,"")))
dff2.m = mutate(dff2.m, se = as.numeric(se), value = as.numeric(value))

ggplot(dff2.m, aes(x=Bin*100, y=value, group=variable, color=variable)) + 
  geom_line() +
  geom_pointrange(aes(ymin=value-(se), ymax=value+(se))) +
  # annotate('text', x = 20, y = 0.9520198-.02, color = 'black', label = 'Median') +
  geom_line(size = 1) +
  #geom_errorbar(,aes(ymin = value - ci, ymax = value + ci)) +
  scale_color_jama() +
  scale_y_continuous(expand = c(0,0), n.breaks=10, limits = c(0,1.15)) +
  scale_x_discrete(limits = c(15,20,30,40,50,60,70,80,90,100)) +
  theme_classic() +
  theme(#aspect.ratio = 1,
    legend.position = c(1,0), legend.justification = c(1,0),legend.background=element_blank()) +
  labs(x = 'Purity estimates\n(bins of 10%)', y = 'Recall') +
  guides(color = guide_legend(title = '', keywidth = unit(2,'lines'))) +
  #geom_vline(xintercept = c(5,10), color = 'darkred', linetype = 'dashed') +
  #annotate('text', x = c(7,12), y = .05, color = 'darkred', label = c('5%', '10%'), size = 4) +
  geom_hline(yintercept = 0.95, linetype = 'dashed', color = 'black')
