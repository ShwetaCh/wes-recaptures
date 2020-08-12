df = fread('/juno/work/ccs/bandlamc/other_stuff/impact_facets_cohort_annotation/impact_facets_annotation_08_05_2020.txt.cohort.txt')
head(df)

im_clin = fread('~/tempo-cohort-level/IM_metadata_040620.txt')
dim(im_clin)  #1636

df_im = left_join(im_clin, df, by = c(DMP = 'tumor_sample')) 
dim(df_im)
head(df_im)

df_im_f = df_im %>% select(DMP, purity_run_Purity, purity_run_Ploidy, has_qc, facets_qc) %>% distinct()
dim(df_im_f)
length(unique(df_im_f$DMP))
head(df_im_f)
table(df_im_f$has_qc)
table(df_im_f$facets_qc)


write.table(distinct(df_im_f),'~/tempo-cohort-level/FACETS-QC-Preview-IMPACT_081220.txt', sep = "\t", quote = F, row.names = F, append = F)
