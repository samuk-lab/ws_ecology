# calculate fst (and outliers?) using OutFLANK
# KMS Nov 2018

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("qvalue", version = "3.8")
# devtools::install_github("whitlock/OutFLANK")

library("tidyverse")
library("OutFLANK")

#################################
# read in raw data
#################################

# the pruned snp table
snp_pruned <- read.table("data/snp_tables/whtstbk_2020_LD_0.2_MAF_0.01.txt", h = T)

# the UNpruned snp table
snp_tab <- read.table("data/snp_tables/whtstbk_2020_MAF_0.01.txt", h = T)

# remove individuals from poorly sampled pops
snp_tab <- snp_tab [, !grepl("LD|SP|QR", names(snp_tab))]

# inject 9s for NAs
snp_tab[is.na(snp_tab)] <- 9

# the meta data
meta_df <- read.csv("metadata/mega_meta.csv") %>%
  mutate(Indiv = gsub("whtstbk_gbs_2012_brds_", "", id)) %>%
  mutate(Indiv = gsub("_", ".", Indiv)) %>%
  mutate(Indiv = case_when(year == "2012" ~ paste0(Indiv, "_2013"),
                           year == "2014" ~ paste0(Indiv, "_2015"),
                           TRUE ~ Indiv)) %>%
  mutate(Indiv = ifelse(Indiv == "CL64.2_2015", "CL64-2_2015", Indiv)) %>%
  select(Indiv, pop, year, cluster, sex)

# population labels
pop_labels <- data.frame(Indiv = names(snp_tab)) %>%
  left_join(meta_df) %>%
  pull(cluster)


####################################
# fst: white vs. common
####################################

# target inds
inds_to_include <- meta_df %>%
  filter(cluster %in% c("cmn", "wht")) %>%
  filter(!is.na(cluster)) %>%
  pull(Indiv)

# subset the master snp table
snp_tab_sub <- snp_tab[,(names(snp_tab) %in% inds_to_include)]

# grab the population labels from the meta data table
pop_labels <- data.frame(Indiv = names(snp_tab_sub)) %>%
  left_join(meta_df) %>%
  pull(cluster) %>%
  droplevels

my_fst <- MakeDiploidFSTMat(t(snp_tab_sub), locusNames = row.names(snp_tab_sub), popNames = pop_labels)

# check He vs fst
plot(my_fst$He, my_fst$FST)

#### Evaluating OutFLANK with trimmed SNPs ####
my_fst_pruned <- my_fst %>% filter(LocusName %in% row.names(snp_pruned))
out_trim <- OutFLANK(my_fst_pruned, NumberOfSamples = 2, 
                     RightTrimFraction = 0.1, LeftTrimFraction = 0.05,
                     qthreshold = 0.01, Hmin = 0.05)

# calculate FST for the unpruned data
# and find outliers

fst_wht_cmn <- pOutlierFinderChiSqNoCorr(my_fst, Fstbar = out_trim$FSTNoCorrbar, 
                                dfInferred = out_trim$dfInferred, qthreshold = 0.01, Hmin=0.05)

####################################
# fst: white vs. cbr
####################################

# target inds
inds_to_include <- meta_df %>%
  filter(cluster %in% c("cbr", "wht")) %>%
  filter(!is.na(cluster)) %>%
  pull(Indiv)

# subset the master snp table
snp_tab_sub <- snp_tab[,(names(snp_tab) %in% inds_to_include)]

# grab the population labels from the meta data table
pop_labels <- data.frame(Indiv = names(snp_tab_sub)) %>%
  left_join(meta_df) %>%
  pull(cluster) %>%
  droplevels

my_fst <- MakeDiploidFSTMat(t(snp_tab_sub), locusNames = row.names(snp_tab_sub), popNames = pop_labels)

# check He vs fst
plot(my_fst$He, my_fst$FST)

#### Evaluating OutFLANK with trimmed SNPs ####
my_fst_pruned <- my_fst %>% filter(LocusName %in% row.names(snp_pruned))
out_trim <- OutFLANK(my_fst_pruned, NumberOfSamples = 2, 
                     RightTrimFraction = 0.1, LeftTrimFraction = 0.05,
                     qthreshold = 0.01, Hmin = 0.05)

# calculate FST for the unpruned data
# and find outliers

fst_wht_cbr <- pOutlierFinderChiSqNoCorr(my_fst, Fstbar = out_trim$FSTNoCorrbar, 
                                         dfInferred = out_trim$dfInferred, qthreshold = 0.01, Hmin=0.05)

####################################
# fst: cmn vs. cbr
####################################

# target inds
inds_to_include <- meta_df %>%
  filter(cluster %in% c("cmn", "cbr")) %>%
  filter(!is.na(cluster)) %>%
  pull(Indiv)

# subset the master snp table
snp_tab_sub <- snp_tab[,(names(snp_tab) %in% inds_to_include)]

# grab the population labels from the meta data table
pop_labels <- data.frame(Indiv = names(snp_tab_sub)) %>%
  left_join(meta_df) %>%
  pull(cluster) %>%
  droplevels

my_fst <- MakeDiploidFSTMat(t(snp_tab_sub), locusNames = row.names(snp_tab_sub), popNames = pop_labels)

# check He vs fst
plot(my_fst$He, my_fst$FST)

#### Evaluating OutFLANK with trimmed SNPs ####
my_fst_pruned <- my_fst %>% filter(LocusName %in% row.names(snp_pruned))
out_trim <- OutFLANK(my_fst_pruned, NumberOfSamples = 2, 
                     RightTrimFraction = 0.1, LeftTrimFraction = 0.05,
                     qthreshold = 0.01, Hmin = 0.05)

# calculate FST for the unpruned data
# and find outliers

fst_cmn_cbr <- pOutlierFinderChiSqNoCorr(my_fst, Fstbar = out_trim$FSTNoCorrbar, 
                                         dfInferred = out_trim$dfInferred, qthreshold = 0.01, Hmin=0.05)

####################################
# split locus names
####################################

# split locus names to chr pos in outflank objects
# DRY fail but hey

fst_cmn_cbr <- fst_cmn_cbr %>%
  mutate(chr = gsub("_.*", "", LocusName)) %>%
  mutate(pos = as.numeric(gsub(".*_", "", LocusName))) %>%
  mutate(comparison = "cmn_cbr")

fst_wht_cmn <- fst_wht_cmn %>%
  mutate(chr = gsub("_.*", "", LocusName))%>%
  mutate(pos = as.numeric(gsub(".*_", "", LocusName)))%>%
  mutate(comparison = "wht_cmn")

fst_wht_cbr <- fst_wht_cbr %>%
  mutate(chr = gsub("_.*", "", LocusName))%>%
  mutate(pos = as.numeric(gsub(".*_", "", LocusName)))%>%
  mutate(comparison = "wht_cbr")

# join fst data sets:
fst_df <- bind_rows(fst_cmn_cbr, fst_wht_cmn, fst_wht_cbr)

# write out
saveRDS(fst_df, "data/stats/fst_outflank_snp.rds")

####################################
# compute PBS from FST values
####################################

# Yi et al. 2010
# PBS for tibetan population
# T = - log(1 – FST)
# PBS = (TTH + TTD − THD)/2 

# Vijay et al. 2017
# PBS_Pop1=(-log(1-FST(Pop1_Pop2)+
# (-log(1-FST(Pop1_Pop3)))-log(1-FST(Pop2_Pop3)/2

# goal here is to get PBS for the white stickleback
# pop1 = wht, pop2 = cmn, pop3 = cbr
# AND SO:

pbs_df <- fst_df %>% 
  select(comparison, chr, pos, FST) %>%
  spread(key = comparison, value = FST) %>%
  mutate(wht_cmn = ifelse(wht_cmn > 0, -log(1-wht_cmn), -log(1))) %>%
  mutate(wht_cbr = ifelse(wht_cbr > 0, -log(1-wht_cbr), -log(1))) %>%
  mutate(cmn_cbr = ifelse(cmn_cbr > 0, -log(1-cmn_cbr), -log(1))) %>%
  mutate(wht_pbs = (wht_cmn + wht_cbr - cmn_cbr)/2) 

# write out
saveRDS(pbs_df, "data/stats/pbs_df_snp.rds")

####################################
# 100k window summaries
####################################

# create a list of windows from all the possible positions 

window_size <- 100000

# 100k windows
fst_100k <- fst_df %>%
  mutate(w_pos2 = (((pos / window_size) %>% floor) + 1) * window_size)%>%
  mutate(w_pos1 = w_pos2 - (window_size - 1)) %>%
  group_by(comparison, chr, w_pos1, w_pos2) %>%
  summarise(mean_fst = sum(T1)/sum(T2), mean_he = mean(He), n_outliers = sum(OutlierFlag, na.rm = TRUE))

# write out
saveRDS(fst_100k, "data/stats/fst_100k.rds")

