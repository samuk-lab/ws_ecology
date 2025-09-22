# filter joined vcf files

library("tidyverse")
library("tidyr")

chr_files <- list.files("data/vcf/2020_pipeline/tidy_snps", full.names = TRUE)


build_site_stat_dfs <- function(chr_file){
  
  print(chr_file)
  
  # read in the joined chromosome data
  chr_df <- readRDS(chr_file) %>%
    filter(!grepl(".*,.*,.*", gt_AD)) %>%
    separate(gt_AD, into = c("allele1", "allele2")) %>%
    mutate(allele1 = as.numeric(allele1), allele2 = as.numeric(allele2)) %>%
    mutate(depth = allele1 + allele2)
  
  # apply heterozygostity/allele ratio filter
  # operates at site level
  # don't apply to sex chrom!
  # after McKinney et al. 2017
  
  # polyal_hets <- chr_df %>% 
  #   select(CHROM, POS, gt_AD, gt_GT) %>%
  #   filter(gt_GT == "0/1") %>% filter(grepl(".*,.*,.*", gt_AD)) 
  
  # compute biallelic alelle ratios
  bial_ratios <- chr_df %>% 
    select(CHROM, POS, gt_GT, allele1, allele2, depth) %>%
    filter(gt_GT == "0/1") %>%
    mutate(allele_ratio = allele1 / (allele1 + allele2)) %>%
    group_by(CHROM, POS) %>%
    summarise(mean_AR = mean(allele_ratio, na.rm = TRUE), sum_allele1 = sum(allele1, na.rm = TRUE), sum_allele2 = sum(allele2, na.rm = TRUE))
  
  hobs_df <- chr_df %>% 
    group_by(CHROM, POS) %>%
    summarise(n_00 = mean(gt_GT == "0/0", na.rm = TRUE), 
              n_01 = mean(gt_GT == "0/1", na.rm = TRUE), 
              n_11 = mean(gt_GT == "1/1", na.rm = TRUE), 
              n_genotyped = sum(!is.na(gt_GT)), mean_DP = mean(depth, na.rm = TRUE))
  
  left_join(hobs_df, bial_ratios, by = c("CHROM", "POS"))

}

stat_df <- lapply(chr_files, build_site_stat_dfs)
stat_df <- bind_rows(stat_df)

saveRDS(stat_df, "metadata/site_stats.rds")

build_ind_stat_dfs <- function(chr_file){
  
  print(chr_file)
  
  # read in the joined chromosome data
  chr_df <- readRDS(chr_file) %>%
    filter(!grepl(".*,.*,.*", gt_AD)) %>%
    separate(gt_AD, into = c("allele1", "allele2")) %>%
    mutate(allele1 = as.numeric(allele1), allele2 = as.numeric(allele2)) %>%
    mutate(depth = allele1 + allele2)
  
  # apply heterozygostity/allele ratio filter
  # operates at site level
  # don't apply to sex chrom!
  # after McKinney et al. 2017
  
  bial_hets <- chr_df %>% 
    select(CHROM, POS, Indiv, gt_GT, allele1, allele2, depth) %>%
    filter(gt_GT == "0/1")
  
  # polyal_hets <- chr_df %>% 
  #   select(CHROM, POS, gt_AD, gt_GT) %>%
  #   filter(gt_GT == "0/1") %>% filter(grepl(".*,.*,.*", gt_AD)) 
  
  # compute biallelic alelle ratios
  bial_ratios <- bial_hets %>% 
    mutate(allele_ratio = allele1 / (allele1 + allele2)) %>%
    group_by(CHROM, Indiv) %>% summarise(mean_AR = mean(allele_ratio, na.rm = TRUE))
  
  hobs_df <- chr_df %>% 
    group_by(CHROM, Indiv) %>%
    summarise(n_genotyped_sites = sum(!is.na(gt_GT)), 
              hobs = sum(gt_GT == "0/1", na.rm = TRUE),
              mean_depth = mean(depth, na.rm = TRUE), 
              sd_depth = sd(depth, na.rm = TRUE)) %>%
    arrange(CHROM, Indiv)
  
  left_join(hobs_df, bial_ratios, by = c("CHROM", "Indiv"))
  
}

ind_df <- lapply(chr_files, build_ind_stat_dfs)
ind_df <- bind_rows(ind_df )

saveRDS(ind_df, "metadata/ind_stats.rds")

# plotting results

ind_df <- readRDS("metadata/ind_stats.rds")
stat_df <- readRDS("metadata/site_stats.rds")


# individual stats
ind_df %>% 
  filter(CHROM != "chrY") %>%
  group_by(Indiv) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) %>% 
  filter(n_genotyped_sites > 500) %>%
  filter(mean_depth > 10) %>% 
  #filter(n_genotyped_sites > 2000) %>%
  gather(key = stat, value = value, -Indiv) %>%
  ggplot(aes(x = value))+
  geom_histogram()+
  facet_wrap(~stat, scales = "free")
  #facet_grid(CHROM~stat, scales = "free")
    
# site stats
stat_df %>%
  filter(CHROM != "chrY") %>%
  gather(key = stat, value = value, -CHROM, -POS) %>%
  ggplot(aes(x = value))+
  geom_histogram()+
  facet_wrap(~stat, scales = "free")
  #facet_grid(CHROM~stat, scales = "free")

# why is n_genotyped buimodal?
stat_df %>%
  #filter(CHROM != "chrY") %>%
  ggplot(aes(x= POS, y = n_genotyped))+
  geom_point(size = 0.5)+
  facet_wrap(~CHROM, scales = "free")


# why is n_genotyped buimodal?
stat_df %>%
  #filter(CHROM != "chrY") %>%
  ggplot(aes(x= n_genotyped))+
  geom_point(size = 0.5)+
  facet_wrap(~CHROM, scales = "free")

