# apply site-level filtration
# KMS oct 2018

library("tidyverse")
library("ggridges")

##################################  
# determine and apply filters
################################## 

# read in site level stats
stat_df <- readRDS("metadata/site_stats.rds")

ind_df <- readRDS("metadata/ind_stats.rds")
filt_inds <- ind_df %>% 
  filter(CHROM != "chrY" & CHROM != "chrXIX") %>%
  group_by(Indiv) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) %>% 
  filter(n_genotyped_sites > 500) %>%
  filter(mean_depth > 10) %>%
  pull(Indiv)

# compute mafs, allele ratios, representation etc.
stat_df <- stat_df %>%
  ungroup %>%
  distinct %>%
  filter(CHROM != "chrXIX") %>%
  filter(CHROM != "chrY") %>%
  mutate(maf = case_when((n_01/2) + n_11 < (n_00 + n_01/2) ~ (n_01/2) + n_11,
                   (n_01/2) + n_11 > (n_00 + n_01/2) ~ (n_00 + n_01/2),
                   TRUE ~ (n_01/2) + n_11)) %>%
  mutate(excess_het = n_01 - ((1-maf) * maf * 2)) %>%
  mutate(prop_genotyped = n_genotyped / max(n_genotyped)) %>%
  mutate(summed_AR = (sum_allele1 / (sum_allele2 + sum_allele1))) %>%
  mutate(summed_AR_z_score = scale(summed_AR, scale = TRUE, center = TRUE)) %>%
  mutate(marker_id = 1:length(maf))

# quick plots for filtering
stat_df %>%
  ggplot(aes(x = summed_AR_z_score))+
  geom_histogram()

stat_df %>%
  ggplot(aes(x = prop_genotyped))+
  geom_histogram()+
  facet_wrap(~CHROM)

  

# the primary site-level filters:
# depth, representation, maf, heterozygostiy, allele ratio
# also excess heterozygosity (HWE violations)
stat_df_filt <- stat_df %>% 
  filter(maf > 0.01) %>%
  filter(mean_DP < 100) %>%
  filter(abs(summed_AR_z_score) <= 1.0) %>%
  filter(n_01 < 0.75) %>%
  filter(prop_genotyped > 0.50)

nrow(stat_df_filt)
  
# 151k SNPs pre filter
# 12k SNPs post filter 

# the mckinney et al type plot
# comparing filtered/unfiltered 

#stat_df %>%
stat_df %>%
  mutate(retained = stat_df$marker_id %in% stat_df_filt$marker_id) %>%
  #sample_frac(0.2) %>%
  ggplot(aes(x = n_01, y = summed_AR_z_score, color = retained)) + 
    geom_point(alpha = 0.2, shape = 16) +
    facet_wrap(~CHROM) +
    scale_color_manual(values = c("grey", "red"))

# histograms of all the site stats
# filtering was major!
stat_df %>%
  mutate(filtered = stat_df$marker_id %in% stat_df_filt$marker_id) %>%
  gather(key = stat, value = value, -CHROM, -POS,  -filtered) %>%
  ggplot(aes(x = value, fill = filtered))+
  geom_histogram()+
    facet_wrap(~stat, scales = "free")+
  scale_fill_manual(values = c("grey", "red"))

# reclaim memory
rm(stat_df)

##################################  
# apply filters
################################## 

# function for per-chromosome filtering

chr_files <- list.files("data/vcf/2020_pipeline/tidy_snps/", pattern = ".rds", full.names = TRUE)
chr_files <- chr_files[!grepl("_chrXIX_", chr_files)] # no sex chromosome!
chr_files <- chr_files[!grepl("_chrY_", chr_files)] # no sex chromosome!
chr_files <- chr_files[!grepl("_chrM_", chr_files)] # no sex chromosome!

apply_site_filters <- function(chr_file){
  
  print(chr_file)
  
  chr_df <- readRDS(chr_file)
  #chrom <- gsub(".*/", "", chr_file) %>% gsub("\\..*", "", .)
  chrom <- gsub(".*/", "", chr_file) %>% gsub("_GATK_raw_snps.rds", "", .) %>% gsub("whtstbk_", "", .)
  
  print(chrom)
  
  # determine target sites
  chr_sites <- stat_df_filt %>%
    filter(CHROM == chrom) %>%
    pull(POS)
  
  print(paste0("filtering for ", length(chr_sites), " sites..."))
  
  # apply filter
  chr_df <- chr_df %>%
    filter(POS %in% chr_sites) %>%
    filter(Indiv %in% filt_inds)

  print(paste0("resulted in ", length(chr_df %>% pull(POS) %>% unique), " sites..."))
  
  chr_df
  
}

filtered_sites_df <- lapply(chr_files, apply_site_filters)
filtered_sites_df <- bind_rows(filtered_sites_df)

# read in meta data
meta_df <- read.csv("metadata/mega_meta.csv")
meta_df <- meta_df %>%
  mutate(Indiv = gsub("whtstbk_gbs_2012_brds_", "", id)) %>%
  mutate(Indiv = gsub("_", ".", Indiv)) %>%
  mutate(Indiv = case_when(year == "2012" ~ paste0(Indiv, "_2013"),
                           year == "2014" ~ paste0(Indiv, "_2015"),
                           TRUE ~ Indiv)) %>%
  mutate(Indiv = ifelse(Indiv == "CL64.2_2015", "CL64-2_2015", Indiv))

write.csv(meta_df, "metadata/mega_meta_2020.csv")

meta_df %>% filter(year == 2014 ,!is.na(cluster), cluster != "dk") %>%
  arrange(cluster, Indiv) %>%
  select(Indiv, cluster) %>%
  write.table("metadata/pixy_popfile.txt", row.names = FALSE, quote = FALSE)

# should return "integer(0)"
which(!((filtered_sites_df$Indiv %>% unique) %in% meta_df$Indiv))

filtered_sites_df <- left_join(filtered_sites_df, meta_df, by = "Indiv") %>%
  select(Indiv, pop, year, cluster, sex, everything()) %>%
  distinct()

filtered_sites_df %>%
  select(CHROM, POS) %>%
  distinct %>% 
  group_by(CHROM) %>%
  tally %>%
  View

saveRDS(filtered_sites_df , "data/geno_df_2020_filtered.rds")

nrow(filtered_sites_df)

# add in X and Y data


