# sex chromosome recaluclation of site statistics

library("tidyverse")
library("tidyr")

ind_df <- readRDS("metadata/ind_stats.rds")

filt_inds <- ind_df %>% 
  filter(CHROM != "chrY") %>%
  group_by(Indiv) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) %>% 
  filter(n_genotyped_sites > 500) %>%
  filter(mean_depth > 10) %>%
  pull(Indiv)

# chrXIX (the X chromosome)

chr_file <- list.files("data/vcf/2020_pipeline/tidy_snps", full.names = TRUE, pattern = "chrXIX")

print(chr_file)

# read in the joined chromosome data
chr_df <- readRDS(chr_file) %>%
  filter(!grepl(".*,.*,.*", gt_AD)) %>%
  separate(gt_AD, into = c("allele1", "allele2")) %>%
  mutate(allele1 = as.numeric(allele1), allele2 = as.numeric(allele2)) %>%
  mutate(depth = allele1 + allele2)

head(chr_df)

# get sex data
meta_df <- read.csv("metadata/mega_meta.csv", h = T)

# check indivs are all in the meta list
meta_df <- meta_df %>%
  mutate(Indiv = gsub("whtstbk_gbs_2012_brds_", "", id)) %>%
  mutate(Indiv = gsub("_", ".", Indiv)) %>%
  mutate(Indiv = case_when(year == "2012" ~ paste0(Indiv, "_2013"),
                           year == "2014" ~ paste0(Indiv, "_2015"),
                           TRUE ~ Indiv)) %>%
  mutate(Indiv = ifelse(Indiv == "CL64.2_2015", "CL64-2_2015", Indiv))

# should return "integer(0)"
chr_df$Indiv[which(!((chr_df$Indiv %>% unique) %in% meta_df$Indiv))]

# join in sex data
# remove males!
chr_df <- meta_df %>%
  select(Indiv, pop, year, cluster, sex) %>%
  left_join(chr_df, .) %>%
  filter(sex == "F") %>%
  filter(!is.na(sex))
  
# compute biallelic alelle ratios for hets
bial_ratios <- chr_df %>% 
  select(CHROM, POS, gt_GT, allele1, allele2, depth) %>%
  filter(gt_GT == "0/1") %>%
  mutate(allele_ratio = allele1 / (allele1 + allele2)) %>%
  group_by(CHROM, POS) %>%
  summarise(mean_AR = mean(allele_ratio, na.rm = TRUE), 
            sum_allele1 = sum(allele1, na.rm = TRUE), sum_allele2 = sum(allele2, na.rm = TRUE))

hobs_df <- chr_df %>% 
  group_by(CHROM, POS) %>%
  summarise(n_00 = mean(gt_GT == "0/0", na.rm = TRUE), 
            n_01 = mean(gt_GT == "0/1", na.rm = TRUE), 
            n_11 = mean(gt_GT == "1/1", na.rm = TRUE), 
            n_genotyped = sum(!is.na(gt_GT)), mean_DP = mean(depth, na.rm = TRUE))


chr19_df <- left_join(hobs_df, bial_ratios, by = c("CHROM", "POS")) %>% ungroup

chr19_df <- chr19_df %>%
  mutate(maf = case_when((n_01/2) + n_11 < (n_00 + n_01/2) ~ (n_01/2) + n_11,
                         (n_01/2) + n_11 > (n_00 + n_01/2) ~ (n_00 + n_01/2),
                         TRUE ~ (n_01/2) + n_11)) %>%
  mutate(summed_AR = (sum_allele1 / (sum_allele2 + sum_allele1))) %>%
  mutate(summed_AR_z_score = scale(summed_AR, scale = TRUE, center = TRUE))


chr19_df %>%
  filter(n_genotyped > 100) %>% 
  #filter(maf > 0.01) %>%
  ggplot(aes(x = n_01, y = summed_AR_z_score, color = n_genotyped))+
  geom_point(size = 2) +
  scale_color_viridis_c()

chr19_df %>%
  ggplot(aes(x = POS, y = n_11))+
  geom_point()


# build list of sites to retain for females
chr19_sites_F <- chr19_df %>% 
  ungroup %>%
  distinct %>%
  mutate(maf = case_when((n_01/2) + n_11 < (n_00 + n_01/2) ~ (n_01/2) + n_11,
                         (n_01/2) + n_11 > (n_00 + n_01/2) ~ (n_00 + n_01/2),
                         TRUE ~ (n_01/2) + n_11)) %>%
  mutate(excess_het = n_01 - ((1-maf) * maf * 2)) %>%
  mutate(prop_genotyped = n_genotyped / max(n_genotyped)) %>%
  mutate(summed_AR = (sum_allele1 / (sum_allele2 + sum_allele1))) %>%
  mutate(summed_AR_z_score = scale(summed_AR, scale = TRUE, center = TRUE)) %>%
  mutate(marker_id = 1:length(maf)) %>% 
  filter(maf > 0.01) %>%
  filter(mean_DP < 100) %>%
  filter(abs(summed_AR_z_score) <= 1.0) %>%
  filter(n_01 < 0.6) %>%
  filter(prop_genotyped > 0.50)

  

# filter chr19 sites for females only
chr19_df_F <- chr_df %>%
  filter(sex  == "F" & POS %in% chr19_sites_F$POS) %>%
  select(Indiv, pop, year, cluster, sex, everything())


saveRDS(chr19_df_F, "data/geno_df_XIX_2020_filtered.rds")


# chrY
chr_file <- list.files("data/vcf/2020_pipeline/tidy_snps", full.names = TRUE, pattern = "chrY")

print(chr_file)

# read in the joined chromosome data
chrY_df <- readRDS(chr_file) %>%
  filter(!grepl(".*,.*,.*", gt_AD)) %>%
  separate(gt_AD, into = c("allele1", "allele2")) %>%
  mutate(allele1 = as.numeric(allele1), allele2 = as.numeric(allele2)) %>%
  mutate(depth = allele1 + allele2)

head(chrY_df)

# get sex data
meta_df <- read.csv("metadata/mega_meta.csv", h = T)

# check indivs are all in the meta list
meta_df <- meta_df %>%
  mutate(Indiv = gsub("whtstbk_gbs_2012_brds_", "", id)) %>%
  mutate(Indiv = gsub("_", ".", Indiv)) %>%
  mutate(Indiv = case_when(year == "2012" ~ paste0(Indiv, "_2013"),
                           year == "2014" ~ paste0(Indiv, "_2015"),
                           TRUE ~ Indiv)) %>%
  mutate(Indiv = ifelse(Indiv == "CL64.2_2015", "CL64-2_2015", Indiv))

# should return "character(0)"
chrY_df$Indiv[which(!((chrY_df$Indiv %>% unique) %in% meta_df$Indiv))]

# join in sex data
# remove females!
chrY_df <- meta_df %>%
  select(Indiv, pop, year, cluster, sex) %>%
  left_join(chrY_df, .) %>%
  filter(!is.na(sex)) %>%
  filter(sex == "M")


# compute biallelic alelle ratios for hets
bial_ratios <- chrY_df %>% 
  select(CHROM, POS, gt_GT, allele1, allele2, depth) %>%
  filter(gt_GT == "0/1") %>%
  mutate(allele_ratio = allele1 / (allele1 + allele2)) %>%
  group_by(CHROM, POS) %>%
  summarise(mean_AR = mean(allele_ratio, na.rm = TRUE), 
            sum_allele1 = sum(allele1, na.rm = TRUE), sum_allele2 = sum(allele2, na.rm = TRUE))

hobs_df <- chrY_df %>% 
  group_by(CHROM, POS) %>%
  summarise(n_00 = mean(gt_GT == "0/0", na.rm = TRUE), 
            n_01 = mean(gt_GT == "0/1", na.rm = TRUE), 
            n_11 = mean(gt_GT == "1/1", na.rm = TRUE), 
            n_genotyped = sum(!is.na(gt_GT)), mean_DP = mean(depth, na.rm = TRUE))


chrY_df_stats <- left_join(hobs_df, bial_ratios, by = c("CHROM", "POS")) %>% ungroup

chrY_df_stats <- chrY_df_stats %>%
  mutate(maf = case_when((n_01/2) + n_11 < (n_00 + n_01/2) ~ (n_01/2) + n_11,
                         (n_01/2) + n_11 > (n_00 + n_01/2) ~ (n_00 + n_01/2),
                         TRUE ~ (n_01/2) + n_11)) %>%
  mutate(summed_AR = (sum_allele1 / (sum_allele2 + sum_allele1))) %>%
  mutate(summed_AR_z_score = scale(summed_AR, scale = TRUE, center = TRUE))


chrY_df_stats %>%
  filter(n_genotyped > 20) %>%
  #filter(maf > 0.01) %>%
  ggplot(aes(x = n_01, y = summed_AR_z_score, color = n_genotyped))+
  geom_point(size = 2) +
  scale_color_viridis_c()


chrY_df_stats %>%
  filter(maf > 0.01) %>%
  gather(key = stat, value = value, -CHROM, -POS) %>%
  ggplot(aes(x = value))+
  geom_histogram()+
  facet_wrap(~stat, scales = "free")

chrY_df_stats %>%
  ggplot(aes(x = POS, y = mean_AR, color = n_01))+
  geom_point()+
  scale_color_viridis_c()

# build list of sites to retain for males
chrY_sites_M <- chrY_df_stats %>% 
  ungroup %>%
  distinct %>%
  mutate(maf = case_when((n_01/2) + n_11 < (n_00 + n_01/2) ~ (n_01/2) + n_11,
                         (n_01/2) + n_11 > (n_00 + n_01/2) ~ (n_00 + n_01/2),
                         TRUE ~ (n_01/2) + n_11)) %>%
  mutate(excess_het = n_01 - ((1-maf) * maf * 2)) %>%
  mutate(prop_genotyped = n_genotyped / max(n_genotyped)) %>%
  mutate(summed_AR = (sum_allele1 / (sum_allele2 + sum_allele1))) %>%
  mutate(summed_AR_z_score = scale(summed_AR, scale = TRUE, center = TRUE)) %>%
  mutate(marker_id = 1:length(maf)) %>% 
  filter(maf > 0.01) %>% 
  filter(mean_DP < 100)



# filter chrY sites for males only
chrY_df_M <- chrY_df %>%
  filter(POS %in% chrY_sites_M$POS) %>%
  select(Indiv, pop, year, cluster, sex, everything())


saveRDS(chrY_df_M, "data/geno_df_Y_2020_filtered.rds")

  
