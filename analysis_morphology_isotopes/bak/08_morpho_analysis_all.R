######################################################################
# Combine all phenotypic data into master table :o
# Kieran Samuk - Apr 2016
######################################################################


######################################################################
# Libraries
######################################################################

library("dplyr")
library("tidyr")
library("broom")
library("ggplot2")
library("ggthemes")

list.files("functions", full.names = TRUE) %>% sapply(source) %>% invisible

select <- dplyr::select #allows dplyr to use select when MASS package is loaded

# the white stickleback palatte
whtstbk_palatte <- c("#008FD5", "#BDBDBD")

setwd("analysis_morphology")

######################################################################
# input data
######################################################################

# get file names for files to be combined
pheno_files <- list.files("data/collated", pattern = "raw", full.names = TRUE)

# meta data file from genotypic analysis
meta_df <- read.csv("metadata/mega_meta.csv")

meta_df <- meta_df %>%
  rename(geno_sex = sex) %>%
  filter(!(pop %in% c("DK", "LC")))

meta_df$sequenced <- 1

######################################################################
# joining data frames
######################################################################

# all but raker data are individuals from 2014 (so can join on ids directly)
files_2014 <- grep("raker", pheno_files, invert = TRUE, value = TRUE)

pheno_dfs <- lapply(files_2014, read.table, header = TRUE, stringsAsFactors = FALSE)

# fix PCA data for joining
pheno_dfs[[3]] <- pheno_dfs[[3]] %>%
  mutate(year = gsub(".*_", "", id) %>% as.numeric) %>%
  mutate(id = gsub("_.*", "", id))

pca_dat <- pheno_dfs[[3]]
pheno_dfs <-pheno_dfs[-3]

# initialize pheno_df with genotype/pop metadata
pheno_df <- meta_df

for (i in pheno_dfs){
  
  pheno_df <- full_join(pheno_df, i, by = c("id"))
}

# fix missing and year info and remove duplicated columns
pheno_df <- pheno_df %>%
  mutate(geno_sex = as.character(geno_sex)) %>%
  mutate(sex = as.character(sex)) %>%
  mutate(year = ifelse(is.na(year), 2014, year))

pheno_df <- left_join(pheno_df, pca_dat) 
# fix sex info
pheno_df$joint_sex <- paste0(pheno_df$geno_sex, pheno_df$sex) %>% 
  gsub("NA*", "", .) %>%
  substr(1,1) %>%
  ifelse(. == "", NA, .)

# scrub 2012 ids of weird slug
pheno_df <- pheno_df %>%
  mutate(id = gsub("whtstbk_gbs_2012_brds_", "", id))

# join in 2012 raker data (carefully)
raker_file <- grep("raker", pheno_files, value = TRUE)
raker_df <- read.table(raker_file, header = TRUE, stringsAsFactors = FALSE) 

raker_df$species <- ifelse(raker_df$species == "common", "cmn", "wht")
raker_df$year <- 2012

pheno_df <- full_join(pheno_df, raker_df)

# add in population codes
pheno_df <- pheno_df %>%
  mutate(pop = gsub("[^A-Z]*", "", id) %>% substr(1,2))

write.table(pheno_df, "data/pheno_df_master.txt", quote = FALSE, row.names = FALSE)

######################################################################
# subset data and plot
######################################################################

target_traits <- c("luminance_mean", "egg_number", "egg_diameter_mean", "testis_weight", "body_depth2", 
                   "std_length", "spine_dorsal1", "spine_dorsal2", "spine_dorsal3", 
                   "spine_pelvic", "ventral_rakers", "dorsal_rakers", "raker_long_count", 
                   "raker_short_count", "plate_count")

trait_type <- c("body", "gonad", "gonad", "gonad", "body", "body", "armor", "armor", "armor", "armor", "rakers", "rakers", "rakers", "rakers", "armor")

targ_df <- data.frame(trait = target_traits, trait_type)

pheno_long <- pheno_df %>%
  select(-geno_sex, -sex) %>%
  rename(sex = joint_sex) %>%
  gather(key = trait, value = value, -id, -sex, -pop, -year, -region, -species, -cluster, -sequenced) %>%
  mutate(species = ifelse(!is.na(cluster), as.character(cluster), as.character(species))) %>%
  mutate(species = ifelse(species == "wht", "white", "common")) %>%
  filter(trait %in% target_traits) %>%
  left_join(targ_df)



######################################################################
# basic stats
######################################################################

# compute basic summary stats for each trait
sum_stats <- pheno_long %>%
  select(-id, -sex,  -pop, -year, -region, -cluster, -sequenced, -trait_type) %>%
  .[complete.cases(.),] %>%
  group_by(trait, species) %>%
  summarise_each(funs(mean, median, sd, length))

######################################################################
# size corrected fits
######################################################################

pheno_long_sl <- pheno_df %>%
  gather(key = trait, value = value, -id, -sex, -geno_sex, -joint_sex, -pop, -year, -region, -species, -cluster, -sequenced, -std_length) %>%
  mutate(species = ifelse(!is.na(cluster), as.character(cluster), as.character(species))) %>%
  mutate(species = ifelse(species == "wht", "white", "common")) %>%
  filter(trait %in% target_traits)

# size-corrected fits
# gaussian and poisson fits (worked out later)
sl_fits <- pheno_long_sl %>% 
  left_join(targ_df) %>%
  filter(!(trait_type %in% c("rakers", "gonad"))) %>%
  select(-id, -year, -region, -cluster, -sequenced) %>%
  .[complete.cases(.),] %>%  
  group_by(trait) %>%
  do(fit_lm = anova(lm(value ~ std_length + species * sex, .)), fit_glm = anova(glm(abs(value) ~ std_length + species * sex, ., family = poisson), test = "Chisq"))

# gaussian fits
sl_tab_gauss <- sl_fits %>%
  tidy(fit_lm)  %>%
  mutate(statistic = round(statistic, 2)) %>%
  ungroup %>%
  select(term, trait, p.value, df, statistic) %>%
  unite(stats, statistic, df, p.value) %>%
  spread(key = term, value = stats)  %>%
  mutate(resid_df = gsub("NA_|_NA", "", Residuals)) %>%
  select(-Residuals)

#poisson fits
sl_tab_poiss <- sl_fits %>%
  tidy(fit_glm) %>%
  filter(term != "NULL") %>%
  mutate(Resid..Dev = round(Resid..Dev, 4), Deviance= round(Deviance, 4)) %>%
  ungroup %>% 
  select(term, trait, Deviance, Resid..Dev, df, Resid..Df, p.value) %>%
  unite(stats, Deviance, Resid..Dev, df, Resid..Df, p.value) %>%
  spread(key = term, value = stats)

# gonadal traits

sl_fits_gon <- pheno_long_sl %>% 
  left_join(targ_df) %>%
  filter(trait_type %in% c("gonad")) %>%
  select(-id, -year, -region, -cluster, -sequenced) %>%
  .[complete.cases(.),] %>%  
  group_by(trait) %>%
  do(fit_lm = anova(lm(value ~ std_length + species, .)), fit_glm = anova(glm(abs(value) ~ std_length + species, ., family = poisson), test = "Chisq"))

# gaussian fits
sl_tab_gon <- sl_fits_gon %>%
  tidy(fit_lm)  %>%
  mutate(statistic = round(statistic, 2)) %>%
  ungroup %>%
  select(term, trait, p.value, df, statistic) %>%
  unite(stats, statistic, df, p.value) %>%
  spread(key = term, value = stats)  %>%
  mutate(resid_df = gsub("NA_|_NA", "", Residuals)) %>%
  select(-Residuals)

# raker
sl_fits_rak <- pheno_long_sl %>% 
  left_join(targ_df) %>%
  filter(grepl("raker", trait_type)) %>%
  select(-id, -year, -sex, -region, -cluster, -sequenced) %>%
  filter(!is.na(value)) %>%
  group_by(trait) %>%
  do(fit_glm = anova(glm(abs(value) ~ std_length + species, ., family = poisson), test = "Chisq")) %>%
  unnest %>% tidy

# poisson fits
sl_tab_rak <- sl_fits_rak %>%
  tidy(fit_glm)  %>%
  filter(term != "NULL") %>%
  mutate(Resid..Dev = round(Resid..Dev, 4), Deviance= round(Deviance, 4)) %>%
  ungroup %>% 
  select(term, trait, Deviance, Resid..Dev, df, Resid..Df, p.value) %>%
  unite(stats, Deviance, Resid..Dev, df, Resid..Df, p.value) %>%
  spread(key = term, value = stats)

######################################################################
# output nice table
######################################################################

# format summary stats and calculate cohen'sd
table_stats <- sum_stats %>%
  ungroup %>%
  unite(summary_stats, mean, sd, median, length) %>%
  spread(key = species, value = summary_stats) %>%
  separate(common, c("common_mean", "common_sd", "common_median", "common_n"), sep = "_", convert = TRUE) %>%
  separate(white, c("white_mean", "white_sd", "white_median", "white_n"), sep = "_", convert = TRUE) %>%
  mutate(pooled_sd = sqrt((((white_n - 1) * white_sd^2) + ((common_n - 1)*common_sd^2)) / (white_n + common_n - 2))) %>%
  mutate(cohens_d = (white_mean - common_mean) / pooled_sd) %>% 
  select(-pooled_sd)

# prettify and filter unwanted variables
table_stats <- table_stats[,-1] %>% 
  lapply(round, digits = 2) %>% 
  data.frame(trait = table_stats$trait, .) %>%
  filter(!grepl("mode|resid|csize|rgb", trait)) %>%
  select(-matches("median")) %>%
  arrange(desc(abs(cohens_d)))

######### GAUSSIAN MODELS

# grab corrected gaussian anova results
gauss_anova_table_sl <- sl_tab_gauss %>%
  separate(species, c("species_f", "species_num_df", "species_p_value"), sep = "_", convert = TRUE) %>%
  separate(`species:sex`, c("sxs_f", "sxs_num_df", "sxs_p_value"), sep = "_", convert = TRUE) %>%
  separate(std_length, c("sl_f", "sl_num_df", "sl_p_value"), sep = "_", convert = TRUE) %>%
  mutate(species_f = paste0("F",species_num_df,",",resid_df," = ",species_f)) %>%
  mutate(sl_f = paste0("F",sl_num_df,",",resid_df," = ",sl_f)) %>%
  mutate(sxs_f = paste0("F",sxs_num_df,",",resid_df," = ",sxs_f)) %>%
  select(-species_num_df, -resid_df, -sl_num_df, -sxs_num_df) %>%
  select(trait, sl_f, sl_p_value, species_f, species_p_value, species_f, sxs_f, sxs_p_value)

# join in gauss anova and prettify f statistic info
table_stats_norm <- left_join(table_stats %>% filter(trait %in% sl_tab_gauss$trait), gauss_anova_table_sl, by = "trait")

######### GONDAL TRAITS

# grab corrected gaussian anova results
gon_anova_table_sl <- sl_tab_gon %>%
  separate(species, c("species_f", "species_num_df", "species_p_value"), sep = "_", convert = TRUE) %>%
  separate(std_length, c("sl_f", "sl_num_df", "sl_p_value"), sep = "_", convert = TRUE) %>%
  mutate(species_f = paste0("F",species_num_df,",",resid_df," = ",species_f)) %>%
  mutate(sl_f = paste0("F",sl_num_df,",",resid_df," = ",sl_f)) %>%
  select(-species_num_df, -resid_df, -sl_num_df) %>%
  select(trait, sl_f, sl_p_value, species_f, species_p_value)

# join in gauss anova and prettify f statistic info
table_stats_gon <- left_join(table_stats %>% filter(trait %in% sl_tab_gon$trait), gon_anova_table_sl, by = "trait")

######### RAKERS

# grab corrected gaussian anova results
poiss_anova_table_rak <- sl_tab_rak%>%
  separate(species, c("species_poiss_deviance", "species_poiss_resid_deviance", 
                      "species_df", "species_resid_df", "species_p_value"), sep = "_", convert = TRUE) %>%
  separate(std_length, c("std_poiss_deviance", "std_poiss_resid_deviance", 
                         "std_df", "std_resid_df", "sl_p_value"), sep = "_", convert = TRUE) %>%
  mutate(species_f = paste0("D", species_df, ",",species_resid_df," = ", species_poiss_deviance)) %>% 
  mutate(sl_f = paste0("D", std_df, ",",std_resid_df," = ", std_poiss_deviance)) %>% 
  select(trait, sl_f, sl_p_value, species_f, species_p_value)

# join in gauss anova and prettify f statistic info
table_stats_rak <- left_join(table_stats %>% filter(trait %in% sl_tab_rak$trait), poiss_anova_table_rak)

# join all
table_stats <- bind_rows(table_stats_norm, table_stats_gon) %>% bind_rows(table_stats_rak)

table_stats$test_type <- ifelse(grepl("D1", table_stats$sl_f), "Poisson, ANODE", "Gaussian, ANOVA")

table_stats <- table_stats %>%
  select(-matches("^species|^sl|_std_"))

pvals <- c(table_stats$species_p_value, table_stats$sxs_p_value)
#pvals <- pvals[!is.na(pvals)]

#pvals %>% round (digits = 4)
p_adj <- p.adjust(pvals, method = "BH") %>% round (digits = 8)

table_stats$pval_species_cor <- p_adj[1:length(table_stats$pval_species_cor)]
table_stats$pval_species_uncor <- p_adj[(length(p_adj) - length(table_stats$pval_species_uncor) +1) : length(p_adj)]

write.table(table_stats, "figures/Table1.txt", row.names = FALSE, quote = FALSE, sep = "\t")

table_stats %>%
  mutate(common_se = common_sd/sqrt(common_n), white_se = white_sd/sqrt(white_n)) %>%
  select(trait, common_mean, common_se, white_mean, white_se) %>%
  gather(key = mean_sp, value = mean_value, -trait, -white_se, -common_se) %>%
  gather(key = se_sp, value = se_value, -trait, -mean_sp, -mean_value) %>%
  mutate(species = ifelse(grepl("white", mean_sp), "white", "common")) %>%
  select(-mean_sp, -se_sp) %>%
  ggplot(aes(x = trait, y = mean_value, color = species, group = species))+
  geom_point()

pheno_long$trait %>% unique

pheno_reorder <- c("std_length", "body_depth2", "luminance_mean", "plate_count",
                   "spine_dorsal1", "spine_dorsal2","spine_dorsal3","spine_pelvic",
                   "raker_short_count", "raker_long_count", "ventral_rakers", "dorsal_rakers",
                   "egg_number", "egg_diameter_mean", "testis_weight")

pheno_long %>%
  filter(!is.na(species)) %>%
  mutate(trait = factor(trait, levels = pheno_reorder)) %>%
  ggplot(aes(x = species, y = value, color = species))+
  geom_jitter(width = 0.25, size = 1)+
  stat_summary(fun.data = "mean_cl_normal", 
               geom = "pointrange", color = "black", pch = 19, size = 0.5)+
  facet_wrap(~trait, scales = "free")+
  theme_bw()+
  xlab("Species")+
  ylab("Value")+
  theme(legend.position = "none",
        strip.background = element_blank(),
        panel.grid = element_blank())+
  scale_color_manual(values = whtstbk_palatte)+
  scale_fill_manual(values = whtstbk_palatte)

