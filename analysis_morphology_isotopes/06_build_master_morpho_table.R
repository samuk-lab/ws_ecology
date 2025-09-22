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

# initial joint of meta_df to 2014 data

# initialize pheno_df
pheno_df <- meta_df

for (i in pheno_dfs){
  
  pheno_df <- full_join(pheno_df, i, by = c("id"))
}

# fix missing and year info and remove duplicated columns
pheno_df <- pheno_df %>%
  mutate(geno_sex = as.character(geno_sex)) %>%
  mutate(sex = as.character(sex)) %>%
  mutate(year = ifelse(is.na(year), 2014, year))

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
# make master table
######################################################################
pheno_long <- pheno_df %>%
  gather(key = trait, value = value, -id, -sex, -geno_sex, -joint_sex, -pop, -year, -region, -species, -cluster, -sequenced) %>%
  mutate(species = ifelse(!is.na(cluster), as.character(cluster), as.character(species))) %>%
  mutate(species = ifelse(species == "wht", "white", "common"))

pheno_long_sl <- pheno_df %>%
  gather(key = trait, value = value, -id, -sex, -geno_sex, -joint_sex, -pop, -year, -region, -species, -cluster, -sequenced, -std_length) %>%
  mutate(species = ifelse(!is.na(cluster), as.character(cluster), as.character(species))) %>%
  mutate(species = ifelse(species == "wht", "white", "common"))

pheno_long_sl %>%
  filter(grepl("PC", trait)) %>%
  ggplot(aes(x = species, y =value))+
  geom_jitter()+
  stat_summary(fun.y=median, fun.ymin=median, fun.ymax=median, 
               geom="crossbar", width = 0.9, color = "black", fatten = 2)+
  facet_wrap(~trait)

######################################################################
# basic stats
######################################################################

# compute basic summary stats for each trait
sum_stats <- pheno_long %>%
  select(-id, -sex, -geno_sex, -joint_sex, -pop, -year, -region, -cluster, -sequenced) %>%
  .[complete.cases(.),] %>%
  group_by(trait, species) %>%
  summarise_each(funs(mean, median, sd, length))

######################################################################
# size corrected fits
######################################################################

# size-corrected fits
# gaussian and poisson fits (worked out later)
sl_fits <- pheno_long_sl %>% 
  select(-id, -sex, -geno_sex, -year, -joint_sex, -region, -cluster, -sequenced) %>%
  .[complete.cases(.),] %>%  
  group_by(trait) %>%
  do(fit_lm = anova(lm(value ~ std_length + species, .)), fit_glm = anova(glm(abs(value) ~ std_length + species, ., family = poisson), test = "Chisq"))

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

######################################################################
# unsize corrected fits
######################################################################

# unsize-corrected fits
# gaussian and poisson fits (worked out later)
unsl_fits <- pheno_long %>% 
  select(-id, -sex, -geno_sex, -joint_sex, -year, -region, -cluster, -sequenced) %>%
  .[complete.cases(.),] %>%
  group_by(trait) %>%
  do(fit_lm = anova(lm(value ~ species, .)), fit_glm = anova(glm(abs(value) ~ species, ., family = poisson), test = "Chisq"))

# gaussian fits
unsl_tab_gauss <- unsl_fits %>%
  tidy(fit_lm)  %>%
  mutate(statistic = round(statistic, 2)) %>%
  ungroup %>%
  select(term, trait, p.value, df, statistic) %>%
  unite(stats, statistic, df, p.value) %>%
  spread(key = term, value = stats)  %>%
  mutate(resid_df = gsub("NA_|_NA", "", Residuals)) %>%
  select(-Residuals)

#poisson fits
unsl_tab_poiss <- unsl_fits %>%
  tidy(fit_glm) %>%
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

# grab uncorrected gaussian anova results
gauss_anova_table <- unsl_tab_gauss %>%
  separate(species, c("species_f", "species_num_df", "species_p_value"), sep = "_", convert = TRUE) %>%
  rename(species_resid_df = resid_df)

# join in gauss anova and prettify f statistic info
table_stats <- left_join(table_stats, gauss_anova_table) %>%
  mutate(species_f = paste0("F",species_num_df,",",species_resid_df," = ",species_f)) %>%
  select(-species_num_df, -species_resid_df)

# grab corrected gaussian anova results
gauss_anova_table_sl <- sl_tab_gauss %>%
  separate(species, c("species_f_sl", "species_num_df", "species_p_value_sl"), sep = "_", convert = TRUE) %>%
  separate(std_length, c("sl_f", "sl_num_df", "sl_p_value"), sep = "_", convert = TRUE) %>%
  mutate(species_f_sl = paste0("F",species_num_df,",",resid_df," = ",species_f_sl)) %>%
  mutate(sl_f = paste0("F",sl_num_df,",",resid_df," = ",sl_f)) %>%
  select(-species_num_df, -resid_df, -sl_num_df) %>%
  select(trait, sl_f, sl_p_value, species_f_sl, species_p_value_sl)

# join in gauss anova and prettify f statistic info
table_stats <- left_join(table_stats, gauss_anova_table_sl)

######### POISSON MODELS

#Deviance, Resid..Dev, df, Resid..Df, Pr..Chi.

# grab uncorrected poisson anode results
poiss_anova_table <- unsl_tab_poiss %>%
  separate(species, c("species_poiss_deviance", "species_poiss_resid_deviance", 
                      "species_poiss_df", "species_poiss_resid_df", "species_poiss_unsl_pvalue"), sep = "_", convert = TRUE) %>%
  mutate(species_poiss_unsl_d = paste0("D", species_poiss_df, ",",species_poiss_resid_df," = ",species_poiss_deviance)) %>%
  select(trait, matches("unsl"), species_poiss_unsl_pvalue)

# join in poisson and prettify d statistic info
table_stats <- left_join(table_stats, poiss_anova_table)

# grab corrected poisson anova results
poiss_anova_table_sl <- sl_tab_poiss %>%
  separate(species, c("species_poiss_deviance_sl", "species_poiss_resid_deviance_sl", 
                      "species_poiss_df_sl", "species_poiss_resid_df_sl", "species_poiss_pvalue_sl"), sep = "_", convert = TRUE) %>%
  separate(std_length, c("sl_poiss_deviance", "sl_poiss_resid_deviance", 
                      "sl_poiss_df", "sl_poiss_resid_df", "sl_poiss_pvalue"), sep = "_", convert = TRUE) %>%
  mutate(species_poiss_sl_d = paste0("D", species_poiss_df_sl, ",",species_poiss_resid_df_sl," = ",species_poiss_deviance_sl)) %>%
  mutate(sl_poiss_d = paste0("D", sl_poiss_df, ",",sl_poiss_resid_df," = ",sl_poiss_deviance)) %>%
  select(trait, matches("_d$"), species_poiss_pvalue_sl, sl_poiss_pvalue)
  

# join in corrected poisson and prettify f statistic info
table_stats <- left_join(table_stats, poiss_anova_table_sl)

table_stats$test_type <- table_stats$trait %in% c("egg_number", "raker_long_count", "plate_count", "raker_short_count")

table_stats$test_species_uncor <- ifelse(table_stats$test_type, table_stats$species_poiss_unsl_d, table_stats$species_f)
table_stats$pval_species_uncor <- ifelse(table_stats$test_type, table_stats$species_poiss_unsl_pvalue, table_stats$species_p_value)
table_stats$test_std_length <- ifelse(table_stats$test_type, table_stats$sl_poiss_d, table_stats$sl_f) 
table_stats$pval_std_length <- ifelse(table_stats$test_type, table_stats$sl_poiss_pvalue, table_stats$sl_p_value) 
table_stats$test_species_cor <- ifelse(table_stats$test_type, table_stats$species_poiss_sl_d, table_stats$species_f_sl)  
table_stats$pval_species_cor <- ifelse(table_stats$test_type, table_stats$species_poiss_pvalue_sl, table_stats$species_p_value_sl) 
  
table_stats$test_type <- ifelse(table_stats$test_type == TRUE, "Poisson, ANODE", "Gaussian, ANOVA")

table_stats <- table_stats %>%
  select(-matches("^species|^sl|_std_"))

pvals <- c(table_stats$pval_species_cor, table_stats$pval_species_uncor)
#pvals <- pvals[!is.na(pvals)]

#pvals %>% round (digits = 4)
p_adj <- p.adjust(pvals, method = "BH") %>% round (digits = 8)

table_stats$pval_species_cor <- p_adj[1:length(table_stats$pval_species_cor)]
table_stats$pval_species_uncor <- p_adj[(length(p_adj) - length(table_stats$pval_species_uncor) +1) : length(p_adj)]

write.table(table_stats, "figures/Table1.txt", row.names = FALSE, quote = FALSE, sep = "\t")
