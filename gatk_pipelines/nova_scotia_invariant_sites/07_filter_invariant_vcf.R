library("tidyverse")
library("parallel")

# find the tidy invariant data 
tidy_files <-  list.files("by_chromo", pattern = ".rds", full.names = TRUE)

########################################
# vcf 2: the filtering
########################################
filter_tidy_rds <- function(tidy_file){
  
  message(tidy_file)
  
  vcf_all <- readRDS(tidy_file)
  
  vcf_all <- vcf_all %>%
    filter(!is.na(gt_GQ)) %>%
    filter(gt_DP >= 10) %>%
    select(-QUAL, -gt_AD, -gt_DP, -gt_GQ)
  
  saveRDS(vcf_all, file = gsub(".rds", "_filt.rds", tidy_file) %>% gsub("by_chromo", "by_chromo/filt", .))
  
}
# remove nas and apply depth filter


# filter for representation

mclapply(tidy_files, filter_tidy_rds, mc.cores = 2)