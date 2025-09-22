# invariant sites
# for calculating pi/dxy

library("vcfR")
library("tidyverse")
library("parallel")

# unzip vcfs

unzip_tidy_vcf <- function(vcf_gz){
  
  message(vcf_gz)
  
  # unzip vcf
  vcf_file <- gsub(".gz", "", vcf_gz)
  system(paste0("gunzip -c ", vcf_gz, " > ", vcf_file))
  
  # read in file
  vcf_df <- read.vcfR(vcf_file)
  
  file.remove(vcf_file)
  
  # convert to tidy vcf 
  vcf_df <- vcfR2tidy(vcf_df)
  
  # join fixed and genotype fields
  vcf_fix <- vcf_df$fix %>% 
    select(ChromKey, CHROM, POS, REF, ALT, QUAL) 
  
  vcf_df <- vcf_df$gt %>%
    select(ChromKey, POS, Indiv, gt_AD, gt_DP, gt_GQ, gt_GT) %>%
    left_join(vcf_fix, .) %>%
    select(-ChromKey)
    
  saveRDS(vcf_df, file = gsub(".vcf", "_tidy.rds", vcf_file))
          
}

# process each vcf file
vcf_files <- list.files("by_chromo", pattern = ".gz", full.names = TRUE)

vcf_files <- vcf_files[1:10]

#lapply(vcf_files, unzip_tidy_vcf)

mclapply(vcf_files, unzip_tidy_vcf, mc.cores = 2, mc.silent = FALSE, mc.preschedule = FALSE)



