# process a raw vcf from GATK 
# into a tidy formatted genotype data frame
# KMS Aug 2018

library("vcfR")
library("tidyverse")

dir.create("data/vcf/2020_pipeline/tidy_invar")

# get the full y sequence
vcf_file <- list.files("data/vcf/2020_pipeline/invar", pattern = ".*chrY.*.vcf.gz$", full.names = TRUE)

# read in GATK vcf as a tidy vcf
wht_vcf <- read.vcfR(file = vcf_file, nrow = 10000) %>% vcfR2tidy(single_frame = TRUE)
wht_vcf <- wht_vcf$dat %>%
  select(-ID, -FILTER) %>%
  filter(!(is.na(ALT)))

wht_vcf_stats <- wht_vcf %>%
  select(CHROM, POS, REF, ALT, MQ, FS, DP, BaseQRankSum, QUAL) %>%
  gather(key = stat, value = value, -CHROM, -POS, -REF, -ALT) 

wht_vcf_stats <- wht_vcf_stats %>%
  group_by(POS,stat) %>%
  summarise(mean_stat = mean(value, na.rm = TRUE))

# apply site-level filters to SNPs and indels
# GATK hard filter suggestions SNPs "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5"
# GATK hard filter suggestions for INDELs "QD < 2.0 || FS > 200.0"
# (these the criteria for *exclusion*, so are inverted below)
# position-based filters excluded, as GBS/RAD results in extreme position bias

# apply the basic depth/quality filter

wht_vcf <- wht_vcf %>%
  filter(gt_DP >= 7) 

# force garbage collection
gc()

# apply the variant-class specific filtering
# (invar sites only get DP and RGQ filters)
snps <- wht_vcf %>%
  filter(nchar(ALT) == 1) %>% 
  filter(QD > 2.0, FS < 60.0, MQ > 40.0, MQRankSum > -12.5)%>%
  filter(gt_GQ >= 20)

#indels <- wht_vcf %>% 
#  filter(grepl("[A-Z\\*]{2,}",ALT)) %>% 
#  filter(QD > 2.0, FS < 200.0) %>%
# filter(gt_GQ >= 20)

#invar <- wht_vcf %>%
#  filter(is.na(ALT)) %>%
#  filter(gt_RGQ >= 20)

# rejoin data for output
# (drops all site-level stats)
wht_vcf <- snps %>%
  select(CHROM, POS, REF, ALT, Indiv, gt_AD, gt_DP, gt_GQ, gt_GT, 
         gt_MIN_DP, gt_PGT, gt_PID, gt_PL, gt_RGQ, gt_SB, gt_GT_alleles) %>%
  arrange(CHROM, POS, Indiv)

file_name <- vcf_file %>% gsub(".*\\/|.vcf.gz", "", .) %>% paste0("data/vcf/2020_pipeline/tidy_snps/", ., ".rds")

saveRDS(wht_vcf, file = file_name)

