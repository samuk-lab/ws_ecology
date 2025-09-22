# convert the tidy genotypes 
# into a SNP table and a GDS file (SNPRelate)
# KMS Aug 2018

library("tidyverse")
library("SNPRelate")

# create a numeric snp table

geno_df <- readRDS("data/geno_df_2020_filtered.rds")%>%
  filter(CHROM != "chrXIX" ,CHROM != "chrY")

snp_tab <- geno_df  %>%
  select(Indiv, CHROM, POS, REF, ALT, gt_GT) %>%
  extract(gt_GT, c("gt1", "gt2"), "([[:alnum:]]+)/([[:alnum:]]+)", convert = TRUE) %>%
  mutate(genotype = gt1 + gt2, chrom_pos = paste0(CHROM, "_", POS)) %>%
  filter(genotype %in% c(0,1,2)) %>%
  select(Indiv, chrom_pos, genotype) 

snp_tab <- snp_tab %>%
  spread(Indiv, genotype)

write.table(snp_tab, file = gzfile("data/snp_tables/wht_cmn_2020_snp_table.txt.gz"))

snp_tab[,-1][which(snp_tab[,-1] > 2)] <- NA

write.table(snp_tab, file = gzfile("data/snp_tables/wht_cmn_2020_snp_table_bial.txt.gz"))

##########################################
# parse (numeric) snp table to GDS format
##########################################

snp_tab <- read.table(file = gzfile("data/snp_tables/wht_cmn_2020_snp_table_bial.txt.gz"))

sample_id <- names(snp_tab[,-1])
snp_id <- snp_tab$chrom_pos %>% as.character
snp_position <- snp_tab$chrom_pos %>% gsub(".*_", "", .)
snp_chr <- snp_tab$chrom_pos %>% gsub("_.*", "", .)

snp_tab <- data.matrix(snp_tab[,-1])
row.names(snp_tab) <- snp_id

snpgdsCreateGeno("data/snp_relate/whtstbk_2020.gds", genmat = snp_tab, 
                 sample.id = sample_id, snpfirstdim = TRUE, 
                 snp.id = snp_id, snp.chromosome = snp_chr, snp.position = snp_position)

##########################################
# chrY
##########################################

geno_df <- readRDS("data/geno_df_Y_2020_filtered.rds")

snp_tab <- geno_df  %>%
  select(Indiv, CHROM, POS, REF, ALT, gt_GT) %>%
  extract(gt_GT, c("gt1", "gt2"), "([[:alnum:]]+)/([[:alnum:]]+)", convert = TRUE) %>%
  mutate(genotype = gt1 + gt2, chrom_pos = paste0(CHROM, "_", POS)) %>%
  filter(genotype %in% c(0,1,2)) %>%
  select(Indiv, chrom_pos, genotype) 

snp_tab <- snp_tab %>%
  spread(Indiv, genotype)

sample_id <- names(snp_tab[,-1])
snp_id <- snp_tab$chrom_pos %>% as.character
snp_position <- snp_tab$chrom_pos %>% gsub(".*_", "", .)
snp_chr <- snp_tab$chrom_pos %>% gsub("_.*", "", .)

snp_tab <- data.matrix(snp_tab[,-1])
row.names(snp_tab) <- snp_id

snpgdsCreateGeno("data/snp_relate/whtstbk_2020_chrY.gds", genmat = snp_tab, 
                 sample.id = sample_id, snpfirstdim = TRUE, 
                 snp.id = snp_id, snp.chromosome = snp_chr, snp.position = snp_position)

##########################################
# chrXIX
##########################################

geno_df <- readRDS("data/geno_df_XIX_2020_filtered.rds")

snp_tab <- geno_df  %>%
  select(Indiv, CHROM, POS, REF, ALT, gt_GT) %>%
  extract(gt_GT, c("gt1", "gt2"), "([[:alnum:]]+)/([[:alnum:]]+)", convert = TRUE) %>%
  mutate(genotype = gt1 + gt2, chrom_pos = paste0(CHROM, "_", POS)) %>%
  filter(genotype %in% c(0,1,2)) %>%
  select(Indiv, chrom_pos, genotype) 

snp_tab <- snp_tab %>%
  spread(Indiv, genotype)

sample_id <- names(snp_tab[,-1])
snp_id <- snp_tab$chrom_pos %>% as.character
snp_position <- snp_tab$chrom_pos %>% gsub(".*_", "", .)
snp_chr <- snp_tab$chrom_pos %>% gsub("_.*", "", .)

snp_tab <- data.matrix(snp_tab[,-1])
row.names(snp_tab) <- snp_id

snpgdsCreateGeno("data/snp_relate/whtstbk_2020_chrXIX.gds", genmat = snp_tab, 
                 sample.id = sample_id, snpfirstdim = TRUE, 
                 snp.id = snp_id, snp.chromosome = snp_chr, snp.position = snp_position)
