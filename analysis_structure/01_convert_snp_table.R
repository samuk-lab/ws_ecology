# covnert snp table to fast structure

library("tidyverse")
library("radiator")

snp_tab <- read.table("data/snp_tables/whtstbk_2020_LD_0.2_MAF_0.01.txt")

snp_tab <- snp_tab %>% t %>% data.frame
ind_id <- rownames(snp_tab)
snp_id <- colnames(snp_tab)
snp_tab <- data.frame(ind_id, snp_tab)

snp_tab <- gather(snp_tab, key = locus, value = genotype, -ind_id)

# radiator format

rad_df <- snp_tab %>%
  mutate(POP_ID = gsub("[0-9]+_", "", ind_id)) %>%
  select(ind_id, POP_ID, locus, genotype) %>%
  mutate(genotype = case_when(genotype == 0 ~ "001/001",
                             genotype == 1 ~ "001/002",
                             genotype == 2 ~ "002/002"))
names(rad_df) <- c("INDIVIDUALS", "POP_ID", "MARKERS", "GT")
rad_tw <- radiator::tidy_wide(rad_df, import.metadata = TRUE)

rad_tw2 <- calibrate_alleles(rad_tw, biallelic = TRUE)
strata <- generate_strata(rad_tw, pop.id = TRUE)

tmp <- filter_rad(rad_tw, strata = strata )

genomic_converter(rad_tw2, output = "plink", strata = strata)
genomic_converter(rad_tw, output = "faststructure")

# convert from genotype to "double haplotype" format
snp_tab <- snp_tab %>%
  mutate(hap1 = case_when(genotype == 0 ~ 0,
                          genotype == 1 ~ 1,
                          genotype == 2 ~ 1)) %>%
  mutate(hap2 = case_when(genotype == 0 ~ 0,
                          genotype == 1 ~ 0,
                          genotype == 2 ~ 1))

# format back to wide format
tmp <- snp_tab %>%
  select(-genotype) %>%
  gather(key = hap, value = genotype, -ind_id, -locus) %>%
  mutate(hap_id = paste0(ind_id, "_", hap)) %>%
  select(hap_id, locus, genotype) %>%
  spread(key = locus, value = genotype)

# -9 codes for missing in .str
tmp[is.na(tmp)] <- -9

# remove haplotype codes
tmp$hap_id <- gsub("_hap.*", "", tmp$hap_id)

# preview
tmp[1:10,1:10]

# add in empty meta columns

tmp <- data.frame(meta1 = "META1", meta2 = "META2", meta3 = "META3", 
           meta4 = "META4", meta5 = "META5", tmp) 

# write to file
write.table(tmp, file = "analysis_faststructure/data/whtstbk_pruned_2020.str", col.names = FALSE, row.names = FALSE, quote = FALSE)


