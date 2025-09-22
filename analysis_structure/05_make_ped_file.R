library("SNPRelate")
library("tidyverse")
library("plotly")
library("GGally")

meta_df <- read.csv("metadata/mega_meta.csv", stringsAsFactors = FALSE, header = TRUE)

meta_df <- meta_df %>%
  mutate(Indiv = gsub("whtstbk_gbs_2012_brds_", "", id)) %>%
  mutate(Indiv = gsub("_", ".", Indiv)) %>%
  mutate(Indiv = case_when(year == "2012" ~ paste0(Indiv, "_2013"),
                           year == "2014" ~ paste0(Indiv, "_2015"),
                           TRUE ~ Indiv)) %>%
  mutate(Indiv = ifelse(Indiv == "CL64.2_2015", "CL64-2_2015", Indiv))


genofile <- snpgdsOpen("data/snp_relate/whtstbk_2020.gds", allow.duplicate = TRUE)

# ld pruning 
# also filters for MAF and missingness
# exclude XIX and Y

snpset <- snpgdsLDpruning(genofile, autosome.only = FALSE, maf = 0.01, 
                          missing.rate = 0.50, method = "corr", slide.max.bp = 75000,
                          ld.threshold = 0.20)
# extract the snp set
snpset_vec <- unlist(snpset)
names(snpset_vec) <- NULL

genmat <- snpgdsGetGeno(genofile, snp.id = snpset_vec, with.id = TRUE)
snpgdsClose(genofile)

snpgdsCreateGeno("analysis_faststructure/data/whtstbk_2020_nochr.gds", genmat$genotype, snp.chromosome=rep("1", 8116), sample.id = genmat$sample_id, snp.id = genmat$snp.id)


#genofile <- snpgdsOpen("data/snp_relate/whtstbk_2020.gds", allow.duplicate = TRUE)
genofile <- snpgdsOpen("analysis_faststructure/data/whtstbk_2020_nochr.gds", allow.duplicate = TRUE)

snpgdsGDS2BED(genofile, bed.fn = "analysis_faststructure/data/output_plink/whtstbk_2020_nochr")
