library("tidyverse")

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

samples <- genmat$sample.id

meta_df <- read.csv("metadata/mega_meta.csv", stringsAsFactors = FALSE, header = TRUE)

meta_df <- meta_df %>%
  mutate(Indiv = gsub("whtstbk_gbs_2012_brds_", "", id)) %>%
  mutate(Indiv = gsub("_", ".", Indiv)) %>%
  mutate(Indiv = case_when(year == "2012" ~ paste0(Indiv, "_2013"),
                           year == "2014" ~ paste0(Indiv, "_2015"),
                           TRUE ~ Indiv)) %>%
  filter(Indiv %in% genmat$sample.id)

props <- read.table("analysis_faststructure/data/admixture/whtstbk_2020_nochr.2.P", h = F)
pops <- read.table("analysis_faststructure/data/admixture/whtstbk_2020_nochr.2.Q", h = F) 

names(pops) <- c("p1", "p2")

data.frame(meta_df %>% select(Indiv, pop, region, cluster), pops) %>% 
  gather(key = pop, value = prop, -Indiv, -region, -pop, -cluster)%>%
  ggplot(aes(fill=pop, y=prop, x=Indiv)) + 
  geom_bar(position="stack", stat="identity") +
  facet_wrap(~cluster, scales = "free_x")+
  theme(axis.text.x = element_text(angle = 90))
