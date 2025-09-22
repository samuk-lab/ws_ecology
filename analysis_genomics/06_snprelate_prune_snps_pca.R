# prune snps based on ld and write to file
# also, perform a PCA on pruned data
# KMD Aug 2020

# source("https://bioconductor.org/biocLite.R")
# biocLite("SNPRelate")

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


#################################
# prepare gds for PCA (e.g. prune, maf filter, etc.)
#################################

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

# no ld pruning
#snpset_vec <- snpgdsSelectSNP(genofile, autosome.only = FALSE, missing.rate = 0.5, maf = 0.01)

# invesgiate chr IX ld block
#snpset_ix <- snpset_vec[grepl("chrIX_", snpset_vec)]

#snpset_vec <- snpgdsSelectSNP(genofile, missing.rate = 0.25, maf = 0.01, autosome.only = "chrY")
pca <- snpgdsPCA(genofile, snp.id = snpset_vec, num.thread = 3, eigen.cnt = 16, autosome.only = FALSE, bayesian = TRUE)

# create PCA data frame
tab <- data.frame(Indiv = pca$sample.id,
                  PC1 = pca$eigenvect[,1],  
                  PC2 = pca$eigenvect[,2],  
                  PC3 = pca$eigenvect[,3],  
                  PC4 = pca$eigenvect[,4],  
                  PC5 = pca$eigenvect[,5],  
                  PC6 = pca$eigenvect[,6],  
                  PC7 = pca$eigenvect[,7],  
                  PC8 = pca$eigenvect[,8],  
                  PC9 = pca$eigenvect[,9],  
                  PC10 = pca$eigenvect[,10],
                  PC11 = pca$eigenvect[,11],
                  stringsAsFactors = FALSE)



#convert ids to match metadata
#tab$id <- tab$id %>% gsub("whtstbk_gbs_|brds_", "", .)
pca_df <- left_join(tab, meta_df) %>%
  filter(!is.na(PC1))

#pca_df$id <- pca_df$id %>% gsub("whtstbk_gbs_|brds_", "", .)
#pca_df$id <- pca_df$id %>% gsub("NG-5241_[0-9]*_STD_", "", .)

whtstbk_palatte <- c("#008FD5", "#FFFFFF")

(pca_plot <- pca_df %>% 
    filter(!is.na(pop)) %>%
    filter(!is.na(year)) %>%
    filter(!is.na(sex)) %>%
    ggplot(aes(x = PC1, y = PC2, color = cluster, label = Indiv, shape = sex)) +
    xlab("PC1") +
    ylab("PC2")+
    geom_point(size = 3)+
    theme_bw()+
    scale_color_brewer(type = "qual", palette = "Set1"))

p <- ggplotly(pca_plot, dynamicTicks = TRUE)
p

# inspect pairs of EVs 
ggpairs(pca_df %>% filter(!is.na(sex)), columns = 2:6, aes(color = cluster, shape = sex))

# inspect pairs of EVs for explanatory power
ggpairs(pca_df, columns = 2:6, aes(color = region, shape = sex))

# looks good
# PC 1 = white/non white
# PC 2 = ~sex, but also geo (bras d'or non-bras d'or)
p <- ggplotly(pca_plot, dynamicTicks = TRUE)
p

var_exp <- data.frame(pc_num = 1:15, var_exp = pca$varprop[1:15])

# output PCA data
write.table(pca_df, "data/pca/pca_df.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(var_exp, "data/pca/pca_var_df.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)

# output pruned snp table
out_snps_tmp <- snpgdsGetGeno(genofile, with.id = TRUE, snp.id = snpset_vec)
out_snps <- out_snps_tmp$genotype

out_snps <- out_snps_tmp$genotype
row.names(out_snps) <- out_snps_tmp$snp.id
colnames(out_snps) <- out_snps_tmp$sample.id

write.table(out_snps, file = "data/snp_tables/whtstbk_2020_LD_0.2_MAF_0.01.txt")

####################################
# exploratory PCA: chr XIX 
####################################

genofile_XIX <- snpgdsOpen("data/snp_relate/whtstbk_2020_chrXIX.gds", allow.duplicate = TRUE)

snpset <- snpgdsLDpruning(genofile_XIX, autosome.only = FALSE, maf = 0.01, 
                          missing.rate = 0.50, method = "corr", slide.max.bp = 75000,
                          ld.threshold = 0.20)
# extract the snp set
snpset_vec <- unlist(snpset)
names(snpset_vec) <- NULL

# do the actual PCA
pca <- snpgdsPCA(genofile_XIX, snp.id = snpset_vec, num.thread = 3, eigen.cnt = 16, autosome.only = FALSE, bayesian = TRUE)


# create PCA data frame
tab <- data.frame(Indiv = pca$sample.id,
                  PC1 = pca$eigenvect[,1],  
                  PC2 = pca$eigenvect[,2],  
                  PC3 = pca$eigenvect[,3],  
                  PC4 = pca$eigenvect[,4],  
                  PC5 = pca$eigenvect[,5],  
                  PC6 = pca$eigenvect[,6],  
                  PC7 = pca$eigenvect[,7],  
                  PC8 = pca$eigenvect[,8],  
                  PC9 = pca$eigenvect[,9],  
                  PC10 = pca$eigenvect[,10],
                  PC11 = pca$eigenvect[,11],
                  stringsAsFactors = FALSE)

#convert ids to match metadata
#tab$id <- tab$id %>% gsub("whtstbk_gbs_|brds_", "", .)
pca_df <- left_join(tab, meta_df) %>%
  filter(!is.na(PC1))

#pca_df$id <- pca_df$id %>% gsub("whtstbk_gbs_|brds_", "", .)
#pca_df$id <- pca_df$id %>% gsub("NG-5241_[0-9]*_STD_", "", .)

whtstbk_palatte <- c("#008FD5", "#FFFFFF")

(pca_plot <- pca_df %>% 
  filter(!is.na(pop)) %>%
  filter(!is.na(year)) %>%
    filter(!is.na(sex)) %>%
  ggplot(aes(x = PC1, y = PC2, color = cluster, label = Indiv, shape = sex)) +
  xlab("PC1") +
  ylab("PC2")+
  geom_point(size = 3)+
  theme_bw()+
  scale_color_brewer(type = "qual", palette = "Set1"))

p <- ggplotly(pca_plot, dynamicTicks = TRUE)
p

####################################
# exploratory PCA: chr Y 
####################################

genofile_Y <- snpgdsOpen("data/snp_relate/whtstbk_2020_chrY.gds", allow.duplicate = TRUE)

snpset <- snpgdsLDpruning(genofile_Y, autosome.only = FALSE, maf = 0.01, 
                          missing.rate = 0.50, method = "corr", slide.max.bp = 75000,
                          ld.threshold = 0.0)
# extract the snp set
snpset_vec <- unlist(snpset)
names(snpset_vec) <- NULL

# do the actual PCA
pca <- snpgdsPCA(genofile_Y, snp.id = snpset_vec, num.thread = 3, eigen.cnt = 16, autosome.only = FALSE, bayesian = TRUE)


# create PCA data frame
tab <- data.frame(Indiv = pca$sample.id,
                  PC1 = pca$eigenvect[,1],  
                  PC2 = pca$eigenvect[,2],  
                  PC3 = pca$eigenvect[,3],  
                  PC4 = pca$eigenvect[,4],  
                  PC5 = pca$eigenvect[,5],  
                  PC6 = pca$eigenvect[,6],  
                  PC7 = pca$eigenvect[,7],  
                  PC8 = pca$eigenvect[,8],  
                  PC9 = pca$eigenvect[,9],  
                  PC10 = pca$eigenvect[,10],
                  PC11 = pca$eigenvect[,11],
                  stringsAsFactors = FALSE)

#convert ids to match metadata
#tab$id <- tab$id %>% gsub("whtstbk_gbs_|brds_", "", .)
pca_df <- left_join(tab, meta_df) %>%
  filter(!is.na(PC1))

#pca_df$id <- pca_df$id %>% gsub("whtstbk_gbs_|brds_", "", .)
#pca_df$id <- pca_df$id %>% gsub("NG-5241_[0-9]*_STD_", "", .)

whtstbk_palatte <- c("#008FD5", "#FFFFFF")

(pca_plot <- pca_df %>% 
    filter(!is.na(pop)) %>%
    filter(!is.na(year)) %>%
    filter(!is.na(sex)) %>%
    ggplot(aes(x = PC3, y = PC2, color = cluster, label = Indiv, shape = sex)) +
    #xlab("PC1") +
    #ylab("PC2")+
    geom_point(size = 3)+
    theme_bw()+
    scale_color_brewer(type = "qual", palette = "Set1"))

p <- ggplotly(pca_plot, dynamicTicks = TRUE)
p


########################################
# NO LD filter (for joining to big set + OUTFLANK)
########################################

genofile2 <- snpgdsOpen("data/snp_relate/whtstbk_2020.gds", allow.duplicate = TRUE)

# ld pruning 
snpset <- snpgdsSelectSNP(genofile2, autosome.only = FALSE, maf = 0.01, missing.rate = 0.50)

snpset_vec <- unlist(snpset)
names(snpset_vec) <- NULL

out_snps_tmp <- snpgdsGetGeno(genofile2, with.id = TRUE, snp.id = snpset_vec, sample.id = pca_samples)
out_snps <- out_snps_tmp$genotype

out_snps <- out_snps_tmp$genotype
row.names(out_snps) <- out_snps_tmp$snp.id
colnames(out_snps) <- out_snps_tmp$sample.id

write.table(out_snps, file = "data/snp_tables/whtstbk_2020_MAF_0.01.txt")

# chrXIX: pruned

snpset_vec_xix_pruned <- snpgdsLDpruning(genofile_XIX , autosome.only = FALSE, maf = 0.01, 
                          missing.rate = 0.70, method = "corr", slide.max.bp = 75000,
                          ld.threshold = 0.20)
snpset_vec_xix_pruned <- unlist(snpset_vec_xix_pruned)
names(snpset_vec_xix_pruned) <- NULL
snpset_vec_xix_pruned <- snpset_vec_xix_pruned[grepl("chrXIX", snpset_vec_xix_pruned)]
out_snps_tmp <- snpgdsGetGeno(genofile_XIX, with.id = TRUE, snp.id = snpset_vec_xix_pruned)
out_snps <- out_snps_tmp$genotype

out_snps <- out_snps_tmp$genotype
row.names(out_snps) <- out_snps_tmp$snp.id
colnames(out_snps) <- out_snps_tmp$sample.id
write.table(out_snps, file = "data/snp_tables/whtstbk_2020_MAF_0.01_XIX_pruned.txt")


snpset_vec_xix_unpruned <- snpgdsSelectSNP(genofile_XIX, autosome.only = FALSE, missing.rate = 0.7, maf = 0.01)
snpset_vec_xix_unpruned <- unlist(snpset_vec_xix_unpruned)
names(snpset_vec_xix_unpruned) <- NULL
snpset_vec_xix_unpruned <- snpset_vec_xix_unpruned[grepl("chrXIX", snpset_vec_xix_unpruned)]

out_snps_tmp <- snpgdsGetGeno(genofile_XIX, with.id = TRUE, snp.id = snpset_vec_xix_unpruned)
out_snps <- out_snps_tmp$genotype

out_snps <- out_snps_tmp$genotype
row.names(out_snps) <- out_snps_tmp$snp.id
colnames(out_snps) <- out_snps_tmp$sample.id
write.table(out_snps, file = "data/snp_tables/whtstbk_2020_MAF_0.01_XIX_unpruned.txt")

# chrY: no pruning
genofile_Y <- snpgdsOpen("data/snp_relate/whtstbk_2020_chrY.gds", allow.duplicate = TRUE)

snpset_vec_y_unpruned <- snpgdsSelectSNP(genofile_Y, autosome.only = FALSE)
snpset_vec_y_unpruned <- unlist(snpset_vec_y_unpruned)
names(snpset_vec_y_unpruned) <- NULL
snpset_vec_y_unpruned <- snpset_vec_y_unpruned[grepl("chrY", snpset_vec_y_unpruned)]

out_snps_tmp <- snpgdsGetGeno(genofile_Y, with.id = TRUE, snp.id = snpset_vec_y_unpruned)
out_snps <- out_snps_tmp$genotype

out_snps <- out_snps_tmp$genotype
row.names(out_snps) <- out_snps_tmp$snp.id
colnames(out_snps) <- out_snps_tmp$sample.id
write.table(out_snps, file = "data/snp_tables/whtstbk_2020_MAF_0.01_Y_unpruned.txt")

sample.id <- read.gdsn(index.gdsn(genofile_Y, "sample.id"))
pops <- data.frame(Indiv = sample.id) %>% 
  left_join(meta_df) %>% 
  filter(cluster %in% c("wht","cmn")) %>%
  select(Indiv, cluster) %>%
  mutate(cluster = droplevels(cluster))
 
fst_dat <- snpgdsFst(genofile_Y, method = c("W&C84"), population = pops$cluster, 
                     sample.id = pops$Indiv, 
                     autosome.only = FALSE,
                     snp.id = snpset_vec_y_unpruned, with.id = TRUE)

fst_y_dat <- data.frame(snp_id = fst_dat$snp.id, FST = fst_dat$FstSNP) %>%
  separate(snp_id, into = c("chrom", "pos")) %>%
  mutate(OutlierFlag = FALSE)

write_rds(fst_y_dat ,"data/stats/fst_snp_Y.rds")
  

