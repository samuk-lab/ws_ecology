install.packages("Momocs")
library("Momocs")

######################################################################
# Analysis of morphometric landmark data
# Kieran Samuk - Apr 2016
######################################################################


######################################################################
# Libraries
######################################################################

library("tidyverse")
library("car")
library("broom")
library("geomorph")
library("ggthemes")

list.files("functions", full.names = TRUE, pattern = "sherr") %>% sapply(source) %>% invisible

select<- dplyr::select #allows dplyr to use select when MASS package is loaded

# the white stickleback palatte
whtstbk_palatte <- c("#008FD5", "#FFFFFF")

setwd("analysis_morphology")

options(device="quartz")

######################################################################
# input data
######################################################################

# raw morphometrics data for 2014 (std. length, spines, landmarks)
morpho_df <- read.csv("data/corrected_landmark_data.csv")

# harmonize names


# meta data file from genotypic analysis
#meta_df <- read.csv("metadata/mega_meta.csv")

# meta data file from genotypic analysis
meta_df <- read.table("metadata/meta_data_clean.txt", h = T)

# harmonize ids
morpho_df <- morpho_df %>%
  select(-X, -sex) %>%
  mutate(population = as.character(population)) %>%
  mutate(population = ifelse(grepl("CP", population), "CP", population)) %>%
  mutate(population = ifelse(grepl("AL", population), "AL", population)) %>%
  mutate(population = ifelse(grepl("GC", population), "GC", population)) %>%
  mutate(id = paste0(population, individual, "_2014")) %>%
  select(id, ID, everything())

# join meta data to morpho data
morpho_df <- left_join(morpho_df, meta_df)

######################################################################
# transform data + generalized procrustes
######################################################################

# Create a 3D array with just the x/y coordinates 

outlier_individuals <- c("GC15_2014", "MH9_2014")

morpho_sub <- morpho_df %>%
  filter(!is.na(x.1)) %>%
  filter(!is.na(cluster)) %>%
  filter(!grepl("AL-S", ID)) %>%
  filter(!(id %in% outlier_individuals)) %>%
  filter(!grepl("GC", id))  %>% 
  distinct %>%
  select(-individual, -species,  -pop_cluster, -pop, -ID) %>%
  gather(key = landmark, value = position, -id, -population, -cluster, -sex) %>%
  mutate(land_axis = ifelse(grepl("x", landmark), "x", "y")) %>%
  mutate(landmark = gsub("[xyXY].", "", landmark)) %>% 
  filter(!(landmark %in% c(17:19, 8, 10))) %>%
  spread(key = land_axis, value = position) %>%
  arrange(id, as.numeric(landmark))

morpho_list<-list()
morpho_list$coo <- split(morpho_sub %>% select(x,y), morpho_sub$id) 
morpho_list$coo <- lapply(morpho_list$coo, function(x) matrix(c(as.numeric(x$x), as.numeric(x$y)), nrow = 14))

fac <- morpho_sub %>% select(id, cluster, population, sex) %>% distinct

morpho_coo <- Ldk(morpho_list$coo, fac = fac)

coo_scale(morpho_coo)

morph_pg <- morpho_coo %>%
  fgProcrustes

morpho_pg_scale <- morph_pg %>%
  coo_alignxax() %>%
  coo_scale()

coo_slide(morpho_pg_scale, id = 12) %>%
  pile() %>%
  efourier(6, norm=FALSE) %>%
  plot

pca_out <- morpho_pg_scale %>%
  PCA

plot(pca_out, ~cluster, xax = 1, yax = 2, cex = 1, pch = 21,
     nc.shp = 3, nr.shp = 3, pos.shp = "full", size.shp = 2,
     border.shp = "black", col.shp = "black", rug = FALSE, 
     center.origin = TRUE)
  
lda_out <- LDA(pca_out, fac = "cluster")
plot(lda_out)
plot_table(lda_out, fac1 = "cluster")
classification_metrics(lda_out)

KMEANS(lda_out , centers = 3)

manova_test <- MANOVA_PW(pca_out, ~cluster)
manova_test$manovas

MSHAPES(morpho_pg_scale, fac = "cluster") 

# select landmarks
# ignoring eye/pelvis (17-19) for now
landmark_df <- morpho_sub %>%
  select(matches("^[x,y]{1}\\.")) %>%
  select(-x.19, -y.19, -x.18, -y.18, -x.17, -y.17, -x.10, -y.10, -x.8, -y.8)
