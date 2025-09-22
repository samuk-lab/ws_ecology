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

setwd("analysis_morphology")
list.files("functions", full.names = TRUE, pattern = "sherr") %>% sapply(source) %>% invisible

select<- dplyr::select #allows dplyr to use select when MASS package is loaded

# the white stickleback palatte
whtstbk_palatte <- c("#008FD5", "#FFFFFF")

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
  select(-X, -ID, -sex) %>%
  mutate(population = as.character(population)) %>%
  mutate(population = ifelse(grepl("CP", population), "CP", population)) %>%
  mutate(population = ifelse(grepl("AL", population), "AL", population)) %>%
  mutate(population = ifelse(grepl("GC", population), "GC", population)) %>%
  mutate(id = paste0(population, individual, "_2014")) %>%
  select(id, everything())

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
  filter(!(id %in% outlier_individuals)) %>%
  filter(!grepl("GC", id))
  
# select landmarks
# ignoring eye/pelvis (17-19) for now
landmark_df <- morpho_sub %>%
  select(matches("^[x,y]{1}\\.")) %>%
  select(-x.19, -y.19, -x.18, -y.18, -x.17, -y.17, -x.10, -y.10, -x.8, -y.8)

# create landmark array (14 landmarks in 2 columns each)
landmark_array <- arrayspecs(landmark_df, 14, 2)

# generalized procrustes analysis (scale/align landmarks)
landmark_gpa <- gpagen(landmark_array) 

# confirm gpa is not insane
plot(landmark_gpa)

#landmark_gpa_2d <- two.d.array(landmark_gpa$coords) #2D Data frame of procrustes coordinates

# add in scaling factor into morpho_sub
a <- landmark_gpa$Csize

morpho_sub <- data.frame(morpho_sub, csize = landmark_gpa$Csize)

######################################################################
# pca of landmark data
######################################################################

plot_groups <- ifelse(grepl("cbr", morpho_sub$pop_cluster), "cbr", as.character(morpho_sub$cluster))
plot_groups <- factor(plot_groups)
morpho_pca <- plotTangentSpace(landmark_gpa$coords, axis1 = 1, axis2 = 2, groups = factor(morpho_sub$cluster), verbose = TRUE, label = morpho_sub$id)

plotTangentSpace(landmark_gpa$coords, axis1 = 1, axis2 = 2,  groups = factor(morpho_sub$cluster), verbose = TRUE)
plotTangentSpace(landmark_gpa$coords, axis1 = 3, axis2 = 4,  groups = factor(morpho_sub$cluster), verbose = TRUE)
plotTangentSpace(landmark_gpa$coords, axis1 = 4, axis2 = 5,  groups = factor(morpho_sub$cluster), verbose = TRUE)

species <- ifelse(morpho_sub$cluster == "cmn", "common", "white")

morpho_sub$cluster <- plot_groups

sherrat_pca_plot(pc1 = 1, pc2 = 2, landmark_gpa = landmark_gpa, file_name = "figures/Figure5A.pdf",
                 group_factor = plot_groups, label = morpho_sub$id, correct_allometry = TRUE)


sherrat_pca_plot(pc1 = 3, pc2 = 4, landmark_gpa = landmark_gpa, file_name = "figures/Figure5B.pdf",
                 groups =  plot_groups, label = morpho_sub$id, correct_allometry = TRUE)

# screen plot of pcs
pca_scree <- morpho_pca$pc.summary$importance[2,]
names(pca_scree) <- 1:length(pca_scree)
barplot(pca_scree, cex.lab = 0.1)

MASS::parcoord(morpho_pca$pc.scores[,1:16], col=as.numeric(as.factor(morpho_sub$cluster)))

######################################################################
# comparison of models
######################################################################

# create a data frame for procD operations
gdf <- geomorph.data.frame(landmark_gpa, cluster = morpho_sub$cluster, sex = morpho_sub$sex) # geomorph data frame

# simple test of allometry
whtstbk_allometry1 <- procD.lm(coords~log(Csize), data=gdf) 
whtstbk_allometry2 <- procD.lm(coords~log(Csize)+sex, data=gdf)
whtstbk_allometry3 <- procD.lm(coords~log(Csize)+sex+cluster, data=gdf)
whtstbk_allometry4 <- procD.lm(coords~log(Csize)+sex*cluster, data=gdf)

anova(whtstbk_allometry1, whtstbk_allometry2, whtstbk_allometry3, whtstbk_allometry4)
anova(whtstbk_allometry4)
######################################################################
# allometry
######################################################################

# create a data frame for procD operations
gdf <- geomorph.data.frame(landmark_gpa, cluster = morpho_sub$cluster, sex = morpho_sub$sex) # geomorph data frame

# simple test of allometry
whtstbk_allometry1 <- procD.lm(coords~log(Csize), data=gdf) 

whtstbk_allometry2 <- procD.lm(coords~(Csize)+sex*cluster, data=gdf)

# anova table
summary(whtstbk_allometry1)
summary(whtstbk_allometry2)

plotAllometry(whtstbk_allometry2, size = gdf$Csize, logsz = TRUE, method = "RegScore")

PLS <- two.b.pls(log(gdf$Csize), gdf$coords, print.progress = FALSE)
PLS
plot(PLS)

plotAllometry(whtstbk_allometry2, size = gdf$Csize, logsz = TRUE, method = "CAC")

fit.unique <- procD.lm(coords ~ log(Csize) * sex/cluster, 
                       data = gdf, print.progress = FALSE) # unique allometries
fit.common <- procD.lm(coords ~ log(Csize) + sex/cluster, 
                       data = gdf, print.progress = FALSE) 

anova(fit.common, fit.unique, print.progress = FALSE)

plotAllometry(fit.unique, size = gdf$Csize, logsz = TRUE, method = "RegScore",
              pch = 19, col = as.factor(gdf$cluster))

# a more sane plot
pca_df <- data.frame(id = morpho_sub$id, sex = morpho_sub$sex, csize = landmark_gpa$Csize,
                     cluster = morpho_sub$cluster, morpho_pca$pc.scores[,1:10])

pca_df_long <- gather(pca_df, key = pc, value = score, -id, -sex, -csize, -cluster)

# how do all the pcs scale with body size?
pca_df_long %>%
  ggplot(aes(y = score, x = csize, color = cluster))+
  geom_point()+
  geom_smooth(method = "lm", color = "black", fill = "black", size = 1, se = FALSE)+
  facet_wrap(~pc)+
  theme_base()+
  scale_color_manual(values = c("#77AB43", "#008FD5", "#BDBDBD"))

# how do all the pcs compare between species?
pca_df_long %>%
  ggplot(aes(y = score, x = cluster))+
  geom_jitter(color = "grey", width = 0.5)+
  stat_summary(fun.y=median, fun.ymin=median, fun.ymax=median, 
               geom="crossbar", width = 0.9, color = "black", fatten = 2)+
  facet_wrap(~pc)+
  theme_base()


# create residual pcs (regress on csize)

resid <- pca_df_long %>%
  group_by(cluster, pc) %>%
  do(augment(lm(score ~ csize, data=.))) %>%
  ungroup %>%
  select(.resid) %>%
  unlist %>% as.numeric

pca_df_long$resid_score <- resid

# plot residuals vs. size (sanity check)
pca_df_long %>%
  ggplot(aes(y = resid_score, x = csize, color = cluster))+
  geom_point()+
  geom_smooth(method = "lm", color = "black", fill = "black", size = 1, se = FALSE)+
  facet_wrap(~pc)+
  theme_base()+
  scale_color_manual(values = c("#77AB43", "#008FD5", "#BDBDBD"))

# how do the residuals pcs compare between species?
pca_df_long %>%
  #mutate(pc = factor(pc, levels = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))) %>%
  mutate(species = ifelse(!is.na(cluster), as.character(cluster), as.character(species))) %>%
  mutate(species = ifelse(species == "white", "white", "common")) %>%
  ggplot(aes(y = resid_score, x = cluster, color = cluster))+
  geom_jitter(width = 0.5, size = 2)+
  stat_summary(fun.y=median, fun.ymin=median, fun.ymax=median, 
               geom="crossbar", width = 0.9, color = "black", fatten = 4)+
  facet_wrap(~pc)+
  theme_base()+
  scale_color_manual(values = c("#77AB43", "#008FD5", "#BDBDBD"))

# write collated output
pca_df_long %>% 
  select(-score) %>% 
  spread(key = pc, value = resid_score) %>%
  select(-sex, -cluster) %>%
  write.table(., file = "data/collated/raw_pc_scores.txt", quote = FALSE, row.names = FALSE)

#pca_df_long %>% 
#  select(-resid_score) %>% 
#  spread(key = pc, value = score) %>%
#  select(-sex, -cluster) %>%
#  write.table(., file = "data/collated/raw_pc_scores.txt", quote = FALSE, row.names = FALSE)

######################################################################
# Anova: do pcs differ between groups?
######################################################################

pca_df <- pca_df_long %>%
  mutate(species = ifelse(!is.na(cluster), as.character(cluster), as.character(species))) %>%
  mutate(species = ifelse(species == "white", "white", "common")) %>%
  select(-score) %>%
  spread(key = pc, value = resid_score)


pca_df %>% 
  dplyr::select(-id, -sex, -csize, -cluster) %>%
  select(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10) %>%
  GGally::ggpairs(columns = 1:10)

lm(PC2 ~ sex * cluster, data = pca_df) %>% anova
lm(PC3 ~ sex * cluster, data = pca_df) %>% anova
lm(PC4 ~ sex * cluster, data = pca_df) %>% anova
lm(PC5 ~ sex * cluster, data = pca_df) %>% anova
lm(PC6 ~ sex * cluster, data = pca_df) %>% anova

mano <- manova(cbind(PC2, PC3, PC4, PC5, PC6, PC7) ~ sex * cluster, data = pca_df)
summary(mano)
Anova(mano)

# nope

######################################################################
# Anova: do pcs differ between groups?
######################################################################

landmark_gpa_2d <- two.d.array(landmark_gpa$coords) #2D Data frame of procrustes coordinates

# create a data frame for procD operations
gdf <- geomorph.data.frame(landmark_gpa, species = as.factor(morpho_sub$cluster), sex = morpho_sub$sex, pop = morpho_sub$pop_cluster) # geomorph data frame
gdf$Csize <- gdf$Csize %>% as.numeric

manova_fit <- procD.lm(coords ~ log(Csize)*sex*species, iter=499, data = gdf)
anova(manova_fit)

gdf_anova <- procD.lm(landmark_gpa_2d~log(Csize), ~log(Csize) + species, groups = ~species, iter = 999, data = gdf) 
summary(gdf_anova)

pfit1 <- procD.lm(coords~log(Csize)+pop, data = gdf)

######################################################################
# LDA
######################################################################
library("MASS")

landmark_gpa_2d <- two.d.array(landmark_gpa$coords) #2D Data frame of procrustes coordinates

# create a data frame for procD operations
gdf <- geomorph.data.frame(landmark_gpa, species = as.factor(morpho_sub$cluster), sex = morpho_sub$sex, pop = morpho_sub$pop_cluster) # geomorph data frame
gdf$Csize <- gdf$Csize %>% as.numeric

lda_df <- data.frame(species = gdf$species, sex = gdf$sex, landmark_gpa_2d)

lda.species <- lda(gdf$species ~ landmark_gpa_2d)

land_nam <- colnames(lda.species$means)

sp_mean_lda <- data.frame(species = rownames(lda.species$means), data.frame(lda.species$means)) %>%
  gather(key = landmark, value = position, -species) %>%
  mutate(land_axis = ifelse(grepl("X", landmark), "x", "y")) %>%
  distinct %>%
  mutate(landmark = gsub("landmark_gpa_2d", "", landmark) %>% gsub(".[XY]", "", .)) %>%
  spread(key = land_axis, value = species_mean) %>%
  arrange(species, as.numeric(landmark))

sp_mean_lda %>%
  ggplot(aes(x = x, y = y, color = species))+
  geom_point()+
  theme_bw()

#Linear discriminant function analysis 

lda.species <- lda(GPA.2D, landmarks$species) #what does this mean? 

# the default plot, looks meh
plot(lda.species)

# project the original values in lda space
plda <- predict(object = lda.species,
                newdata = lda_df[,-c(1:2)])

# the percent of variance explained by the LD funcitons
prop.lda <- lda.species$svd^2/sum(lda.species$svd^2) 
prop.lda <- round(prop.lda*100)

# data frame of projected data with species names
lda.project <- data.frame(sex = gdf$population, 
                          species = gdf$species, 
                          ID = landmarks$ID,
                          ld1 = plda$x[,1], 
                          ld2 = plda$x[,2])
########
##use this data frame to look at lda values versus original x/y

lda.xy <- merge(lda.project, landmarks, by=c("population", "individual", "species", "ID", "X"))
lda.xy <- subset(lda.xy, select=-c(X)) 

#create a linear regression coefficient for each x:ld1 and y:ld2
lda.lm<-data.frame(matrix(NA, nrow = nrow(lda.xy) , ncol = length(PC.scores)))
colnames(lda.lm)<-c("LMx.1","LMx.2","LMx.3","LMx.4","LMx.5","LMx.6","LMx.7","LMx.8",
                    "LMx.9","LMx.10","LMx.11","LMx.12","LMx.13","LMx.14","LMx.15",
                    "LMx.16","LMx.17","LMx.18","LMx.19",
                    "LMy.1","LMy.2","LMy.3","LMy.4","LMy.5","LMy.6","LMy.7","LMy.8",
                    "LMy.9","LMy.10","LMy.11","LMy.12","LMy.13","LMy.14","LMy.15",
                    "LMy.16","LMy.17","LMy.18","LMy.19")
ld1 <- lda.xy$ld1
ld2 <- lda.xy$ld2

#separate loops for X and Y cause my looping skillz are weak

for (l in 1:19){
  pie <- as.matrix(lda.xy[paste("x", l, sep=".")])
  for(i in 1:length(pie)){ 
    lm.1<-lm(ld1 ~ pie)
    lda.lm[,l] <- lm.1$fitted.values
  }
}

for (l in 1:19){
  cake <- as.matrix(lda.xy[paste("y", l, sep=".")])
  for(i in 1:length(cake)){ 
    lm.2<-lm(ld2 ~ cake)
    lda.lm[,(l+19)] <- lm.2$fitted.values
  }
}


#######

# plot with ggplot
lda.project %>%
  ggplot(aes(color = species, x = ld1,y = ld2))+
  geom_point(size = 3) +
  labs(x = paste0("LD1 (", prop.lda[1], "%)"),
       y = paste0("LD2 (", prop.lda[2], "%)"))

# no "Both" locations
# ooo pretty
lda.project %>%
  filter(species != "B") %>%
  ggplot(aes(color = species, x = ld1,y = ld2))+
  geom_point(size = 3) +
  labs(x = paste0("LD1 (", prop.lda[1], "%)"),
       y = paste0("LD2 (", prop.lda[2], "%)")) +
  theme_classic() +
  guides(color=guide_legend(title="Species")) +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))+
  theme(axis.text=element_text(size=20), axis.title=element_text(size=22), 
        legend.text=element_text(size=20), legend.title=element_text(size=20))

