######################################################################
# Libraries
######################################################################
library("vegan")
library("dplyr")
library("tidyr")
library("broom")
library("ggplot2")
library("ggthemes")
library("visreg")
library("visdat")
library("broom")
library("wesanderson")
library("smatr")

list.files("functions", full.names = TRUE) %>% sapply(source) %>% invisible
select <- dplyr::select #allows dplyr to use select when MASS package is loaded

# the white stickleback palatte
whtstbk_palatte <- c("#008FD5", "#BDBDBD")

setwd("analysis_morphology")

pheno_df <- read.table("data/pheno_df_master.txt", h = T)

pheno_sub <- pheno_df %>% 
  select(-id, -pop, -year, -region, -cluster, -sequenced, -geno_sex, -sex) %>%
  select(species, joint_sex, everything()) %>%
  filter(!is.na(species)) %>%
  mutate(species = as.factor(species))

pheno_sub <- pheno_sub %>%
  select(-matches("mode|rgb|ratio|depth2|PC7|PC8|PC9|PC10|isotope_d15N_resid|isotope_d13C_resid|egg_diameter_mean|egg_diameter_sd|testis_width_mean|testis_length_mean")) %>%
  select(-isotope_d13C,-isotope_d15N) %>%
  filter(!is.na(luminance_mean) | !is.na(spine_dorsal1) )
  #filter(complete.cases(.))

pheno_long <-pheno_sub %>%
  gather(key = trait, value = value, -species, -joint_sex, -csize, -std_length) 

# reject standard length

std_len_df <- pheno_sub %>%
  gather(key = trait, value = value, -species, -joint_sex, -csize) %>%
  filter(trait == "std_length") %>%
  mutate(std_length = NA)

pheno_long <- bind_rows(pheno_long, std_len_df)

##############################################
# Plot of raw phenotypic differences
##############################################

trait_ord <- c("std_length", "luminance_mean","body_depth1","PC1","PC2","PC3","PC4","PC5","PC6",
                   "spine_dorsal1","spine_dorsal2","spine_dorsal3","spine_pelvic", "plate_count",
                   "raker_long_count", "raker_short_count", "ventral_rakers", "dorsal_rakers",
                   "testis_weight","egg_number","egg_weight")


trait_recode <- c("Standard Length", "Body Lightness","Body Depth", "Shape PC1", "Shape PC2", "Shape PC3", "Shape PC4", "Shape PC5", "Shape PC6",
                      "1st Dorsal Spine", "2nd Dorsal Spine","3rd Dorsal Spine","Pelvic Spine","No. Lateral Plate",
                      "No. Long GR", "No. Short GR", "Long GR Length", "Short GR Length",
                      "Testis Weight","Egg Number","Egg Weight")

figureS1 <- pheno_sub %>%
  gather(key = trait, value = value, -species, -joint_sex, -csize) %>%
  filter(!is.na(species)) %>%
  filter(trait != "csize") %>%
  mutate(trait_recode = factor(trait, levels = trait_ord, labels = trait_recode)) %>%
  #mutate(trait = factor(trait, levels = pheno_reorder)) %>%
  ggplot(aes(x = species, y = value, color = species, fill = species))+
  geom_jitter(width = 0.2, size = 1)+
  stat_summary(fun.data = "mean_cl_normal", 
               geom = "pointrange", color = "black", fill = "white", 
               pch = 21)+
  #stat_summary(aes(y=value_corr), fun.data = "mean_cl_normal", 
  #             geom = "pointrange", color = "black", fill = "grey",
  #             pch = 21, size = 0.65,  position=position_nudge(x = 0.1, y = 0))+
  facet_wrap(~trait_recode, scales = "free", nrow = 6, ncol = 4)+
  theme_bw()+
  xlab("")+
  ylab("Value")+
  theme(legend.position=c(0.90, 0.05),
        legend.justification = 0.5,
        legend.text.align = 0.5,
        legend.title.align = 0.5,
        strip.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_blank())+
  scale_color_manual(values = whtstbk_palatte)+
  scale_fill_manual(values = whtstbk_palatte)

ggsave(figureS1, filename = "figures/figureS1_raw.pdf", device = "pdf", height = 7, width = 7, useDingbats=FALSE)

##############################################
# models
##############################################

fit_model <- function(trait_df){
  
  trait_df <- trait_df %>%
    mutate(value = scale(value))
  
  print(unique(trait_df$trait))
  
  if(grepl("std_length", unique(trait_df$trait))){
    
    lm(data = trait_df, value ~ joint_sex + species) %>%
      tidy%>%
      data.frame(., trait = unique(trait_df$trait))
    
  }
  
  if(grepl("test|egg", unique(trait_df$trait))){
    
    lm(data = trait_df, value ~ csize + species) %>%
      tidy%>%
      data.frame(., trait = unique(trait_df$trait))
    
  } else if(grepl("raker|plate", unique(trait_df$trait))){
    
    trait_df <- trait_df %>%
      filter(value > 0)
    
    glm(data = trait_df, value ~ std_length + species, family = "poisson") %>%
      tidy%>%
      data.frame(., trait = unique(trait_df$trait))
    
  } else if(grepl("luminance_mean", unique(trait_df$trait))){
    
    lm(data = trait_df, value ~ joint_sex + species) %>%
      tidy %>%
      data.frame(., trait = unique(trait_df$trait))
    
  } else{
    
    lm(data = trait_df, value ~ csize + joint_sex + species) %>%
      tidy %>%
      data.frame(., trait = unique(trait_df$trait))
    
    
  }
  
  
}

pheno_long_corr <- pheno_long %>%
  split(x = ., f = pheno_long$trait)

pheno_mod <- lapply(pheno_long_corr, fit_model) %>%
  bind_rows

# Plot of all trait regression coefficients
# Figure 4

trait_ord <- rev(c("luminance_mean","std_length", "body_depth1","PC1","PC2","PC3","PC4","PC5","PC6",
               "spine_dorsal1","spine_dorsal2","spine_dorsal3","spine_pelvic", "plate_count",
               "raker_long_count", "raker_short_count", "ventral_rakers", "dorsal_rakers",
               "testis_weight","egg_number","egg_weight"))


trait_recode <- rev(c("Body Lightness","Standard Length","Body Depth", "Shape PC1", "Shape PC2", "Shape PC3", "Shape PC4", "Shape PC5", "Shape PC6",
                  "1st Dorsal Spine Length", "2nd Dorsal Spine Length","3rd Dorsal Spine Length","Pelvic Spine Length","No. Lateral Plates",
                  "No. Long Gill Rakers", "No. Short Gill Rakers", "Long Gill Raker Length", "Short Gill Raker Length",
                   "Testis Weight","Egg Number","Egg Weight"))

overlay_cols <- c(wes_palette("Darjeeling1", 5, type = "d") )
wes_palettes

overlay_cols <- c("#CDCDCD", "#FFFFFF","#CDCDCD", "#FFFFFF#")
y_extents <- c(-1.75, 1.75)

fig4 <-pheno_mod %>%
  filter(term == "specieswht") %>%
  mutate(trait_recode = factor(trait, levels = trait_ord, labels = trait_recode)) %>%
  mutate(estimate_max = estimate + (std.error*2)) %>%
  mutate(estimate_min = estimate - (std.error*2)) %>%
  mutate(p_adjust = p.adjust(p.value, method = "holm")) %>%
  ggplot(aes(y= estimate, x = trait_recode, 
             ymax= estimate_max, ymin = estimate_min, fill = (p_adjust < 0.1)))+
  geom_errorbar(color = "black", width = 0)+
  geom_point(size = 3, pch = 21)+
  geom_hline(yintercept = 0, linetype = 2)+
  annotate("rect", xmin = 12.5, xmax =  21.5, ymin = y_extents[1], ymax = y_extents[2], alpha = .35, fill = overlay_cols[1] )+
  annotate("rect", xmin = 7.5, xmax =  12.5, ymin = y_extents[1], ymax = y_extents[2],alpha = .35, fill = overlay_cols[2])+
  annotate("rect", xmin = 3.5, xmax =  7.5, ymin = y_extents[1], ymax = y_extents[2],alpha = .35, fill = overlay_cols[3])+
  annotate("rect", xmin = 0.5, xmax =  3.5, ymin = y_extents[1], ymax = y_extents[2],alpha = .35, fill = overlay_cols[5])+
  geom_errorbar(color = "black", width = 0)+
  geom_point(size = 3, pch = 21)+
  geom_hline(yintercept = 0, linetype = 2)+
  coord_flip()+
  theme_bw(base_size = 14)+
  ylab("Scaled Regression Coefficient")+
  xlab("Trait")+
  theme(legend.position = "none",
        legend.justification = c(0, 1),
        #legend.position = c(0, 1),
        #legend.box = "horizontal",
        #legend.direction = "horizontal",
        legend.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA)
        #plot.margin = unit(c(0,0,0,0), "cm"),
        #panel.spacing = unit(c(0,0,0,0), "cm"))+
  )+
  scale_fill_manual(values = c("black", "red"))+
  scale_y_continuous(limits = y_extents, breaks = scales::pretty_breaks(6))

ggsave(fig4, filename = "figures/Figure4_raw.pdf", height = 6, width = 6, useDingbats = FALSE)

# table of all regression results

pheno_mod %>%
  filter(term == "specieswht") %>%
  mutate(trait_recode = factor(trait, levels = trait_ord, labels = trait_recode)) %>%
  mutate(estimate_max = estimate + (std.error*2)) %>%
  mutate(estimate_min = estimate - (std.error*2)) %>%
  mutate(p_adjust = p.adjust(p.value, method = "BH")) %>%
  select(trait_recode, estimate, std.error, statistic, p.value)

##############################################
# sma models
##############################################

fit_model <- function(trait_df){
  
  trait_df <- trait_df %>%
    mutate(value = scale(value))
  
  print(unique(trait_df$trait))
  
  if(grepl("std_length", unique(trait_df$trait))){
    
    lm(data = trait_df, value ~ joint_sex + species) %>%
      tidy%>%
      data.frame(., trait = unique(trait_df$trait))
    
  }
  
  if(grepl("test|egg", unique(trait_df$trait))){
    
    sma(data = trait_df, value ~ csize + species) %>%
      tidy %>%
      data.frame(., trait = unique(trait_df$trait))
    
  } else if(grepl("raker|plate", unique(trait_df$trait))){
    
    trait_df <- trait_df %>%
      filter(value > 0)
    
    glm(data = trait_df, value ~ std_length + species, family = "poisson") %>%
      tidy%>%
      data.frame(., trait = unique(trait_df$trait))
    
  } else if(grepl("luminance_mean", unique(trait_df$trait))){
    
    lm(data = trait_df, value ~ joint_sex + species) %>%
      tidy %>%
      data.frame(., trait = unique(trait_df$trait))
    
  } else{
    
    lm(data = trait_df, value ~ csize + joint_sex + species) %>%
      tidy %>%
      data.frame(., trait = unique(trait_df$trait))
    
    
  }
  
  
}

pheno_long_corr <- pheno_long %>%
  split(x = ., f = pheno_long$trait)

pheno_mod <- lapply(pheno_long_corr, fit_model) %>%
  bind_rows

pheno_long$trait %>% unique

# perform a standardized major axis regression

# traits are as such:
#[1] "luminance_mean"    "egg_number"        "egg_weight"        "testis_weight"     "body_depth1"      
#[6] "spine_dorsal1"     "spine_dorsal2"     "spine_dorsal3"     "spine_pelvic"      "PC1"              
#[11] "PC2"               "PC3"               "PC4"               "PC5"               "PC6"              
#[16] "ventral_rakers"    "dorsal_rakers"     "raker_long_count"  "raker_short_count" "plate_count"      
#[21] "std_length"        
trait_df <- pheno_long %>%
  filter(trait == "dorsal_rakers")

sma_mod <- sma(data= trait_df, value ~ std_length + species, log = "xy")
summary(sma_mod)

# plot the sma and the model I regression together
par(mfrow=c(1,2))
plot(sma_mod, col = c("red", "blue"), pch = 16, cex = 1.0, lwd = 3)

mod1 <- trait_df %>%
  mutate(log_std_length = log(std_length)) %>%
  lm(data = ., log(value) ~ log_std_length + species)

visreg::visreg(mod1, "log_std_length", by = "species", overlay = TRUE, 
                 points=list(cex=1.0, pch=16, col = c("red", "blue")), 
                 line=list(lwd=3.0, col = c("red", "blue")),
                 ) 

trait_df %>%
  ggplot(aes(y = value, x = std_length, color = species))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE)

lm(data= trait_df, value ~ std_length + species ) %>% car::Anova(type = 3 )

trait_df %>%
  ggplot(aes(x = std_length, y= value, color = "species"))+
  

