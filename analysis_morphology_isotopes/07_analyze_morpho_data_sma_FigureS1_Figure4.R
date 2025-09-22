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

#options(device = "cairo")

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
# sma models
##############################################

trait_df <- pheno_long %>%
  filter(trait == "egg_number")

trait_df %>%
  ggplot(aes(x = body_size_scale, y = value_scale)) +
  geom_point(aes(color = species))+
  geom_smooth(method = "lm", formula = I(y ~ x))

trait_df %>%
  ggplot(aes(x = body_size_scale, y = value_scale, color = species)) +
  geom_boxplot()

sma(data = trait_df, value_scale ~ body_size_scale + species) %>% plot
lm(data = trait_df, value_scale ~ species) %>% anova

fit_model_sma <- function(trait_df){
  
  print(unique(trait_df$trait))
  
  # check which size measures are available/more complete 
  std_len_n <- trait_df %>% filter(!is.na(std_length) & !is.na(value)) %>% nrow
  csize_n <- trait_df %>% filter(!is.na(csize) & !is.na(value)) %>% nrow
  
  if(std_len_n >= csize_n){
    
    trait_df$body_size <- trait_df$std_length
    
  } else {
    
    trait_df$body_size <- trait_df$csize
    
  }
  
  trait_df <- trait_df %>%
    filter(!is.na(value)) %>%
    filter(!is.na(body_size)) %>%
    mutate(value_log = log(value + 1)) %>%
    mutate(value_scale = scale(value)) %>%
    mutate(body_size_log = log(body_size + 1)) %>%
    mutate(body_size_scale = scale(body_size)) %>%
    filter(value_scale > -4) # remove extreme outliters (4 SDs)
  
  # fit the preliminary sma model
  # also control for sex here, if applicable (not supported by sma for some reason)
  if(grepl(".*raker.*", unique(trait_df$trait)) == TRUE| grepl(".*plate.*", unique(trait_df$trait)) == TRUE){
    
    sma_mod1 <- sma(data = trait_df, value_scale ~ body_size_scale + species)
    
  } else if(!(unique(trait_df$trait) %in% c("egg_number", "egg_weight", "testis_weight", "std_length"))){
    
    trait_df <- trait_df %>% filter(!is.na(joint_sex))
    sex_resid <- lm(data = trait_df, value_scale ~ joint_sex)
    trait_df$value_resid_sex <- residuals(sex_resid)
    sma_mod1 <- sma(data = trait_df, value_resid_sex ~ body_size_scale)
    
  } else{
    
    sma_mod1 <- sma(data = trait_df, value_scale ~ body_size_scale)
    
  }
  
  # check if trait and body size are significantly associated
  if(sma_mod1$pval[[1]] <= 0.05 & sma_mod1$r2 >= 0.4 & unique(trait_df$trait) != "std_length" ){
    
    
    trait_df$value_resid <- residuals(sma_mod1)
    
    #[1] "luminance_mean"    "egg_number"        "egg_weight"        "testis_weight"     "body_depth1"      
    #[6] "spine_dorsal1"     "spine_dorsal2"     "spine_dorsal3"     "spine_pelvic"      "PC1"              
    #[11] "PC2"               "PC3"               "PC4"               "PC5"               "PC6"              
    #[16] "ventral_rakers"    "dorsal_rakers"     "raker_long_count"  "raker_short_count" "plate_count"      
    #[21] "std_length"
    
      # refit model with species id as a factor
      
      if((unique(trait_df$trait) %in% c("egg_number", "egg_weight", "testis_weight", "std_length"))|
         grepl(".*raker|plate.*", unique(trait_df$trait))){
        
        sma_mod2 <- sma(data = trait_df, value_scale ~ body_size_scale + species)
        
      } else{
        
        sma_mod2 <- sma(data = trait_df, value_resid_sex ~ body_size_scale + species)
        
      }
      
      # collect sma results into a dataframe
    
      # compute the estimated difference in elevations + approximate SE
      trait_est <- sma_mod2$coef$wht$`coef(SMA)`[1] - sma_mod2$coef$cmn$`coef(SMA)`[1]
      
      # this is the average width of the two CIs produced by sma
      # (so an approximation of the total CI, which is not computed by smatr directly)
      trait_se <- (((sma_mod2$coef$wht$`upper limit`[1] - sma_mod2$coef$wht$`lower limit`[1]) +
        (sma_mod2$coef$cmn$`upper limit`[1] - sma_mod2$coef$cmn$`lower limit`[1]))/2)/2
      
      sma_df <- data.frame(trait = unique(trait_df$trait), body_size_r2 = sma_mod1$r2[[1]], body_size_pval = sma_mod1$pval[[1]],
                           trait_test_type = "wald_chi", trait_est = trait_est, trait_se = trait_se , trait_pval = sma_mod2$gtr$p)
      
      
      # # fit model I regressions (for comparison)
      # if(unique(trait_df$trait) == c("std_length")){
      #   
      #   trait_lm <- lm(data = trait_df, value_scale ~ species) %>%
      #     tidy %>%
      #     data.frame(., trait = unique(trait_df$trait)) %>%
      #     filter(term != "(Intercept)")
      #   
      # } else if((unique(trait_df$trait) %in% c("egg_number", "egg_weight", "testis_weight", "std_length"))|
      #           grepl(".*raker|plate.*", unique(trait_df$trait))){
      #   
      #   trait_lm <- lm(data = trait_df, value_resid ~ species) %>%
      #     tidy %>%
      #     data.frame(., trait = unique(trait_df$trait)) %>%
      #     filter(term != "(Intercept)")
      #   
      # } else{
      #   
      #   trait_lm <- lm(data = trait_df, value_resid_sex ~ species) %>%
      #     tidy %>%
      #     data.frame(., trait = unique(trait_df$trait)) %>%
      #     filter(term != "(Intercept)")
      #   
      # }
      
      
      #lm_df <- data.frame(trait = trait_lm$trait, body_size_r2 = sma_mod1$r2[[1]], body_size_pval = sma_mod1$pval[[1]],
      #                    trait_test_type = "anova", trait_pval =  trait_lm$p.value)
      
    return(sma_df)
      
    } else{
      
    if(grepl(".*raker.*", unique(trait_df$trait)) == TRUE | grepl(".*plate.*", unique(trait_df$trait)) == TRUE){
      
      trait_lm <- glm(data = trait_df, (value_scale + 5) ~ species, family = "poisson") %>%
        tidy %>%
        data.frame(., trait = unique(trait_df$trait)) %>%
        filter(term != "(Intercept)")
      
      trait_pval <- trait_lm$p.value[1]
      trait_se <- trait_lm$std.error[1]
      trait_est <- trait_lm$estimate[1]
      
    
    } else if((unique(trait_df$trait) %in% c("egg_number", "egg_weight", "testis_weight", "std_length"))){
    
    trait_lm <- lm(data = trait_df, value_scale ~ species) %>%
      tidy %>%
      data.frame(., trait = unique(trait_df$trait)) %>%
      filter(term != "(Intercept)")
    
    trait_pval <- trait_lm$p.value[1]
    trait_se <- trait_lm$std.error[1]
    trait_est <- trait_lm$estimate[1]
  
    } else{
      
      trait_lm <- lm(data = trait_df, value_scale ~ joint_sex + species) %>%
        tidy %>%
        data.frame(., trait = unique(trait_df$trait)) %>%
        filter(term != "(Intercept)")
      
      trait_pval <- trait_lm$p.value[2]
      trait_se <- trait_lm$std.error[2]
      trait_est <- trait_lm$estimate[2]
      
      
    }
      
    data.frame(trait = trait_lm$trait[1], body_size_r2 = sma_mod1$r2[[1]], body_size_pval = sma_mod1$pval[[1]],
               trait_test_type = "anova", trait_est = trait_est, trait_se = trait_se, trait_pval =  trait_pval)
    
}

}



pheno_long_corr <- pheno_long %>%
  split(x = ., f = pheno_long$trait)

options(scipen = 10)

pheno_mod <- lapply(pheno_long_corr, fit_model_sma) %>%
  bind_rows %>%
  mutate(trait_pval = signif(trait_pval))

pheno_mod$p_adjust <- p.adjust(pheno_mod$trait_pval)

##############################################
# figure 4
##############################################

# Plot of all trait regression coefficients

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

fig4 <- pheno_mod %>%
  mutate(trait_recode = factor(trait, levels = trait_ord, labels = trait_recode)) %>%
  mutate(estimate_max = trait_est + (trait_se*1.96)) %>%
  mutate(estimate_min = trait_est - (trait_se*1.96)) %>%
  mutate(p_adjust = p.adjust(trait_pval, method = "BH")) %>%
  ggplot(aes(y = trait_est, x = trait_recode, 
             ymax = estimate_max, ymin = estimate_min, fill = (p_adjust < 0.1)))+
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

changeSciNot <- function(n) {
  output <- format(n, scientific = TRUE) #Transforms the number into scientific notation even if small
  output <- sub("e", "*10^", output) #Replace e with 10^
  output <- sub("\\+0?", "", output) #Remove + symbol and leading zeros on expoent, if > 1
  output <- sub("-0?", "-", output) #Leaves - symbol but removes leading zeros on expoent, if < 1
  output
}

pheno_mod %>%
  mutate(trait_recode = factor(trait, levels = trait_ord, labels = trait_recode)) %>%
  mutate(p_adjust = p.adjust(trait_pval, method = "BH")) %>%
  select(-trait) %>%
  select(trait_recode, everything()) %>% 
  mutate(body_size_pval = changeSciNot(body_size_pval)) %>%
  mutate(trait_pval = changeSciNot(trait_pval)) %>%
  mutate(p_adjust = changeSciNot(p_adjust)) %>%
  mutate(body_size_r2 = changeSciNot(body_size_r2)) %>%
  setNames(nm = c("Trait", "Body Size r2", "Body Size p-value", "Test Type", "Species Difference (Average)", "Species Difference (SE)", "Test p-value", "Test Adjusted p-value")) %>%
  xlsx::write.xlsx("figures_good/TableS1_power10.xlsx")
  

