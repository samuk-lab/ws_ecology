######################################################################
# Analysis of isotopic abundance data
# Kieran Samuk - Apr 2016
######################################################################


######################################################################
# Libraries
######################################################################

library("ggplot2")
library("dplyr")
library("tidyr")
library("ggplot2")
library("ggthemes")
library("cowplot")

#list.files("functions", full.names = TRUE) %>% sapply(source) %>% invisible

select <- dplyr::select #allows dplyr to use select when MASS package is loaded

# the white stickleback palatte
whtstbk_palatte <- c("#008FD5", "#BDBDBD")

setwd("analysis_morphology")

######################################################################
# input data
######################################################################

# raw morphometrics data for 2014 (std. length, spines, landmarks)
iso_df <- read.csv("data/isotope_data.csv")

# meta data file from genotypic analysis
meta_df <- read.csv("metadata/mega_meta.csv")

# harmonize ids
meta_df <- meta_df %>%
  select(id, cluster, sex) %>%
  rename(geno_sex = sex)

iso_df  <- left_join(iso_df, meta_df, by = "id")

iso_df <- iso_df %>%
  mutate(group = ifelse(cluster == "wht", "white", "common")) %>%
  mutate(pop = gsub("[^A-Z]*", "", as.character(id))) %>%
  mutate(cn.ratio = C.amount / N.amount)


iso_df$d15N.resid <- lm(data = iso_df, d15N + cn.ratio ~ geno_sex, na.action = "na.exclude") %>% residuals
iso_df$d13C.resid <- lm(data = iso_df, d13C + cn.ratio ~ geno_sex, na.action = "na.exclude") %>% residuals

iso_df %>%
  select(id, d13C, d15N, d13C.resid, d15N.resid, cn.ratio) %>%
  rename(isotope_d13C = d13C, isotope_d15N = d15N, 
         isotope_d13C_resid = d13C.resid, isotope_d15N_resid = d15N.resid, isotope_CN_ratio = cn.ratio) %>%
  write.table(file = "data/collated/raw_isotope_data.txt", quote = FALSE, row.names = FALSE)

figure5 <- iso_df %>%
  mutate(group = ifelse(pop == "LN", "cbr", group)) %>%
  filter(!is.na(group)) %>%
  ggplot(aes(x = d13C, y = d15N, label = pop, color = group, fill = group)) +
  #geom_vline(xintercept = 0, color = "grey")+
  #geom_hline(yintercept = 0, color = "grey")+
  geom_point(size = 0, pch = 21, color = "black")+
  stat_ellipse(geom = "polygon", alpha = 0.1) +
  geom_point(size = 3.5, pch = 21, color = "black")+
  #geom_text(size = 4) +
  theme_bw(base_size = 16)+
  xlab(expression({delta}^13*C~'\u2030'))+
  ylab(expression({delta}^15*N~'\u2030'))+
  theme(legend.position = "none")+
  scale_color_manual(values = c("#77AB43", "#008FD5", "#FFFFFF")) +
  scale_fill_manual(values = c("#77AB43", "#008FD5", "#BDBDBD"))
  
ggsave(figure5, filename = "figures/figure5_raw.pdf", device = "pdf", height = 7, width = 7, useDingbats=FALSE)


######################################################################
# manova of isotope abundances
######################################################################

#group <- as.numeric(as.factor(iso_df$group))

# manova controlling for pop structure
manova_sir <- manova(cbind(d15N, d13C) ~ pop + group , data = iso_df)

anova(manova_sir)

# Analysis of Variance Table
# 
# Df  Pillai approx F num Df den Df  Pr(>F)    
#   (Intercept)   1 0.99917    64171      2    106 < 2e-16 ***
#   pop           4 1.18969       39      8    214 < 2e-16 ***
#   group         1 0.04481        2      2    106 0.08808 .  
#   Residuals   107                                           
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


