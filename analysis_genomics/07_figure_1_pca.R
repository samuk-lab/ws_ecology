#install.packages("scatterpie")
#install.packages("ggmap")
#devtools::install_github("thomasp85/patchwork")
#install.packages("mapdata")
#install.packages("egg")

library("tidyverse")
library("ggmap")
library("scatterpie")
library("patchwork")
library("mapdata")
library("egg")
library("wesanderson")


################################################################################
# set up data for plotting 
################################################################################

pop_dat <- read.csv("metadata/whtstbk_site_coordinates.csv", header = TRUE, stringsAsFactors = FALSE)
meta_df <- read.csv("metadata/mega_meta.csv", stringsAsFactors = FALSE, header = TRUE)
pca_df <- read.table("data/pca/pca_df.txt", stringsAsFactors = FALSE, header = TRUE)
pca_var <- read.table("data/pca/pca_var_df.txt", stringsAsFactors = FALSE, header = TRUE)


#harmonize names
meta_df$pop[meta_df$pop=="WR"] <- "RR"

# calculate cluster proportions

prop_df <- meta_df %>%
  filter(!is.na(cluster)) %>%
  filter(cluster != "dk") %>%
  group_by(pop, cluster) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  select(pop, cluster, freq) %>%
  spread(key = cluster, value = freq)

prop_df[is.na(prop_df)] <- 0 

pop_dat <- left_join(pop_dat, prop_df)

world <- map_data("worldHires", region = "Canada")

################################################################################
# Figure 1A: PCA
################################################################################

var_pc1 <- pca_var %>%
  filter(pc_num == 1) %>%
  pull(var_exp)

var_pc2 <- pca_var %>%
  filter(pc_num == 2) %>%
  pull(var_exp)

plot_pal <- wes_palette("Zissou1", n = 4, type = "discrete")

fig1_A <- pca_df %>% 
   filter(!is.na(pop)) %>%
   filter(!is.na(year)) %>%
   filter(!is.na(sex)) %>%
   filter(Indiv != "CL63_2015") %>% # outlier
   ggplot(aes(x = PC1, y = PC2, fill = cluster, label = Indiv, shape = sex)) +
   xlab(paste0("PC1 ", "(", round(var_pc1*100, digits = 2), "%)")) +
   ylab(paste0("PC2 ", "(", round(var_pc2*100, digits = 2), "%)")) +
   geom_point(size = 3.5, color = "grey10", alpha = 1.0)+
   #stat_ellipse(geom = "polygon", alpha = 0.4, color = "grey50") + 
   theme_bw()+
   theme(legend.position = c(1.0, 0.00),
         legend.justification = c(1, 0),
         #legend.box = "horizontal",
         #legend.direction = "horizontal",
         legend.background = element_blank(),
         panel.border = element_rect(color = "black", fill = NA),
         plot.margin = unit(c(0,0,0,0), "cm"),
         panel.spacing = unit(c(0,0,0,0), "cm"))+
   annotate("text", x = -0.045, y = 0.07, label = "Common\n (Bras d'Or)", fill = "#008FD5")+
   annotate("text", x = -0.03, y =-0.07, label = "Common\n (Mainland)", fill = "#77AB43")+
   annotate("text", x = 0.045, y = 0.075, label = "White", fill = "#FFFFFF")+
   #scale_color_brewer(type = "qual", palette = "Set1", guide = "none")+
   #scale_fill_manual(values = c("#f08a5d" ,"#3fc1c9", "#FFFFFF"), guide = "none")+
   scale_fill_manual(values = c("#77AB43", "#008FD5", "#FFFFFF"), guide = "none")+
   scale_shape_manual(values = c(21,24), name = "Sex") +
   guides(shape = guide_legend(
    label.position = "left",
    label.hjust = 1))

# scree_inset <- pca_var %>%
#   ggplot(aes(x = pc_num, y = var_exp))+
#   geom_bar(stat = "identity", fill = "white", color = "black")+
#   theme_minimal()+
#   theme(panel.grid = element_blank())+
#   ylab("PVE")+
#   xlab("PC Number")
# 
# fig1_A <- fig1_A + 
#   annotation_custom(
#     ggplotGrob(scree_inset), 
#     xmin = 0.05, xmax = 0.10, ymin = 0.03, ymax = 0.09
#   )

################################################################################
# Figure 1B: population pie charts for nova scotia
################################################################################

fig1_B <- world %>%
  filter(lat > 42.25, lat < 47.25) %>%
  filter(long > -66.5, long < -59.75) %>%
  ggplot(aes(long, lat)) +
  geom_map(map=world, aes(map_id=region), fill="grey75", color="grey55") +
  coord_fixed()+
  geom_scatterpie(aes(x=long, y=lat, group = ecotypes, r = 0.2),
                    data=pop_dat, cols=c("cbr", "cmn", "wht"), 
                    color=1, alpha=1.0) +
  
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw()+
  theme(legend.position = c(1.01, -0.02),
        legend.justification = c(1, 0),
        legend.background = element_blank(),
        #legend.direction = "horizontal",
        panel.border = element_rect(color = "black", fill = NA),
        plot.margin = unit(c(0,0,0,0), "cm"),
        panel.spacing = unit(c(0,0,0,0), "cm"),
        axis.text.x = element_text(margin=unit(c(0,0,0,0), "cm"))
        )+
  scale_fill_manual(values = c("#77AB43", "#008FD5", "#FFFFFF"), 
                    labels = c("Common (BD)", "Common (ML)", "White"), name = NULL )+
  guides(fill = guide_legend(
    label.position = "left",
    label.hjust = 1))


################################################################################
# Figure 1C: population pie charts for nova scotia (zoomed)
################################################################################

# zoom on straight of canso
# xlim=c(-62.125,-60.5), ylim=c(45.25,46.2)

pop_df_canso <- pop_dat %>%
  filter(!(pop %in% c("CL", "SH")))

fig1_C <- world %>%
    filter(lat > 45.25, lat < 46.2) %>%
    filter(long > -62.125, long < -60.5) %>%
    ggplot(aes(long, lat)) +
    geom_map(map=world, aes(map_id=region), fill="grey75", color="grey55") +
    coord_fixed() + 
    geom_scatterpie(aes(x=long, y=lat, group = ecotypes, r = 0.07),
                         data=pop_df_canso, cols=c("cbr", "cmn", "wht"), 
                         color=1, alpha=1.0) +
    #geom_point(aes(x=long, y=lat),data=pop_df_canso)+
    xlab("Longitude")+
    ylab("")+
    theme_bw()+
    theme(legend.position = c(1.01, -0.02),
          legend.justification = c(1, 0),
          #legend.justification = 'right',
          #legend.direction = "horizontal",
          legend.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          plot.margin = unit(c(0,0,0,0), "cm"),
          panel.spacing = unit(c(0,0,0,0), "cm")
    )+
  scale_fill_manual(values = c("#77AB43", "#008FD5", "#FFFFFF"), 
                    labels = c("Common (BD)", "Common (ML)", "White"), name = NULL )+
  guides(fill = guide_legend(
    label.position = "left",
    label.hjust = 1))


################################################################################
# compose figure 1
################################################################################

fig1 <- fig1_A / (fig1_B | fig1_C) + plot_layout(ncol = 1, heights = c(1.2, 1), guides = "auto")
ggsave("figures/Figure1_raw.pdf", height = 8, width = 8, plot = fig1, useDingbats=FALSE)
