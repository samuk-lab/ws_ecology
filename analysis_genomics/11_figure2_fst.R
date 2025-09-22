library("tidyverse")
library("patchwork")

window_size <- 100000



# outflank output 

outflank <- readRDS("../data/stats/fst_outflank_snp.rds") %>%
  filter(comparison == "wht_cmn") %>%
  select(chr, pos, OutlierFlag, FST)

outflank_xix <- readRDS("../data/stats/fst_outflank_snp_XIX.rds") %>%
  filter(comparison == "wht_cmn") %>%
  select(chr, pos, OutlierFlag, FST)

outflank_y <- readRDS("../data/stats/fst_outflank_snp_Y.rds") %>%
  filter(comparison == "wht_cmn") %>%
  select(chr, pos, OutlierFlag, FST)

outflank <- bind_rows(outflank_xix, outflank_y, outflank)

outflank %>% lm(data = ., FST~1) %>% confint(level = 0.90)

site_stats <- outflank %>%
  mutate(OutlierFlag = ifelse(is.na(OutlierFlag), FALSE, OutlierFlag)) 

# hack for adjusting Y chromosome outlier cutoff
site_stats <- site_stats %>%
  mutate(OutlierFlag = ifelse(chr == "chrY" & FST > 0.5, TRUE, OutlierFlag))

#plot theme

genome_plot_theme <- theme_bw()+
  theme(rect = element_rect(fill = "black"),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0.0, "cm"),
        plot.background = element_rect(fill = "black", color = "black"),
        panel.background = element_rect(fill = "black", color = "black"),
        strip.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line.x = element_blank(),
        strip.placement = "inside",
        legend.position ="none",
        text = element_text(color = "white", size = 16),
        axis.text = element_text(color = "white", size = 16),
        strip.text = element_text(color = "white", size = 10))

############################################
# figure 2: fst outliers
############################################

site_stats$chr <- site_stats$chr %>% gsub("chr", "", .)
site_stats$chr <- factor(site_stats$chr, levels = c(as.roman(1:18) %>% as.character  %>% c(., "XX", "XXI", "XIX", "Y", "Un")))

fst_plot <- site_stats %>%
  mutate(FST = ifelse(FST < 0, 0, FST)) %>%
  mutate(chrom_col = (as.numeric(chr) %% 2) == 0 ) %>% 
  ggplot(aes(x = pos, y = FST, color = OutlierFlag))+
  geom_rect(aes(fill = chrom_col), xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf,alpha = 0.3, color = NA) +
  geom_point(size = 1.5, alpha = 0.8)+
  #geom_point(size = 1.5, alpha = 0.5, shape = 16)+
  facet_grid(.~chr, scales = "free_x", space = "free", switch = "x")+
  genome_plot_theme+
  scale_fill_manual(values=c("grey10", "grey15")) +
  scale_color_manual(values = c("grey75","red"))+
  ylab(expression("Weir and Cockerham's F"[ST]))+
  xlab("Chromosome") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.0))

fst_plot 

ggsave(fst_plot , file = "../figures/Figure2_raw_talk.png", height = 5, width = 8)

############################################
# figure s3: pi and dxy
############################################  

figureS3 <- stats_df %>%
  mutate(chrom_col = (as.numeric(chr) %% 2) == 0 ) %>% 
  filter(grepl("pi|dxy", stat)) %>%
  ggplot(aes(x = w_pos1, y = value, color = stat))+
  geom_rect(aes(fill = chrom_col),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf,alpha = 0.3, color = NA) +
  geom_line()+
  #geom_point(size = 1.5, alpha = 0.5, shape = 16)+
  facet_grid(stat_type~chr, scales = "free", space = "free_x", switch = "x")+
  genome_plot_theme+
  scale_fill_manual(values=c("grey93", "grey98")) +
  #scale_color_manual(values=c("grey", "red")) +
  #ylab(expression("Weir and Cockerham's F"[ST]))+
  xlab("Chromosome")+
  ylab("Statistic")


#fig2 <- fst_plot / window_stats + plot_layout(ncol = 1, heights = c(1, 1), guides = "auto")

