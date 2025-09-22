library("tidyverse")
library("patchwork")

window_size <- 100000


# pixy outpuit
pi_files <- list.files("analysis_pixy/output", pattern = "pi.txt", full.names = TRUE)
pi_dat <- lapply(pi_files, read.table, h = T)
pi_dat <- bind_rows(pi_dat)


dxy_files <- list.files("analysis_pixy/output", pattern = "dxy.txt", full.names = TRUE)
dxy_dat <- lapply(dxy_files, read.table, h = T)
dxy_dat <- bind_rows(dxy_dat)

# outflank output (fst + pbs)
fst_100k <- readRDS("data/stats/fst_100k.rds")
pbs <- readRDS("data/stats/pbs_df_snp.rds")

outflank <- readRDS("data/stats/fst_outflank_snp.rds") %>%
  filter(comparison == "wht_cmn") %>%
  select(chr, pos, OutlierFlag, FST)

outflank_xix <- readRDS("data/stats/fst_outflank_snp_XIX.rds") %>%
  filter(comparison == "wht_cmn") %>%
  select(chr, pos, OutlierFlag, FST)

outflank <- bind_rows(outflank_xix, outflank)

outflank %>% lm(data = ., FST~1) %>% confint(level = 0.90)

  
site_stats <- left_join(outflank, pbs) %>%
  mutate(OutlierFlag = ifelse(is.na(OutlierFlag), FALSE, OutlierFlag))

# bind all windowed stats into a single dataframe

dxy_dat <- dxy_dat %>% 
  mutate(stat_type = "dxy") %>%
  mutate(stat = paste0(pop1, "_", pop2, "_dxy")) %>%
  select(stat, stat_type, chromosome, window_pos_1, avg_dxy) %>%
  rename(value = avg_dxy) %>%
  rename(w_pos1 = window_pos_1) %>%
  rename(chr = chromosome)

pi_dat <- pi_dat  %>%
  mutate(stat_type = "pi") %>%
  mutate(stat = paste0("pi_", pop)) %>%
  rename(value = avg_pi) %>%
  select(stat, stat_type, chromosome, window_pos_1, value) %>%
  rename(w_pos1 = window_pos_1) %>%
  rename(chr = chromosome)

fst_100k <- fst_100k %>%
  ungroup %>%
  mutate(stat = paste0("fst_", comparison)) %>%
  mutate(stat_type = "fst") %>%
  rename(value = mean_fst) %>%
  select(stat, stat_type,chr, w_pos1, value)

pbs <- pbs %>%
  mutate(w_pos2 = (((pos / window_size) %>% floor) + 1) * window_size)%>%
  mutate(w_pos1 = w_pos2 - (window_size - 1)) %>%
  mutate(value = wht_pbs) %>%
  mutate(stat = "pbs_wht") %>%
  mutate(stat_type = "pbs") %>%
  filter(value > 0) %>%
  select(stat, stat_type, chr, w_pos1, value)

chr_levels <- c(as.roman(1:21) %>% as.character, "Un", "Y")

stats_df <- bind_rows(dxy_dat, pi_dat, fst_100k, pbs) %>%
  ungroup %>%
  mutate(chr = gsub("chr", "", chr)) %>%
  mutate(chr = factor(chr, levels = chr_levels))

#plots

genome_plot_theme <- theme_bw()+
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0.0, "cm"),
        strip.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line.x = element_blank(),
        strip.placement = "inside",
        legend.position ="none")

############################################
# figure 2: fst outliers
############################################

site_stats$chr <- site_stats$chr %>% gsub("chr", "", .)
site_stats$chr <- factor(site_stats$chr, levels = c(as.roman(1:21) %>% as.character  %>% c(., "Un")))

fst_plot <- site_stats %>%
  mutate(FST = ifelse(FST < 0, 0, FST)) %>%
  mutate(chrom_col = (as.numeric(chr) %% 2) == 0 ) %>% 
  ggplot(aes(x = pos, y = FST, color = OutlierFlag))+
  geom_rect(aes(fill = chrom_col), xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf,alpha = 0.3, color = NA) +
  geom_point(size = 0.75, alpha = 0.8)+
  #geom_point(size = 1.5, alpha = 0.5, shape = 16)+
  facet_grid(.~chr, scales = "free_x", space = "free", switch = "x")+
  genome_plot_theme+
  #ylab(expression("Weir and Cockerham's F"[ST]))+
  scale_fill_manual(values=c("grey93", "grey98")) +
  scale_color_manual(values = c("grey75","red"))+
  ylab(expression("Weir and Cockerham's F"[ST]))+
  xlab("Chromosome") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.5))

fig2 <- fst_plot 

ggsave(fig2, "figures/Figure2_raw.pdf", height = 4, width = 8, plot = fig2, useDingbats=FALSE)

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

