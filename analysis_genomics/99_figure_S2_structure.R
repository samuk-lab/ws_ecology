# plot structure & harvester results

library("tidyverse")
library("tidyr")
library("patchwork")


meta_df <- read.csv("metadata/mega_meta.csv", stringsAsFactors = FALSE, header = TRUE)

meta_df <- meta_df %>%
  mutate(Indiv = gsub("whtstbk_gbs_2012_brds_", "", id)) %>%
  mutate(Indiv = gsub("_", ".", Indiv)) %>%
  mutate(Indiv = case_when(year == "2012" ~ paste0(Indiv, "_2013"),
                           year == "2014" ~ paste0(Indiv, "_2015"),
                           TRUE ~ Indiv)) %>%
  mutate(Indiv = ifelse(Indiv == "CL64.2_2015", "CL64-2_2015", Indiv))

parse_structure <- function(str_file){
  
  str_lines <- readLines(str_file)
  
  start_line <- which(grepl("*Label (%Miss)*", str_lines)) + 1
  end_line <- which(grepl("*Estimated Allele Frequencies*", str_lines)) - 3
  
  str_lines <- str_lines[start_line:end_line] 
  str_lines <- gsub("^[ 0-9]+ ", "", str_lines)
  str_lines <- strsplit(str_lines, split = " +") 
  
  nk <- length(str_lines[[1]]) - 3
  
  k_names <- paste0("k", 1:nk)
  
  str_df <- do.call(rbind, str_lines) %>% data.frame
  names(str_df) <- c("Indiv", "missing", "sep", k_names)
  
  str_df %>%
    mutate(missing = gsub("\\(|\\)", "", missing)) %>%
    select(-sep)
    
  
}


k1 <- list.files("data/structure_raw/k1/", full.names = TRUE)[1]
k2 <- list.files("data/structure_raw/k2/", full.names = TRUE)[1]
k3 <- list.files("data/structure_raw/k3/", full.names = TRUE)[1]
k4 <- list.files("data/structure_raw/k4/", full.names = TRUE)[1]
k5 <- list.files("data/structure_raw/k5/", full.names = TRUE)[1]

plot_theme <- theme_bw() +
    theme(axis.text.x  = element_blank(), 
      axis.ticks       = element_blank(), 
      axis.line        = element_blank(),
      axis.title       = element_blank(),
      strip.background = element_blank())

k2_plot <- parse_structure(k2) %>%
  gather(key = structure_pop, value = anc_prop, -Indiv, -missing) %>%
  mutate(anc_prop = as.numeric(anc_prop)) %>%
  left_join(meta_df) %>%
  filter(!pop=="MM", !pop=="NHR", !pop=="QR", !pop=="MM", !pop=="LD") %>%
  filter(!is.na(cluster)) %>%
  mutate(Indiv = fct_reorder2(Indiv, region , cluster)) %>%
  ggplot(aes(x = Indiv, y = anc_prop, fill = structure_pop))+
  geom_bar(position = "fill", stat = "identity", width = 1)+
  plot_theme +
  facet_grid(~region, space = "free", scales = "free_x")+
  scale_fill_manual(values = c("#008FD5", "#FFFFFF"))+
  scale_y_continuous(expand = c(0,0))+
  labs(fill = NULL) 


k3_plot <- parse_structure(k3) %>%
  gather(key = structure_pop, value = anc_prop, -Indiv, -missing) %>%
  mutate(anc_prop = as.numeric(anc_prop)) %>%
  left_join(meta_df) %>%
  filter(!pop=="MM", !pop=="NHR", !pop=="QR", !pop=="MM", !pop=="LD") %>%
  filter(!is.na(cluster)) %>%
  mutate(Indiv = fct_reorder2(Indiv, region , cluster)) %>%
  ggplot(aes(x = Indiv, y = anc_prop, fill = structure_pop))+
  geom_bar(position = "fill", stat = "identity", width = 1)+
  plot_theme +
  facet_grid(~region, space = "free", scales = "free_x")+
  scale_fill_manual(values = c("#FFFFFF", "#77AB43", "#008FD5"))+
  scale_y_continuous(expand = c(0,0))+
  labs(fill = NULL) 


k4_plot <- parse_structure(k4) %>%
  gather(key = structure_pop, value = anc_prop, -Indiv, -missing) %>%
  mutate(anc_prop = as.numeric(anc_prop)) %>%
  left_join(meta_df) %>%
  filter(!pop=="MM", !pop=="NHR", !pop=="QR", !pop=="MM", !pop=="LD") %>%
  filter(!is.na(cluster)) %>%
  mutate(Indiv = fct_reorder2(Indiv, region , cluster)) %>%
  ggplot(aes(x = Indiv, y = anc_prop, fill = structure_pop))+
  geom_bar(position = "fill", stat = "identity", width = 1)+
  plot_theme +
  facet_grid(~region, space = "free", scales = "free_x")+
  scale_fill_manual(values = c("#008FD5", "#01C3FE", "#FFFFFF", "#77AB43"))+
  scale_y_continuous(expand = c(0,0))+
  labs(fill = NULL) 

k5_plot <- parse_structure(k5) %>%
  gather(key = structure_pop, value = anc_prop, -Indiv, -missing) %>%
  mutate(anc_prop = as.numeric(anc_prop)) %>%
  left_join(meta_df) %>%
  filter(!pop=="MM", !pop=="NHR", !pop=="QR", !pop=="MM", !pop=="LD") %>%
  filter(!is.na(cluster)) %>%
  mutate(Indiv = fct_reorder2(Indiv, region , cluster)) %>%
  ggplot(aes(x = Indiv, y = anc_prop, fill = structure_pop))+
  geom_bar(position = "fill", stat = "identity", width = 1)+
  plot_theme +
  facet_grid(~region, space = "free", scales = "free_x")+
  scale_fill_manual(values = c("#77AB43", "#FEB801", "#008FD5", "#01C3FE", "#FFFFFF"))+
  scale_y_continuous(expand = c(0,0))+
  labs(fill = NULL) 

struc_plot <- k2_plot/k3_plot/k4_plot/k5_plot

ggsave(filename = "figures/FigureS2_raw.pdf", struc_plot, width = 8, height = 6)  
