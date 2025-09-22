library("tidyverse")
library("reshape2")
library("patchwork")

fs_file <- "dadi_files/results/IM_data_spectra.fs"

plot_fs_file <- function(fs_file, ndims1 = 41, ndims2 = 41, vmin = 1){
  
  im_data_fs <- read.table(fs_file, skip = 1, header = F) %>% t %>%
    data.frame()
  
  im_data_df <- matrix(data = im_data_fs$X1, nrow = ndims1, ncol = ndims2) %>% melt
  im_mask_df <- matrix(data = im_data_fs$X2, nrow = ndims1, ncol = ndims2) %>% melt
  
  names(im_data_df) <- c("count_1", "count_2", "density")
  names(im_mask_df) <- c("count_1", "count_2", "mask")
  
  data_plot <- left_join(im_data_df, im_mask_df) %>%
    filter(mask == 0) %>%
    filter(density > vmin) %>% 
    mutate(density = log(log(density+1))) %>% 
    ggplot(aes(x = count_1, y = count_2, fill = density))+
    geom_raster(interpolate = FALSE)+
    geom_abline(slope = 1, intercept = 0, size = 0.5)+
    theme_bw()+
    xlab("Common Allele Count")+
    ylab("White Allele Count")+
    labs(fill = "log(Frequency)")+
    theme(panel.grid = element_blank(),
          legend.position = c(1.0, 0.00),
          legend.justification = c(1, 0),
          legend.background = element_blank())+
    scale_fill_viridis_c()+
    scale_x_continuous(limits = c(0, NA), expand = c(0, 0)) + 
    scale_y_continuous(limits = c(0, NA), expand = c(0, 0))+
    guides(fill = guide_colorbar(
      label.position = "left",
      label.hjust = 1,
      title.hjust = 1,
      title.position = "left",
      title.theme = element_text(angle = 90)))
  
  return(data_plot)
  
}


im_data <- plot_fs_file("dadi_files/results/IM_data_spectra.fs")
im_fit <- plot_fs_file("dadi_files/results/IM_fit_spectra.fs")

snm_data <- plot_fs_file("dadi_files/results/SNM_data_spectra.fs")
snm_fit <- plot_fs_file("dadi_files/results/SNM_fit_spectra.fs")

split_data <- plot_fs_file("dadi_files/results/split_no_mig_data_spectra.fs")
split_fit <- plot_fs_file("dadi_files/results/split_no_mig_fit_spectra.fs")

im_plot <- im_data + im_fit
snm_plot <- snm_data + snm_fit
split_plot <- split_data + split_fit

snm_plot / split_plot / im_plot

ggsave(im_plot, filename = "figures/dadi_raw/im_dadi_raw.png", width = 8, height = 4, device = "png")
ggsave(snm_plot, filename = "figures/dadi_raw/snm_dadi_raw.png", width = 8, height = 4, device = "png")
ggsave(split_plot, filename = "figures/dadi_raw/split_dadi_raw.png", width = 8, height = 4, device = "png")

ggsave(im_plot, filename = "figures/dadi_raw/im_dadi_raw.pdf", width = 8, height = 4, device = "pdf")
ggsave(snm_plot, filename = "figures/dadi_raw/snm_dadi_raw.pdf", width = 8, height = 4, device = "pdf")
ggsave(split_plot, filename = "figures/dadi_raw/split_dadi_raw.pdf", width = 8, height = 4, device = "pdf")
