build_ld_df <- function(pop_name, stats_files, window = FALSE){
  
  # determine files
  stats_files_pop <- stats_files[grep(paste0(stats_folder,"/",pop_name), stats_files)]  
  stats_files_pop <- stats_files_pop[grep("r2", stats_files_pop)] 
  
  # read into a list df
  stats_dfs <- read.table(stats_files_pop[1], header = TRUE, stringsAsFactors = FALSE)
  
  # the magic: calculate site-wise ld (r2)
  # requires >3 pairwise r2's, each with >30 genotypes 
  calc_sitewise_ld <- . %>%
    filter(N_INDV > 30) %>%
    mutate(dist = abs(POS1 - POS2)) %>%
    group_by(CHR, POS1) %>%
    summarise(r2 = mean(R.2), n_sites = n()) %>%
    filter(n_sites > 3) %>%
    ungroup
  
  ld_df <- calc_sitewise_ld(stats_dfs)
  
  if (window == TRUE){
    ld_df$w_pos2 <- (((ld_df$POS1 / 50000) %>% floor) + 1)*50000
    ld_df$w_pos1 <- ld_df$w_pos2 - 49999
    
    ld_df <- ld_df %>%
      group_by(CHR, w_pos1, w_pos2) %>%
      summarise(r2 = mean(r2)) %>%
      ungroup
    
    names(ld_df)[4] <- paste0("r2_", pop_name)
    names(ld_df)[1] <- "chr"
    names(ld_df)[2] <- "pos1"
    names(ld_df)[3] <- "pos2"
    ld_df 
    
  } else{
    names(ld_df)[3] <- paste0("r2_", pop_name)
    names(ld_df)[1] <- "chr"
    names(ld_df)[2] <- "pos1"
    
    ld_df <- ld_df %>%
      select(-n_sites) %>% 
      ungroup()
    ld_df
  }
  
}