find_gene_by_chr_pos <- function(index, chr_target, glazer_genes, slack = 0, window_size = 75000){
  
  if(window_size > 0){
    
    chr_target <- chr_target
    pos_target <- index
    
    targ1 <- pos_target
    targ2 <- pos_target + (window_size - 2)
    
    sites <- targ1 : targ2
    
    gene_df1 <- glazer_genes %>%
      filter(chr == chr_target) %>%
      filter(pos1 > targ1) %>%
      filter(pos2 < targ2) %>%
      group_by(gene_id) %>%
      summarise(overlap = (any(sites %in% pos1:pos2))) %>%
      filter(overlap == TRUE)
    #filter(pos1 - slack <= targ1, pos2 + slack >= targ2)
    
    if(nrow(gene_df1) > 0){
      data.frame(chr = chr_target, pos1 = targ1, pos2 = targ2, gene = unique(gene_df1$gene_id))
    }else{
      data.frame(chr = chr_target, pos1 = targ1, pos2 = targ2, gene = NA)
    }
    
    
  } else{
    
    chr_target <- chr_target
    pos_target <- index
    
    gene_df <- glazer_genes %>%
      filter(chr == chr_target) %>%
      filter(pos1 - slack <= pos_target, pos2 + slack >= pos_target)
    
    gene <- unique(gene_df$gene_id)
    if(length(gene) == 0){
      gene <- NA
    }
    gene_df <- data.frame(chr = chr_target, pos = pos_target, gene)
    gene_df
    
  }
  
}