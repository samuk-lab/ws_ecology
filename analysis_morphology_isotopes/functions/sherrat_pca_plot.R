
# generate a pca plot like the one described here:
# http://www.r-bloggers.com/tips-tricks-7-plotting-pca-with-tps-grids/

sherrat_pca_plot <- function(pc1 = 1, pc2 = 2, landmark_gpa = landmark_gpa, group_factor = factor(morpho_sub$cluster), ..., correct_allometry = TRUE, file_name = ""){
  
  if(correct_allometry){
    
    #perform pca
    morpho_pca <- plotTangentSpace(landmark_gpa$coords, axis1 = pc1, axis2 = pc2,
                                   groups = group_factor, verbose = TRUE)
    
    #harvest pca data (1:6)
    pca_df <- data.frame(id = morpho_sub$id, sex = morpho_sub$sex, csize = landmark_gpa$Csize,
                         cluster = morpho_sub$cluster, morpho_pca$pc.scores[,1:6])
    
    # wide to long
    pca_df_long <- gather(pca_df, key = pc, value = score, -id, -sex, -csize, -cluster)
    
    # correct for body size ("csize" i.e. the gpa scaling factor)
    resid <- pca_df_long %>%
      group_by(cluster, pc) %>%
      do(augment(lm(score ~ csize, data=.))) %>%
      ungroup %>%
      select(.resid) %>%
      unlist %>% as.numeric()
    
    pca_df_long$resid_score <- resid
    
    pca_df_resid <- pca_df_long %>%
      select(-score)
    
    # remake the pca data frame using the residual scores
    pca_df_resid <- spread(pca_df_resid, key = pc, value = resid_score)
    
    # apply the whtstbk palatte
    
    #77AB43 = green
    #008FD5 = blue 
    #FFFFFF = white 
    
    whtstbk_palatte <- c("#77AB43", "#008FD5", "#FFFFFF")
    gp <- group_factor
    col.gp <- whtstbk_palatte
    names(col.gp) <- levels(gp)
    col.gp <- col.gp[match(gp, names(col.gp))] # col.gp must not be a factor
    
    pdf(file = file_name, height = 7, width = 8.5, useDingbats = FALSE)
    
    # layout!
    mat <- matrix(c(4,5,0,1,1,2,1,1,3), 3)
    layout(mat, widths=c(1,1,1), heights=c(1,1,0.6))# set the size of the rows and columns
    
    # Item 1 to plot, the graph of PC1 vs PC2
    par(mar=c(4, 4, 1, 1)) # sets the margins
    
    # x and y labels
    xlab <- paste("Size Corrected Principal Component ", pc1, " (", round(morpho_pca$pc.summary$importance[2,pc1]*100, 1), "%)", sep="")
    ylab <- paste("Size Corrected Principal Component ", pc2, " (", round(morpho_pca$pc.summary$importance[2,pc2]*100, 1), "%)", sep="")
    
    
    # plot pca scores
    plot(pca_df_resid[,4+pc1], pca_df_resid[,4+pc2], pch = 21, 
         cex = 2, bg = col.gp, xlab = xlab, ylab = ylab)
    #legend(-0.09, 0.07, legend= unique(gp), pch=19,  col=unique(col.gp))

    legend("topleft", legend = levels(gp), pt.cex = 2,
           pch = 21, pt.bg = whtstbk_palatte, bty = "n")
    
    ref <- mshape(landmark_gpa$coords)# assign mean shape for use with plotRefToTarget below
    
    # Item 2 to plot, the first TPS grid; here we use the outline option to add to the visualisation
    grid_par <- gridPar(grid.col = "grey", tar.pt.size = 2)
    par(mar = c(0,0,0,0)) # sets the margins
    
    refs_to_plot <- c((pc1):(pc1+1), (pc2):(pc2+1))
    nope <- lapply(refs_to_plot, function(x) plotRefToTarget(ref, morpho_pca$pc.shapes[[x]], useRefPts = FALSE, gridPars = grid_par, method = "TPS"))
    
    dev.off()
    
  } else{
    
    # create pca object
    morpho_pca <- plotTangentSpace(landmark_gpa$coords, axis1 = pc1, axis2 = pc2,
                                   groups = group_factor, verbose = TRUE)
    
    # x and y labels
    xlab <- paste("Principal Component ", pc1, " (", round(morpho_pca$pc.summary$importance[2,pc1]*100, 1), "%)", sep="")
    ylab <- paste("Principal Component ", pc2, " (", round(morpho_pca$pc.summary$importance[2,pc2]*100, 1), "%)", sep="")
    
    # apply the whtstbk palatte
    
    #77AB43 = green
    #008FD5 = blue 
    #FFFFFF = white 
    
    whtstbk_palatte <- c("#008FD5", "#FFFFFF")
    gp <- factor(morpho_sub$cluster)
    col.gp <- whtstbk_palatte
    names(col.gp) <- levels(gp)
    col.gp <- col.gp[match(gp, names(col.gp))] # col.gp must not be a factor
    
    # plot to pdf
    pdf(file = "figures/Figure1.pdf", height = 8.5, width = 8.5, useDingbats = FALSE)
    
    # layout!
    mat <- matrix(c(4,5,0,1,1,2,1,1,3), 3)
    layout(mat, widths=c(1,1,1), heights=c(1,1,0.6))# set the size of the rows and columns
    
    # Item 1 to plot, the graph of PC1 vs PC2
    par(mar=c(4, 4, 1, 1)) # sets the margins
    
    # plot pca scores
    plot(morpho_pca$pc.scores[,pc1], morpho_pca$pc.scores[,pc2], pch = 21, 
         cex = 2, bg = col.gp, xlab = xlab, ylab = ylab, asp = TRUE)
    #legend(-0.09, 0.07, legend= unique(gp), pch=19,  col=unique(col.gp))
    
    legend("topleft", legend = levels(gp), pt.cex = 2,
           pch = 21, pt.bg = whtstbk_palatte, bty = "n")
    
    ref <- mshape(landmark_gpa$coords)# assign mean shape for use with plotRefToTarget below
    
    # Item 2 to plot, the first TPS grid; here we use the outline option to add to the visualisation
    grid_par <- gridPar(grid.col = "grey", tar.pt.size = 2)
    par(mar = c(0,0,0,0)) # sets the margins
    
    nope <- lapply(1:4, function(x) plotRefToTarget(ref, useRefPts = FALSE, morpho_pca$pc.shapes[[x]], gridPars = grid_par, method = "TPS"))
    dev.off()
  }
  
  
}


