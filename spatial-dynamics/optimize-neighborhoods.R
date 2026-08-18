optimize_neighborhood=function(neighborhoods,thresholds_= seq(0.5,0.95,0.05),out_dir,make_gif=F,min_cells=100){
  
  require(dplyr)
  require(ggplot2)
  require(glue)
  
  ## thresholding on delta smoooth:
  target_celltypes=unique(neighborhoods$targets)[unique(neighborhoods$targets)!='else']
  
  delta_z_spreading=NULL
  target_densities=NULL
  thresholded_neighborhoods=list()
  objective=NULL
  
  dir.create(glue('{out_dir}/refinement/'),showWarnings = F)
  dir.create(glue('{out_dir}/refinement-objective/'),showWarnings = F)
  
  for(i in thresholds_){
  
    tmp.circuit_region=neighborhoods %>% subset(delta_z_smooth > quantile(neighborhoods$delta_z_smooth,i))
    if(dim(tmp.circuit_region)[1] >= min_cells & sum(!is.na(table(tmp.circuit_region$targets)[target_celltypes])) == length(target_celltypes) ){
      message(glue('N = {dim(tmp.circuit_region)[1]}'))
      message(glue('Min target N = {min(table(tmp.circuit_region$targets)[target_celltypes])}'))
      celltypes.sorted <- names(sort(table(neighborhoods$targets)[target_celltypes], decreasing = T))
      colors_ <- sapply(target_celltypes, function(x) ggsci::pal_igv()(match(x, target_celltypes))[match(x, target_celltypes)])
      
      overlay_points <- do.call("rbind", lapply(celltypes.sorted, function(celltype) {
        celltype.tmp <- gsub("[+]", "[+]", celltype)
        overlay_points.tmp <- data.frame(
          x = tmp.circuit_region$x[grep(celltype.tmp, as.character(tmp.circuit_region$targets))], # X-coordinates
          y = tmp.circuit_region$y[grep(celltype.tmp, as.character(tmp.circuit_region$targets))], # Y-coordinates
          color = colors_[[celltype]] # Color (can be a vector for multiple colors)
        )
        return(overlay_points.tmp)
      }))
      
      p2 <- ggplot(tmp.circuit_region) +
        geom_point(aes(x, y),size = 1.5, alpha = 0.2) +
        geom_point(data = overlay_points, aes(x = x, y = y), color = overlay_points$color, size = 3, shape = 16) +
        theme_minimal() +
        coord_fixed(ratio = max(abs(tmp.circuit_region$y)) / max(tmp.circuit_region$x))
      p2 = p2+ggtitle(glue('delta_z_smooth > {round(quantile(neighborhoods$delta_z_smooth,i),3)}'))
      print(p2)
      
      target_celltype_density = sum(table(tmp.circuit_region$targets)[unlist(target_celltypes)])/dim(neighborhoods)[1]
      
      
      target_densities[length(target_densities)+1]=target_celltype_density
      delta_z_spreading[length(delta_z_spreading)+1]=(exp( 0.5*quantile(tmp.circuit_region$delta_z,0.5) + 0.25*(quantile(tmp.circuit_region$delta_z,0.25)+quantile(tmp.circuit_region$delta_z,0.75)) )*exp( 0.5*quantile(tmp.circuit_region$delta_z_smooth,0.5) + 0.25*(quantile(tmp.circuit_region$delta_z_smooth,0.25)+quantile(tmp.circuit_region$delta_z_smooth,0.75)) ))
      
      objective.tmp=target_densities*log(delta_z_spreading)
      objective[[length(objective)+1]]=objective.tmp
      png(glue('{out_dir}/refinement-objective/{i}.png'))
      print(plot(objective.tmp,type='l'))
      dev.off()
      
      thresholded_neighborhoods[[length(thresholded_neighborhoods)+1]]=tmp.circuit_region
      
      ggsave(glue('{out_dir}/refinement/{i}.png'),p2)
    }
    
    
  }
  par(mfrow=c(1,2))
  plot(delta_z_spreading,type='l')
  plot(target_densities,type='l')
  objective_global=objective[[length(objective)]]
  plot(objective_global,type='l')
  max_idx=match(max(objective_global),objective_global)
  optimal_threshold=thresholds_[max_idx]
  optimally_refined_neighborhood=thresholded_neighborhoods[[max_idx]]
  
  message(glue('optimal neighborhood found at idx = {max_idx} | threshold: {thresholds_[max_idx]}'))
  abline(h=max(objective_global))
  abline(v=max_idx,lty=4)
  
  png(glue('{out_dir}/final-objective.png'))
  par(mfrow=c(1,1))
  plot(objective_global,type='l')
  abline(h=max(objective_global))
  abline(v=max_idx,lty=4)
  dev.off()
  
  library(ggplot2)
  
  # Function to create the appropriate color scale
  create_min_zero_max_scale <- function(data_values,name="Cluster-wise ES Z") {
    min_val <- min(data_values, na.rm = TRUE)
    max_val <- max(data_values, na.rm = TRUE)
    
    # Define the colors we'll use
    colors <- c("blue", "white", "red")
    
    # Define where those colors will appear on the scale (as normalized positions)
    if (min_val >= 0) {
      # All values are positive - use only white to red portion
      positions <- c(0, 0, 1)  # Blue and white both at 0, red at max
      breaks <- c(0, max_val)
      labels <- c("0", paste("Max:", round(max_val, 2)))
    } else if (max_val <= 0) {
      # All values are negative - use only blue to white portion
      positions <- c(0, 1, 1)  # Blue at min, white and red both at 0
      breaks <- c(min_val, 0)
      labels <- c(paste("Min:", round(min_val, 2)), "0")
    } else {
      # Values span across zero - use full blue-white-red scale
      # Calculate proportions for proper spacing
      total_range <- max_val - min_val
      zero_pos <- abs(min_val) / total_range
      
      positions <- c(0, zero_pos, 1)
      breaks <- c(min_val, 0, max_val)
      labels <- c(paste("Min:", round(min_val, 2)), "0", paste("Max:", round(max_val, 2)))
    }
    
    # Return the scale
    scale_color_gradientn(
      colors = colors,
      values = positions,
      breaks = breaks,
      labels = labels,
      limits = c(min_val, max_val),
      name = name
    )
  }
  
  p=ggplot(optimally_refined_neighborhood, aes(x = x, y = y, col = delta_z)) +
    geom_point() +
    create_min_zero_max_scale(optimally_refined_neighborhood$delta_z) +
    theme_minimal() +
    coord_fixed(ratio = max(abs(optimally_refined_neighborhood$y)) / max(optimally_refined_neighborhood$x))
 
  ggsave(glue('{out_dir}/optimal_neighborhood.pdf'),p,h=20,w=20)
  ggsave(glue('{out_dir}/optimal_neighborhood.png'),p)
  
  overlay_points <- do.call("rbind", lapply(celltypes.sorted, function(celltype) {
    celltype.tmp <- gsub("[+]", "[+]", celltype)
    overlay_points.tmp <- data.frame(
      x = optimally_refined_neighborhood$x[grep(celltype.tmp, as.character(optimally_refined_neighborhood$targets))], # X-coordinates
      y = optimally_refined_neighborhood$y[grep(celltype.tmp, as.character(optimally_refined_neighborhood$targets))], # Y-coordinates
      color = colors_[[celltype]], # Color (can be a vector for multiple colors),
      celltype = optimally_refined_neighborhood$targets[grep(celltype.tmp, as.character(optimally_refined_neighborhood$targets))]
    )
    return(overlay_points.tmp)
  }))
  
  p2 <- ggplot(optimally_refined_neighborhood) +
    geom_point(aes(x, y),size = 1.5, alpha = 0.2) +
    geom_point(data = overlay_points, aes(x = x, y = y,color=celltype), size = 3, shape = 16) +
    theme_minimal() +
    coord_fixed(ratio = max(abs(optimally_refined_neighborhood$y)) / max(optimally_refined_neighborhood$x))
  p2 = p2+ggtitle(glue('delta_z_smooth > {round(quantile(neighborhoods$delta_z_smooth,thresholds_[max_idx]),3)}'))+scale_color_manual(values = colors_, name = "Celltype")
  
  ggsave(glue('{out_dir}/optimal_neighborhood-targets.pdf'),p2,h=20,w=20)
  ggsave(glue('{out_dir}/optimal_neighborhood-targets.png'),p2)
  
  if(make_gif){
    ## make refinement GIF
    library(magick)
    library(magrittr)
    
    list.files(path=glue('{out_dir}refinement/'), pattern = '*.png', full.names = TRUE) %>% 
      image_read() %>% # reads each path file
      image_join() %>% # joins image
      image_animate(fps=4) %>% # animates, can opt for number of loops
      image_write(glue("{out_dir}/refinement.gif")) # write to current dir
    
    list.files(path=glue('{out_dir}/refinement-objective/'), pattern = '*.png', full.names = TRUE) %>% 
      image_read() %>% # reads each path file
      image_join() %>% # joins image
      image_animate(fps=4) %>% # animates, can opt for number of loops
      image_write(glue("{out_dir}/refinement-objective.gif")) # write to current dir
    
  }
  return(optimally_refined_neighborhood)
}
