# Center (for PCA) and spherical normalize (amplitude-blindness)
normalize <- function(X) {
  mean_vec <- colMeans(X)
  X_centered <- sweep(X, 2, mean_vec, FUN = "-")
  row_norms <- sqrt(rowSums(X_centered^2))
  X_normalized <- X_centered / (row_norms + 1e-8)
  return(X_normalized)
}

# Do tde and optionally reduce the dimension of the matrix
tde <- function(timeseries, dim, lag, nc=NULL){
  tde_matrix = buildTakens(time.series = timeseries, embedding.dim = dim, time.lag = lag)
  if (is.null(nc)) return(tde_matrix) # No dimension reduction
  tde_matrix = normalize(tde_matrix)
  tde_matrix = prcomp(tde_matrix, center=FALSE)$x[, 1:nc]
  return(tde_matrix)
}


plot_band <- function(persistence_diagram, c_hat, dim, x_lim = 4, y_lim = 4, rb_scale = 5, sublevel=FALSE, col="#619CFF", shp = 17 ){
  pers_hom = persistence_diagram %>% data.frame() %>% filter(dimension==dim)
  if(!sublevel){
    # Switch the birth and death times
    death_ = pers_hom$death; birth_ = pers_hom$birth; 
    pers_hom$death = birth_; pers_hom$birth = death_
    y_label = "Birth"; x_label = "Death"
  }else{
    y_label = "Death"; x_label = "Birth"
  }
  pers_hom$ind = as.factor(5*pers_hom$dimension)
  pers_H1 = pers_hom %>% filter(dimension==1) %>% mutate(pers = death - birth)
  bxy = pers_H1[which.max(pers_H1$pers), c(2, 3)]
  # Compute the breaks
  breaks <- pretty(pers_hom$death)  # Use the pretty function to get usual tick marks
  breaks <- c(breaks, max(breaks) +  breaks[2])  # Ensure the maximum value is included
  # Compute the labels
  labels <- as.character(breaks)
  labels[labels == as.character(max(breaks))] <- expression(infinity)  

  ribbondata <- data.frame(
    x=c(-rb_scale, pers_hom$birth, rb_scale),
    ymin=c(-rb_scale, pers_hom$birth, rb_scale),
    ymax=c(-rb_scale + c_hat, pers_hom$birth + c_hat, rb_scale + c_hat)
    )
  ggplot(data = pers_hom, aes(x = birth, y = death)) + 
    geom_point( aes(shape = ind, color = ind), size = 2.5) + 
    geom_abline(intercept = 0, slope = 1) +
    geom_abline(intercept = c_hat, slope = 1, linewidth = 0, color=NA) +
    geom_ribbon(data = ribbondata, 
                aes(x = x, ymin = ymin, ymax = ymax, fill = "Band"), 
                outline.type = 'lower', color = NA, inherit.aes = FALSE, alpha = 0.15) +
    ylab(y_label) + xlab(x_label) + coord_cartesian(xlim = c(0, x_lim), ylim = c(0, y_lim)) +
    scale_shape_manual(name = "", labels = c(bquote(H[.(dim)])), values = shp) + 
    scale_colour_manual(name="", labels= c(bquote(H[.(dim)])), values=col) +
    theme( text = element_text(size=18), 
           legend.box.background = element_rect(color = "black"), 
           legend.background = element_blank(),
           legend.title=element_blank(), 
           legend.position = c(0.8, 0.2), 
           legend.margin = margin(1, 1, 1, 1),
           legend.spacing.x = unit(0, "mm"), 
           legend.spacing.y = unit(0, "mm"))
}


plot_diagram <- function(pers_hom, flip=FALSE, x_lim=c(0, 1), y_lim=c(0, 1)){
  pers_hom$ind <- as.factor(5*pers_hom$dimension)
  breaks <- pretty(y_lim)#pretty(pers_hom$death)  # Use the pretty function to get usual tick marks
  breaks <- c(breaks) #c(breaks, max(breaks) +  breaks[2])  # Ensure the maximum value is included
  # Compute the labels
  labels <- as.character(breaks)
  labels[labels == as.character(max(breaks))] <- expression(infinity) 
  if(!flip){
    ggplot(pers_hom) + 
      geom_point(aes(x = birth, y = death, shape = ind, color = ind), size = 2.0) + 
      #geom_hline(yintercept = max(breaks), linetype = "dashed", color = "red") +
      geom_point(aes(x = 0, y = max(breaks)), color = "red", size=0.5) + xlim(c(0, max(breaks))) +
      scale_y_continuous(breaks = breaks, labels =  labels) +
      geom_abline(intercept = 0, slope = 1) +
      xlab("Birth") + ylab("Death") + xlim(x_lim) + #ylim(y_lim) +
      scale_shape_manual(name = "", labels = c(expression(H[0]), expression(H[1])), values = c(16, 17)) + 
      scale_colour_manual(name="", labels=c(expression(H[0]), expression(H[1])), values=c("Red", "#619CFF")) +
      theme(text = element_text(size=14),
            legend.box.background = element_rect(color = "black"),
            legend.title=element_blank(), 
            legend.position = c(0.8, 0.2), 
            legend.margin = margin(0, 0, 0, 0),
            legend.spacing.x = unit(0, "mm"),
            legend.spacing.y = unit(0, "mm"))
  }
  else{
    ggplot(pers_hom) + 
      geom_point(aes(x = death, y = birth, shape = ind, color = ind), size = 2.0) + 
      geom_abline(intercept = 0, slope = 1) + ylab("Birth") + xlab("Death") + xlim(x_lim) + ylim(y_lim) +
      scale_shape_manual(name = "", labels = c(expression(H[0]), expression(H[1])), values = c(16, 17)) + 
      scale_colour_manual(name="", labels=c(expression(H[0]), expression(H[1])), values=c("Red", "#619CFF")) +
      theme(text = element_text(size=14),
            legend.box.background = element_rect(color = "black"),
            legend.title=element_blank(), 
            legend.position = c(0.8, 0.2), 
            legend.margin = margin(0, 0, 0, 0),
            legend.spacing.x = unit(0, "mm"),
            legend.spacing.y = unit(0, "mm"))
  }
}






