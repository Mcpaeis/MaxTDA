rm(list = ls())
library(TDA)
library(ggplot2)
library(tidyverse)
theme_set(theme_bw())
img.dir = '../images/'

# Load the data
noisy_manifold = read.csv('../data/graphicalabstract/noisy_manifold.csv') %>% select(X0, X1)
kde_samples_thresh = read.csv('../data/graphicalabstract/thresholded_data.csv') %>% select(X0, X1)
kde_samples_nothresh = read.csv('../data/graphicalabstract/non_thresholded_data.csv') %>% select(X0, X1)

# Add the circular geometry and plot
circle_unif <- function(n, r=1){
  theta <- runif(n, 0, 2*pi)  # uniform distribution of angles
  x <- cos(theta)*r
  y <- sin(theta)*r
  circle_data <- data.frame(x1=x, x2=y) %>% arrange(atan2(x2, x1))
  return(circle_data)
}
sparse_geometry = circle_unif(1200) %>% rename(sparse_x0=x1, sparse_x1=x2)
dense_geometry = (circle_unif(1200, 0.5) + 1.25) %>% rename(dense_x0=x1, dense_x1=x2)
noisy_manifold_ = cbind(noisy_manifold, sparse_geometry, dense_geometry)
noisy_manifold_ %>% ggplot(aes(x = X0, y=X1)) + 
  xlim(c(-1.4, 2.1)) + ylim(c(-1.4, 2.1)) +
  geom_point(alpha=0.75, size=0.5) + xlab('x') + ylab('y') +
  geom_path(aes(x=sparse_x0, y=sparse_x1), color='blue', alpha=0.4) +
  geom_path(aes(x=dense_x0, y=dense_x1), color='red', alpha=0.4) +
  theme(axis.title.y = element_text(angle = 0, vjust=0.5), axis.title = element_text(size = 18), axis.text = element_text(size = 16))
#ggsave(paste(img.dir, 'noisy_manifold_ga.png', sep=''), dpi=300, height = 4, width = 4)

# Make the density surface plots
ggplot(data=kde_samples_thresh,  aes(x=X0, y=X1)) + xlab('x') + ylab('y') +
  stat_density2d(geom="tile", aes(fill=..density..^0.25), contour=FALSE) +
  xlim(c(-1.4, 2.1)) + ylim(c(-1.4, 2.1)) +
  scale_fill_gradientn(colours = colorRampPalette(c("white", blues9))(256), 
                       guide = guide_colorbar(title = 'Density', 
                                              barheight = unit(3, "mm"),
                                              barwidth = unit(30, "mm"),
                                              direction = "horizontal", 
                                              title.position = "top")) +
  theme(text = element_text(size=14),
        legend.box.background = element_rect(color = "black"),
        legend.title=element_blank(), 
        axis.title.y = element_text(angle = 0, vjust=0.5),
        legend.position = c(0.74, 0.05), 
        legend.margin = margin(0, 0, 0, 0),
        legend.text = element_text(size = 7),
        axis.title = element_text(size = 18), 
        axis.text = element_text(size = 16),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"))

#ggsave(paste(img.dir, 'thresholded_samples_ga.png', sep=''), dpi=300, height = 4, width = 4)

ggplot(data=kde_samples_nothresh,  aes(x=X0, y=X1)) + xlab('x') + ylab('y') +
  stat_density2d(geom="tile", aes(fill=..density..^0.25), contour=FALSE) +
  xlim(c(-1.4, 2.1)) + ylim(c(-1.4, 2.1)) +
  scale_fill_gradientn(colours = colorRampPalette(c("white", blues9))(256), 
                       guide = guide_colorbar(title = 'Density', 
                                              barheight = unit(3, "mm"),
                                              barwidth = unit(30, "mm"),
                                              direction = "horizontal", 
                                              title.position = "top")) +
  theme(text = element_text(size=14),
        legend.box.background = element_rect(color = "black"),
        legend.title=element_blank(), 
        axis.title.y = element_text(angle = 0, vjust=0.5),
        legend.position = c(0.74, 0.05), 
        legend.margin = margin(0, 0, 0, 0),
        legend.text = element_text(size = 7),
        axis.title = element_text(size = 18), 
        axis.text = element_text(size = 16),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"))

#ggsave(paste(img.dir, 'nonthresholded_samples_ga.png', sep=''), dpi=300, height = 4, width = 4)
