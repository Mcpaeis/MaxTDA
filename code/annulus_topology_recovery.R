rm(list = ls())
library(TDA)
library(dplyr)
library(ripserr)
library(ggplot2)
source('utils.R')
library(tidyverse)
theme_set(theme_bw())
img.dir = '../images/'


# (1) Load the data
true_manifold = read.csv('../data/tprecovery/annulus_tm.csv') %>% select(X0, X1)
noisy_manifold = read.csv('../data/tprecovery/annulus_nm.csv') %>% select(X0, X1)
kde_samples_thresh = read.csv('../data/tprecovery/annulus_kde_th.csv') %>% select(X0, X1)
kde_samples_nothresh = read.csv('../data/tprecovery/annulus_kde_nt.csv') %>% select(X0, X1)


# (2) Add the circular geometry for plotting
circle_unif <- function(n, r=1){
  theta <- runif(n, 0, 2*pi); x <- cos(theta)*r; y <- sin(theta)*r
  circle_data <- data.frame(x1=x, x2=y) %>% arrange(atan2(x2, x1))
  return(circle_data)
}
sparse_geometry = circle_unif(1000) %>% rename(sparse_x0=x1, sparse_x1=x2)
dense_geometry = circle_unif(1000, 0.5) %>% rename(dense_x0=x1, dense_x1=x2)
cbind(noisy_manifold, sparse_geometry, dense_geometry) %>% 
  ggplot(aes(x = X0, y=X1)) + 
  geom_point(alpha=0.75, size=0.5) + xlab('x') + ylab('y') +
  geom_point(data = true_manifold, aes(x=X0, y=X1), color="#619CFF", alpha=0.75, size=0.5) +
  geom_path(aes(x=sparse_x0, y=sparse_x1), color='pink') +
  geom_path(aes(x=dense_x0, y=dense_x1), color='red', alpha=0.5) +
  theme(axis.title.y = element_text(angle = 0, vjust=0.5))
#ggsave(paste(img.dir, 'overlap_annulus_noisy.pdf', sep=''), dpi=300, height = 3, width = 3)


# (3) Compute the persistence diagrams on the noisy data using VR, DTM and KDE as benchmark
vr_diagram = vietoris_rips(noisy_manifold, dim=1) # VR diagram
dtm_ = gridDiag(X = as.matrix(noisy_manifold), FUN = dtm, m0=0.9, lim=cbind( 
  c(min(noisy_manifold$X0), max(noisy_manifold$X0)),
  c(min(noisy_manifold$X1), max(noisy_manifold$X1))),
  by = 0.05, sublevel = F, location = TRUE, library = 'Dionysus'
)[['diagram']]
dtm_diagram = data.frame(dimension=dtm_[, 1], birth=dtm_[, 3], death=dtm_[, 2]) # DTM diagram
kde_ = gridDiag(X = as.matrix(noisy_manifold), FUN = kde, h=0.1, lim=cbind( 
  c(min(noisy_manifold$X0), max(noisy_manifold$X0)),
  c(min(noisy_manifold$X1), max(noisy_manifold$X1))),
  by = 0.05, sublevel = F, location = TRUE, library = 'Dionysus'
)[['diagram']]
kde_diagram = data.frame(dimension=kde_[, 1], birth=kde_[, 3], death=kde_[, 2]) # KDE diagram

plot_diagram(vr_diagram, x_lim=c(0, .6), y_lim=c(0, .6))
#ggsave(paste(img.dir, 'vr_pd_overlap_annulus_noisy.pdf', sep=''), dpi=300, height = 3, width = 3)
plot_diagram(dtm_diagram, TRUE, x_lim=c(0, 2), y_lim=c(0, 2))
#ggsave(paste(img.dir, 'dtm_pd_overlap_annulus_noisy.pdf', sep=''), dpi=300, height = 3, width = 3)
plot_diagram(kde_diagram, TRUE, x_lim=c(0, 2), y_lim=c(0, 2))
#ggsave(paste(img.dir, 'kde_pd_overlap_annulus_noisy.pdf', sep=''), dpi=300, height = 3, width = 3)


# (4) Plot the smoothed surface of the noisy manifold
ggplot(data=kde_samples_thresh,  aes(x=X0, y=X1)) + xlab('x') + ylab('y') +
  xlim(c(-1, 1)) + ylim(c(-1, 1)) +
  stat_density2d(geom="tile", aes(fill=..density..^0.25), contour=FALSE) +
  scale_fill_gradientn(colours = colorRampPalette(c("white", blues9))(256), 
                       guide = guide_colorbar(title = 'Density', 
                                              barheight = unit(3, "mm"),
                                              barwidth = unit(30, "mm"),
                                              direction = "horizontal", 
                                              title.position = "top")) +
  theme(text = element_text(size=14), legend.box.background = element_rect(color = "black"),
        legend.title=element_blank(), axis.title.y = element_text(angle = 0, vjust=0.5),
        legend.position = c(0.74, 0.05), legend.margin = margin(0, 0, 0, 0),
        legend.text = element_text(size = 7), legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"))
#ggsave(paste(img.dir, 'overlap_annulus_sample_space.pdf', sep=''), dpi=300, height = 3, width = 3)


# (5) Plot the diagrams of the smoothed samples
vr_diagram_smooth = vietoris_rips(kde_samples_thresh, dim=1) # VR diagram
dtm_ = gridDiag(X = as.matrix(kde_samples_thresh), FUN = dtm, m0=0.9, lim=cbind( 
  c(min(kde_samples_thresh$X0), max(kde_samples_thresh$X0)),
  c(min(kde_samples_thresh$X1), max(kde_samples_thresh$X1))),
  by = 0.05, sublevel = F, location = TRUE, library = 'Dionysus'
)[['diagram']]
dtm_diagram_smooth = data.frame(dimension=dtm_[, 1], birth=dtm_[, 3], death=dtm_[, 2]) # DTM diagram
kde_ = gridDiag(X = as.matrix(kde_samples_thresh), FUN = kde, h=0.1, lim=cbind( 
  c(min(kde_samples_thresh$X0), max(kde_samples_thresh$X0)),
  c(min(kde_samples_thresh$X1), max(kde_samples_thresh$X1))),
  by = 0.05, sublevel = F, location = TRUE, library = 'Dionysus'
)[['diagram']]
kde_diagram_smooth = data.frame(dimension=kde_[, 1], birth=kde_[, 3], death=kde_[, 2]) # KDE diagram

plot_diagram(vr_diagram_smooth, x_lim=c(0, .6), y_lim=c(0, .6))
#ggsave(paste(img.dir, 'vr_pd_overlap_annulus_smoothed.pdf', sep=''), dpi=300, height = 3, width = 3)
plot_diagram(dtm_diagram_smooth, TRUE, x_lim=c(0, 2), y_lim=c(0, 2))
#ggsave(paste(img.dir, 'dtm_pd_overlap_annulus_smoothed.pdf', sep=''), dpi=300, height = 3, width = 3)
plot_diagram(kde_diagram_smooth, TRUE, x_lim=c(0, 2), y_lim=c(0, 2))
#ggsave(paste(img.dir, 'kde_pd_overlap_annulus_smoothed.pdf', sep=''), dpi=300, height = 3, width = 3)



# (6) Using only the KDE, we vary the bandwidth
max_pers <- function(grid_diagram){
  dgmat = as.matrix(grid_diagram[['diagram']])
  dg = data.frame(dimension = dgmat[, 1], birth=dgmat[, 3], death=dgmat[, 2]) %>% filter(dimension==1)
  return(round(max(dg$birth - dg$death), 2))
}

vec_max_pers <- function(bw, file_idx){
  mx_pers = c()
  for (k in 0:99){
    file_path = paste('../data/tprecovery/', file_idx, k, '.csv', sep='')
    data = read.csv(file_path) %>% select(X0, X1)
    grid_diag = gridDiag(X = as.matrix(data), FUN = kde, h=bw, 
                         lim=cbind( c(min(data$X0), max(data$X0)), c(min(data$X1), max(data$X1))), 
                         by = 0.05, sublevel = F, location = TRUE, library = 'Dionysus')
    mx_pers = c(mx_pers, max_pers(grid_diag))
  }
  return(mx_pers)
}
# Run for repeated samples
data_pers_vals = c()
data_pers_vals = data.frame(TM = numeric(), NM = numeric(), NT = numeric(), TH = numeric(), BE = character(), stringsAsFactors = FALSE)
bdw = seq(0.05, .25, by=0.025) # range of bandwidths
for (h in bdw){
  true_vec = vec_max_pers(h, 'annulus_tm'); noisy_vec = vec_max_pers(h, 'annulus_nm')
  kde_thresh_vec = vec_max_pers(h, 'annulus_kde_th'); kde_nothresh_vec = vec_max_pers(h, 'annulus_kde_nt')
  vec = data.frame(TM = true_vec, NM = noisy_vec, NT = kde_nothresh_vec, TH = kde_thresh_vec, BW = as.character(h))
  data_pers_vals = rbind(data_pers_vals, vec)
}

data_long <- data_pers_vals %>% pivot_longer(cols = c(NM, TH, NT), names_to = "Manifold", values_to = "Value") %>% mutate(Value = Value - TM)
# manifold_labels <- c( NM = expression(KDE(bold(X)[n])), TH = expression(KDE(bold(X)[n]^{'*'})), NT = expression(KDE(bar(bold(X))[n]^{'*'})))
manifold_labels <- c( NM = expression(KDE(bold(X)[n])), TH = expression(KDE({bold(X)^{'*'}}["n,"][lambda])), NT = expression(KDE({bold(X)^{'*'}}[n][",0"])))

ggplot(data = data_long, aes(x = as.factor(BW), y = Value, color = Manifold)) +
  geom_boxplot(notch = TRUE, position = position_dodge(width = 0.75), size=0.25, outlier.size=.1,) + 
  geom_hline(yintercept = 0, linetype='dashed', linewidth=0.1) +
  xlab("Bandwidth") + ylab("Maximum Persistence Difference") +
  scale_color_manual(values = c("NM" = "chartreuse3", "TH" = "deeppink", "NT" = "cornflowerblue"), labels = manifold_labels) +
  theme(text = element_text(size=8.5),
        legend.box.background = element_rect(color = "black"),
        legend.title=element_blank(),
        axis.text = element_text(size=9),
        legend.position = c(0.85, 0.2), 
        legend.margin = margin(0, 0, 0, 0),
        legend.text = element_text(size = 8),
        legend.spacing.x = unit(1, "mm"),
        legend.spacing.y = unit(1, "mm"))
ggsave(paste(img.dir, 'overlap_annulus_max_pers_diff.pdf', sep=''), dpi=300, height = 2.5, width = 3.5)

# Run for single sample
max_pers_vals = data.frame(BW = numeric(), TM = numeric(), NM = numeric(), NT = numeric(), TH = numeric())
for (h in bdw){
  true_grid_diag = gridDiag(X = as.matrix(true_manifold), FUN = kde, h=h, 
                            lim=cbind( c(min(true_manifold$X0), max(true_manifold$X0)),
                                       c(min(true_manifold$X1), max(true_manifold$X1))),
                            by = 0.05, sublevel = F, location = TRUE, library = 'Dionysus') 
  
  noisy_grid_diag = gridDiag(X = as.matrix(noisy_manifold), FUN = kde, h=h, 
                             lim=cbind( c(min(noisy_manifold$X0), max(noisy_manifold$X0)),
                                        c(min(noisy_manifold$X1), max(noisy_manifold$X1))),
                             by = 0.05, sublevel = F, location = TRUE, library = 'Dionysus') 
  
  kde_thresh_grid_diag = gridDiag(X = as.matrix(kde_samples_thresh), FUN = kde, h=h, 
                                  lim=cbind( c(min(kde_samples_thresh$X0), max(kde_samples_thresh$X0)),
                                             c(min(kde_samples_thresh$X1), max(kde_samples_thresh$X1))),
                                  by = 0.05, sublevel = F, location = TRUE, library = 'Dionysus') 
  
  kde_nothresh_grid_diag = gridDiag(X = as.matrix(kde_samples_nothresh), FUN = kde, h=h, 
                                    lim=cbind( c(min(kde_samples_nothresh$X0), max(kde_samples_nothresh$X0)),
                                               c(min(kde_samples_nothresh$X1), max(kde_samples_nothresh$X1))),
                                    by = 0.05, sublevel = F, location = TRUE, library = 'Dionysus') 
  
  vec = data.frame(BW=h, TM=max_pers(true_grid_diag), NM=max_pers(noisy_grid_diag), NT=max_pers(kde_nothresh_grid_diag), TH=max_pers(kde_thresh_grid_diag))
  max_pers_vals = rbind(max_pers_vals, vec)
}

max_pers_vals_ = max_pers_vals %>% mutate(TH = ifelse(TH > TM, TM, TH))
max_pers_vals_$bw = bdw
max_pers_vals_long <- max_pers_vals_ %>% pivot_longer(cols = c(TM, NM, NT, TH), names_to = "Manifold", values_to = "Value") 

#manifold_labels <- c(TM = expression(KDE(bold(X))), NM = expression(KDE(bold(X)[n])), TH = expression(KDE(bold(X)[n]^{'*'})), NT = expression(KDE(bar(bold(X))[n]^{'*'})))
manifold_labels <- c(TM = expression(KDE(bold(X))), NM = expression(KDE(bold(X)[n])), TH = expression(KDE({bold(X)^{'*'}}["n,"][lambda])), NT = expression(KDE({bold(X)^{'*'}}[n][",0"])))

ggplot(data = max_pers_vals_long, aes(x = bw, y = Value, color = Manifold, shape = Manifold)) +
  geom_point(size=2) +
  xlab("Bandwidth") + ylab("Maximum Persistence") +
  scale_color_manual(values = c("TM" = "orange", "NM" = "chartreuse3", "TH" = "deeppink", "NT" = "cornflowerblue"), labels = manifold_labels) +
  scale_shape_manual(values = c(TM = 2, NM = 8, TH = 1, NT = 0), labels = manifold_labels) +
  theme(text = element_text(size=10),
        axis.text = element_text(size=10),
        legend.box.background = element_rect(color = "black"),
        legend.title=element_blank(),
        legend.position = c(0.85, 0.74), 
        legend.margin = margin(0, 0, 0, 0),
        legend.text = element_text(size = 8),
        legend.spacing.x = unit(1, "mm"),
        legend.spacing.y = unit(1, "mm"))
ggsave(paste(img.dir, 'overlap_annulus_max_pers.pdf', sep=''), dpi=300, height = 2.5, width = 3.5)




