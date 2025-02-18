rm(list = ls())
library(TDA)
library(dplyr)
library(foreach)
source('utils.R')
library(doParallel)
theme_set(theme_bw())
library(scatterplot3d)
img.dir = '../images/'

# Read in the data
true_manifold = read.csv('../data/varyingdatadistribution/threshold_relaxation_true_manifold.csv') %>% select(V1, V2, V3)
noisy_manifold = read.csv('../data/varyingdatadistribution/threshold_relaxation_noisy_manifold.csv') %>% select(V1, V2, V3)
kde_samples_thresh = read.csv('../data/varyingdatadistribution/threshold_relaxation_kde_samples_thresh.csv') %>% select(X0, X1, X2)
kde_samples_nothresh = read.csv('../data/varyingdatadistribution/threshold_relaxation_kde_samples_nothresh.csv') %>% select(X0, X1, X2)


noisy_manifold_plt  = prcomp(noisy_manifold)$x
par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
scatterplot3d(noisy_manifold_plt[,1:3], pch = 16,grid=TRUE, box=F, cex.symbols = .75, lab.z = 5,
              angle = -153, xlab = '', ylab = 'y', zlab = 'z', cex.lab = 1, cex.axis = .7)
mtext("x", side = 1, line = -1.25, at = -8, cex = 1)

# We vary the kde bandwidth and compute the max persistence
max_pers <- function(grid_diagram){
  dgmat = as.matrix(grid_diagram)
  dg = data.frame(dimension = dgmat[, 1], birth=dgmat[, 3], death=dgmat[, 2]) %>% filter(dimension==1)
  return(round(max(dg$birth - dg$death), 2))
}

vec_max_pers <- function(bw, file_idx, grd = 0.05){
  mx_pers = c()
  for (k in 0:10){
    file_path = paste('../data/varyingdatadistribution/', file_idx, '.csv', sep='')
    data = read.csv(file_path)[, 2:4]
    grid_diag = gridDiag(X = as.matrix(data), FUN = kde, h=bw, lim=cbind( c(min(data[, 1]), max(data[, 1])), 
                                              c(min(data[, 2]), max(data[, 2])), c(min(data[, 3]), max(data[, 3]))), 
                         by = grd, sublevel = F, location = TRUE, library = 'Dionysus')[['diagram']]
    mx_pers = c(mx_pers, max_pers(grid_diag))
  }
  return(mx_pers)
}
#Uncomment the below block of code to generate mps at various bandwidths...it takes a while to run
###-------
# data_pers_vals = data.frame(TM = numeric(), NM = numeric(), NT = numeric(), TH = numeric(), BW = character(), GRD = numeric(), stringsAsFactors = FALSE)
# bdw = seq(0.005, 0.15, by=0.005)
# for (h in bdw){
#     true_vec = vec_max_pers(h, 'threshold_relaxation_true_manifold', 0.05)
#     noisy_vec = vec_max_pers(h, 'threshold_relaxation_noisy_manifold', 0.05)
#     kde_thresh_vec = vec_max_pers(h, 'threshold_relaxation_kde_samples_thresh', 0.05)
#     kde_nothresh_vec = vec_max_pers(h, 'threshold_relaxation_kde_samples_nothresh', 0.05)
#     vec = data.frame(TM = true_vec, NM = noisy_vec, NT = kde_nothresh_vec, TH = kde_thresh_vec, BW = as.character(h), GRD = 0.05)
#     data_pers_vals = rbind(data_pers_vals, vec)
# }
# write.csv(data_pers_vals, '../data/varyingdatadistribution/threshold_relaxation_grid_values.csv')
###------
data_pers_vals = read.csv("../data/varyingdatadistribution/threshold_relaxation_grid_values.csv") %>% select(-X)

data_long <- data_pers_vals %>% pivot_longer(cols = c(NM, TH, NT), names_to = "Manifold", values_to = "Value") %>% mutate(Value = Value - TM)
data_long_ <- data_long %>% group_by(BW, Manifold) %>% summarize(Value = mean(Value), .groups = 'drop')
manifold_labels <- c( NM = expression(KDE(bold(X)[n])), TH = expression(KDE({bold(X)^{'*'}}["n,"][lambda])), NT = expression(KDE({bold(X)^{'*'}}[n][",0"])))

ggplot(data = data_long_, aes(x = as.factor(BW), y = Value, color = Manifold)) +
  geom_point(position = position_dodge(width = 0.75), size=1) + 
  geom_line(aes(group = Manifold, linetype=Manifold), position = position_dodge(width = 0.75),) + 
  geom_hline(yintercept = 0, linetype='dashed', linewidth=0.1) +
  xlab("Bandwidth") + ylab("Maximum Persistence Difference") +
  scale_color_manual(values = c("NM" = "chartreuse3", "TH" = "deeppink", "NT" = "cornflowerblue"), labels = manifold_labels) +
  scale_linetype_manual(values = c("NM" = 1, "TH" = 2, "NT" = 3), labels = manifold_labels) +
  scale_x_discrete(labels = function(x) sprintf("%.4f", as.numeric(as.character(x)))) +
  theme(text = element_text(size=9),
        legend.box.background = element_rect(color = "black"),
        legend.title=element_blank(),
        axis.text = element_text(size=11),
        axis.text.x = element_text(angle = 90),
        legend.position = c(0.85, 0.22), 
        legend.margin = margin(0, 0, 0, 0),
        legend.text = element_text(size = 10),
        legend.spacing.x = unit(1, "mm"),
        legend.spacing.y = unit(1, "mm"))
ggsave(paste(img.dir, 'threshold_relxation_max_pers_diff.pdf', sep=''), dpi=300, height = 2.7, width = 6)



# We now construct the persistence diagrams with the optimal parameters
by = 0.05
kde_thresh_lim = cbind( 
  c(min(kde_samples_thresh$X0), max(kde_samples_thresh$X0)),
  c(min(kde_samples_thresh$X1), max(kde_samples_thresh$X1)),
  c(min(kde_samples_thresh$X2), max(kde_samples_thresh$X2)) )
kde_thresh = gridDiag(X = as.matrix(kde_samples_thresh), FUN = kde, h=0.02, maxdimension = 1, lim=kde_thresh_lim,
                      by = 0.05, sublevel = F, location = TRUE, library = 'Dionysus', printProgress=TRUE
)[['diagram']]
kde_diagram_smooth_thresh = data.frame(dimension=kde_thresh[, 1], birth=kde_thresh[, 3], death=kde_thresh[, 2]) 

kde_noisy_lim = cbind( 
  c(min(noisy_manifold$V1), max(noisy_manifold$V1)),
  c(min(noisy_manifold$V2), max(noisy_manifold$V2)),
  c(min(noisy_manifold$V3), max(noisy_manifold$V3)) )
kde_noisy = gridDiag(X = as.matrix(noisy_manifold), FUN = kde, h=0.015, maxdimension = 1, lim=kde_noisy_lim,
                     by = 0.05, sublevel = F, location = TRUE, library = 'Dionysus', printProgress=TRUE
)[['diagram']]
kde_diagram_smooth_noisy = data.frame(dimension=kde_noisy[, 1], birth=kde_noisy[, 3], death=kde_noisy[, 2])


# Finally, we construct confidence bands --- we used parallel computations to speed up the process

# Uncomment this block of code to compute the bottleneck distances
###----------------------------
# num_cores <- detectCores()  # Leave one core free for system processes
# cl <- makeCluster(num_cores)
# registerDoParallel(cl)
# 
# # Create empty vectors to store results
# n_iterations <- 1000
# results <- foreach(k = 0:(n_iterations-1), .verbose = TRUE, .combine = rbind, .packages = c('dplyr', 'TDA')) %dopar% {
#   # Read files
#   noisy_manifold_ <- read.csv(paste("../data/varyingdistribution/temp_bt_mfld/bootstrapped_noisy_manifold", k, '.csv', sep = '')) %>% select(X0, X1, X2)
#   kde_samples_thresh_ <- read.csv(paste("../data/varyingdatadistribution/temp_bt_kde_samples/bootstrapped_thresh_kde_samples", k, '.csv', sep = '')) %>% select(X0, X1, X2)
#   # KDE persistence diagram threshold
#   kde_thresh_lim_ <- cbind(c(min(kde_samples_thresh_$X0), max(kde_samples_thresh_$X0)), c(min(kde_samples_thresh_$X1), max(kde_samples_thresh_$X1)),
#                            c(min(kde_samples_thresh_$X2), max(kde_samples_thresh_$X2)) )
#   kde_thresh_ <- gridDiag(X = as.matrix(kde_samples_thresh_), FUN = kde, h = 0.02, maxdimension = 1, lim = kde_thresh_lim_,
#                           by = by, sublevel = F, location = TRUE, library = 'Dionysus', printProgress = F )[['diagram']]
#   btl_dist_kde <- bottleneck(kde_thresh_, kde_thresh)
#   # Noisy manifold persistence diagram threshold
#   kde_noisy_lim_ <- cbind( c(min(noisy_manifold_$X0), max(noisy_manifold_$X0)), c(min(noisy_manifold_$X1), max(noisy_manifold_$X1)),
#                            c(min(noisy_manifold_$X2), max(noisy_manifold_$X2)))
#   kde_noisy_ <- gridDiag(X = as.matrix(noisy_manifold_), FUN = kde, h = 0.015, maxdimension = 1, lim = kde_noisy_lim_,
#                          by = by, sublevel = F, location = TRUE, library = 'Dionysus', printProgress = F )[['diagram']]
#   btl_dist_noisy <- bottleneck(kde_noisy_, kde_noisy)
#   # Return results for this iteration
#   c(k=k, btl_dist_kde = btl_dist_kde, btl_dist_noisy = btl_dist_noisy)
# }
# # Stop the cluster
# stopCluster(cl)
# # Extract results into separate vectors
# results_df <- as.data.frame(results)
# colnames(results_df) <- c("k", "btl_dist_kde", "btl_dist_noisy")
# write.csv(results_df, '../data/varyingdatadistribution/cc_btl_dist_kde_noisy_1000.csv')
###----------------------------

cc_df = read.csv('../data/varyingdatadistribution/cc_btl_dist_kde_noisy_1000.csv')
btl_dist_kde = cc_df$btl_dist_kde
btl_dist_noisy = cc_df$btl_dist_noisy

(cc_kde_thresh_ = quantile(btl_dist_kde, 0.95)*(2))
(cc_kde_noisy_ = quantile(btl_dist_noisy, 0.95)*(2))

plot_band(kde_diagram_smooth_thresh, c_hat=cc_kde_thresh_, dim=1, rb_scale=50, x_lim = 15.5, y_lim = 15.5, sublevel=F)
#ggsave(paste(img.dir, 'kde_pd_conf_threshold_relaxation_thresholded2.pdf', sep=''), dpi=300, height = 3, width = 3)
plot_band(kde_diagram_smooth_noisy, c_hat=cc_kde_noisy_, dim=1, rb_scale=50, x_lim = 15.5, y_lim = 15.5, sublevel=F)
#ggsave(paste(img.dir, 'kde_pd_conf_threshold_relaxation_noisy2.pdf', sep=''), dpi=300, height = 3, width = 3)


