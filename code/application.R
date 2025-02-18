rm(list = ls())
library(TDA)
library(ripserr)
library(tidyverse)
library(scatterplot3d)
library(nonlinearTseries)
theme_set(theme_bw())
source('utils.R')
img.dir = '../images/'
data.dir = '../data/application/'

# Read in the exoplanet data
exo_data = read.csv(paste(data.dir, "exoplanet_time_series.csv", sep = ''))

# Add random noise to the planet+spot rv signal
set.seed(100)
n_obs     = nrow(exo_data %>% filter(type=="Planet+Spot"))
rnd_noise = rnorm(n=n_obs, sd=.3)
epsilon   = runif(n_obs, min = -2, max = 2)
ps_clean  = exo_data[exo_data$type=='Planet+Spot', 'rv']
ps_noisy  = approx(1:n_obs, ps_clean, xout = (1:n_obs) + epsilon, method = "linear", rule = 2)$y
exo_data[exo_data$type=='Planet+Spot', 'rv'] = ps_noisy + rnd_noise

# Individually extract the various components
spot        = exo_data %>% filter(type=="Spot")
planet      = exo_data %>% filter(type=="Planet")
planet_spot = exo_data %>% filter(type=="Planet+Spot")
# Re-level for plotting
exo_data$type = factor(exo_data$type, levels = c("Planet", "Spot", "Planet+Spot"))

exo_data %>% ggplot(aes(x = time, y=rv)) + xlab('Time (Days)') + ylab('Radial Velocity (m/sec)') +
  geom_point(aes(color=type, shape = type)) + 
  theme(legend.position = "top", legend.margin = margin(b = -11),  
        legend.background = element_rect(fill = "transparent", color = NA)) + 
  scale_shape_manual(values = c(19, 17, 0)) +
  labs(color = NULL, shape=NULL) + geom_line(aes(color = type)) #
ggsave(paste(img.dir,  'exo-planet-signals.pdf', sep = ''), width = 6, height = 3, dpi = 300)

# Estimate the delay step using average mutual information
tau.ami_spot       = timeLag(spot$rv, selection.method = 'first.minimum', technique = "ami", lag.max = 100, do.plot = T)
tau.ami_planet     = timeLag(planet$rv, selection.method = 'first.minimum', technique = "ami", lag.max = 100, do.plot = T)
tau.ami_planetspot = timeLag(planet_spot$rv, selection.method = 'first.minimum', technique = "ami", lag.max = 100, do.plot = T)

# Estimate the minimum embedding dimensions
emb.dim_spot        = estimateEmbeddingDim(spot$rv, time.lag = tau.ami_spot, max.embedding.dim = 8, ylim = c(0, 2))
emb.dim_spot        = 7 # max consistent embedding dimension
emb.dim_planet      = estimateEmbeddingDim(planet$rv, time.lag = tau.ami_planet, max.embedding.dim = 23)
emb.dim_planet      = 16 # max consistent embedding dimension
emb.dim_planetspot  = estimateEmbeddingDim(planet_spot$rv, time.lag = tau.ami_planetspot, max.embedding.dim = 19)
emb.dim_planetspot  = 16 # max consistent embedding dimension close the desired signal-planet

# Construct the TDE embedding and do the PCA
n_comps = 2
TDE_spot        = tde(spot$rv, emb.dim_spot, tau.ami_spot, n_comps)
TDE_planet      = tde(planet$rv, emb.dim_planet, tau.ami_planet, n_comps)
TDE_planetspot1 = tde(planet_spot$rv, emb.dim_planet, tau.ami_planet, n_comps)
TDE_planetspot2 = tde(planet_spot$rv, emb.dim_spot, tau.ami_spot, n_comps)

# Merge for plotting
cd <- data.frame(rbind(TDE_planet, TDE_spot, TDE_planetspot1, TDE_planetspot2)) %>% 
  mutate(source = c(
    rep(c("Planet"), nrow(TDE_planet)),rep(c("Spot"), nrow(TDE_spot)),
    rep(c("Planet+Spot[Planet Parameters]"), nrow(TDE_planetspot1)),
    rep(c("Planet+Spot[Spot Parameters]"), nrow(TDE_planetspot2))
    ))
cd$source = factor(cd$source, levels = c("Planet", "Spot", "Planet+Spot[Planet Parameters]", "Planet+Spot[Spot Parameters]"))

ggplot(cd, aes(x = PC1, y = PC2)) +
  geom_point(color="darkorange") + labs(x = 'x', y = 'y') + 
  theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.1), "cm"), axis.title.y = element_text(angle = 360, vjust = 0.5)) +
  facet_wrap(~source, scale='free', ncol=4)
ggsave(paste(img.dir, 'four_rv_embeddings.pdf', sep=''), dpi=300, height = 2.75, width = 9)

# Write the three matrices to file
write.csv(TDE_spot, paste(data.dir, 'tde_spot.csv', sep = ''), row.names = FALSE)
write.csv(TDE_planet, paste(data.dir, 'tde_planet.csv', sep = ''), row.names = FALSE)
write.csv(TDE_planetspot1, paste(data.dir, 'tde_planetspot1.csv', sep = ''), row.names = FALSE)
write.csv(TDE_planetspot2, paste(data.dir, 'tde_planetspot2.csv', sep = ''), row.names = FALSE)
# </-- Run the application.ipynb file to get smoothed data and the dtm diagrams --/>

# Read in the pre-computed dtm diagrams
dgm_spot        = read.csv(paste(data.dir, 'spot_dgm.csv', sep = '')) %>% select(-X)
dgm_planet      = read.csv(paste(data.dir, 'planet_dgm.csv', sep = '')) %>% select(-X)
dgm_smooth1     = read.csv(paste(data.dir, 'smooth1_dgm.csv', sep = '')) %>% select(-X)
dgm_smooth2     = read.csv(paste(data.dir, 'smooth2_dgm.csv', sep = '')) %>% select(-X)
dgm_planetspot1 = read.csv(paste(data.dir, 'planetspot1_dgm.csv', sep = '')) %>% select(-X)
dgm_planetspot2 = read.csv(paste(data.dir, 'planetspot2_dgm.csv', sep = '')) %>% select(-X)

# Compute the quantile
alpha = 0.05
btlnk_dist       = read.csv(paste(data.dir, 'bottleneck_distances.csv', sep = '')) %>% select(-X)
cc_1_spot        = quantile(btlnk_dist$spot, 1-alpha)[[1]]
cc_1_planet      = quantile(btlnk_dist$planet, 1-alpha)[[1]]
cc_1_smooth1     = quantile(btlnk_dist$smooth1, 1-alpha)[[1]]
cc_1_smooth2     = quantile(btlnk_dist$smooth2, 1-alpha)[[1]]
cc_1_planetspot1 = quantile(btlnk_dist$planetspot1, 1-alpha)[[1]]
cc_1_planetspot2 = quantile(btlnk_dist$planetspot2, 1-alpha)[[1]]


# Plot the regular persistence diagram
h = 5; w=5
plot_band(dgm_planet, 2*cc_1_planet, dim=1, x_lim = 1, y_lim = 1, rb_scale = 20, sublevel=T)
ggsave(paste(img.dir, 'ssp_planet_h1.pdf', sep=''), dpi=300, height = h, width = w)

plot_band(dgm_planetspot1, 2*cc_1_planetspot1, dim=1, x_lim = 1, y_lim = 1, rb_scale = 20, sublevel=T)
ggsave(paste(img.dir, 'ssp_planetspot1_h1.pdf', sep=''), dpi=300, height = h, width = w)

plot_band(dgm_planetspot2, 2*cc_1_planetspot2, dim=1, x_lim = 1, y_lim = 1, rb_scale = 20, sublevel=T)
ggsave(paste(img.dir, 'ssp_planetspot2_h1.pdf', sep=''), dpi=300, height = h, width = w)

plot_band(dgm_spot, 2*cc_1_spot, dim=1, x_lim = 1, y_lim = 1, rb_scale = 15, sublevel=T)
ggsave(paste(img.dir, 'ssp_spot_h1.pdf', sep=''), dpi=300, height = h, width = w)

plot_band(dgm_smooth1, 2*cc_1_smooth1, dim=1, x_lim = 1, y_lim = 1, rb_scale = 35, sublevel=T)
ggsave(paste(img.dir, 'ssp_smooth1_h1.pdf', sep=''), dpi=300, height = h, width = w)

plot_band(dgm_smooth2, 2*cc_1_smooth2, dim=1, x_lim = 1, y_lim = 1, rb_scale = 35, sublevel=T)
ggsave(paste(img.dir, 'ssp_smooth2_h1.pdf', sep=''), dpi=300, height = h, width = w)


