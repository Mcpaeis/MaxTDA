rm(list = ls())
library(ggtda)
library(ripserr)
library(patchwork)
library(tidyverse)
theme_set(theme_bw())
img.dir = '../images/'

# Read in the data
DMatrix_ <- read.csv('../data/graphicalabstract/toy-data.csv') %>% select(-X)
df <- DMatrix_[12:21, ]

# Plot the point cloud
ggplot(data =df,  aes(x = V1, y = V2, linewidth = 1.5)) + 
  geom_point() + 
  xlim(c(.75, 4.25)) + ylim(c(.75, 4.25)) +
  coord_fixed(ratio = 1)+
  theme_bw() +
  xlab("X") + ylab("Y") + 
  theme(legend.position  = "none", text = element_text(size=20), axis.title.y = element_text(angle = 0, vjust = 0.5))
#ggsave(paste(img.dir, 'toy-point-cloud.pdf', sep=''), width = 3, height = 3, dpi = 300)

# Make the two complexes
prox = 0.8
ggplot(data =df,  aes(x = V1, y = V2, linewidth = 1.5)) + 
  geom_point() + 
  xlim(c(.75, 4.25)) + ylim(c(.75, 4.25)) +
  coord_fixed(ratio = 1)+
  theme_bw() +
  stat_disk(radius = prox/2, fill = "aquamarine3") +
  stat_vietoris2(diameter = prox, fill = "darkgoldenrod") +
  stat_vietoris1(diameter = prox, alpha = .3) +
  stat_vietoris0() +
  xlab("X") + ylab("Y") + 
  theme(legend.position  = "none", text = element_text(size=20), axis.title.y = element_text(angle = 0, vjust = 0.5))
#ggsave(paste(img.dir, 'point-cloud-delta08.pdf', sep=''), width = 3, height = 3, dpi = 300)

prox = 1.5
ggplot(data =df,  aes(x = V1, y = V2, linewidth = 1.5)) + 
  geom_point() + 
  xlim(c(.75, 4.25)) + ylim(c(.75, 4.25)) +
  coord_fixed(ratio = 1) +
  theme_bw() +
  stat_disk(radius = prox/2, fill = "aquamarine3") +
  stat_vietoris2(diameter = prox, fill = "darkgoldenrod") +
  stat_vietoris1(diameter = prox, alpha = .3) +
  stat_vietoris0() +
  xlab("X") + ylab("Y") + 
  theme(legend.position  = "none", text = element_text(size=20), axis.title.y = element_text(angle = 0, vjust = 0.5))
#ggsave(paste(img.dir, 'point-cloud-delta15.pdf', sep=''), width = 3, height = 3, dpi = 300)

# Plot the persistence diagram
pers_hom <- vietoris_rips(df, dim=1)
pers_hom$ind = as.factor(5*pers_hom$dimension)

breaks <- pretty(pers_hom$death)  # Use the pretty function to get usual tick marks
breaks <- c(breaks, max(breaks) +  breaks[2])  # Ensure the maximum value is included
# Compute the labels
labels <- as.character(breaks)
labels[labels == as.character(max(breaks))] <- expression(infinity) 

ggplot(pers_hom) + 
  geom_point(aes(x = birth, y = death, shape = ind, color = ind), size = 2.0) + 
  geom_hline(yintercept = max(breaks), linetype = "dashed", color = "red") +
  geom_point(aes(x = 0, y = max(breaks)), color = "red") +
  xlim(c(0, max(breaks))) +
  scale_y_continuous(breaks = breaks, labels =  labels) +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Birth") + ylab("Death") + #xlim(c(0, 1.5)) + ylim(c(0, 2)) +
  scale_shape_manual(name = "",
                     labels = c(expression(H[0]), expression(H[1])),
                     values = c(16, 17)) + 
  scale_colour_manual(name="",
                      labels=c(expression(H[0]), expression(H[1])),
                      values=c("Red", "#619CFF")) +
  theme(text = element_text(size=20),
        legend.box.background = element_rect(color = "black"),
        legend.title=element_blank(), 
        legend.position = c(0.8, 0.2), 
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"))
#ggsave(paste(img.dir, 'toy-persistence-diagram.pdf', sep = ''), width = 3, height = 3, dpi = 300)




