####### Script Information ########################
# Brandon P.M. Edwards
# BBS Point Level
# 12-figure-detectability-change.R
# Created August 2024
# Last Updated August 2024

####### Import Libraries and External Files #######

library(terra)
library(ggpubr)
library(tidyterra)
library(napops)
library(ggforce)
theme_set(theme_pubclean())

####### Set Constants #############################


####### Read Data #################################


####### Main Code #################################

# Create the maps
map_1 <- rast(ncol = 800, nrow = 800,
                 xmin = 0, xmax = 800,
                 ymin = 0, ymax = 800)
values(map_1) <- rbinom(ncell(map_1), size = 1, p = 0.1)

m1_plot <- ggplot() +
  geom_spatraster(data = map_1) +
  geom_circle(aes(x0 = 400, y0 = 400, r = 400), color = "darkred", inherit.aes = TRUE) +
  geom_point(aes(x = 400, y = 400), color = "darkred") +
  theme(legend.position = "none") +
  NULL


map_2 <- rast(ncol = 800, nrow = 800,
              xmin = 0, xmax = 800,
              ymin = 0, ymax = 800)
values(map_2) <- rbinom(ncell(map_2), size = 1, p = 0.5)

m2_plot <- ggplot() +
  geom_spatraster(data = map_2) +
  geom_circle(aes(x0 = 400, y0 = 400, r = 400), color = "darkred", inherit.aes = TRUE) +
  geom_point(aes(x = 400, y = 400), color = "darkred") +
  theme(legend.position = "none") +
  NULL


map_3 <- rast(ncol = 800, nrow = 800,
              xmin = 0, xmax = 800,
              ymin = 0, ymax = 800)
values(map_3) <- rbinom(ncell(map_3), size = 1, p = 0.9) 

m3_plot <- ggplot() +
  geom_spatraster(data = map_3) +
  geom_circle(aes(x0 = 400, y0 = 400, r = 400), color = "darkred", inherit.aes = TRUE) +
  geom_point(aes(x = 400, y = 400), color = "darkred") +
  theme(legend.position = "none") +
  NULL

# Now create the detectability plots
fc_values <- seq(0,1,by = 0.05)
distance_values <- c(rep(50, length(fc_values)),
                     rep(100, length(fc_values)),
                     rep(200, length(fc_values)),
                     rep(400, length(fc_values)))

q <- percept(species = "OVEN",
             model = "best",
             forest = rep(fc_values, 4),
             road = rep(TRUE, length(fc_values) * 4),
             pairwise = TRUE,
             distance = distance_values)
q$Distance <- factor(distance_values,
                     levels = c("50",
                                "100",
                                "200",
                                "400"))

detect_1 <- ggplot() +
  geom_line(data = q, aes(x = Forest, y = q, group = Distance, color = Distance)) +
  geom_vline(xintercept = 0.90) +
  ylim(0,1) +
  xlab("Forest Coverage") + 
  ylab("Perceptibility") +
  theme(legend.position = "none") +
  NULL

detect_2 <- ggplot() +
  geom_line(data = q, aes(x = Forest, y = q, group = Distance, color = Distance)) +
  geom_vline(xintercept = 0.50) +
  ylim(0,1) +
  xlab("Forest Coverage") + 
  ylab("Perceptibility") +
  theme(legend.position = "none") +
  NULL

detect_3 <- ggplot() +
  geom_line(data = q, aes(x = Forest, y = q, group = Distance, color = Distance)) +
  geom_vline(xintercept = 0.10) +
  ylim(0,1) +
  xlab("Forest Coverage") + 
  ylab("Perceptibility") +
  theme(legend.position = "none") +
  NULL

ggarrange(m1_plot, m2_plot, m3_plot,
          detect_1, detect_2, detect_3,
          ncol = 3, nrow = 2)
