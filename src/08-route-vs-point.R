####### Script Information ########################
# Brandon P.M. Edwards
# BBS Point Level
# <08-route-vs-point.R>
# Created April 2024
# Last Updated August 2024

####### Import Libraries and External Files #######

library(bbsBayes2)
library(ggpubr)
library(cmdstanr)
library(bayesplot)
theme_set(theme_pubclean())
bayesplot::color_scheme_set("red")

####### Set Constants #############################

sp <- "OVEN"
mid_year <- 2016

####### Read Data #################################

route_indices <- readRDS(file = paste0("output/indices/", sp, "_route.RDS"))
point_indices <- readRDS(file = paste0("output/indices/", sp, "_point.RDS"))

####### Compare Indices ###########################

regions_both <- intersect(route_indices$indices$region, point_indices$indices$region)

route_df <- route_indices$indices[which(route_indices$indices$region %in% regions_both &
                                          route_indices$indices$region_type == "stratum"), ]
route_df$yr_region <- paste0(route_df$year, "-", route_df$region)

point_df <- point_indices$indices[which(point_indices$indices$region %in% regions_both &
                                          point_indices$indices$region_type == "stratum"), ]
point_df$yr_region <- paste0(point_df$year, "-", point_df$region)

# Make sure order of indices align
region_order <- match(route_df$yr_region, point_df$yr_region)
point_df <- point_df[region_order, ]

route_df_mid <- route_df[which(route_df$year == mid_year), ]
point_df_mid <- point_df[which(point_df$year == mid_year), ]

index_model_data <- list(N = nrow(route_df_mid),
                   y = route_df_mid$index,
                   x = point_df_mid$index,
                   beta_mean_prior = 50)

index_comp_model <- cmdstan_model(stan_file = "models/slr-model.stan")

index_comp_model_run <- index_comp_model$sample(
  data = index_model_data,
  iter_warmup = 1000,
  iter_sampling = 2000,
  chains = 4,
  parallel_chains = 4,
  refresh = 10
)

index_mod_summary <- index_comp_model_run$summary()

index_comp_model_draws <- index_comp_model_run$draws(variables = c("intercept", "beta"), format = "df")

to_plot <- data.frame(point = index_model_data$x,
                      route = index_model_data$y)
indices_comp_plot <- ggplot(data = to_plot, aes(x = point, y = route)) + 
  geom_point(alpha = 0.3) +
  geom_abline(intercept = index_comp_model_draws$intercept, slope = index_comp_model_draws$beta, color = "grey", alpha = 0.1) +
  geom_abline(intercept = mean(index_comp_model_draws$intercept),
              slope = mean(index_comp_model_draws$beta),
              color = "black", size = 1) +
  geom_abline(slope = 50, intercept = 0, color = "red", linetype = "dashed") +
  xlab("Index of Abundance (POINT)") +
  ylab("Index of Abundance (ROUTE)") +
  NULL

####### Compare Trends ############################

route_indices$indices <- route_indices$indices[which(route_indices$indices$region %in% regions_both), ]
point_indices$indices <- point_indices$indices[which(point_indices$indices$region %in% regions_both), ]

route_trends <- generate_trends(route_indices)
point_trends <- generate_trends(point_indices)

n_trends <- nrow(route_trends$trends)
# n_trends - 1 so that we don't model the continental trend, only stratum-level
trend_model_data <- list(N = n_trends - 1,
                         y = route_trends$trends$trend[2:n_trends],
                         x = point_trends$trends$trend[2:n_trends],
                         beta_mean_prior = 1)

trend_comp_model <- cmdstan_model(stan_file = "models/slr-model.stan")

trend_comp_model_run <- trend_comp_model$sample(
  data = trend_model_data,
  iter_warmup = 1000,
  iter_sampling = 2000,
  chains = 4,
  parallel_chains = 4,
  refresh = 10
)

trend_mod_summary <- trend_comp_model_run$summary()

trend_comp_model_draws <- trend_comp_model_run$draws(variables = c("intercept", "beta"), format = "df")

to_plot <- data.frame(point = trend_model_data$x,
                      route = trend_model_data$y)
trend_comp_plot <- ggplot(data = to_plot, aes(x = point, y = route)) + 
  geom_point(alpha = 0.3) +
  geom_abline(intercept = trend_comp_model_draws$intercept, slope = trend_comp_model_draws$beta, color = "grey", alpha = 0.1) +
  geom_abline(intercept = mean(trend_comp_model_draws$intercept),
              slope = mean(trend_comp_model_draws$beta),
              color = "black", size = 1) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  xlab("Trend (POINT)") +
  ylab("Trend (ROUTE)") +
  NULL

####### Output ####################################

write.table(file = "output/index_comp_model.csv", x = index_mod_summary, sep = ",", row.names = FALSE)
write.table(file = "output/trend_comp_model.csv", x = trend_mod_summary, sep = ",", row.names = FALSE)

png(filename = "output/plots/route-vs-point.png",
    width = 6, height = 3, res = 300, units = "in")
ggarrange(indices_comp_plot, trend_comp_plot, nrow = 1,
          labels = c("A", "B"))
dev.off()
