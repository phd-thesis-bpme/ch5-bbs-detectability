####### Script Information ########################
# Brandon P.M. Edwards
# BBS Point Level
# 09-varprop-diagnostics.R
# Created May 2024
# Last Updated August 2024

####### Import Libraries and External Files #######

library(cmdstanr)
library(bayesplot)
library(ggpubr)
theme_set(theme_pubclean())
bayesplot::color_scheme_set("gray")

####### Set Constants #############################

sp <- "OVEN"
mid_year <- 2016

####### Read Data #################################

detect_model <- readRDS(paste0("output/model_runs/", sp, "-detectability.rds"))
varprop_model <- readRDS(paste0("output/model_runs/", sp, "-varprop.rds"))

detect_indices <- readRDS(file = paste0("output/indices/", sp, "_detectability.RDS"))
varprop_indices <- readRDS(file = paste0("output/indices/", sp, "_varprop.RDS"))

####### Detect vs Varprop #########################

regions_both <- intersect(detect_indices$indices$region, varprop_indices$indices$region)

detect_df <- detect_indices$indices[which(detect_indices$indices$region %in% regions_both &
                                          detect_indices$indices$region_type == "stratum"), ]
detect_df$yr_region <- paste0(detect_df$year, "-", detect_df$region)

varprop_df <- varprop_indices$indices[which(varprop_indices$indices$region %in% regions_both &
                                              varprop_indices$indices$region_type == "stratum"), ]
varprop_df$yr_region <- paste0(varprop_df$year, "-", varprop_df$region)

# Make sure order of indices align
region_order <- match(detect_df$yr_region, varprop_df$yr_region)
varprop_df <- varprop_df[region_order, ]

detect_df_mid <- detect_df[which(detect_df$year == mid_year), ]
varprop_df_mid <- varprop_df[which(varprop_df$year == mid_year), ]

model_data <- list(N = nrow(detect_df_mid),
                   n_strata = length(unique(detect_df_mid$region)),
                   detect_index = detect_df_mid$index ,
                   varprop_index = varprop_df_mid$index,
                   stratum = as.numeric(factor(detect_df_mid$region)))

model <- cmdstan_model(stan_file = "models/detect-vs-varprop.stan")

model_run <- model$sample(
  data = model_data,
  iter_warmup = 1000,
  iter_sampling = 2000,
  chains = 4,
  parallel_chains = 4,
  refresh = 10
)

mod_summary <- model_run$summary()

model_draws <- model_run$draws(variables = c("intercept", "BETA"), format = "df")

to_plot <- data.frame(detect = model_data$detect_index,
                      varprop = model_data$varprop_index)
comp_plot <- ggplot(data = to_plot, aes(x = detect, y = varprop)) + 
  geom_abline(intercept = model_draws$intercept, slope = model_draws$BETA, color = "grey", alpha = 0.1) +
  geom_abline(intercept = mean(model_draws$intercept),
              slope = mean(model_draws$BETA),
              color = "black", size = 1) +
  geom_point(alpha = 0.1) +
  xlab("Index of Abundance (DETECT)") +
  ylab("Index of Abundance (VARPROP)") +
  NULL

####### Varprop Parameters ########################

varprop_summary <- varprop_model$model_fit$summary(variables = c("zeta", "xi"))

avail_plot <- bayesplot::mcmc_areas(varprop_model$model_fit$draws("zeta")) +
  scale_y_discrete(labels = c("Intercept", "Ordinal Day", "(Ordinal Day)^2"))

percept_plot <- bayesplot::mcmc_areas(varprop_model$model_fit$draws("xi")) +
  scale_y_discrete(labels = c("Intercept","Forest","Road","RoadForest"))

####### Output ####################################

write.table(x = mod_summary, file = "output/varprop_comparison_mod_mid.csv", sep = ",", row.names = FALSE)
write.table(x = varprop_summary, file = "output/varprop_summary.csv", sep = ",", row.names = FALSE)

png(filename = paste0("output/plots/varprop-diagnostics_mid.png"),
    width = 6, height = 4, res = 600, units = "in")
ggarrange(comp_plot, ggarrange(avail_plot, percept_plot, nrow = 2,
                               labels = c("B", "C")), nrow = 1,
          labels = "A")
dev.off()
