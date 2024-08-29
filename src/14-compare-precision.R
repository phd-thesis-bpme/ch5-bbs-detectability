####### Script Information ########################
# Brandon P.M. Edwards
# BBS Point Level
# <14-compare-precision.R>
# Created August 2024
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

####### Read Data #################################

route_indices <- readRDS(file = paste0("output/indices/", sp, "_route.RDS"))
point_indices <- readRDS(file = paste0("output/indices/", sp, "_point.RDS"))
detect_indices <- readRDS(file = paste0("output/indices/", sp, "_detectability.RDS"))
varprop_indices <- readRDS(file = paste0("output/indices/", sp, "_varprop.RDS"))

####### Compare Indices ###########################

route_trends <- generate_trends(indices = route_indices)$trends
point_trends <- generate_trends(indices = point_indices)$trends
detect_trends <- generate_trends(indices = detect_indices)$trends
varprop_trends <- generate_trends(indices = varprop_indices)$trends

route_trends$ci_width <- route_trends$trend_q_0.95 - route_trends$trend_q_0.05
point_trends$ci_width <- point_trends$trend_q_0.95 - point_trends$trend_q_0.05
detect_trends$ci_width <- detect_trends$trend_q_0.95 - detect_trends$trend_q_0.05
varprop_trends$ci_width <- varprop_trends$trend_q_0.95 - varprop_trends$trend_q_0.05

ci_width_df <- data.frame(Model = c(rep("ROUTE", times = nrow(route_trends)),
                                    rep("POINT", times = nrow(point_trends)),
                                    rep("DETECT", times = nrow(detect_trends)),
                                    rep("VARPROP", times = nrow(varprop_trends))),
                          CI_Width = c(as.vector(route_trends$ci_width),
                                       as.vector(point_trends$ci_width),
                                       as.vector(detect_trends$ci_width),
                                       as.vector(varprop_trends$ci_width)))
ci_width_df$Model <- factor(ci_width_df$Model,
                               levels = c("ROUTE", "POINT", "DETECT", "VARPROP"))

ci_width_boxplot <- ggplot(data = ci_width_df) +
  geom_boxplot(aes(x = Model, y = CI_Width, group = Model)) +
  ylab("90% CI Width (Trend)") +
  ylim(0, max(ci_width_df$CI_Width + 10)) +
  NULL

route_summary <- as.vector(unname(summary(route_trends$ci_width))) 
point_summary <- as.vector(unname(summary(point_trends$ci_width)))
detect_summary <- as.vector(unname(summary(detect_trends$ci_width)))
varprop_summary <- as.vector(unname(summary(varprop_trends$ci_width)))

model <- c("ROUTE", "POINT", "DETECT", "VARPROP")
summary_table <- rbind(route_summary, point_summary, detect_summary, varprop_summary)
summary_table <- as.data.frame(cbind(model, summary_table))
names(summary_table) <- c("model", "min", "q1", "median", "mean", "q3", "max")

####### Output ####################################

write.table(file = "output/precision_summary.csv", x = summary_table, sep = ",", row.names = FALSE)

png(filename = "output/plots/precision_comparison.png",
    width = 5, height = 3, res = 300, units = "in")
print(ci_width_boxplot)
dev.off()
