####### Script Information ########################
# Brandon P.M. Edwards
# BBS Point Level
# <03-run-point-model.R>
# Created June 2023
# Last Updated April 2024

####### Import Libraries and External Files #######

library(bbsBayes2)

####### Set Constants #############################

species_list <- c("Ovenbird", "Wood Thrush", "White-throated Sparrow", 
                  "Black-throated Green Warbler", "Winter Wren",
                  "Pileated Woodpecker", "Blue-headed Vireo", 
                  "American Goldfinch", "Eastern Wood-pewee")
sp_code <- c("OVEN", "WOTH", "WTSP", "BTNW", "WIWR", "PIWO",
             "BHVI", "AMGO", "EAWP")

st <- "latlong"

####### Read Data #################################

bbs_counts <- readRDS(file = "data/generated/bbs_counts.RDS")
bbs_sites <- readRDS(file = "data/generated/bbs_sites_fc.RDS"); bbs_sites$rt_st <- NULL
bbs_species <- load_bbs_data(level = "stop")$species

####### Main Code #################################

bbs_data <- list(birds = bbs_counts,
                 routes = bbs_sites,
                 species = bbs_species)

for (i in 1:length(species_list))
{
  sp <- species_list[i]
  sp_out <- sp_code[i]
  
  print(sp)
  bbs_stratified <- stratify(by = st, level = "stop", species = sp, data_custom = bbs_data)
  
  # Limit analysis to only Canada
  bbs_stratified$routes_strata <- 
    bbs_stratified$routes_strata[which(bbs_stratified$routes_strata$country == "CA"), ]
  
  mod_prepped <- prepare_data(strata_data = bbs_stratified,
                              min_year = 2010,
                              min_n_routes = 1) %>%
    prepare_spatial(strata_map = load_map(st)) %>%
    prepare_model(model = "gamye",
                  model_file = "models/gamye_spatial_route_hyperparameter.stan", 
                  model_variant = "spatial")
  
  route_original_list <- strsplit(mod_prepped$raw_data$route,
                             split = "-",
                             fixed = TRUE)
  mod_prepped$raw_data$route_original <- paste0(sapply(route_original_list, "[[", 2),
                           "-",
                           sapply(route_original_list, "[[", 3))
  mod_prepped$model_data$route_original <- as.numeric(as.factor(mod_prepped$raw_data$route_original))
  mod_prepped$model_data$n_route_original <- length(unique(mod_prepped$model_data$route_original))
  
  model_run <- run_model(model_data = mod_prepped,
                         output_basename = paste0(sp_out, "-point_RH"),
                         output_dir = "output/model_runs",
                         overwrite = TRUE)
}
