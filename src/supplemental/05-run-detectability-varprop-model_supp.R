####### Script Information ########################
# Brandon P.M. Edwards
# BBS Point Level
# <05-run-detectability-varprop-model.R>
# Created June 2024
# Last Updated August 2024

####### Import Libraries and External Files #######

library(bbsBayes2)
library(napops)
library(suntools)
library(lutz)
library(terra)

####### Set Constants #############################

species_list <- c("Ovenbird", "Wood Thrush","Pileated Woodpecker", 
                  "Blue-headed Vireo", "Black-throated Green Warbler", 
                  "Winter Wren", "American Goldfinch", "Eastern Wood-pewee", 
                  "White-throated Sparrow")
sp_code_list <- c("OVEN", "WOTH", "PIWO", "BHVI", "BTNW", "WIWR", "AMGO", 
                  "EAWP", "WTSP")
st <- "latlong"

####### Read Data #################################

bbs_counts <- readRDS(file = "data/generated/bbs_counts.RDS")
bbs_sites <- readRDS(file = "data/generated/bbs_sites_fc.RDS"); bbs_sites$rt_st <- NULL
bbs_species <- load_bbs_data(level = "stop")$species

nalcms_2010 <- rast("data/raw/spatial/nalcms/CAN_Land_cover_2010v2_30m_TIF/CAN_NALCMS_landcover_2010v2_30m/data/CAN_NALCMS_landcover_2010v2_30m.tif")
nalcms_2015 <- rast("data/raw/spatial/nalcms/can_land_cover_2015v3_30m_tif/CAN_NALCMS_landcover_2015v3_30m/data/CAN_NALCMS_landcover_2015v3_30m.tif")
nalcms_2020 <- rast("data/raw/spatial/nalcms/can_land_cover_2020_30m_tif/CAN_NALCMS_landcover_2020_30m/data/CAN_NALCMS_landcover_2020_30m.tif")


####### Main Code #################################

bbs_data <- list(birds = bbs_counts,
                 routes = bbs_sites,
                 species = bbs_species)
rm(bbs_counts); rm(bbs_sites); rm(bbs_species); gc()

for (i in 1:length(species_list))
{
  sp <- species_list[i]
  sp_code <- sp_code_list[i]
  
  bbs_stratified <- stratify(by = st, level = "stop", species = sp, data_custom = bbs_data)
  
  # Limit analysis to only Canada
  bbs_stratified$routes_strata <- 
    bbs_stratified$routes_strata[which(bbs_stratified$routes_strata$country == "CA"), ]
  
  prepared_data <- prepare_data(strata_data = bbs_stratified,
                                min_year = 2010,
                                min_n_routes = 1) %>%
    prepare_spatial(strata_map = load_map(st))
  
  mod_prepped <- prepare_model(prepared_data = prepared_data,
                               model_file = "models/gamye_spatial_detectability_varprop_RH.stan", 
                               model = "gamye", 
                               model_variant = "spatial")
  
  route_original_list <- strsplit(mod_prepped$raw_data$route,
                                  split = "-",
                                  fixed = TRUE)
  route_original <- paste0(sapply(route_original_list, "[[", 2),
                           "-",
                           sapply(route_original_list, "[[", 3))
  route_original <- as.numeric(as.factor(route_original))
  mod_prepped$model_data$n_route_original <- length(unique(route_original))
  
  route_lookup <- data.frame(Site = mod_prepped$model_data$site,
                             Route = route_original)
  route_lookup <- route_lookup[!duplicated(route_lookup$Site), ]
  route_lookup <- route_lookup[order(route_lookup$Site), ]
  mod_prepped$model_data$route_lookup <- as.vector(route_lookup$Route)
  

  #' Consider doing this stuff below (i.e. extracting the ordinal day and other
  #' detectability variables) within bbsBayes. There might be opportunity to add
  #' that functionality into bbsBayes so that others can have detectability offsets
  #' as well. Plus I think that might be the easiest way to then subset the offsets
  #' properly based on the BBS data that are going to be included in the model.
  
  subsetted_data <- mod_prepped$raw_data
  
  # Ordinal Day
  od <- as.POSIXlt(sprintf("%04d-%02d-%02d", 
                           subsetted_data$year,
                           subsetted_data$month,
                           subsetted_data$day),
                   format = "%Y-%m-%d")$yday
  
  # Time since local sunrise
  coords <- matrix(c(subsetted_data$longitude, subsetted_data$latitude),
                   nrow = nrow(subsetted_data))
  
  date_time <- as.POSIXct(sprintf("%04d-%02d-%02d %s", 
                                  subsetted_data$year,
                                  subsetted_data$month,
                                  subsetted_data$day,
                                  subsetted_data$start_time),
                          format = "%Y-%m-%d %H%M")
  
  time_zone <- tz_lookup_coords(lat = subsetted_data$latitude,
                                lon = subsetted_data$longitude)
  
  utc_offset <- numeric(length = length(time_zone))
  for (i in seq_along(time_zone))
  {
    utc_offset[i] <- lutz::tz_offset(date_time[i], 
                                     tz = time_zone[i])$utc_offset_h
  }
  
  utc_time <- lubridate::force_tz(date_time  + (utc_offset*-1 * 60 * 60),
                                  tzone = "UTC")
  
  sunrise_times <- sunriset(crds = coords,
                            dateTime = utc_time,
                            direction = "sunrise",
                            POSIXct.out = TRUE)
  
  tssr <- as.numeric(difftime(utc_time, sunrise_times$time, units = "hours"))
  
  #' Determine which detectability model to use for each species.
  #' I may eventually decide to choose the "fullest best model", where if there is
  #' a tied (or within 5 delta AIC) model, I hcoose the fuller one.
  p_best_model_df <- coef_removal(species = sp_code)[,c("Model", "AIC")]
  p_best_model <- p_best_model_df[which(p_best_model_df$AIC == min(p_best_model_df$AIC)), "Model"]
  
  q_best_model_df <- coef_distance(species = sp_code)[, c("Model", "AIC")]
  q_best_model <- q_best_model_df[which(q_best_model_df$AIC == min(q_best_model_df$AIC)), "Model"]
  
  kappa_p <- avail_fd(species = sp_code,
                      model = p_best_model,
                      od = od,
                      tssr = tssr,
                      pairwise = TRUE,
                      time = rep(3, times = length(od)))
  
  kappa_q <- percept_fd(species = sp_code,
                        model = q_best_model,
                        road = rep(TRUE, times = nrow(subsetted_data)),
                        forest = subsetted_data$forest_coverage,
                        pairwise = TRUE,
                        distance = rep(400, times = length(od)))
  
  p <- avail(species = sp_code,
             model = p_best_model,
             od = od,
             tssr = tssr,
             pairwise = TRUE,
             time = rep(3, times = length(od)))$p
  
  q <- percept(species = sp_code,
               model = q_best_model,
               road = rep(TRUE, times = nrow(subsetted_data)),
               forest = subsetted_data$forest_coverage,
               pairwise = TRUE,
               distance = rep(400, times = length(od)))$q
  
  p_vcv <- get_vcv(species = sp_code,
                   model_type = "rem",
                   model_num = p_best_model)
  q_vcv <- get_vcv(species = sp_code, 
                   model_type = "dis", 
                   model_num = q_best_model)
  
  #' In order to complete the next step (that is, making an array for kappa_ that
  #' contain indices which correspond to specific stratum-year-site combos), we need
  #' to append a matrix of kappa_p and kappa_q onto the end of the current kappa_p
  #' and kappa_q, that correspond to stratum-year-site combos that were NOT run
  #' but would be considered in index calculation. In essense, if a site
  #' was not sampled during a year (or more) of the time series we were interested
  #' in, we still would make an estimate for it with the index of abundance calculations.
  #' thus we need information about the detectability at that site, but don't yet have that
  #' because the BBS data only show a route IF it was sampled. okay here we go.
  
  raw_data <- mod_prepped$raw_data
  model_data <- mod_prepped$model_data
  
  raw_data$strat_prep <- model_data$strat
  raw_data$year_prep <- model_data$year
  raw_data$site_prep <- model_data$site
  raw_data$yr_ste <- paste0(raw_data$year_prep, "-",
                                raw_data$site_prep)
  
  #' Now need to generate a list of ALL potential strat-year-site combos,
  #' if all sites were sampled every year of the time series.
  year_all <- rep(seq(1, model_data$n_years), times = model_data$n_sites)
  yr_ste_all <- paste0(year_all, "-",
                       rep(seq(1, model_data$n_sites), each = model_data$n_years))
  
  # Which site-year combos were not run?
  not_run_df <- data.frame(yr_ste = setdiff(yr_ste_all, raw_data$yr_ste))
  yr_and_st <- strsplit(not_run_df$yr_ste, split = "-")
  not_run_df$year <- unlist(lapply(yr_and_st, `[[`, 1))
  not_run_df$site <- unlist(lapply(yr_and_st, `[[`, 2))
  
  # get latitude and longitude of these sites
  site_lookup <- raw_data[!duplicated(raw_data$site_prep), ]
  not_run_df <- merge(not_run_df, site_lookup[, c("site_prep", "longitude", "latitude", "strat_prep")],
                      by.x = "site", by.y = "site_prep")

  # get average od and tssr per site
  raw_data$od <- od
  raw_data$tssr <- tssr
  mean_od_lookup <- aggregate(od ~ site_prep, data = raw_data, FUN = mean)
  mean_tssr_lookup <- aggregate(tssr ~ site_prep, data = raw_data, FUN = mean)
  not_run_df <- merge(not_run_df, mean_od_lookup,
                      by.x = "site", by.y = "site_prep")
  not_run_df <- merge(not_run_df, mean_tssr_lookup,
                      by.x = "site", by.y = "site_prep")
  
  # get forest coverage...again
  not_run_df$forest_coverage <- NA
  not_run_df$year_actual <- as.numeric(not_run_df$year) + 2009
  for (s in unique(not_run_df$site))
  {
    temp <- not_run_df[which(not_run_df$site == s), ]
    
    forest_df <- data.frame(year = seq(2010, max(2020, max(temp$year_actual))))
    
    forest_df <- merge(forest_df, temp[, c("year_actual", "latitude", "longitude")],
                       by.x = "year", by.y = "year_actual", all = TRUE)
    
    #' Check that all missing latitude and longitudes can safely just be copied 
    #' based on latitudes and longitudes that are availabile. I.e., if the route
    #' location had changed overtime, we have to account for this
    
    #if (var(forest_df$latitude, na.rm = TRUE) == 0)
    # {
    forest_df[which(is.na(forest_df$latitude)), "latitude"] <- median(forest_df$latitude, na.rm = TRUE)
    # }
    
    #  if (var(forest_df$longitude, na.rm = TRUE) == 0)
    #{
    forest_df[which(is.na(forest_df$longitude)), "longitude"] <- median(forest_df$longitude, na.rm = TRUE)
    #}
    
    # landcover in 2010
    lat <- forest_df[which(forest_df$year == 2010), "latitude"]
    lon <- forest_df[which(forest_df$year == 2010), "longitude"]
    
    bbs_point <- project(vect(matrix(c(lon, lat), ncol = 2),
                              crs = "+proj=longlat"),
                         crs(nalcms_2010))
    bbs_radius <- buffer(bbs_point, width = 400)
    landcover_types <- extract(nalcms_2010, bbs_radius)
    n_forest <- nrow(landcover_types[which(grepl("forest", 
                                                 landcover_types$Class_Name)),])
    fc_2010 <- n_forest / nrow(landcover_types)
    
    # landcover in 2015
    lat <- forest_df[which(forest_df$year == 2015), "latitude"]
    lon <- forest_df[which(forest_df$year == 2015), "longitude"]
    
    bbs_point <- project(vect(matrix(c(lon, lat), ncol = 2),
                              crs = "+proj=longlat"),
                         crs(nalcms_2015))
    bbs_radius <- buffer(bbs_point, width = 400)
    landcover_types <- extract(nalcms_2015, bbs_radius)
    n_forest <- nrow(landcover_types[which(landcover_types$CAN_NALCMS_landcover_2015v3_30m <= 6),])
    fc_2015 <- n_forest / nrow(landcover_types)
    
    # landcover in 2020
    lat <- forest_df[which(forest_df$year == 2020), "latitude"]
    lon <- forest_df[which(forest_df$year == 2020), "longitude"]
    
    bbs_point <- project(vect(matrix(c(lon, lat), ncol = 2),
                              crs = "+proj=longlat"),
                         crs(nalcms_2020))
    bbs_radius <- buffer(bbs_point, width = 400)
    landcover_types <- extract(nalcms_2020, bbs_radius)
    n_forest <- nrow(landcover_types[which(landcover_types$CAN_NALCMS_landcover_2020_30m <= 6),])
    fc_2020 <- n_forest / nrow(landcover_types)
    
    # Interpolate missing values
    # 2010 - 2015
    values_tp1 <- seq(from = fc_2010,
                      to = fc_2015,
                      length.out = 6)
    
    # 2015 - 2020
    values_tp2 <- seq(from = fc_2015,
                      to = fc_2020,
                      length.out = 6)
    
    #' Since we don't know the current rate of landcover change for any years past 2020,
    #' we will just copy the 2020 values for the remaining years
    if (nrow(forest_df[which(forest_df$year > 2020), ]) > 0)
    {
      values_tp3 <- rep(fc_2020,
                        times = nrow(forest_df[which(forest_df$year > 2020), ])) 
    }else {
      values_tp3 <- NULL
    }
    
    fc_vals <- c(values_tp1, values_tp2[2:length(values_tp2)], values_tp3)
    forest_df$forest_coverage <- fc_vals
    
    forest_df_red <- forest_df[which(forest_df$year %in% temp$year_actual), ]
    
    not_run_df[which(not_run_df$site == s), "forest_coverage"] <- forest_df_red$forest_coverage
  }
  
  not_run_df$str_yr_ste <- paste0(not_run_df$strat_prep, "-", not_run_df$yr_ste)
  
  kappa_p_notrun <- avail_fd(species = sp_code,
                      model = p_best_model,
                      od = not_run_df$od,
                      tssr = not_run_df$tssr,
                      pairwise = TRUE,
                      time = rep(3, times = length(not_run_df$od)))
  
  kappa_q_notrun <- percept_fd(species = sp_code,
                        model = q_best_model,
                        road = rep(TRUE, times = nrow(not_run_df)),
                        forest = not_run_df$forest_coverage,
                        pairwise = TRUE,
                        distance = rep(400, times = length(not_run_df$od)))
  
  #' Need to make an array for both kappa_p and kappa_q that contains indices which
  #' correspond to the specific stratum-year-site combination that each value of
  #' kappa_p and kappa_q refer to. We must do this to match the data setup in the model
  #' for calculating n_t. In particular, this array is going to have the following dimensions:
  #' n_strata
  #' max_n_obs_sites_strata
  #' n_years
  #' Then, we can fill each entry with the index that corresponds to the specific
  #' strata, year, and site combination (site will be obtained from ste_mat)
  ste_mat <- mod_prepped$model_data$ste_mat
  str_yr_ste <- paste0(mod_prepped$model_data$strat, "-",
                       mod_prepped$model_data$year, "-",
                       mod_prepped$model_data$site)
  str_yr_ste_all <- c(str_yr_ste, not_run_df$str_yr_ste)
  kappa_indices <- array(data = 0, dim = c(dim(ste_mat), mod_prepped$model_data$n_years))
  
  for (y in 1:mod_prepped$model_data$n_years)
  {
    for (s in 1:mod_prepped$model_data$n_strata)
    {
      for (t in 1:mod_prepped$model_data$n_obs_sites_strata[s])
      {
        sys_current <- paste0(s, "-", y, "-", ste_mat[s,t])
        #if (sys_current %in% str_yr_ste_all)
      #  {
          kappa_indices[s,t,y] <- which(str_yr_ste_all == sys_current)
        #}else
        #{
       #   kappa_indices[s,t,y] <- NA
        #}
      }
    }
  }
  
  
  
  detectability_data_list <- list(n_avail_covs = ncol(kappa_p),
                                  n_percept_covs = ncol(kappa_q),
                                  kappa_p = kappa_p,
                                  kappa_p_notrun = kappa_p_notrun,
                                  kappa_q = kappa_q,
                                  kappa_q_notrun = kappa_q_notrun,
                                  n_notrun = nrow(kappa_p_notrun),
                                  kappa_indices = kappa_indices,
                                  vcv_p = p_vcv,
                                  vcv_q = q_vcv,
                                  p = p,
                                  q = q)
  
  mod_prepped$model_data <- c(mod_prepped$model_data, detectability_data_list)
  
  model_run <- run_model(model_data = mod_prepped,
                         chains = 4,
                         parallel_chains = 4,
                         output_basename = paste0(sp_code, "-varprop_RH"),
                         output_dir = "output/model_runs",
                         overwrite = TRUE)
}
