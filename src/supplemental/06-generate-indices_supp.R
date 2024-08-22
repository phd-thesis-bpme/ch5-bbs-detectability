####### Script Information ########################
# Brandon P.M. Edwards
# BBS Point Level
# <06-generate-indices.R>
# Created December 2023
# Last Updated June 2024

####### Import Libraries and External Files #######

library(bbsBayes2)

####### Set Constants #############################

sp_list <- c("OVEN")

####### Calculate Indices and Output ##############

for (sp in sp_list)
{
  route <- readRDS(paste0("output/model_runs/", sp, "-route.rds"))
  point <- readRDS(paste0("output/model_runs/", sp, "-point_RH.rds"))
  detectability <- readRDS(paste0("output/model_runs/", sp, "-detectability_RH.rds"))
  varprop <- readRDS(paste0("output/model_runs/", sp, "-varprop.rds"))
  
  indices_route <- generate_indices(model_output = route)
  saveRDS(object = indices_route,
          file = paste0("output/indices/", sp, "_route.RDS"))
  
  indices_point <- generate_indices(model_output = point)
  saveRDS(object = indices_point,
          file = paste0("output/indices/", sp, "_point_RH.RDS"))
  
  indices_detectability <- generate_indices(model_output = detectability)
  saveRDS(object = indices_detectability,
          file = paste0("output/indices/", sp, "_detectability_RH.RDS")) 
  
  indices_varprop <- generate_indices(model_output = varprop)
  saveRDS(object = indices_varprop,
          file = paste0("output/indices/", sp, "_varprop.RDS"))
}
