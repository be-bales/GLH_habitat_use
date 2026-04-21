## make decadal averages for 2000-2009 and 2010-2019 for comparison
library(rasterVis)
library(RColorBrewer)

## call prediction file for year and season, make into raster and assign
## est values to the raster
library(terra)
years <- 2000:2019
seasons <- c("spring", "fall")

land_vect <- vect("Data/Land_shp.shp")
land_vect <- terra::project(land_vect, "EPSG:32622")
land_vect <- terra::crop(land_vect, study_area)
plot(land_vect)

env_stack_all <- rast("Data/env_stack_all.tif")

### Climate/predictive rasters sdmTMB
stack_2000_spring <- readRDS("prediction_stacks/stack_2000_spring.rds")
stack_2000_spring <- as.data.frame(stack_2000_spring, na.rm = FALSE)
complete_cases <- complete.cases(stack_2000_spring)


make_decadal_mean_from_preds <- function(years,
                                         seasons,
                                         pred_dir = "predictions",
                                         template_raster = env_stack_all[[1]],
                                         complete_cases,
                                         pred_col = "est_rf",
                                         out_file) {
  
  pred_rasters <- list()
  k <- 1
  
  for (yr in years) {
    for (ss in seasons) {
      
      cat("Loading predictions for:", yr, ss, "\n")
      
      preds_file <- sprintf("%s/preds_%s_%s.rds", pred_dir, yr, ss)
      preds_yr_ss <- readRDS(preds_file)
      
      # make empty raster
      r <- deepcopy(template_raster)
      values(r) <- NA
      
      # fill using your exact workflow
      values(r)[complete_cases] <- preds_yr_ss[[pred_col]]
      
      # store
      pred_rasters[[k]] <- r
      k <- k + 1
      
      rm(preds_yr_ss, r)
      gc()
    }
  }
  
  # stack all rasters
  pred_stack <- rast(pred_rasters)
  
  # decadal mean
  decadal_mean <- app(pred_stack, mean, na.rm = TRUE)
  
  # save
  writeRaster(decadal_mean, out_file, overwrite = TRUE)
  
  return(decadal_mean)
}

all_random_2010_2019 <- make_decadal_mean_from_preds(
  years = 2010:2019,
  seasons, 
  pred_dir = "predictions",
  template_raster = env_stack_all[[1]],
  complete_cases = complete_cases,   
  pred_col = "est_rf",
  out_file = "decadal_est_rf_2010_2019.tif"
)

decadal_epsilon_2010_2019 <- make_decadal_mean_from_preds(
  years = 2010:2019,
  seasons,
  pred_dir = "predictions",
  template_raster = env_stack_all[[1]],
  complete_cases = complete_cases,   
  pred_col = "epsilon_st",
  out_file = "decadal_epsilon_2010_2019.tif"
)

present_est_2000_2019 <- make_decadal_mean_from_preds(
  years = 2000:2019,
  seasons,
  pred_dir = "predictions",
  template_raster = env_stack_all[[1]],
  complete_cases = complete_cases,   
  pred_col = "est",
  out_file = "present_est_2000_2019.tif"
)

decades_stack <- c(decadal_mean_2000_2009, decadal_mean_2010_2019)
diff_Raster <- decadal_mean_2010_2019 -decadal_mean_2000_2009
plot(diff_Raster)


make_climate_mean_from_preds <- function(scenarios,
                                         years,
                                         seasons,
                                         pred_dir = "maxent_climate_predictions", ## change as needed between models
                                         template_raster = env_stack_all[[1]],
                                         complete_cases,
                                         pred_col = "est",
                                         out_file) {
  
  pred_rasters <- list()
  k <- 1
  
  for (ssp in scenarios){
  for (yr in years) {
    for (ss in seasons) {
      
      cat("Loading predictions for:", ssp, yr, ss, "\n")
      
      preds_file <- sprintf("%s/preds_%s_%s_%s.rds", pred_dir, ssp, yr, ss)
      preds_ssp_yr_ss <- readRDS(preds_file)
      
      # make empty raster
      r <- deepcopy(template_raster)
      values(r) <- NA
      
      # fill using your exact workflow
      values(r)[complete_cases] <- preds_ssp_yr_ss
      
      # store
      pred_rasters[[k]] <- r
      k <- k + 1
      
      rm(preds_ssp_yr_ss, r)
      gc()
    }
  }}
  
  # stack all rasters
  pred_stack <- rast(pred_rasters)
  
  # decadal mean
  decadal_mean <- app(pred_stack, mean, na.rm = TRUE)
  
  # save
  writeRaster(decadal_mean, out_file, overwrite = TRUE)
  
  return(decadal_mean)
}

maxent_raster_20_585 <- make_climate_mean_from_preds(
  scenarios = 585,
  years = 2020:2029,
  seasons,
  pred_dir = "maxent_climate_predictions",
  template_raster = env_stack_all[[1]],
  complete_cases = complete_cases,   
  pred_col,
  out_file = "test_maxent.tif"
)

#plot(maxent_raster_40)

writeRaster(maxent_raster_20_585, "climate_rasters/maxent_raster_20_585.tif")

decadal_mean_2040_2049_245 <- make_climate_mean_from_preds(
  scenarios = 245,
  years = 2040:2049,
  seasons,
  pred_dir = "climate_predictions",
  template_raster = env_stack_all[[1]],
  complete_cases = complete_cases,   
  pred_col = "est",
  out_file = "decadal_mean_2040_2049_245.tif"
)

decadal_mean_2040_2049_585 <- make_climate_mean_from_preds(
  scenarios = 585,
  years = 2040:2049,
  seasons,
  pred_dir = "climate_predictions",
  template_raster = env_stack_all[[1]],
  complete_cases = complete_cases,   
  pred_col = "est",
  out_file = "decadal_mean_2040_2049_585.tif"
)

plot(decadal_mean_2020_2029_585)

minus_Rast <- decadal_mean_2040_2049_585 - decadal_mean_2020_2029_585
plot(minus_Rast)

decadal_mean_2000_2009 <- rast("climate_rasters/decadal_mean_2000_2009.tif")
plot(decadal_mean_2000_2009)
decadal_mean_2010_2019 <- rast("climate_rasters/decadal_mean_2010_2019.tif")
plot(decadal_mean_2010_2019)
decadal_mean_2010_2019 <- rast("climate_rasters/decadal_mean_2010_2019.tif")


### Climate/predictive rasters MaxEnt
maxent_preds <- readRDS("maxent_climate_predictions/preds_119_2020_fall.rds")

stack_2000_spring <- as.data.frame(stack_2000_spring, na.rm = FALSE)

complete_cases <- complete.cases(stack_2000_spring)
