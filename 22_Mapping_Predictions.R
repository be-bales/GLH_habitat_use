## make decadal averages for 2000-2009 and 2010-2019 for comparison

env_stack_all <- rast("Data/env_stack_all.tif")
ncell(static_stack[[1]])

## present day rasters 
preds_2018 <- readRDS("predictions/preds_2018.rds")
preds_2018 <- as.data.frame(preds_2018, na.rm = FALSE)
complete_cases_preesent <- complete.cases(preds_2018)
str(complete_cases_preesent)

## ## select predictions and see raster 
preds_2007_spring <- readRDS("predictions/preds_2007_spring.rds")
prediction_raster <- deepcopy(env_stack_all[[1]])
ncell(env_stack_all[[1]])
values(prediction_raster) <- NA
values(prediction_raster)[complete_cases] <- preds_2018$est
plot(prediction_raster)

##turn into a function making decadal means 
make_decadal_mean_from_preds <- function(years,
                                         #seasons,
                                         pred_dir = "predictions",
                                         template_raster = env_stack_all[[1]],
                                         complete_cases,
                                         pred_col = "est",
                                         out_file) {
  
  pred_rasters <- list()
  k <- 1
  
  for (yr in years) {
    #for (ss in seasons) {
      
      cat("Loading predictions for:", yr, "\n")
      
      preds_file <- sprintf("%s/preds_%s.rds", pred_dir, yr)
      preds_yr <- readRDS(preds_file)
      
      # make empty raster
      r <- deepcopy(template_raster)
      values(r) <- NA
      
      # fill using your exact workflow
      values(r)[complete_cases] <- preds_yr[[pred_col]]
      
      # store
      pred_rasters[[k]] <- r
      k <- k + 1
      
      rm(preds_yr, r)
      gc()
    }
  
  
  # stack all rasters
  pred_stack <- rast(pred_rasters)
  
  # decadal mean
  decadal_mean <- app(pred_stack, mean, na.rm = TRUE)
  
  # save
  writeRaster(decadal_mean, out_file, overwrite = TRUE)
  
  return(decadal_mean)
  }


Pred_2010_2019 <- make_decadal_mean_from_preds(
  years = 2000:2019,
  pred_dir = "predictions",
  template_raster = env_stack_all[[1]],
  complete_cases = complete_cases,   
  pred_col = "est",
  out_file = "preds_est_2000_2019.tif"
)

### Climate/predictive rasters sdmTMB
preds_119_2020 <- readRDS("climate_predictions/preds_119_2030.rds")
preds_119_2020 <- as.data.frame(preds_119_2020, na.rm = FALSE)
complete_cases_climate <- complete.cases(preds_119_2020)
str(complete_cases_climate)

make_climate_mean_from_preds <- function(scenarios, ##OLD
                                         years,
                                         #seasons,
                                         pred_dir = "climate_predictions", ## change as needed between models
                                         template_raster = env_stack_all[[1]],
                                         complete_cases,
                                         pred_col = "est",
                                         out_file) {
  
  pred_rasters <- list()
  k <- 1
  
  for (ssp in scenarios){
  for (yr in years) {
    #for (ss in seasons) {
      
      cat("Loading predictions for:", ssp, yr, "\n")
      
      preds_file <- sprintf("%s/preds_%s_%s.rds", pred_dir, ssp, yr)
      preds_ssp_yr <- readRDS(preds_file)
      
      # make empty raster
      r <- deepcopy(template_raster)
      values(r) <- NA
      
      # fill using your exact workflow
      values(r) <- preds_ssp_yr[[pred_col]]
      
      # store
      pred_rasters[[k]] <- r
      k <- k + 1
      
      rm(preds_ssp_yr, r)
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


stack_file <- sprintf("climate_prediction_stacks/stack_119_2020.rds")
prediction_stack_ssp_year <- readRDS(stack_file)
prediction_frame_ssp_year <- as.data.frame(prediction_stack_ssp_year, na.rm = FALSE)
complete_cases <- complete.cases(prediction_frame_ssp_year)
str(complete_cases)

static_stack <- rast("Data/static_stack.tiff")

make_sdm_climate_rasters <- function(scenarios,
                                     years,
                                     pred_dir = "climate_predictions",
                                     template_raster = env_stack_all[[1]],
                                     complete_cases,
                                     out_file) {
  
  pred_rasters <- list()
  k <- 1
  
  for (ssp in scenarios){
    for (yr in years) {
      
      cat("Loading predictions for:", ssp, yr, "\n")
      
      preds_file <- sprintf("%s/preds_%s_%s.rds", pred_dir, ssp, yr)
      preds_ssp_yr <- as.data.frame(readRDS(preds_file))
      
      # make empty raster
      r <- deepcopy(template_raster)
      values(r) <- NA
      
      # fill using your exact workflow
      values(r)[complete_cases] <- preds_ssp_yr$est
      
      # store
      pred_rasters[[k]] <- r
      k <- k + 1
      
      rm(preds_ssp_yr, r)
      gc()
    }}
  
  
  # stack all rasters
  pred_stack <- rast(pred_rasters)
  
  # decadal mean
  decadal_mean <- app(pred_stack, mean, na.rm = TRUE)
  
  # save
  writeRaster(decadal_mean, out_file, overwrite = TRUE)
  
  return(decadal_mean)
}

make_sdm_climate_rasters(
  scenarios = 119,
  years = 2040:2049,
  pred_dir = "climate_predictions",
  template_raster = env_stack_all[[1]],
  complete_cases = complete_cases,   
  out_file = "sdm_119_2040s.tif"
)


### Climate/predictive rasters MaxEnt
maxent_preds <- readRDS("maxent_climate_predictions/preds_119_2020_fall.rds")

stack_2000_spring <- as.data.frame(stack_2000_spring, na.rm = FALSE)

complete_cases <- complete.cases(stack_2000_spring)

make_maxent_climate_rasters <- function(scenarios,
                                        years,
                                         pred_dir = "maxent_climate_predictions",
                                         template_raster = env_stack_all[[1]],
                                         complete_cases,
                                         out_file) {
  
  pred_rasters <- list()
  k <- 1
  
  for (ssp in scenarios){
  for (yr in years) {

    cat("Loading predictions for:", ssp, yr, "\n")
    
    preds_file <- sprintf("%s/preds_%s_%s.rds", pred_dir, ssp, yr)
    preds_ssp_yr <- readRDS(preds_file)
    
    # make empty raster
    r <- deepcopy(template_raster)
    values(r) <- NA
    
    # fill using your exact workflow
    values(r)[complete_cases] <- preds_ssp_yr
    
    # store
    pred_rasters[[k]] <- r
    k <- k + 1
    
    rm(preds_ssp_yr, r)
    gc()
  }}
  
  
  # stack all rasters
  pred_stack <- rast(pred_rasters)
  
  # decadal mean
  decadal_mean <- app(pred_stack, mean, na.rm = TRUE)
  
  # save
  writeRaster(decadal_mean, out_file, overwrite = TRUE)
  
  return(decadal_mean)
}

make_maxent_climate_rasters(
  scenarios = 585,
  years = 2040:2049,
  pred_dir = "maxent_climate_predictions",
  template_raster = env_stack_all[[1]],
  complete_cases = complete_cases,   
  out_file = "maxent_585_2040s.tif"
)