## refit model with extra years 
extra_years <- c(2020:2050)

## refit with new years
mod_climate <- sdmTMB(
  data= pres_abs_scaled,
  formula = presence ~ s(depth) + s(temp) + s(sal) + s(bpi_coarse) + 
    s(bpi_fine) + tri + smf + season,
  mesh = mesh3,
  family = binomial(link = "logit"),
  spatial = "on", 
  time = "year",
  spatiotemporal = "RW",
  extra_time = extra_years
)

summary(mod_climate)
sanity(mod_climate)

## download bio-ORACLE data for SSP 1-1.9, 2-4.5, and 5-8.5 
## until 2050

scenarios <- c("SSP 1-1.9", "SSP 2-4.5", "SSP 5-8.5")
decades <- list("2020-2030" = 2025,
                "2030-2040" = 2035,
                "2040-2050" = 2045)

 #####
time <- c("2020-01-01T00:00:00Z", "2090-01-01T00:00:00Z")
latitude = c(35,79)
longitude = c(-80,-42)
constraints = list(time, latitude, longitude)
names(constraints) = c("time", "latitude", "longitude")

## 1-1.9 

temp_119 <- "thetao_ssp119_2020_2100_depthmean"
sal_119 <- "so_ssp119_2020_2100_depthmean"
#do_245 <- "o2_ssp119_2020_2100_depthmean"
temp_layer_119 <- download_layers(temp_119, temp_variable, constraints)
sal_layer_119 <-download_layers(sal_119, sal_variable, constraints)
## 2020-2030 
temp_layer_119_2020_2030 <- temp_layer_119[[1]]
sal_layer_119_2020_2030 <- sal_layer_119[[1]]
## 2030-2040
temp_layer_119_2030_2040 <- temp_layer_119[[2]]
sal_layer_119_2030_2040 <- sal_layer_119[[2]]
#do_layer_119_2030_2040 <- do_layer_119[[2]]
## 2040-2050
temp_layer_119_2040_2050 <- temp_layer_119[[3]]
sal_layer_119_2040_2050 <- sal_layer_119[[3]]
#do_layer_119_2040_2050 <- do_layer_119[[3]]

## 2-4.5
temp_245 <- "thetao_ssp245_2020_2100_depthmean"
sal_245 <- "so_ssp245_2020_2100_depthmean"
#do_245 <- "o2_ssp245_2020_2100_depthmean"
temp_layer_245 <- download_layers(temp_245, temp_variable, constraints)
sal_layer_245 <-download_layers(sal_245, sal_variable, constraints)
#do_layer_245 <-download_layers(do_245, do_variable, constraints)

## 2020-2030
temp_layer_245_2020_2030 <- temp_layer_245[[1]]
sal_layer_245_2020_2030 <- sal_layer_245[[1]]
#do_layer_245_2020_2030 <- do_layer_245[[1]]

## 2030-2040
temp_layer_245_2030_2040 <- temp_layer_245[[2]]
sal_layer_245_2030_2040 <- sal_layer_245[[2]]
#do_layer_245_2030_2040 <- do_layer_245[[2]]

## 2040-2050 
temp_layer_245_2040_2050 <- temp_layer_245[[3]]
sal_layer_245_2040_2050 <- sal_layer_245[[3]]
#do_layer_245_2040_2050 <- do_layer_245[[3]]

## 5-8.5

temp_585 <- "thetao_ssp585_2020_2100_depthmean"
sal_585 <- "so_ssp585_2020_2100_depthmean"
do_585 <- "o2_ssp585_2020_2100_depthmean"
temp_layer_585 <- download_layers(temp_585, temp_variable, constraints)
sal_layer_585 <-download_layers(sal_585, sal_variable, constraints)
#do_layer_585 <-download_layers(do_585, do_variable, constraints)

## 2020-2030
temp_layer_585_2020_2030 <- temp_layer_585[[1]]
sal_layer_585_2020_2030 <- sal_layer_585[[1]]
#do_layer_585_2020_2030 <- do_layer_585[[1]]

## 2030-2040
temp_layer_585_2030_2040 <- temp_layer_585[[2]]
sal_layer_585_2030_2040 <- sal_layer_585[[2]]
#do_layer_585_2030_2040 <- do_layer_585[[2]]

## 2040-2050 
temp_layer_585_2040_2050 <- temp_layer_585[[3]]
sal_layer_585_2040_2050 <- sal_layer_585[[3]]
#do_layer_585_2040_2050 <- do_layer_585[[3]]

## match the freaking GEODATA BS
climate_vars <- c(
  temp_layer_119_2020_2030,
  temp_layer_245_2020_2030,
  temp_layer_585_2020_2030,
  sal_layer_119_2020_2030,
  sal_layer_245_2020_2030,
  sal_layer_585_2020_2030,
  temp_layer_119_2030_2040,
  temp_layer_245_2030_2040,
  temp_layer_585_2030_2040,
  sal_layer_119_2030_2040,
  sal_layer_245_2030_2040,
  sal_layer_585_2030_2040,
  temp_layer_119_2040_2050,
  temp_layer_245_2040_2050,
  temp_layer_585_2040_2050,
  sal_layer_119_2040_2050,
  sal_layer_245_2040_2050,
  sal_layer_585_2040_2050
)

## convert -9999 to NA for Bio-ORACLE layers 
climate_vars <- lapply(climate_vars, function(x) {
  x[x == -9999] <- NA
  return(x)
})

# Then resample, crop, and mask
climate_vars <- lapply(climate_vars, function(x) {
  resample(x, depth, method = "bilinear")
})

## match projection and extents
climate_vars <- lapply(climate_vars, function(x) {
  terra::project(x, "EPSG:32622")})

climate_vars <- lapply(climate_vars, function(x) crop(x, study_area))

## create a mask for all
climate_vars <- lapply(climate_vars, function(x) mask(x, NAFO_vect))

## create stack for all variable rasters
climate_vars <- rast(climate_vars)
names(climate_vars) <-c(
  "temp_119_2020_2030", 
  "temp_245_2020_2030", 
  "temp_585_2020_2030", 
  
  "sal_119_2020_2030", 
  "sal_245_2020_2030", 
  "sal_585_2020_2030", 
  
  "temp_119_2030_2040", 
  "temp_245_2030_2040", 
  "temp_585_2030_2040",
  
  "sal_119_2030_2040", 
  "sal_245_2030_2040", 
  "sal_585_2030_2040", 
  
  "temp_119_2040_2050",  
  "temp_245_2040_2050",  
  "temp_585_2040_2050", 
  
  "sal_119_2040_2050", 
  "sal_245_2040_2050",  
  "sal_585_2040_2050") 

writeRaster(climate_vars, "temp_data/climate_vars.tif", overwrite = T)
climate_vars <- rast("temp_data/climate_vars.tif")
names(climate_vars)
plot(climate_vars)

dir.create("climate_rasters", showWarnings = FALSE)
#####
temp_SSP119_2020_2030 <- climate_vars[[1]]
temp_SSP245_2020_2030 <- climate_vars[[2]]
temp_SSP585_2020_2030 <- climate_vars[[3]]

sal_SSP119_2020_2030 <- climate_vars[[4]]
sal_SSP245_2020_2030 <- climate_vars[[5]]
sal_SSP585_2020_2030 <- climate_vars[[6]]

temp_SSP119_2030_2040 <- climate_vars[[7]]
temp_SSP245_2030_2040 <- climate_vars[[8]]
temp_SSP585_2030_2040 <- climate_vars[[9]]

sal_SSP119_2030_2040 <- climate_vars[[10]]
sal_SSP245_2030_2040 <- climate_vars[[11]]
sal_SSP585_2030_2040 <- climate_vars[[12]]

temp_SSP119_2040_2050 <- climate_vars[[13]]
temp_SSP245_2040_2050 <- climate_vars[[14]]
temp_SSP585_2040_2050 <- climate_vars[[15]]

sal_SSP119_2040_2050 <- climate_vars[[16]]
sal_SSP245_2040_2050 <- climate_vars[[17]]
sal_SSP585_2040_2050 <- climate_vars[[18]]

#####

## make nice 
climate_layers <- list(
  "2020-2030" = list(
    "119" = list(
      temp = "climate_rasters/temp_SSP119_2020_2030.tif",
      sal  = "climate_rasters/sal_SSP119_2020_2030.tif"
    ),
    "245" = list(
      temp = "climate_rasters/temp_SSP245_2020_2030.tif",
      sal  = "climate_rasters/sal_SSP245_2020_2030.tif"
    ),
    "585" = list(
      temp = "climate_rasters/temp_SSP585_2020_2030.tif",
      sal  = "climate_rasters/sal_SSP585_2020_2030.tif"
    )
  ),
  "2030-2040" = list(
   "119" = list(
      temp = "climate_rasters/temp_SSP119_2030_2040.tif",
      sal  = "climate_rasters/sal_SSP119_2030_2040.tif"
    ),
    "245" = list(
      temp = "climate_rasters/temp_SSP245_2030_2040.tif",
      sal  = "climate_rasters/sal_SSP245_2030_2040.tif"
    ),
    "585" = list(
      temp = "climate_rasters/temp_SSP585_2030_2040.tif",
      sal  = "climate_rasters/sal_SSP585_2030_2040.tif"
    )
  ),
  "2040-2050" = list(
    "119" = list(
      temp = "climate_rasters/temp_SSP119_2040_2050.tif",
      sal  = "climate_rasters/sal_SSP119_2040_2050.tif"
    ),
    "245" = list(
      temp = "climate_rasters/temp_SSP245_2040_2050.tif",
      sal  = "climate_rasters/sal_SSP245_2040_2050.tif"
    ),
    "585" = list(
      temp = "climate_rasters/temp_SSP585_2040_2050.tif",
      sal  = "climate_rasters/sal_SSP585_2040_2050.tif"
    )
  )
)
saveRDS(climate_layers, "temp_data/climate_layers.rds")

#####



## Create predictions showing presence probability changing 
## decadally by SSP
years <- 2040:2049
scenarios <- c(119, 245, 585)

get_decade_key <- function(year) {
  if (year >= 2020 && year < 2030) return("2020-2030")
  if (year >= 2030 && year < 2040) return("2030-2040")
  if (year >= 2040 && year < 2050) return("2040-2050")
  stop("Year out of range")
}

jobs <- expand.grid(
  ssp = scenarios,
  yr  = years,
  ss  = seasons,
  stringsAsFactors = FALSE
)
jobs$ssp <- as.character(jobs$ssp)

dir.create("climate_prediction_stacks", showWarnings = FALSE)


for(i in 1:nrow(jobs)) {
  
  ssp <- jobs$ssp[i]
  yr  <- jobs$yr[i]
  ss  <- jobs$ss[i]
  
  cat("Creating raster stack:", ssp, yr, ss, "\n")
  
  ## set season_raster values set for either spring, fall (already existing raster)
  season_raster_temp <- deepcopy(static_stack[[1]])
  values(season_raster_temp) <- ss
  season_raster_temp <- mask(season_raster_temp, NAFO_vect)
  
  ## set year_raster values for the specific year (already existing raster)
  year_raster_temp <- deepcopy(static_stack[[1]])
  values(year_raster_temp) <- yr
  year_raster_temp <- mask(year_raster_temp, NAFO_vect)
  
  ## select correct decadal layer
  dec_key <- get_decade_key(yr)
  temp_layer_use <- rast(climate_layers[[dec_key]][[ssp]]$temp)
  sal_layer_use  <- rast(climate_layers[[dec_key]][[ssp]]$sal)
  
  ## loop combining the new temp, sal, season, and year rasters back into the raster stack
  prediction_stack_ssp_year_season <- c(temp_layer_use, 
                                    sal_layer_use,
                                    static_stack, 
                                    season_raster_temp, 
                                    year_raster_temp)
  
  ## these are the names for the raster layers needed
  ## within each stack for the model to predict correctly 
  names(prediction_stack_ssp_year_season) <- c("temp", 
                                           "sal",
                                           "depth", 
                                           "bpi_fine", 
                                           "bpi_coarse", 
                                           "tri", 
                                           "smf", 
                                           "X",
                                           "Y",
                                           "season",
                                           "year") 
  
  ## store with unique key
  saveRDS(
    prediction_stack_ssp_year_season,
    file = sprintf("climate_prediction_stacks/stack_%s_%s_%s.rds", 
                   gsub(" ", "_", ssp), yr, ss)
  )
  
  rm(prediction_stack_ssp_year_season, season_raster_temp, year_raster_temp)
  gc()
}

## now predict off the raster

dir.create("climate_prediction_frames", showWarnings = FALSE)
dir.create("climate_predictions", showWarnings = FALSE)

foreach(i = 1:nrow(jobs)) %dopar% {
                      
  time_finish_processing <- Sys.time()
  
                      ssp <- jobs$ssp[i]
                      yr  <- jobs$yr[i]
                      ss  <- jobs$ss[i]
    
    cat("\nProcessing:", ssp, yr, ss, "\n")
    
    stack_file <- sprintf("climate_prediction_stacks/stack_%s_%s_%s.rds",ssp, yr, ss)
    prediction_stack_ssp_year_season <- readRDS(stack_file)
    
    cat("Changing stack to frame for:", ssp, yr, ss, "\n")
    
    prediction_frame_year_season <- as.data.frame(prediction_stack_ssp_year_season, na.rm = FALSE)
    
    rm(prediction_stack_ssp_year_season)
    gc()
    
    complete_cases <- complete.cases(prediction_stack_ssp_year_season)
    prediction_frame_ssp_year_season <- prediction_stack_ssp_year_season[complete_cases,] 
    
    ## standardize using means and SD of training data 
    prediction_frame_ssp_year_season_scaled <- prediction_frame_ssp_year_season %>%
      mutate(temp = (temp-temp_mean) / temp_sd,
             sal = (sal-sal_mean )/ sal_sd,
             smf = (smf-smf_mean)/smf_sd,
             bpi_coarse = (bpi_coarse-bpi_coarse_mean)/bpi_coarse_sd,
             bpi_fine = (bpi_fine-bpi_fine_mean)/bpi_fine_sd,
             tri = (tri-tri_mean)/tri_sd)
    
    ## make season a factor just in case it didn't change  
    prediction_frame_ssp_year_season_scaled$season <- factor(prediction_frame_ssp_year_season_scaled$season,
                                                         levels = c("spring", "summer", "fall"),
                                                         labels = c(1,2,3))
    
    cat ("Predicting for:", ssp, yr, ss, "\n")
    
    ## now loop every single stack frame through this predict function
    preds_ssp_year_season <- predict(mod_climate, 
                                 newdata=prediction_frame_ssp_year_season_scaled, 
                                 type = "response")
    
    saveRDS(
      preds_ssp_year_season,
      file = sprintf("climate_predictions/preds_%s_%s_%s.rds", ssp, yr, ss)
    )
    
    rm(prediction_frame_ssp_year_season,
       prediction_frame_ssp_year_season_scaled,
       preds_ssp_year_season)
    
    gc()
    
}

stopCluster(cl)
registerDoSEQ()

## see one raster
stack_119_2045_fall <- readRDS(
  "climate_prediction_stacks/stack_119_2045_fall.rds"
)

stack_119_2045_fall
plot(stack_119_2045_fall)
stack_119_2045_fall <- as.data.frame(stack_119_2045_fall)

## get decadal means
decadal_mean_2020_2029 <- make_decadal_mean_from_preds(
  years = 2020:2029,
  seasons,
  pred_dir = "predictions",
  template_raster = env_stack_all[[1]],
  complete_cases = complete_cases,   
  pred_col = "est",
  out_file = "decadal_mean_2020_2029.tif"
)

plot(decadal_mean_2020_2029)

decades_stack <- c(decadal_mean_2000_2009, decadal_mean_2010_2019)
diff_Raster <- decadal_mean_2010_2019 -decadal_mean_2000_2009
plot(diff_raster)

