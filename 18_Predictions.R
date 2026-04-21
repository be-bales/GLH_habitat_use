

##now done with table -- want to create rasters 
temp_1 <- temp_layer[[1]] ##2000-2009
temp_2 <- temp_layer[[2]] #2010-2019

sal_1 <- sal_layer[[1]] ##2000-2009
sal_2 <- sal_layer[[2]] #2010-2019

do_1 <- do_layer[[1]] ##2000-2009
do_2 <- do_layer[[2]] #2010-2019

## give raster stack layers for x and y in lat/long and UTM
## make a copy of one raster 
X_raster <- deepcopy(depth)
Y_raster <- deepcopy(depth)

## assign coords and UTM to new rasters using values function  
coords <- crds(depth, na.rm =F)
coords <- as.data.frame(coords)
coords <- add_utm_columns(coords, 
                          ll_names = c("x", "y"),
                          utm_names = c("X", "Y"),
                          utm_crs = 32622)
values(X_raster) <- coords[,3]
values(Y_raster) <- coords[,4]
names(X_raster) <- "X"
names(Y_raster) <- "Y"

range(Y_raster)

plot(Y_raster)

## match all env data to same extent, resolution, and mask,
## uses 2000-2010 data
vars <- list(
  temp_1,
  temp_2,
  sal_1,
  sal_2,
  do_1,
  do_2,
  depth, 
  bpi_coarse, 
  bpi_fine,
  slope,
  tri_layer,
  smf_comb,
  X_raster,
  Y_raster)

## import NAFO bounds to mask 
NAFO_vect <- vect("Data/NAFO_Projected_ExportFeatures.shp")
NAFO_vect <- terra::project(NAFO_vect, "EPSG:32622")
crs(NAFO_vect)
ext(NAFO_vect)
plot(NAFO_vect)

study_area <- ext(NAFO_vect)

## convert -9999 to NA for Bio-ORACLE layers 
vars_na_fixed <- lapply(vars, function(x) {
  x[x == -9999] <- NA
  return(x)
})

# Then resample, crop, and mask
resampled_vars <- lapply(vars_na_fixed, function(x) {
  resample(x, depth, method = "bilinear")
})

## match projection and extents
extent_vars <- lapply(resampled_vars, function(x) {
  terra::project(x, "EPSG:32622")})

cropped_vars <- lapply(extent_vars, function(x) crop(x, study_area))

## create a mask for all
masked_vars <- lapply(cropped_vars, function(x) mask(x, NAFO_vect))

## create stack for all variable rasters
env_stack_all <- rast(masked_vars)

#writeRaster(env_stack_all, "Data/env_stack_all.tif", overwrite = T)
env_stack_all <- terra::rast("Data/env_stack_all.tif")

plot(env_stack_all)

## make season and years
season_raster <- deepcopy(env_stack_all[[1]])
values(season_raster) <- "spring"
names(season_raster) <- "season"
season_raster <- mask(season_raster, NAFO_vect)

year_raster <- deepcopy(env_stack_all[[1]])
values(year_raster) <- 2000
names(year_raster) <- "year"
year_raster <- mask(year_raster, NAFO_vect)

env_stack_w_season <- c(env_stack_all, season_raster, year_raster)

names(env_stack_w_season) <- c("temp",
                      "temp2",
                      "sal", 
                      "sal2",
                      "do",
                      "do2",
                      "depth",
                      "bpi_fine",
                      "bpi_coarse",
                      "slope",
                      "tri",
                      "smf",
                      "X",
                      "Y",
                      "season",
                      "year")


env_stack_frame <- as.data.frame(env_stack_w_season, na.rm =F)
complete_cases <- complete.cases(env_stack_frame)
env_stack_frame <- env_stack_frame[complete_cases,] 

## standardize using means and SD of training data 
env_stack_frame_scaled <- env_stack_frame %>%
  mutate(temp = (temp-temp_mean) / temp_sd,
         depth = (depth - depth_mean)/ depth_sd,
         sal = (sal-sal_mean )/ sal_sd,
         smf = (smf-smf_mean)/smf_sd,
         bpi_coarse = (bpi_coarse-bpi_coarse_mean)/bpi_coarse_sd,
         bpi_fine = (bpi_fine-bpi_fine_mean)/bpi_fine_sd,
         tri = (tri-tri_mean)/tri_sd)

## season factor 
env_stack_frame_scaled$season <- factor(env_stack_frame_scaled$season,
                                        levels = c("spring", "fall"),
                                        labels = c(1,2))

preds_spring_2000_sdm <- predict(mod_rw_temp, 
                             newdata=env_stack_frame_scaled, 
                             type = "response")

prediction_raster <- deepcopy(env_stack_all[[1]])
values(prediction_raster) <- NA
values(prediction_raster)[complete_cases] <- preds_spring_2000_sdm$epsilon_st
plot(prediction_raster)

# test 
## take one training point 

one_point <- pres_abs_scaled[1,]

one_point$season <- factor(one_point$season,
                                 labels = c("spring","summer", "fall"),
                                 levels  = c(1,2,3))

diff_point <- pres_abs_scaled[20331,]

diff_point$X <- -24.80714 + 5

prediction_point <- predict(mod_rw_temp,
                            newdata = diff_point,
                            type = "response")


## make a static stack
#####
static_vars <- list(
  env_stack_all[[7]], ## depth
  env_stack_all[[8]], ## bpi_coarse
  env_stack_all[[9]], ##bpi_fine 
  env_stack_all[[11]], ## tri
  env_stack_all[[12]], ## smf
  env_stack_all[[13]], ## X
  env_stack_all[[14]]) ## Y

static_stack <- rast(static_vars)
names(static_stack) <- c("depth", 
                         "bpi_coarse", 
                         "bpi_fine", 
                         "tri", 
                         "smf", 
                         "X", 
                         "Y") 
plot(static_stack)
writeRaster(static_stack, "Data/static_stack.tiff")
static_stack <- rast("Data/static_stack.tiff")


layers <- list(
  env_stack_all[[1]], 
  env_stack_all[[3]], 
  env_stack_all[[5]],  
  env_stack_all[[7]], 
  env_stack_all[[8]], 
  env_stack_all[[9]], 
  env_stack_all[[10]],
  env_stack_all[[11]],
  env_stack_all[[12]]) 

layers_stack <- rast(layers)
names(layers_stack) <- c("temperature", 
                         "salinity", 
                         "dissolved oxygen", 
                         "depth", 
                         "BPI coarse", 
                         "BPI fine", 
                         "slope", 
                         "TRI",
                         "SMF") 
plot(layers_stack)

#####

## get decadal layers 
temp_1 <- env_stack_all[[1]]
temp_2 <- env_stack_all[[2]]
sal_1 <- env_stack_all[[3]]
sal_2 <- env_stack_all[[4]]

plot(temp_1)

## make stacks for every season (spring, summer, fall) in every year from 2000-2019
seasons <- c("spring", "fall")
years <- c(2000:2019)

## loop making all stacks for each season, each year 

dir.create("prediction_stacks", showWarnings = FALSE)

for (yr in years) {
  for (ss in seasons) {
    
    cat("Creating raster stack:", yr, ss, "\n")
    
    ## set season_raster values set for either spring, summer, fall (already existing raster)
    season_raster_temp <- deepcopy(static_stack[[1]])
    values(season_raster_temp) <- ss
    season_raster_temp <- mask(season_raster_temp, NAFO_vect)
    
    ## set year_raster values for the specific year (already existing raster)
    year_raster_temp <- deepcopy(static_stack[[1]])
    values(year_raster_temp) <- yr
    year_raster_temp <- mask(year_raster_temp, NAFO_vect)
    
    ## select correct decadal layer
    if (yr <= 2009) {
      temp_layer_use <- temp_1
      sal_layer_use <- sal_1
    } else {
      temp_layer_use <- temp_2
      sal_layer_use <- sal_2
    }
    
    ## loop combining the new temp, sal, season, and year rasters back into the raster stack
    prediction_stack_year_season <- c(temp_layer_use, 
                                      sal_layer_use,
                                      static_stack, 
                                      season_raster_temp, 
                                      year_raster_temp)
    
    ## these are the names for the raster layers needed
    ## within each stack for the model to predict correctly 
    names(prediction_stack_year_season) <- c("temp", 
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
      prediction_stack_year_season,
      file = sprintf("prediction_stacks/stack_%s_%s.rds", yr, ss)
    )
    
    rm(prediction_stack_year_season)
    gc()
  }}

## see one raster
stack_2000_summer <- readRDS(
  "prediction_stacks/stack_2000_summer.rds"
)

stack_2000_summer
plot(stack_2000_summer)
stack_2000_summer <- as.data.frame(stack_2000_summer)

complete_cases <- complete.cases(stack_2000_summer)
saveRDS(complete_cases, "Data/complete_cases.rds")

## loop turning EACH stack into a data frame, subset complete cases,
## standarize using training data, make season a factor, and predict

## make a smaller set
years <- c(2000:2009)

dir.create("prediction_frames", showWarnings = FALSE)
dir.create("predictions", showWarnings = FALSE)

for (yr in years) {
  for (ss in seasons) {
    
    cat("\nProcessing:", yr, ss, "\n")
    
    stack_file <- sprintf("prediction_stacks/stack_%s_%s.rds", yr, ss)
    prediction_stack_year_season <- readRDS(stack_file)
    
    cat("Changing stack to frame for:", yr, ss, "\n")
    
    prediction_frame_year_season <- as.data.frame(prediction_stack_year_season, na.rm = FALSE)
    
    rm(prediction_stack_year_season)
    gc()
    
    complete_cases <- complete.cases(prediction_frame_year_season)
    prediction_frame_year_season <- prediction_frame_year_season[complete_cases,] 
    
    ## standardize using means and SD of training data 
    prediction_frame_year_season_scaled <- prediction_frame_year_season %>%
      mutate(temp = (temp-temp_mean) / temp_sd,
             depth = (depth - depth_mean) / depth_sd,
             sal = (sal-sal_mean )/ sal_sd,
             smf = (smf-smf_mean)/smf_sd,
             bpi_coarse = (bpi_coarse-bpi_coarse_mean)/bpi_coarse_sd,
             bpi_fine = (bpi_fine-bpi_fine_mean)/bpi_fine_sd,
             tri = (tri-tri_mean)/tri_sd)
    
    ## make season a factor just in case it didn't change  
    prediction_frame_year_season_scaled$season <- factor(prediction_frame_year_season_scaled$season,
                                                         levels = c("spring", "fall"),
                                                         labels = c(1,2))
    
    cat ("Predicting for:", yr, ss, "\n")
    
    ## now loop every single stack frame through this predict function
    preds_year_season <- predict(mod_rw_temp, 
                                 newdata=prediction_frame_year_season_scaled, 
                                 type = "response")
    
    saveRDS(
      preds_year_season,
      file = sprintf("predictions/preds_%s_%s.rds", yr, ss)
    )
    
    rm(prediction_frame_year_season,
       prediction_frame_year_season_scaled,
       preds_year_season)
    
    gc()
    
  }
}

 
## select predictions
years <- 2000:2019
seasons <- c("fall", "spring")


preds_2007_spring <- readRDS("predictions/preds_2007_spring.rds")

prediction_raster <- deepcopy(env_stack_all[[1]])

values(prediction_raster) <- NA

values(prediction_raster)[complete_cases] <- preds_2007_spring$est
plot(prediction_raster)

writeRaster(prediction_raster, "prediction_rasters/preds_2007_spring.tif")

## average predictions for spring/summer/fall 
spring_mean_2000_2009 <- make_decadal_mean_from_preds(
  years = 2000:2009,
  seasons = "spring",
  pred_dir = "predictions",
  template_raster = env_stack_all[[1]],
  complete_cases = complete_cases,   
  pred_col = "est",
  out_file = "spring_mean_2000_2009.tif"
)
spring_mean_2010_2019 <- make_decadal_mean_from_preds(
  years = 2010:2019,
  seasons = "spring",
  pred_dir = "predictions",
  template_raster = env_stack_all[[1]],
  complete_cases = complete_cases,   
  pred_col = "est",
  out_file = "spring_mean_2010_2019.tif"
)
fall_mean_2000_2009 <- make_decadal_mean_from_preds(
  years = 2000:2009,
  seasons = "fall",
  pred_dir = "predictions",
  template_raster = env_stack_all[[1]],
  complete_cases = complete_cases,   
  pred_col = "est",
  out_file = "fall_mean_2000_2009.tif"
)
fall_mean_2010_2019 <- make_decadal_mean_from_preds(
  years = 2010:2019,
  seasons = "fall",
  pred_dir = "predictions",
  template_raster = env_stack_all[[1]],
  complete_cases = complete_cases,   
  pred_col = "est",
  out_file = "fall_mean_2010_2019.tif"
)



