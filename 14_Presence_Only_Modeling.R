## MaxEnt model using maxnet package 
library(pROC)
library(maxnet)

## pull in presence points
maxent_presences <- pres_abs %>%
  filter(presence == 1) %>%
  dplyr::select(lon, lat, year, temp, sal, depth, 
                bpi_coarse, bpi_fine, tri, smf, presence) 

nrow(maxent_presences)

## generate background points on uniform grid within study area
study_area_maxent <- extent(-64, -42, 41, 61)
latlon_extent <- ext(study_area_maxent)
background_resolution_degrees <- 0.25 ## roughly 10x presence points

years <- 2000:2019

background_list <- list()

for(yr in years) {
  
  
  # distinguish correct time period for env. covariates from bio-ORACLE
  if(yr<2010) {time_idx <- 1} else {time_idx <- 2}
  
  #create grid in lat/long 
  background_points <- expand.grid(
    lon = seq(from = xmin(latlon_extent), to = xmax(latlon_extent), by = background_resolution_degrees),
    lat = seq(from = ymin(latlon_extent), to = ymax(latlon_extent), by = background_resolution_degrees))
  
  ## convert to spatial data 
  bg_vect <- terra::vect(background_points, geom = c("lon", "lat"), crs = crs(temp_layer))
  
  background_points$temp <- as.numeric(terra::extract(temp_layer, bg_vect, ID = FALSE)[, time_idx])
  background_points$depth <- as.numeric(terra::extract(depth, bg_vect, ID = FALSE)[, 1])
  background_points$sal <- as.numeric(terra::extract(sal_layer, bg_vect, ID = FALSE)[, time_idx])
  background_points$bpi_coarse <- as.numeric(terra::extract(bpi_coarse, bg_vect, ID = FALSE)[, 1])
  background_points$bpi_fine <- as.numeric(terra::extract(bpi_fine, bg_vect, ID = FALSE)[, 1])
  background_points$tri <- as.numeric(terra::extract(tri, bg_vect, ID = FALSE)[, 1])
  background_points$smf <- as.numeric(terra::extract(smf_comb, bg_vect, ID = FALSE)[, 1])
  
  ## add year and presence  
  background_points$year <- as.numeric(yr)
  background_points$presence <- 0
  
  background_list[[as.character(yr)]] <- background_points
  }
  
## combine all years 

background_points_all <- bind_rows(background_list)
background_points_all <- add_utm_columns(background_points_all, c("lon", "lat"))
background_points_all <- background_points_all %>%
  tidyr::drop_na()

table(background_points_all$year) ## should be ~10x presence points
table(maxent_presences$year)

## bring in presence data
## ensure all pres_abs variables are numeric

maxent_data <- bind_rows(maxent_presences, background_points_all)

## define predictors for model formula 

pred_vars <- c( "temp", "sal", "depth", "bpi_coarse", "bpi_fine",
               "tri", "smf")

# Verify all predictors are numeric vectors
print(sapply(maxent_data[, pred_vars], class))

## fit model 

mod_maxent <- maxnet(
  p = maxent_data$presence,
  data = maxent_data[, pred_vars]
)

#saveRDS(mod_maxent, "Data/maxent_full_model.rds")
#mod_maxent <- readRDS("Data/maxent_full_model.rds")

saveRDS(mod_maxent, "Data/maxent_model_april.rds")
mod_maxent <- readRDS("Data/maxent_model_april.rds")

"Data/maxent_model_april.rds"

summary(mod_maxent)

## make predictions
maxent_predictions <- predict(mod_maxent, maxent_data[, pred_vars], type = "cloglog")

#saveRDS(maxent_predictions, "Data/maxent_predictions.rds")
#maxent_predictions <- readRDS("Data/maxent_predictions.rds")

maxent_data$preds <- maxent_predictions[,1]
#write.csv(maxent_data, "Data/maxent_data_with_preds.csv")

## calculate overall auc
?pROC::roc()
roc <- roc(maxent_data$presence, maxent_data$preds, quiet = TRUE)
auc <- as.numeric(auc(roc))
auc
auc_full <- auc

## export prediction raster
prediction_stack <- readRDS("prediction_stacks/stack_2008.rds")

prediction_frame <- as.data.frame(prediction_stack, na.rm = FALSE)

complete_cases <- complete.cases(prediction_frame)
prediction_frame <- prediction_frame[complete_cases,] 

maxent_predictions_new <- predict(mod_maxent, newdata=prediction_frame, 
                                  type = "cloglog")

prediction_raster <- deepcopy(prediction_stack[[1]])
values(prediction_raster) <- NA
values(prediction_raster)[complete_cases] <- maxent_predictions_new
plot(prediction_raster)

writeRaster(prediction_raster, "prediction_rasters/maxent_preds_2000_2009.tiff")

## permutation importance - 200 permutations
pred_vars <- c( "temp", "sal", "depth", "bpi_coarse", "bpi_fine",
                "tri", "smf")

## for perm. imp

n_perms <- 200

perm_importance <- lapply(pred_vars, function(v) {
  
  delta_aucs <- replicate(n_perms, {
    
    # Copy the full dataset
    perm_data <- maxent_data ##needs to be unscaled data
    
    # Permute only the predictor
    perm_data[[v]] <- sample(perm_data[[v]])
    
    # Predict on full dataset
    perm_out <- predict(mod_maxent,
                        newdata = perm_data,
                        type = "cloglog")
    
    # Extract predictions 
    perm_preds <- perm_out[,1]
    
    auc_perm <- as.numeric(auc(maxent_data$presence, perm_preds))
    auc_full - auc_perm
  })
  
  data.frame(
    variable = v,
    delta_auc = mean(delta_aucs),
    sd_delta_auc = sd(delta_aucs)
  )
})

perm_importance <- do.call(rbind, perm_importance) %>%
  arrange(desc(delta_auc))
perm_importance$delta_auc <- round(perm_importance$delta_auc, 4)
perm_importance


### partial response curves

plot(
  mod_maxent,
  vars = names(mod_maxent$samplemeans),
  common.scale = TRUE,
  type = c("cloglog"),
  plot = TRUE
)

plot(
  mod_maxent,
  vars = names(mod_maxent$samplemeans),
  common.scale = TRUE,
  type = c("link"),
  plot = TRUE
)

## now make these with ggplot for publication chapter prettyness

library(dplyr)

pred_vars <- names(mod_maxent$samplemeans)

response_df <- lapply(pred_vars, function(v) {
  
  # Sequence over observed range
  x_seq <- seq(
    min(maxent_presences[[v]], na.rm = TRUE),
    max(maxent_presences[[v]], na.rm = TRUE),
    length.out = 500)
  
  # Start with means of all predictors
  newdata <- as.data.frame(as.list(mod_maxent$samplemeans))
  
  # Repeat for each x value
  newdata <- newdata[rep(1, length(x_seq)), ]
  
  # Replace variable with sequence
  newdata[[v]] <- x_seq
  
  # Predict 
  preds <- predict(mod_maxent, newdata, type = "cloglog")
  
  data.frame(
    var = v,
    value = x_seq,
    suitability = as.numeric(preds))
})

response_df <- bind_rows(response_df)

pred_vars <- c( "temp", "sal", "depth", "bpi_coarse", "bpi_fine",
                "tri", "smf")

response_df <- response_df %>%
  rename("temperature (°C)" = temp,
         "depth (m)" = depth,
         "salinity" = sal,
         "BPI coarse" = bpi_coarse,
         "BPI fine" = bpi_fine,
         "terrain ruggedness index" = tri,
         "sediment mobility frequency" = smf)

ggplot(response_df,
       aes(x = value, y = suitability)) +
  geom_line(size = 1) +
  facet_wrap(~ variable, scales = "free") +
  theme_bw() +
  labs(
    x = "Predictor value",
    y = "Predicted suitability (cloglog)",
    title = "MaxEnt Response Curves"
  )


