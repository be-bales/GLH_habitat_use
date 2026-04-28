## 2020-2025 Model Validation 
## create ground-truth presence data
gt_data <- full_database %>%
  dplyr::filter(year >= 2020, year <= 2025) %>% 
  dplyr::select(
    month,
    day, 
    year,
    trawl_depth = set.depth.mean,
    trawl_temp = bot.temp,
    lon = long.start,
    lat = lat.start,
    number,
    weight, 
    season,
    NAFOdiv)

gt_data <- gt_data %>%
  mutate(
    presence = ifelse(number == 0, 0, 1)
  )

get_crs(gt_data, c("lon", "lat"))
gt_data <- add_utm_columns(gt_data, c("lon", "lat"))

## download bio-ORACLE average values for 2020-2030, SSP 1-1.9
layers <- list_layers("oxygen")
info_layer("thetao_ssp119_2020_2100_depthmean")

temp_119 <- "thetao_ssp119_2020_2100_depthmean"
sal_119 <- "so_ssp119_2020_2100_depthmean"
do_119 <- "o2_ssp119_2020_2100_depthmean"

time <- c("2020-01-01T00:00:00Z", "2090-01-01T00:00:00Z")
latitude = c(35,79)
longitude = c(-80,-42)
constraints = list(time, latitude, longitude)
names(constraints) = c("time", "latitude", "longitude")

## new variables for environment
temp_variable = "thetao_mean"
sal_variable = "so_mean"
do_variable = "o2_mean"

## new future layers
temp_layer_119 <- download_layers(temp_119, temp_variable, constraints)
sal_layer_119 <-download_layers(sal_119, sal_variable, constraints)
do_layer_119 <-download_layers(do_119, do_variable, constraints)

temp_layer_119_2020_2030 <- temp_layer_119[[1]]
sal_layer_119_2020_2030 <- sal_layer_119[[1]]
do_layer_119_2020_2030 <- do_layer_119[[1]]

# Check what years were downloaded
names(temp_layer_119)

## make vector of trawl points
v2 <- vect(gt_data, geom = c("lon", "lat"), crs= crs(temp_layer))

## extract bio-ORACLE data at trawl points
v2$temperature <- terra::extract(temp_layer_119_2020_2030, v2, ID = FALSE)
v2$salinity <- terra::extract(sal_layer_119_2020_2030, v2, ID = FALSE)
v2$depth <- terra::extract(depth, v2, ID = FALSE)
v2$bpi_coarse <- terra::extract(bpi_coarse, v2, ID = FALSE)
v2$bpi_fine <- terra::extract(bpi_fine, v2, ID = FALSE)
v2$tri <- terra::extract(tri_layer, v2, ID = FALSE)
v2$smf_comb <- terra::extract(smf_comb, v2, ID = FALSE)
##append back to dataframe
gt_data$temp <- v2$temperature
gt_data$sal <- v2$salinity
gt_data$depth <- v2$depth
gt_data$bpi_coarse <- v2$bpi_coarse
gt_data$bpi_fine <- v2$bpi_fine 
gt_data$tri <- v2$tri 
gt_data$smf <- v2$smf_comb

## standardize covariates
future_data <- gt_data %>%
  mutate(temp = (temp-temp_mean) / temp_sd,
         depth = (depth-depth_mean) / depth_sd,
         sal = (sal-sal_mean )/ sal_sd,
         smf = (smf-smf_mean)/smf_sd,
         bpi_coarse = (bpi_coarse-bpi_coarse_mean)/bpi_coarse_sd,
         bpi_fine = (bpi_fine-bpi_fine_mean)/bpi_fine_sd,
         tri = (tri-tri_mean)/tri_sd)

future_data <- future_data %>%
  mutate(met_season = case_when(
    month == 1 ~ "winter",
    month == 2 ~ "winter",
    month == 3 ~ "spring",
    month == 4 ~ "spring",
    month == 5 ~ "spring",
    month == 6 ~ "summer",
    month == 7 ~ "summer",
    month == 8 ~ "summer",
    month == 9 ~ "fall",
    month == 10 ~ "fall",
    month == 11 ~ "fall",
    month == 12 ~ "winter",
  ))
future_data$met_season <- factor(future_data$met_season,
                                   levels = c("winter","spring", "summer","fall"),
                                   labels = c(1,2,3,4))

future_data$season <- factor(future_data$season,
                                        levels = c("spring", "fall"),
                                        labels = c(1,2))

## create new model with extra years
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


sanity(mod_climate)

saveRDS(mod_climate, "Data/mod_climate.rds")
mod_climate <- readRDS("Data/mod_climate.rds")

## new model with no season
mod_climate_ices <- sdmTMB(
  data= pres_abs_ices,
  formula = presence ~ s(depth) + s(temp) + s(sal) + s(bpi_coarse) + s(bpi_fine) + tri + smf,
  mesh = mesh3,
  family = binomial(link = "logit"),
  spatial = "on", 
  time = "year",
  spatiotemporal = "RW",
  extra_time = extra_years
)

summary(mod_climate_ices)
sanity(mod_climate_ices)
saveRDS(mod_climate_ices, "Data/mod_climate_ices_ns_n_anis.rds")

## try another model
mod_climate_ices_w_season <- sdmTMB(
  data= pres_abs_ices,
  formula = presence ~ s(depth) + s(temp) + s(sal) + s(bpi_coarse) + s(bpi_fine) + tri + smf + met_season,
  mesh = mesh3,
  family = binomial(link = "logit"),
  spatial = "on", 
  time = "year",
  spatiotemporal = "RW",
  extra_time = extra_years
)

summary(mod_climate_ices_w_season)
sanity(mod_climate_ices_w_season)
saveRDS(mod_climate_ices_w_season, "Data/mod_climate_ices_w_season.rds")
mod_climate_ices <-readRDS("Data/mod_climate_ices_w_season.rds")

## use climate model, trained on original data, to predict presence at real trawl points 
pred_2020_2025_ices <- predict(mod_climate_ices_w_season, 
                          newdata = future_data, 
                          type = "response"
                          )
pred_2020_2025_ices

#saveRDS(pred_2020_2025, "Data/ground_truth_predictions.rds")
#pred_2020_2025 <- readRDS("Data/ground_truth_predictions.rds")

# AUC 
roc_validation <- roc(pred_2020_2025_ices$presence, plogis(pred_2020_2025_ices$est), quiet = TRUE)
auc_validation <- as.numeric(auc(roc_validation))
auc_validation
ci(roc_validation)

plot(roc_validation)

future_data$est <- pred_2020_2025$est


## now try with MaxEnt 
mod_maxent <- readRDS("Data/maxent_model_april.rds")

maxent_pred_2020_2025 <- predict(mod_maxent, gt_data, type = "cloglog")

gt_data$maxent_pred <- maxent_pred_2020_2025[,1]

#saveRDS(maxent_pred_2020_2025, "Data/maxent_pred_2020_2025.rds")
#maxent_pred_2020_2025 <- readRDS("Data/maxent_pred_2020_2025.rds")

# AUC 
roc_validation <- roc(gt_data$presence, gt_data$maxent_pred, quiet = TRUE)
auc_validation <- as.numeric(auc(roc_validation))
print(auc_validation)
ci(roc_validation)

plot(roc_validation)


write.csv(pred_2020_2025, "Data/pred_2020_2025")
