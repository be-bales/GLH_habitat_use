## modeling script analyzing trawl data and environmental variable relationships

setwd("D:/BrennaStuff/Github/Turbot-Thesis-Brenna")

## load libraries
library(sf)
library(dplyr)
library(terra)
library(sp)
library(raster)
library(pROC)
library(MultiscaleDTM)
library(batchtools)
library(biooracler)
library(ggplot2)
library(sdmTMB)
library(tidyr)
library(corrplot)
library(stats)
library(performance)

## load trawl dataset
full_database <- read.csv("Data/multispecies_trawl_surveys_1971_2025.csv")

full_database$long.end <- full_database$long.end* -1
full_database$long.start <- full_database$long.start* -1

## select fall only
## full_database <- full_database %>% filter(season == "fall")

## format for presence-absence modeling, 1996 onwards (spatiotemporal GLMM)
pres_abs <- full_database %>%
  dplyr::filter(year >= 2000, year <= 2019) %>% 
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

pres_abs <- pres_abs %>%
  mutate(
    presence = ifelse(number == 0, 0, 1)
  )

get_crs(pres_abs, c("lon", "lat"))
pres_abs <- add_utm_columns(pres_abs, c("lon", "lat"))

## set study area extent and crs
study_area <- extent(-64, -42, 41, 61)

## download bio-ORACLE values for 2000-2020
layers <- list_layers("temp")
info_layer("tas_baseline_2000_2020_depthsurf")
temp_dataset_id <- "thetao_baseline_2000_2019_depthmean"
sal_dataset_id <- "so_baseline_2000_2019_depthmean"
do_dataset_id <- "o2_baseline_2000_2018_depthmean"
terrain_dataset_id <- "terrain_characteristics"

time <- c('2000-01-01T00:00:00Z', '2010-01-01T00:00:00Z')
latitude = c(35,79)
longitude = c(-80,-42)
constraints = list(time, latitude, longitude)
names(constraints) = c("time", "latitude", "longitude")

temp_variable = "thetao_mean"
sal_variable = "so_mean"
do_variable = "o2_mean"
tri_variable = "terrain_ruggedness_index"
slope_variable =  "slope"
depth_variable  = "bathymetry_mean"

temp_layer <- download_layers(temp_dataset_id, temp_variable, constraints)
sal_layer <-download_layers(sal_dataset_id, sal_variable, constraints)
do_layer <-download_layers(do_dataset_id, do_variable, constraints)

plot(temp_layer)

time <- c('1970-01-01T00:00:00Z', '1970-01-01T00:00:00Z')
latitude = c(35,79)
longitude = c(-80,-42)
constraints = list(time, latitude, longitude)
names(constraints) = c("time", "latitude", "longitude")

tri_layer <- download_layers(terrain_dataset_id, tri_variable, constraints)
slope_layer <- download_layers(terrain_dataset_id, slope_variable, constraints)

target_crs <- crs(temp_layer)

## load in files for terrain attributes
depth <- rast("Data/clipped_depth.tif")
depth <- terra::project(depth, target_crs)
depth <- crop(depth, study_area)

depth_mean <- download_layers(terrain_dataset_id, depth_variable, constraints)

## calculate terrain variables 
##bpi2 <- BPI(depth, w = c(25,150), unit = "cell", stand= "none", na.rm = TRUE)
##writeRaster(bpi, "Data/bpi.tif", overwrite=TRUE)
bpi_fine <- rast("BPI_layers/BPI_3x25.tif")
projected_bpi_fine <- terra::project(bpi_fine, target_crs)
bpi_fine <- crop(projected_bpi_fine, study_area)
bpi_coarse <- rast("BPI_layers/BPI_25x150.tif") 
projected_bpi_coarse <- terra::project(bpi_coarse, target_crs)
bpi_coarse <- crop(projected_bpi_coarse, study_area)

## add sediment mobility frequency
smf_comb <- rast("Data/SMF_comb.tif")
smf_comb <- terra::project(smf_comb, target_crs)
smf_comb <- crop(smf_comb, study_area)

## crop terrain vars
slope <- crop(slope_layer, study_area)

tri <- crop(tri_layer, study_area)

## make pres_abs occurrence data spatial 
## using terra, use function vect to make csv dataframe spatial 
v <- vect(pres_abs, geom = c("lon", "lat"), crs= crs(temp_layer))

## extract environmental values at points
temp_vals <- terra::extract(temp_layer, v, ID = FALSE)
pres_abs$temp_2000_2010 <- temp_vals[,1] ## add to database
pres_abs$temp_2010_2020 <- temp_vals[,2]
pres_abs$temp <- ifelse(pres_abs$year <= 2010, 
                        pres_abs$temp_2000_2010, 
                        pres_abs$temp_2010_2020) ## match correct value 

sal_vals <- terra::extract(sal_layer, v, ID = FALSE)
pres_abs$sal_2000_2010 <- sal_vals[,1]
pres_abs$sal_2010_2020 <- sal_vals[,2]
pres_abs$sal <- ifelse(pres_abs$year <= 2010, 
                        pres_abs$sal_2000_2010, 
                        pres_abs$sal_2010_2020) 

do_vals <- terra::extract(do_layer, v, ID = FALSE)
pres_abs$do_2000_2010 <- do_vals[,1]
pres_abs$do_2010_2020 <- do_vals[,2]
pres_abs$do <- ifelse(pres_abs$year <= 2010, 
                       pres_abs$do_2000_2010, 
                       pres_abs$do_2010_2020)

## clean!
pres_abs <- pres_abs |> 
  dplyr::select(-contains("2000_2010"), -contains("2010_2020"))

## extract benthic variables in dataset
v$depth <- terra::extract(depth, v, ID = FALSE)
pres_abs$depth <- v$depth

v$bpi_coarse <- terra::extract(bpi_coarse, v, ID = FALSE)
pres_abs$bpi_coarse <- v$bpi_coarse

v$bpi_fine <- terra::extract(bpi_fine, v, ID = FALSE)
pres_abs$bpi_fine <- v$bpi_fine

v$slope <- terra::extract(slope_layer, v, ID = FALSE)
pres_abs$slope <- v$slope

v$tri <- terra::extract(tri_layer, v, ID = FALSE)
pres_abs$tri <- v$tri

v$smf_comb <- terra::extract(smf_comb, v, ID = FALSE)
pres_abs$smf <- v$smf_comb

## recategorize summer into spring
pres_abs <- pres_abs %>%
  mutate(season = ifelse(season == "summer", "spring", season))

#saveRDS(pres_abs, "Data/pres_abs_extracted.rds")
pres_abs <- readRDS("Data/pres_abs_extracted.rds")

## exploratory analysis

## visualize variable distributions


pres_abs_long <- pres_abs %>%
  pivot_longer(
    cols = c(temp, sal, do, depth, bpi_coarse, bpi_fine, slope, tri, smf),
    names_to = "variable",
    values_to = "value"
  )

distributions <- ggplot(data = pres_abs_long %>%
         filter(presence == 1), aes(x = value)) +
  geom_histogram(binwidth = 0.3) +
  facet_wrap(~ variable, ncol = 2, scales = "free") +
  theme_bw() +
  labs(title = " Covariate Data Distributions")
  

ggsave("Figures/distributions.tif",
       dpi = 300,
       width = 17,
       height = 15,
       device = "tiff",
       compression = "lzw",
       units = "cm",
       bg = "white")
  
## view sampling distribution 

pres_abs$year <- as.numeric(pres_abs$year)
pres_abs$season <- factor(pres_abs$season)
str(pres_abs$season)

## import NAFO shapefile
## add NAFO Divisions shapefile 
NAFO_areas <- st_read("Data/NAFO_areas_ch3.shp")
plot(NAFO_areas)
NAFO_areas <- st_transform(NAFO_areas, crs = 4326)

## NAFO centroids and NAFO base map
NAFO_centroids <- st_centroid(NAFO_areas) 
NAFO_centroids_coords <- cbind(NAFO_centroids, st_coordinates(NAFO_centroids))

## changing to just visualizing sampling season w/ NAFO boundaries

sampling_dist_4 <- ggplot() +
  geom_sf(data = NAFO_areas, fill = "grey95", colour = "black") +
  geom_point(data = pres_abs, aes(x = lon, y = lat, color = season),
             size = 0.1,
             alpha = 0.5) +
  coord_sf(xlim = c(-62, -44), ylim = c(42, 59.2), expand = FALSE) +
  geom_sf_text(data = NAFO_centroids_coords, aes(label = Label), size = 2) + 
  scale_x_continuous(name = "Longitude", breaks = seq(-60, -48, by = 6)) +
  scale_y_continuous(name = "Latitude", breaks = seq(58, 42, by = -6)) +
  facet_wrap(~year) +
  theme_bw() + 
  theme(panel.grid=element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size =10),
        text = element_text(size = 10),
        axis.title.x  = element_text(size = 10),
        axis.title.y  = element_text(size = 10))+
  guides(color = guide_legend(override.aes = list(size = 2))) 


ggsave("Figures/sampling_dist_4.tif",
       dpi = 300,
       width = 18,
       height = 18,
       device = "tiff",
       compression = "lzw",
       units = "cm",
       bg = "white")

## visualize this by depth
sampling_depth <- ggplot(pres_abs, aes(x = month, y = depth, color = season)) +
  geom_point(alpha = 0.5, size = 1) +
  facet_wrap(~year) +
  scale_x_continuous(breaks = c(1,4, 7, 10), 
                     labels = c("Jan", "Apr", "Jul", "Nov")) +
  labs(x = "Month",
       y = "Depth (m)",
       color = "Season") +
  theme_bw() + 
  theme(panel.grid=element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        legend.text = element_text(size =10),
        text = element_text(size = 10),
        axis.title.x  = element_text(size = 10),
        axis.title.y  = element_text(size = 10))+
  guides(color = guide_legend(override.aes = list(size = 2)))

ggsave("Figures/sampling_depth.tif",
       dpi = 300,
       width = 18,
       height = 18,
       device = "tiff",
       compression = "lzw",
       units = "cm",
       bg = "white")


## compare bio-ORACLE and trawl depth and temp data 
pres_abs$trawl_depth <- pres_abs$trawl_depth * -1

surfaces_depth <- pres_abs %>%
  dplyr::select(
    month,
    day,
    year,
    X,
    Y,
    trawl_depth,
    bo_depth = depth,
    presence
  )
surfaces_temp <- pres_abs %>%
  dplyr::select(
    month,
    day,
    year,
    X,
    Y,
    trawl_temp,
    bo_temp = temp,
    presence
  )

surfaces_long_depth <- pivot_longer(
  surfaces_depth,
  cols = c(trawl_depth, bo_depth),
  names_to = "depth vals",
  values_to = "depth value"
)

surfaces_long_temp <- pivot_longer(
  surfaces_temp,
  cols = c(trawl_temp, bo_temp),
  names_to = "temp vals",
  values_to = "temp value"
)


p_scatter_temp <- ggplot(surfaces_temp, aes(x = trawl_temp, y = bo_temp)) +
  geom_point(alpha = 0.5, color = "coral") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_smooth(method = "lm", se = TRUE, color = "steelblue") +
  labs(title = "Temperature: Trawl vs Bio-ORACLE",
       x = "Trawl Temperature (°C)",
       y = "Bio-ORACLE Temperature (°C)") +
  theme_bw()

p_scatter_depth <- ggplot(surfaces_depth, aes(x = trawl_depth, y = bo_depth)) +
  geom_point(alpha = 0.5, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_smooth(method = "lm", se = TRUE, color = "coral") +
  labs(title = "Depth: Trawl vs Bio-ORACLE",
       x = "Trawl Depth (m)",
       y = "Bio-ORACLE Depth (m)") +
  theme_bw()

scatter <- p_scatter_temp | p_scatter_depth ##need to fix text sizes

ggsave("Figures/scatter_comparison.tif", 
       dpi = 300,
       width = 6,
       height = 4,
       device = "tiff",
       compression = "lzw",
       units = "cm",
       bg = "white")

## data summary table 
sampling_summary <- pres_abs %>%
  group_by(year, season) %>%
  summarise(
    zones_covered = paste(sort(unique(NAFOdiv)), collapse = ", "),
    depth_range = paste(round(min(depth, na.rm = TRUE)), 
                        "to", 
                        round(max(depth, na.rm = TRUE)), "m"),
    months = paste(sort(unique(month)), collapse = ", "),
    n_points = n(),
    .groups = "drop"
  ) %>%
  arrange(year, season) %>%
  rename(
    Year = year,
    "Sampling Season" = season,
    "Zones Covered" = zones_covered,
    "Depths Sampled" = depth_range,
    "Months" = months,
    "n Points" = n_points)

# View the table
print(sampling_summary, n = Inf)

# Save as CSV
write.csv(sampling_summary, "Results/sampling_summary_table.csv", row.names = FALSE)
