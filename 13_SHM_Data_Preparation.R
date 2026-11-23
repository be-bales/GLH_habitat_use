## modeling script analyzing trawl data and environmental variable relationships

setwd("C:/Users/befat/Desktop/GLH")

remotes::install_github("bio-oracle/biooracler")

install.packages("colorspace")

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
library(colorspace)

## load trawl dataset
full_database <- read.csv("Data/multispecies_trawl_surveys_1971_2025.csv")

full_database$long.end <- full_database$long.end* -1
full_database$long.start <- full_database$long.start* -1

#write.csv(full_database, "Data/trawlpoints.csv")

## select fall only
## full_database <- full_database %>% filter(season == "fall")

## format for presence-absence modeling
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
pres_abs <- pres_abs %>%
  tidyr::drop_na()

## exploratory analysis

## visualize variable distributions
pres_abs_long <- pres_abs %>%
  pivot_longer(
    cols = c(temp, sal, do, depth, bpi_coarse, bpi_fine, slope, tri, smf),
    names_to = "variable",
    values_to = "value"
  )

ggplot(data = pres_abs_long %>%
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
NAFO_areas <- st_transform(NAFO_areas, crs = 4326)

## NAFO centroids and NAFO base map
NAFO_centroids <- st_centroid(NAFO_areas) 
NAFO_centroids_coords <- cbind(NAFO_centroids, st_coordinates(NAFO_centroids))

## import land
land_shp <- st_read("Data/Land_shp.shp")
land_shp_km <- st_transform(land_shp, crs = 32622 )
st_geometry(land_shp_km) <- st_geometry(land_shp_km) / 1000
st_crs(land_shp_km) <- st_crs(32622)

## changing to just visualizing sampling season w/ NAFO boundaries

ggplot() +
  geom_sf(data = NAFO_areas, fill = "grey95", colour = "grey30") +
  geom_point(data = pres_abs, aes(x = lon, y = lat, color = season),
             size = 0.3,
             alpha = 0.5) +
  geom_sf(data = land_shp, fill = "grey95", colour = "grey30") +
  coord_sf(xlim = c(-61, -46), ylim = c(42, 58.5), expand = FALSE) +
  scale_color_manual(values = c("#F8766D", "#00BFC4"),labels = c("Fall", "Spring"),
                     name = "Season") +
  facet_wrap(~year) +
  theme_bw() + 
  theme(panel.grid=element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        panel.spacing = unit(0.1, "cm"),
        legend.key.size = unit(0.5, 'cm'),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
        legend.text = element_text(size =10),
        legend.position = "bottom",
        text = element_text(size = 10),
        strip.text.x = element_text(size =10),
        axis.title  = element_blank(),
        axis.ticks  = element_blank(),
        axis.text = element_blank(),
        legend.box.spacing = unit(0.1,'cm')) +
  guides(color = guide_legend(override.aes = list(size = 2)))


ggsave("Figures/sampling_dist.tif",
       dpi = 300,
       width = 13,
       height = 15,
       device = "tiff",
       compression = "lzw",
       units = "cm",
       bg = "white")

## visualize this by depth
ggplot(pres_abs, aes(x = month, y = depth, color = season)) +
  scale_color_manual(values = c("#F8766D", "#00BFC4"),labels = c("Fall", "Spring"),
                     name = "Season") +
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
        legend.text = element_text(size =9),
        text = element_text(size = 9),
        axis.title.x  = element_text(size = 10),
        axis.title.y  = element_text(size = 10))+
  guides(color = guide_legend(override.aes = list(size = 2)))

ggsave("Figures/sampling_depth.tif",
       dpi = 300,
       width = 15,
       height = 13,
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


ggplot(surfaces_temp, aes(x = trawl_temp, y = bo_temp)) +
  geom_point(alpha = 0.5, color = "coral") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_smooth(method = "lm", se = TRUE, color = "steelblue") +
  labs(
       x = "Trawl temperature (°C)",
       y = "Bio-ORACLE temperature (°C)") +
  theme_bw() +
  theme(
        text = element_text(size = 10),
        axis.title.x  = element_text(size = 10),
        axis.title.y  = element_text(size = 10))

ggplot(surfaces_depth, aes(x = trawl_depth, y = bo_depth)) +
  geom_point(alpha = 0.5, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_smooth(method = "lm", se = TRUE, color = "coral") +
  labs(
       x = "Trawl depth (m)",
       y = "Bio-ORACLE depth (m)") +
  theme_bw()+
theme(
  text = element_text(size = 10),
  axis.title.x  = element_text(size = 10),
  axis.title.y  = element_text(size = 10))


ggsave("Figures/scatter_comparison_depth.tif", 
       dpi = 300,
       width = 10,
       height = 9,
       device = "tiff",
       compression = "lzw",
       units = "cm",
       bg = "white")


## differences map
pres_abs <- pres_abs %>%
  mutate(temp_diff = (abs(trawl_temp) - abs(temp))) ##trawl minus bio-oracle

ggplot(pres_abs, aes(temp_diff)) +
  geom_histogram()

temps <- pres_abs %>%
  summarize(percent_above = sum(temp_diff > 2)/n() * 100,
             percent_below = sum(temp_diff < 2)/n() * 100)

ggplot() +
  geom_sf(data = NAFO_areas, fill = "grey90", color = "black") +  
  geom_sf(data = land_shp, fill = "grey75", colour = "grey30") +
  geom_point(data = pres_abs, aes(x = lon, y = lat, color = temp_diff),
            size = 0.3, alpha = 0.9) +
  scale_color_continuous_diverging(name = "Temperature\ndifferential (°C)", 
                                   palette = "Blue-Red 3") +
  coord_sf(xlim = c(-62, -44), ylim = c(42, 59), expand = FALSE) +
  facet_wrap(~year)+
  scale_x_continuous(name = "Longitude", 
                     breaks = c(-59,  -48)) +
  scale_y_continuous(name = "Latitude", 
                     breaks = c(58, 51, 43)) +
  theme_bw() + 
  theme(panel.grid=element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        text = element_text(size = 10),
        legend.key.size = unit(0.3, 'cm'),
        legend.key.height = unit(0.4, "cm"),
        legend.key.spacing.y = unit(0.1, 'cm'),
        legend.text = element_text(size =10),
        strip.text.x = element_text(size =10),
        axis.title.x  = element_text(size = 10),
        axis.title.y  = element_text(size = 10))

ggsave("Figures/spatial_temp_difference.tif", 
       dpi = 300,
       width = 18,
       height = 16,
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

write.csv(sampling_summary, "Data/sampling_summary_table.csv")
