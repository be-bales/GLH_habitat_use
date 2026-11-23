## uncertainty 
library(stats)
preds <- readRDS("Data/sdmTMB_preds.rds")

predictions <- preds

mod_anis <- readRDS("Data/mod_anis.rds")
mod_ices <- readRDS("Data/mod_ices.rds")
mod <- readRDS("Data/mod_climate_ices_ns_n_anis.rds")
mod_rw_temp <- readRDS("Data/mod_rw_temp.rds")

future_data <- read.csv("Data/future_data.csv")

summary(mod_ices)

sim_ices <- predict(mod, newdata = pres_abs_ices, nsim = 100)

#saveRDS(sim_ices, "Data/sim_ices.rds")
sim <- readRDS("Data/sim.rds")
sim_ices <- readRDS("Data/sim_ices.rds")

predictions$sim_mean <- exp(rowMeans(sim))
predictions$lwr <- apply(exp(sim), 1, quantile, probs = 0.025)
predictions$upr <- apply(exp(sim), 1, quantile, probs = 0.975)
predictions$sd <- round(apply(exp(sim), 1, function(x) sd(x)), 2)
predictions$cv <- round(apply(exp(sim), 1, function(x) sd(x) / mean(x)), 2)

head(predictions)

ggplot(predictions, aes(X, Y, color = cv)) +
  geom_point(size = 1.5) +
  coord_equal() +
  scale_color_viridis_c(option = "magma") +
  theme_minimal()

ggplot(predictions, aes(X, Y, color = sim_mean)) +
  geom_point(size = 1.5) +
  coord_equal() +
  scale_color_viridis_c(option = "magma") +
  theme_minimal()

ggplot(predictions, aes(x = est, y = cv)) +
  geom_point(alpha = 0.2) +
  geom_smooth() +
  labs(
    x = "Predicted suitability",
    y = "Coefficient of variation",
    title = "Uncertainty vs predicted habitat suitability"
  )

ggplot(predictions, aes(x = factor(year), y = cv)) +
  geom_boxplot() +
  labs(
    x = "Year",
    y = "Coefficient of variation",
    title = "Temporal variation in uncertainty"
  ) +
  theme(axis.text.x = element_text(angle = 90))

ggplot(predictions, aes(x = factor(season), y = cv)) +
  geom_boxplot() +
  labs(
    x = "Season",
    y = "Coefficient of variation",
    title = "Seasonal variation in uncertainty"
  )

## try with future points  

sim_val <- predict(mod, newdata = future_data, nsim = 100)

saveRDS(sim_val, "Data/sim_val_june.rds")

sim_val <- readRDS("Data/sim_val_june.rds")

future_data$sim_val <- exp(rowMeans(sim_val))
future_data$lwr <- apply(exp(sim_val), 1, quantile, probs = 0.025)
future_data$upr <- apply(exp(sim_val), 1, quantile, probs = 0.975)
future_data$sd <- round(apply(exp(sim_val), 1, function(x) sd(x)), 2)
future_data$cv <- round(apply(exp(sim_val), 1, function(x) sd(x) / mean(x)), 2)

write.csv(future_data, "Data/uncertainty_points.csv")

land_shp <- st_read("Data/Land_shp.shp")
land_shp_km <- st_transform(land_shp, crs = 32622 )
st_geometry(land_shp_km) <- st_geometry(land_shp_km) / 1000
st_crs(land_shp_km) <- st_crs(32622)

ggplot() +
  geom_sf(data = land_shp, fill = "grey95", colour = "black") +
  geom_sf(data = NAFO_areas, fill = "grey95", colour = "black") +
  geom_point(data = future_data, aes(x=lon, y=lat, color = cv), size = 1.5, alpha = 0.8) +
  coord_sf(xlim = c(-62, -44), ylim = c(42, 59.2), expand = FALSE) +
  scale_color_viridis_c(name = "Coefficient of\nVariation", option = "magma",
                        guide = guide_colorbar(
                          barheight = unit(2.5, "cm"),
                          barwidth = unit(0.6, "cm"))) +
  theme_bw() + 
  theme(panel.grid=element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        legend.position = c(0.98, 0.98),
        legend.justification = c(1, 1),
        legend.margin = margin(2, 2, 2, 2),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.spacing.y = unit(0.05, "cm"),
        legend.text = element_text(size =10),
        text = element_text(size = 10),
        axis.title  = element_blank(),
        axis.ticks  = element_blank(),
        axis.text = element_blank())

ggsave("Figures/uncertainty_2020s_predictions.tif",
       dpi = 300,
       width = 13,
       height = 15,
       device = "tiff",
       compression = "lzw",
       units = "cm",
       bg = "white")

ggplot(future_data, aes(X, Y, color = sim_val)) +
  geom_point(size = 1.5) +
  coord_equal() +
  scale_color_viridis_c(option = "magma") +
  theme_minimal()

mean(pres_abs_ices$cv)
