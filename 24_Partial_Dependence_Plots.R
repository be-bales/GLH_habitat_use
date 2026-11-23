library(ggeffects)

## sdmtMB
temp_pred <- ggpredict(mod_ices_ns, terms =c("temp[all]"))
saveRDS(temp_pred,"predictions/temp_pred_ggpredict.rds") 
temp_pred <- readRDS("predictions/temp_pred_ggpredict.rds")


depth_pred <-ggpredict(mod_ices_ns, terms =c("depth[all]"))
saveRDS(depth_pred,"predictions/depth_pred_ggpredict.rds" )
depth_pred <- readRDS("predictions/depth_pred_ggpredict.rds")


sal_pred <-ggpredict(mod_ices_ns, terms =c("sal[all]"))
saveRDS(sal_pred,"predictions/sal_pred_ggpredict.rds" )
sal_pred <- readRDS("predictions/sal_pred_ggpredict.rds")


#season_pred <-ggpredict(mod_ices_ns, terms = "season") 
#saveRDS(season_pred,"predictions/season_pred_ggpredict.rds" )
#season_pred <- readRDS("predictions/season_pred_ggpredict.rds")


bpi_c_pred <-ggpredict(mod_ices_ns, terms = "bpi_coarse[all]") 
saveRDS(bpi_c_pred,"predictions/bpi_c_pred_ggpredict.rds" )
bpi_c_pred <- readRDS("predictions/bpi_c_pred_ggpredict.rds")


bpi_f_pred <-ggpredict(mod_ices_ns, terms = "bpi_fine[all]") 
saveRDS(bpi_f_pred,"predictions/bpi_f_pred_ggpredict.rds" )
bpi_f_pred <- readRDS("predictions/bpi_f_pred_ggpredict.rds")


tri_pred <-ggpredict(mod_ices_ns, terms = "tri") 
saveRDS(tri_pred,"predictions/tri_pred_ggpredict.rds" )
tri_pred <- readRDS("predictions/tri_pred_ggpredict.rds")


smf_pred <-ggpredict(mod_ices_ns, terms = "smf") 
saveRDS(smf_pred,"predictions/smf_pred_ggpredict.rds" )
smf_pred <- readRDS("predictions/smf_pred_ggpredict.rds")

## make nice plotting 

plot(sal_pred)
depth_pred_df <- as.data.frame(depth_pred)
depth_pred_df$x_unscaled <- depth_pred_df$x * depth_sd + depth_mean
depth_plot_df <- depth_pred_df %>%
  dplyr::select(
    predicted,
    x_unscaled,
    conf.low,
    conf.high
  )
depth_plot_df$var <- "depth"

temp_pred_df <- as.data.frame(temp_pred)
temp_pred_df$x_unscaled <- temp_pred_df$x * temp_sd + temp_mean
temp_plot_df <- temp_pred_df%>%
  dplyr::select(
    x_unscaled,
    predicted,
    conf.low,
    conf.high
  )
temp_plot_df$var <- "temp"

sal_pred_df <- as.data.frame(sal_pred)
sal_pred_df$x_unscaled <- sal_pred_df$x * sal_sd + sal_mean
sal_plot_df <- sal_pred_df%>%
  dplyr::select(
    x_unscaled,
    predicted,
    conf.low,
    conf.high
  )
sal_plot_df$var <- "sal"

bpi_c_pred_df <- as.data.frame(bpi_c_pred)
bpi_c_pred_df$x_unscaled <- bpi_c_pred_df$x * bpi_coarse_sd + bpi_coarse_mean
bpi_c_plot_df <- bpi_c_pred_df%>%
  dplyr::select(
    x_unscaled,
    predicted,
    conf.low,
    conf.high
  )
bpi_c_plot_df$var <- "bpi_coarse"

tri_pred_df <- as.data.frame(tri_pred)
tri_pred_df$x_unscaled <- tri_pred_df$x * tri_sd + tri_mean
tri_plot_df <- tri_pred_df%>%
  dplyr::select(
    x_unscaled,
    predicted,
    conf.low,
    conf.high
  )
tri_plot_df$var <- "tri"

bpi_f_pred_df <- as.data.frame(bpi_f_pred)
bpi_f_pred_df$x_unscaled <- bpi_f_pred_df$x * bpi_fine_sd + bpi_fine_mean
bpi_f_plot_df <- bpi_f_pred_df%>%
  dplyr::select(
    x_unscaled,
    predicted,
    conf.low,
    conf.high
  )
bpi_f_plot_df$var <- "bpi_fine"

smf_pred_df <- as.data.frame(smf_pred)
smf_pred_df$x_unscaled <- smf_pred_df$x * smf_sd + smf_mean
smf_plot_df <- smf_pred_df%>%
  dplyr::select(
    x_unscaled,
    predicted,
    conf.low,
    conf.high
  )
smf_plot_df$var <- "smf"

season_frame <- as.data.frame(season_pred)
season_plot_df <- season_frame %>%
  dplyr::select(
    predicted,
    x_unscaled = x,
    conf.low,
    conf.high
  )
season_plot_df$var <- "season"

p1 <- ggplot(cont_preds, aes(x = x_unscaled, y = predicted, group = var)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_line() +
  facet_wrap(~var, scales = "free_x", ncol= 1) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 0.5)) +
  
  labs(y = "Predicted response", x = NULL) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        text = element_text(size = 10),
              axis.title.x  = element_text(size = 10),
              axis.title.y  = element_text(size = 10))
ggsave("Figures/pdp_cont.tif",
       dpi = 300,
       width = 6,
       height = 11,
       device = "tiff",
       compression = "lzw",
       units = "cm",
       bg = "white")

p2 <- ggplot(season_preds, aes(x = factor(x_unscaled), y = predicted, group = 1)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_line() +
  geom_point() +
  labs(x = "season", y = "Predicted response") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        text = element_text(size = 10),
        axis.title.x  = element_text(size = 10),
        axis.title.y  = element_text(size = 10))

ggsave("Figures/pdp_disc.tif",
       dpi = 300,
       width = 6,
       height = 4,
       device = "tiff",
       compression = "lzw",
       units = "cm",
       bg = "white")

print(temp_pred, collapse_table = TRUE, collapse_ci = TRUE)
print(season_pred, collapse_table = TRUE, collapse_ci = TRUE)
print(sal_pred, collapse_table = TRUE, collapse_ci = TRUE)

ggplot(temp_pred, aes(x = x, y = predicted, colour = group)) +
  geom_line() 
ggplot(depth_pred, aes(x = x, y = predicted, colour = group)) +
  geom_line() 
ggplot(sal_pred, aes(x = x, y = predicted, colour = group)) +
  geom_line() 

mod_no_smoothers$formula
## trying with a different model
temp_pred_no_smooth <- ggpredict(mod_no_smoothers, terms =c("temp"))
plot(temp_pred_no_smooth)


## combine maxent and sdmTMB plots on same graphs

## maxent 
library(ggeffects)

maxent_mod <- readRDS("Data/maxent_full_model.rds")
str(maxent_mod)

temp_pred <- ggpredict(maxent_mod, terms =c("temp"))
saveRDS(temp_pred/"predictions/temp_pred_ggpredict_maxent.rds" )

depth_pred <-ggpredict(maxent_mod, terms =c("depth"))
saveRDS(depth_pred/"predictions/depth_pred_ggpredict_maxent.rds" )

sal_pred <-ggpredict(maxent_mod, terms =c("sal"))
saveRDS(sal_pred/"predictions/sal_pred_ggpredict_maxent.rds" )

bpi_c_pred <-ggpredict(maxent_mod, terms = "bpi_coarse") 
saveRDS(bpi_c_pred/"predictions/bpi_c_pred_ggpredict_maxent.rds" )

bpi_f_pred <-ggpredict(maxent_mod, terms = "bpi_fine") 
saveRDS(bpi_f_pred/"predictions/bpi_f_pred_ggpredict_maxent.rds" )

tri_pred <-ggpredict(maxent_mod, terms = "tri") 
saveRDS(tri_pred/"predictions/tri_pred_ggpredict_maxent.rds" )

smf_pred <-ggpredict(maxent_mod, terms = "smf") 
saveRDS(smf_pred/"predictions/smf_pred_ggpredict_maxent.rds" )

temp_plot_df$model <- "sdmTMB"
depth_plot_df$model <- "sdmTMB"
sal_plot_df$model   <- "sdmTMB"
tri_plot_df$model   <- "sdmTMB"
smf_plot_df$model   <- "sdmTMB"
bpi_c_plot_df$model   <- "sdmTMB"
bpi_f_plot_df$model   <- "sdmTMB"
all_preds <- rbind(temp_plot_df, 
                   sal_plot_df, 
                   depth_plot_df, 
                   smf_plot_df,
                   tri_plot_df,
                   bpi_c_plot_df,
                   bpi_f_plot_df)

all_preds <- all_preds %>%
  mutate(x_unscaled = as.numeric(x_unscaled))

response_df$model <- "MaxEnt"

response_df <- response_df %>%
  rename(
    x_unscaled = value,
    predicted = suitability
  ) %>%
  mutate(
    conf.low = NA,
    conf.high = NA
  )

combined_preds <- bind_rows(all_preds, response_df)

var_labeller <- as_labeller(c(
  "temp" = "Temperature (°C)" ,
  "depth" = "Depth (m)" ,
  "sal"= "Salinity" ,
  "bpi_coarse" = "BPI coarse" ,
  "bpi_fine" = "BPI fine" ,
  "tri"= "Terrain ruggedness index" ,
  "smf" = "Sediment mobility frequency" ))

ggplot(combined_preds, aes(x = x_unscaled, y = predicted, colour = model, fill = model)) +
   geom_line(size = 1) +
   geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.2,
    colour = NA) +
  facet_wrap(~var, labeller = labeller(var = var_labeller), scales = "free_x", ncol = 2) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) + 
  theme_bw() +
  theme(legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size =9),
        text = element_text(size = 9),
        legend.margin = margin(t = -10),
        legend.position = c(0.7,0.1),
        axis.title.x  = element_text(size = 9),
        axis.title.y  = element_text(size =9)) +
  labs(x = "Predictor value",
    y = "Predicted response",
    colour = "Model",
    fill = "Model") 

ggsave("Figures/pdps.tif",
       dpi = 300,
       height = 13,
       width = 8.3,
       device = "tiff",
       compression = "lzw",
       units = "cm",
       bg = "white")










