## sdmTMB presence-absence models

## Assess pearson correlation coefs and VIF values for predictors 

predictors <- pres_abs %>%
  dplyr::select(
    temp,
    sal,
    #do,
    depth,
    bpi_coarse,
    bpi_fine,
    #slope, ## removing from model due to high VIF
    tri, 
    smf
  ) 

cor_matrix <- cor(predictors, use = "complete.obs", method="pearson")
cor_matrix
corrplot(cor_matrix, method = "color", type = "upper", t1.cex = 0.8)
vif_results <- usdm::vif(predictors)
print(vif_results)

#write.csv(vif_results, "Data/vif_results_1.csv")
#write.csv(vif_results, "Data/vif_results_2.csv")
#write.csv(cor_matrix, "Data/cor_matrix.csv")

## standardize covariates
## scale data 
means <- apply(pres_abs[,c("temp", "sal", "smf", "tri", "bpi_coarse", "bpi_fine", "depth")], 2, mean, na.rm=T)
sds <- apply(pres_abs[,c("temp", "sal", "smf", "tri", "bpi_coarse", "bpi_fine", "depth")], 2, sd, na.rm=T)
scaled_data <- data.frame(
  names=names(pres_abs[,c("temp", "sal", "smf", "tri", "bpi_coarse", "bpi_fine", "depth")]),
  means = means,
  sds = sds
)
temp_mean <- scaled_data[1,2]
temp_sd <- scaled_data[1,3]

sal_mean <- scaled_data[2,2]
sal_sd <- scaled_data[2,3]

smf_mean <- scaled_data[3,2]
smf_sd <- scaled_data[3,3]

tri_mean <- scaled_data[4,2]
tri_sd <- scaled_data[4,3]

bpi_coarse_mean <- scaled_data[5,2]
bpi_coarse_sd <- scaled_data[5,3]

bpi_fine_mean <- scaled_data[6,2]
bpi_fine_sd <- scaled_data[6,3]

depth_mean <- scaled_data[7,2]
depth_sd <- scaled_data[7,3]

pres_abs_scaled <- pres_abs %>%
  mutate(depth = (depth-depth_mean) / depth_sd,
         temp = (temp-temp_mean) / temp_sd,
         sal = (sal-sal_mean )/ sal_sd,
         smf = (smf-smf_mean)/smf_sd,
         bpi_coarse = (bpi_coarse-bpi_coarse_mean)/bpi_coarse_sd,
         bpi_fine = (bpi_fine-bpi_fine_mean)/bpi_fine_sd,
         tri = (tri-tri_mean)/tri_sd) %>%
  filter(!is.na(temp))

str(pres_abs_scaled$season)
pres_abs_scaled$season <- factor(pres_abs_scaled$season,
                                 levels = c("spring", "fall"),
                                 labels = c(1,2))

pres_abs_scaled$month <- as.character(pres_abs_scaled$month)

## create meterological season variable
pres_abs_ices <- pres_abs_scaled %>%
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
pres_abs_ices$met_season <- factor(pres_abs_ices$met_season,
                                 levels = c("winter","spring", "summer","fall"),
                                 labels = c(1,2,3,4))

## make SPDE mesh 
mesh3 <- make_mesh(pres_abs_ices, 
                   c("X", "Y"), 
                   fmesher_func = fmesher::fm_mesh_2d_inla,
                   cutoff = 25,
                   max.edge = c(200, 100),
                   offset = c(30, 100) 
                   )
plot(mesh3)
mesh3$mesh$n

## model selection ##

## all covariates
mod_rw <- sdmTMB(
  data= pres_abs,
  formula = presence ~ depth + temp + sal + do + bpi_coarse + bpi_fine + rugosity + smf + as.factor(season),
  mesh = mesh3,
  family = binomial(link = "logit"),
  spatial = "on", 
  time = "year",
  spatiotemporal = "RW"
)

AIC(mod_rw)
sanity(mod_rw)
summary(mod_rw)
tidy(mod_rw, conf.int = TRUE)
tidy(mod_rw, "ran_pars", conf.int = TRUE)

## removing highly correlated DO, adding smoothers 
mod_rw_temp <- sdmTMB(
  data= pres_abs_scaled,
  formula = presence ~ s(depth) + s(temp) + s(sal) + s(bpi_coarse) + 
    s(bpi_fine) + tri + smf + season,
  mesh = mesh3,
  family = binomial(link = "logit"),
  spatial = "on", 
  time = "year",
  spatiotemporal = "RW"
)


AIC(mod_rw_temp)
sanity(mod_rw_temp)
summary(mod_rw_temp)
#saveRDS(mod_rw_temp,"Data/mod_rw_temp.rds")
#mod_rw_temp <- readRDS("Data/mod_rw_temp.rds")

## no smoothers on interpretable terms, no DO
mod_no_season <- sdmTMB(
  data= pres_abs,
  formula = presence ~ depth + temp + sal + bpi_coarse + 
    bpi_fine + tri + smf,
  mesh = mesh3,
  family = binomial(link = "logit"),
  spatial = "on", 
  time = "year",
  spatiotemporal = "RW"
)

AIC(mod_no_season)
sanity(mod_no_season)
summary(mod_no_season)

## AR1 model
pres_abs_ar1 <- pres_abs_scaled
pres_abs_ar1$season <- as.integer(pres_abs_ar1$season
)
mod_ar1 <- sdmTMB(
  data= pres_abs_ar1,
  formula = presence ~ s(depth) + s(temp) + s(sal) + s(bpi_coarse) + 
    s(bpi_fine) + tri + smf + season,
  mesh = mesh3,
  family = binomial(link = "logit"),
  spatial = "on", 
  time = "year",
  spatiotemporal = "AR1"
)

AIC(mod_ar1)
sanity(mod_ar1)
summary(mod_ar1)

## ICES paper model w/ meterological season
mod_ices <- sdmTMB(
  data= pres_abs_ices,
  formula = presence ~ s(depth) + s(temp) + s(sal) + s(bpi_coarse) + 
    s(bpi_fine) + tri + smf + met_season,
  mesh = mesh3,
  family = binomial(link = "logit"),
  spatial = "on", 
  time = "year",
  spatiotemporal = "RW"
)

AIC(mod_ices) 
sanity(mod_ices)
summary(mod_ices)
#saveRDS(mod_ices,"Data/mod_ices.rds")
mod_ices <- readRDS("Data/mod_ices.rds")
mod1 <- tidy(mod_ices, conf.int = TRUE)
write.csv(mod1, "Data/mod_export.csv")
mod2 <- tidy(mod_ices, "ran_pars", conf.int = TRUE)
write.csv(mod2, "Data/mod_export_2.csv")

## ICES model w/o season
mod_ices_no_season <- sdmTMB(
  data= pres_abs_ices,
  formula = presence ~ s(depth) + s(temp) + s(sal) + s(bpi_coarse) + 
    s(bpi_fine) + tri + smf,
  mesh = mesh3,
  family = binomial(link = "logit"),
  spatial = "on", 
  time = "year",
  spatiotemporal = "RW"
)
AIC(mod_ices_no_season)
sanity(mod_ices_no_season)
summary(mod_ices_no_season)
#saveRDS(mod_ices_no_season,"Data/mod_ices_no_season.rds")
mod_ices_ns <- readRDS("Data/mod_ices_no_season.rds")
cAIC(mod_ices_ns, what = "cAIC")

## ICES model no smoothers
mod_ices_simple <- sdmTMB(
  data= pres_abs_ices,
  formula = presence ~ depth + temp + sal + met_season,
  mesh = mesh3,
  family = binomial(link = "logit"),
  spatial = "on", 
  time = "year",
  spatiotemporal = "RW"
)
AIC(mod_ices_simple)
sanity(mod_ices_simple)
summary(mod_ices_simple)
#saveRDS(mod_ices_simple,"Data/mod_ices_simple.rds")

## Model w/ anisotropy 
mod_anis <- sdmTMB(
  data= pres_abs_ices,
  formula = presence ~ s(depth) + s(temp) + s(sal) + met_season,
  mesh = mesh3,
  family =binomial(link = "logit"),
  spatial = "on", 
  time = "year",
  spatiotemporal = "rw",
  anisotropy = TRUE)

AIC(mod_anis)
sanity(mod_anis)
summary(mod_anis)
#saveRDS(mod_anis,"Data/mod_anis.rds")

## bivariate temp/depth and anisotropy
mod_anis_bi <- sdmTMB(
  data= pres_abs_ices,
  formula = presence ~ s(depth, temp),
  mesh = mesh3,
  family =binomial(link = "logit"),
  spatial = "on", 
  time = "year",
  spatiotemporal = "rw",
  anisotropy = TRUE,
  silent = FALSE
)
sdmTMB::cAIC(mod_anis_bi)
sanity(mod_anis_bi)
summary(mod_anis_bi)
#saveRDS(mod_anis_bi,"Data/mod_anis_bi.rds")

## now experiment with season as the time and spatiotemporal as ar1
mod_ar1_season <- sdmTMB(
  data= pres_abs_ices,
  formula = presence ~ s(depth) + s(temp) + s(sal),
  mesh = mesh3,
  family =binomial(link = "logit"),
  spatial = "on", 
  time = "season",
  spatiotemporal = "ar1",
  anisotropy = TRUE,
  silent = FALSE
)

AIC(mod_ar1_season)
sanity(mod_ar1_season)
summary(mod_ar1_season)
#saveRDS(mod_ar1_season,"Data/mod_ar1_season.rds")

## bivariate smoothers + BPI
mod_bi_2 <- sdmTMB(
  data= pres_abs_ices,
  formula = presence ~ s(depth, temp) + bpi_coarse,
  mesh = mesh3,
  family =binomial(link = "logit"),
  spatial = "on", 
  time = "year",
  spatiotemporal = "rw"
)
AIC(mod_bi)
sanity(mod_bi_2)
summary(mod_bi_2)
#saveRDS(mod_bi,"Data/mod_bi.rds")

## season as a random intercept  
mod_int <- sdmTMB(
  data= pres_abs_ices,
  formula = presence ~ s(depth) + s(temp) + s(sal) + s(bpi_coarse) + 
    s(bpi_fine) + tri + smf + (1|season),
  mesh = mesh3,
  family = binomial(link = "logit"),
  spatial = "on", 
  time = "year",
  spatiotemporal = "RW"
)

AIC(mod_int)
sanity(mod_int)
summary(mod_int)
cAIC(mod_int)

## use conditional AIC instead of marginal AIC for models with smoothers
sdmTMB::cAIC(mod_rw_temp)
sdmTMB::cAIC(mod_ices)
sdmTMB::cAIC(mod_bi)
sdmTMB::cAIC(mod_anis)

plot_anisotropy(mod_ar_2)

# Get predictions and residuals 
preds <- predict(mod_ices_ns, type = "response")

roc_full <- roc(pres_abs_scaled$presence, plogis(preds$est), quiet = TRUE)
auc_full <- as.numeric(auc(roc_full))
print(auc_full)
auc_validation <- 0.963864 ## mod_iceS_ns

saveRDS(preds, "Data/sdmTMB_preds_ices.rds")
resids_mle <- residuals(mod_ices, type = "mle-mvn")
qqnorm(resids_mle)
abline(0, 1)
qqline(resids_mle)

tiff(file="QQ_plot.tif",
     width=3, height=2, units="in", res=100)
dev.off()

#preds <- readRDS("Data/sdmTMB_preds.rds")

## view confidence intervals on fixed effects, random effect, and 
## variance parameters

model_intercepts <- tidy(mod_ices, conf.int = TRUE)
write.csv(model_intercepts, "Data/model_intercepts_ices.csv")
ran_pars <- tidy(mod_ices, "ran_pars", conf.int = TRUE)
write.csv(ran_pars, "Data/ran_pars_ices.csv")
                                                                                                                                                                                                                                                                                                                                           
## plot residuals on map
mod_resids <- cbind(mod_ices$data, resids_mle)

## import land
land_shp <- st_read("Data/Land_shp.shp")
land_shp_km <- st_transform(land_shp, crs = 32622 )
st_geometry(land_shp_km) <- st_geometry(land_shp_km) / 1000
st_crs(land_shp_km) <- st_crs(32622)

## convert to spatial
mod_resids_df <- st_as_sf(mod_resids, coords = c("X", "Y"), crs = 32622)
plot(mod_resids_df)

model_residuals_plot_alltime <- ggplot(mod_resids, aes(X, Y, col = resids_mle)) + 
  scale_colour_gradient2(name ="MLE Residuals") +
  geom_point(size = 0.4, alpha = 0.8) +
  geom_sf(data = land_shp_km, fill = "grey90", color = "black", inherit.aes = F) +
  coord_sf(crs =st_crs(32622), ylim = c(4700,6500), xlim = c(-250, 1000)) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid=element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size =8),
        text = element_text(size = 8),
        axis.title.x  = element_text(size = 8),
        axis.title.y  = element_text(size = 8)) +
  labs(title = "Model Residuals 2000-2019")

ggsave("Figures/model_residuals_ices.tif",
       plot = model_residuals_plot_alltime,
       width = 13,
       height = 15,
       units = "cm",
       dpi = 300,
       device = "tiff",
       compression = "lzw",
       bg = "white")


##TSS
preds$pred_01 <- ifelse(preds$est < 0.5, 0, 1)

conmat <- table(preds$pred_01, preds$presence)
true_neg <- conmat[1, 1]
false_neg <- conmat[1, 2]
false_pos <- conmat[2, 1]
true_pos <- conmat[2, 2]

# Calculate TSS:
true_pos_rate <- true_pos / (true_pos + false_neg)
true_neg_rate <- true_neg / (true_neg + false_pos)
TSS <- true_pos_rate + true_neg_rate - 1
TSS 

pred_vars <- c( "temp", "sal", "depth", "bpi_coarse", "bpi_fine",
                "tri", "smf")

n_perms <- 200

perm_importance <- lapply(pred_vars, function(v) {
  
  delta_aucs <- replicate(n_perms, {
    
    # Copy the full dataset including X/Y
    perm_data <- pres_abs_scaled
    
    # Permute only the predictor
    perm_data[[v]] <- sample(perm_data[[v]])
    
    # Predict on full dataset
    perm_out <- predict(mod_ices_ns,
                        newdata = perm_data,
                        type = "response")
    
    # Extract predictions (numeric vector)
    perm_preds <- perm_out$est

    auc_perm <- as.numeric(auc(pres_abs_scaled$presence, perm_preds))
    auc_validation - auc_perm
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

write.csv(perm_importance, "Data/perm_imp_april_sdm.csv")

saveRDS(perm_importance, "Data/perm_importance_ices.rds")