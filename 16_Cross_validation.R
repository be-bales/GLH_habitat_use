## 10-fold cross validation by time period 

## maxent cross validation
maxent_cv_data <- maxent_data %>%
  mutate(
    fold_id = case_when(
      year <= 2001 ~ 1,
      year <= 2003 ~ 2,
      year <= 2005 ~ 3,
      year <= 2007 ~ 4,
      year <= 2009 ~ 5,
      year <= 2011 ~ 6,
      year <= 2013 ~ 7,
      year <= 2015 ~ 8,
      year <= 2017 ~ 9,
      year <= 2019 ~ 10
    )
  )

## check folds by year 

print(table(maxent_cv_data$fold_id, maxent_cv_data$year))

cv_results <- data.frame(
  fold = 1:10,
  n_train = NA_integer_,
  n_test = NA_integer_,
  n_test_presence = NA_integer_,
  auc = NA_real_
)

## store predictions 

all_predictions <- rep(NA_real_, nrow(maxent_cv_data))

## loop CV function over each testing year 

for(i in 1:10) {
  
  cat("Processing fold: ", i , "\n")
  
  ## split data -- train on all years not i and test on year i 
  training_data <- maxent_cv_data %>% filter(fold_id != i)
  testing_data <- maxent_cv_data %>% filter(fold_id == i)
  
  cv_results$n_train[i] <- nrow(training_data)
  cv_results$n_test[i] <- nrow(testing_data)
  cv_results$n_test_presence[i] <- sum(testing_data$presence == 1)
  
  ## refit model on only training data
  mod_maxent_i <- maxnet(
    p = training_data$presence,
    data = training_data[, pred_vars]
  )
  
  ## predict on testing data 
  pred_i <- predict(mod_maxent_i, testing_data[, pred_vars], type = "cloglog")
  
  ## store predictions 
  all_predictions[maxent_cv_data$fold_id == i] <- pred_i
  
  ## calculate auc for the fold 
  roc_i <- roc(testing_data$presence, pred_i, quiet = TRUE)
  cv_results$auc[i] <- as.numeric(auc(roc_i))
  
  cat("Fold AUC:", round(cv_results$auc[i], 4), "\n")
   
}

print(cv_results)

maxent_cv_data$cv_predicted <- all_predictions

saveRDS(maxent_cv_data, "Data/maxent_data_with_cv_predictions.rds")
saveRDS(cv_results, "Data/maxent_cv_results.rds")
maxent_cv_results <- readRDS("Data/maxent_cv_results.rds")
write.csv(maxent_cv_results, "Data/maxent_cv_results.csv")

## sdmTMB cv 
#####
## cross validation data
sdm_cv_data <- pres_abs_ices %>%
  mutate(
    fold_id = case_when(
      year <= 2001 ~ 1,
      year <= 2003 ~ 2,
      year <= 2005 ~ 3,
      year <= 2007 ~ 4,
      year <= 2009 ~ 5,
      year <= 2011 ~ 6,
      year <= 2013 ~ 7,
      year <= 2015 ~ 8,
      year <= 2017 ~ 9,
      year <= 2019 ~ 10
    )
  )

sdm_cv_data <- na.omit(sdm_cv_data)
str(sdm_cv_data)

mesh4 <- make_mesh(sdm_cv_data, 
                   c("X", "Y"), 
                   fmesher_func = fmesher::fm_mesh_2d_inla,
                   cutoff = 25,
                   max.edge = c(200, 100),
                   offset = c(30, 100))
plot(mesh4)
mesh4$mesh$n

## check folds by year 

print(table(sdm_cv_data$fold_id, sdm_cv_data$year))

## use sdmTMB built in cross validation

sdm_cv_simple <- sdmTMB_cv(
  data= sdm_cv_data,
  formula =presence ~ s(depth) + s(temp) + s(sal) + s(bpi_coarse) + s(bpi_fine) + tri + smf,
  mesh = mesh4,
  family =binomial(link = "logit"),
  spatial = "on", 
  time = "year",
  spatiotemporal = "off",
  fold_ids = sdm_cv_data$fold_id,
)

sdm_cv_simple$fold_loglik
sdm_cv_simple$sum_loglik

roc_sdm_simple <- roc(
  sdm_cv_simple$data$presence,
  plogis(sdm_cv_simple$data$cv_predicted),
  quiet = TRUE
)

## overall AUC
auc_sdm_simple <- as.numeric(auc(roc_sdm_simple))
auc_sdm_simple

sqrt(mean((sdm_cv_simple$data$presence - sdm_cv_simple$data$cv_predicted)^2)) 
mean(abs(sdm_cv_simple$data$presence - sdm_cv_simple$data$cv_predicted))

## AUC per fold
fold_aucs_sdm_simple <- sdm_cv_data %>%
  mutate(
    cv_predicted_prob = plogis(sdm_cv_simple$data$cv_predicted)
  ) %>%
  group_by(fold_id) %>%
  summarize(
    n = n(),
    n_presence = sum(presence),
    auc = as.numeric(auc(roc(presence, cv_predicted_prob, quiet = TRUE))))

fold_aucs_sdm_simple

saveRDS(fold_aucs_sdm_ices, "Data/fold_aucs_sdm_ices.rds")
write.csv(fold_aucs_sdm_ices, "Data/fold_aucs_sdm_ices.csv")

sdm_cv_results <- readRDS("Data/fold_aucs_sdm.rds")

### SDM hypothesis testing

sdm_cv_mod_2 <- sdmTMB_cv(
  formula = presence ~ s(depth) + s(temp) + s(sal) + met_season,
  data = sdm_cv_data,
  mesh = mesh4,
  family = binomial(link = "logit"),
  spatial = "on",
  time = "year",
  spatiotemporal = "off",
  fold_ids = sdm_cv_data$fold_id
)

roc_sdm_2 <- roc(
  sdm_cv_mod_2$data$presence,
  plogis(sdm_cv_mod_2$data$cv_predicted),
  quiet = TRUE
)

## overall AUC
auc_sdm_temporal_2 <- as.numeric(auc(roc_sdm_2))
auc_sdm_temporal_2

sdm_cv_mod_2$fold_loglik
sdm_cv_mod_2$sum_loglik

sqrt(mean((sdm_cv_mod_2$data$presence - sdm_cv_mod_2$data$cv_predicted)^2)) 
mean(abs(sdm_cv_mod_2$data$presence - sdm_cv_mod_2$data$cv_predicted))


## AUC per fold
fold_aucs_2 <- sdm_cv_data %>%
  mutate(
    cv_predicted_prob = plogis(sdm_cv_mod_2$data$cv_predicted)
  ) %>%
  group_by(fold_id) %>%
  summarize(
    n = n(),
    n_presence = sum(presence),
    auc = as.numeric(auc(roc(presence, cv_predicted_prob, quiet = TRUE))))

fold_aucs_2

saveRDS(fold_aucs_2, "Data/fold_aucs_2_sdm.rds")
write.csv(fold_aucs_2, "Data/fold_aucs_2_sdm.csv")

