##generate climate predictions (standalone script)

library(terra)
library(dplyr)
library(sdmTMB)

setwd("D:/BrennaStuff/Github/Turbot-Thesis-Brenna")

mod_climate <- readRDS("Data/mod_climate_ices_ns_n_anis.rds")
mod_climate$formula
pres_abs <- readRDS("Data/pres_abs_extracted.rds")

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

years <- 2040:2049
scenarios <- c(585)
#seasons <- c("spring", "fall")

## loop turning EACH stack into a data frame, subset complete cases,
## standarize using training data, make season a factor, and predict

## make a smaller set
if(!dir.exists("prediction_frames")){
  dir.create("prediction_frames", showWarnings = FALSE)
}
if(!dir.exists("climate_predictions")){
  dir.create("climate_predictions", showWarnings = FALSE)
}

#how many chunks to use?
n <- 10

for(ssp in scenarios) {
  for(yr in years) {
   # for(ss in seasons) {
      cat("\nProcessing:", ssp, yr, "\n")
      
      #check if this has been run and skip if it has
      if (
        !file.exists(
          sprintf("climate_predictions/preds_%s_%s.rds", ssp, yr)
        )
      ){
        
        stack_file <- sprintf(
          "climate_prediction_stacks/stack_%s_%s.rds",
          ssp, yr)
        
        prediction_stack_ssp_year <- readRDS(stack_file)
        #prediction_stack_ssp_year <- unwrap(prediction_stack_ssp_year)
        
        cat("Changing stack to frame for:", ssp, yr, "\n")
        
        prediction_frame_ssp_year <- as.data.frame(prediction_stack_ssp_year, na.rm = FALSE)
        
        complete_cases <- complete.cases(prediction_frame_ssp_year)
        prediction_frame_ssp_year <- prediction_frame_ssp_year[complete_cases,] 
        
        ## standardize using means and SD of training data 
        prediction_frame_ssp_year_scaled <- prediction_frame_ssp_year %>%
          mutate(temp = (temp-temp_mean) / temp_sd,
                 depth = (depth -depth_mean)/depth_sd,
                 sal = (sal-sal_mean )/ sal_sd,
                 smf = (smf-smf_mean)/smf_sd,
                 bpi_coarse = (bpi_coarse-bpi_coarse_mean)/bpi_coarse_sd,
                 bpi_fine = (bpi_fine-bpi_fine_mean)/bpi_fine_sd,
                 tri = (tri-tri_mean)/tri_sd)
        
        rm(prediction_frame_ssp_year)
        gc()
        
        ## make season a factor just in case it didn't change  
       # prediction_frame_ssp_year_season_scaled$met_season <- factor(prediction_frame_ssp_year_season_scaled$met_season,
                                                                # levels = c("winter", "spring", "summer" ," fall"),
                                                               #  labels = c(1,2, 3, 4))
        
        cat ("Predicting for:", ssp, yr, "\n")
        
        ## now loop every single stack frame through this predict function
        preds_year <- list()
        chunk <- ceiling(nrow(prediction_frame_ssp_year_scaled)/n)
        for(i in 1:n){
          print(i)
          if(i<n){
            preds_year[[i]] <- predict(mod_climate, 
                                              newdata=prediction_frame_ssp_year_scaled[((i-1)*chunk+1):(i*chunk),], 
                                              type = "response")
          } else{
            preds_year[[i]] <- predict(mod_climate, 
                                              newdata=prediction_frame_ssp_year_scaled[((i-1)*chunk+1):nrow(prediction_frame_ssp_year_scaled),], 
                                              type = "response")
          }
        }
        preds_year <- do.call(rbind, preds_year)
        
        rm(prediction_frame_ssp_year_scaled)
        gc()
        
        saveRDS(
          preds_year,
          file = sprintf("climate_predictions/preds_%s_%s.rds", ssp, yr)
        )
        
        rm(preds_year)
        gc()
      }
    }
  }

