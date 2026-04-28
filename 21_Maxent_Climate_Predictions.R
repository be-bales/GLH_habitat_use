## maxent climate predictions
library(maxnet)
years <- 2020:2049
scenarios <- c(119, 245, 585)

mod_maxent <- readRDS("Data/maxent_model_april.rds")

## loop turning EACH stack into a data frame, subset complete cases,
## standarize using training data, make season a factor, and predict

## make a smaller set
if(!dir.exists("prediction_frames")){
  dir.create("prediction_frames", showWarnings = FALSE)
}
if(!dir.exists("maxent_climate_predictions")){
  dir.create("maxent_climate_predictions", showWarnings = FALSE)
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
          sprintf("maxent_climate_predictions/preds_%s_%s.rds", ssp, yr)
        )
      ){
        
        stack_file <- sprintf(
          "climate_prediction_stacks/stack_%s_%s.rds",
          ssp, yr )
        prediction_stack_ssp_year <- readRDS(stack_file)
      #  prediction_stack_ssp_year <- unwrap(prediction_stack_ssp_year)
        
        cat("Changing stack to frame for:", ssp, yr,"\n")
        
        prediction_frame_ssp_year <- as.data.frame(prediction_stack_ssp_year, na.rm = FALSE)
        
        complete_cases <- complete.cases(prediction_frame_ssp_year)
        prediction_frame_ssp_year <- prediction_frame_ssp_year[complete_cases,] 
        
        prediction_frame_ssp_year_scaled <- prediction_frame_ssp_year
        
        rm(prediction_frame_ssp_year)
        gc()
      
        
        cat ("Predicting for:", ssp, yr, "\n")
        
        ## now loop every single stack frame through this predict function
        preds_year <- list()
        chunk <- ceiling(nrow(prediction_frame_ssp_year_scaled)/n)
        for(i in 1:n){
          print(i)
          if(i<n){
            preds_year[[i]] <- predict(mod_maxent, 
                                              newdata=prediction_frame_ssp_year_scaled[((i-1)*chunk+1):(i*chunk),], 
                                              type = "cloglog")
          } else{
            preds_year[[i]] <- predict(mod_maxent, 
                                              newdata=prediction_frame_ssp_year_scaled[((i-1)*chunk+1):nrow(prediction_frame_ssp_year_scaled),], 
                                              type = "cloglog")
          }
        }
        preds_year <- do.call(rbind, preds_year)
        
        rm(prediction_frame_ssp_year_scaled)
        gc()
        
        saveRDS(
          preds_year,
          file = sprintf("maxent_climate_predictions/preds_%s_%s.rds", ssp, yr)
        )
        
        rm(preds_year)
        gc()
      }
    }
  }

