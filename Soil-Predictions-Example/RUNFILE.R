#Title: RUNFILE- Soil-Predictions-Example
#Authors: Charlotte Rivard & Shree Dangal
#Date: 5/8/20
#Summary: #The following script predicts values for several soil properties
          #from spectral data, using the following machine learning models:
          #1- Partial Least Squares Regression
          #2- Memory Based Learner

#----------------------------------------------#
      # Reference Set Preprocessing #
#----------------------------------------------#
source("Functions/preprocess_functions.R")

# Process Reference Set Spectra 
spectra <- opus_to_dataset("/Data_Raw/SPECTRA")
spectra$spc <- subset_spectral_range(spectra$spc)
spectra$spc <- base_offset(spectra$spc)

# Merge with Reference Set Lab Data
library(readr)
lab <- data.frame(read_csv("Data_Raw/LAB_DATA.csv"))
all_data <- merge(lab, spectra, all.y=TRUE)

# Save all_data file after preprocessing
if(file.exists("./Predictions")==FALSE){dir.create("./Data_Processed")}
save(all_data, file="Data_Processed/ref.ALL.RData")
write.csv(all_data, "Data_Processed/ref.ALL.csv", row.names=FALSE)

# Remove rows with poor lab data
properties <- c("OC","SAND","SILT", "CLAY") #Column names of lab data

for(property in properties){
  prop.data <- all_data
  prop.data <- noNA(prop.data, property) # Remove NAs
  prop.data <- noNeg(prop.data, property) # Remove Negative
  prop.data <- noOut(prop.data, property) # Remove Outliers*
  savename <- paste("Data_Processed/ref",property,"RData", sep=".")
  save(prop.data, file=savename) # Save, Ex: all.OC.RData
  
  #split <- calValSplit(prop.data, property) # Split Calibration & Validation Sets
  #calib <- split[1]; valid <- split[2]
  #save(calib, file=paste("Data_Processed/calib", property,"RData", sep="."))
  #save(valid, file=paste("Data_Processed/valid", property,"RData", sep="."))
}

#----------------------------------------------#
    # Partial Least Squares Regression #
#----------------------------------------------#
library(pls)
source("Functions/plsr_functions.R")

load("./Data_Processed/ref.ALL.RData") #named all_data
all_predictions <- all_data[,-ncol(all_data)] #remove spectra last column

# 1- Create the Model

properties <- c("OC","SAND","SILT", "CLAY")

for(property in properties){
  # Load Data
  refSet = paste("./Data_Processed/ref",property,"RData", sep=".") #Ex: ref.OC.RData
  load(refSet) # variable 'prop.data'
  
  # Create Model
  validType <- "CV" # "CV", "LOO", or "none"
  plsr.model <- plsr(sqrt(get(property))~spc, ncomp=20, data = prop.data , valid=validType) 
  
  # Save Model
  if(file.exists("./Models")==FALSE){dir.create("./Models")}
  save(plsr.model, file = paste("./Models/plsr", property,"RData", sep=".")) #Ex: plsr.OC.RData
  
  # Save Model Coefficients
  saveModCoefs(plsr.model, property)
  
  # Show Model Summary 
  showSummary(plsr.model, property, ncomp_onesigma)
  
  #________________________________________#
  
  # Find Optimal Number of Components
  ncomp_onesigma <- selectNcomp(plsr.model, method = "onesigma")
  
  # Get Predictions
  predType <- "predict" # "fitted", "valid", "predict"
  predSet <- "all_data"
  predictions <- getPredictions(plsr.model, ncomp_onesigma, predType, get(predSet))
  sample_id <- getSampleID(prop.data, get(predSet))
  predTable <- data.frame(sample_id, predictions)
  names(predTable) <- c("sample_id", paste(property, predType, sep="."))
  
  # Save Predictions
  all_predictions <- merge(all_predictions, predTable, all.X=TRUE)
  
  # Get Lab Data
  lab_data <- getLabData(plsr.model, predType, get(predSet), property)
  
  # Save Model Performance
  saveModStats(predictions, lab_data, property, ncomp_onesigma, "PLSR", predType, validType, predSet)
  
}
if(file.exists("./Predictions")==FALSE){dir.create("./Predictions")}
save(all_predictions, "all_predictions.RData")
write.csv(all_predictions, "./Predictions/all_predictions.csv")


 #outlier flagging
output.pred <- data.frame(output.pred)
output.pred$outlier <- 0
oc <- read.csv("/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/output/indigo_output/fratio/bd.anal4.csv")
outs <- c(oc$x) 
outs <- outs[outs>10363]
outs <- outs-10363
output.pred$outlier[outs] <- 1
output.pred <- cbind(TERR_ID,output.pred)

write.csv(output.pred, file ="/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/output/indigo_output/all-predictions/pred.bd.anal4.csv")

#----------------------------------------------#
        # Memory Based Learner Model #
#----------------------------------------------#


#----------------------------------------------#
        # Prediction Set Preprocessing #
#----------------------------------------------#

# Spectral Processing: Prediction Set
#spectra <- opus_to_dataset("/Data_Raw/PRED-SPECTRA")
#{where pds would happen}
#spectra <- subset_spectral_range(spectra)
#spectra <- base_offset(spectra)


