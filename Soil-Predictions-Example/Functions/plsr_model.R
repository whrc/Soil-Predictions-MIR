### PLSR Model 

#---Packages---#
library(pls)

#---Set Properties---#
properties <- c("OC", "SAND", "SILT", "CLAY")

#---{Get Properties}--#
if (FALSE){
  files <- list.files("./Data_Processed", pattern = "all.*.RData", full.names =TRUE)
  filenames_split <- strsplit(files, "/")
  splits <- length(filenames_split[[1]])
  properties <- sapply(filenames_split, function(x) x[[splits]])
  properties <- substr(properties, 5, nchar(properties)-6)
}

#For EACH property:

for(property in properties){
  #---Load Data---#
  data <- load(paste("./Data_Processed/all",property,"RData", sep=".")) #Ex: all.OC.RData
  
  #---Create Model---#
  plsr.model <- plsr(sqrt(get(property))~spc, ncomp=20, data = all_data , valid="LOO")
  
  #---Save Model---#
  save(plsr.model, file = paste("./Models/plsr", property,"RData", sep=".")) #Ex: plsr.OC.RData
}

#---Evaluating the Model---#
summary(plsr.model)
ncomp_onesigma <- selectNcomp(plsr.model, method = "onesigma", plot = TRUE)
plot(plsr.model, ncomp = ncomp_onesigma, asp = 1, line = TRUE)

#---Saving Results---#

#Extracting lab_data
lab_data <- c(plsr.model$model$`sqrt(get(property))`^2)

#Selecting predictions with optimal number of components
sqrt_pred <- data.frame(plsr.model$fitted.values)[ncomp_onesigma] 
pred <- c(sqrt_pred^2)

results <- data.frame(lab_data, pred, colnames=c("lab_data","predictions"))


#---Applying Model---#
predVals <- c(predict(plsr.model, newdata = all_data$spc, ncomp=ncomp_onesigma))^2

test <- data.frame(valid_pred, predVals)
names(test) <- c("valid","predict")
plot(valid ~ predict, data=test)


