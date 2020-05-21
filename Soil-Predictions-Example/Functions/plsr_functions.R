###--- PLSR Functions ---###

showSummary <- function(model, property, ncomp_onesigma){
  print(paste(property,"PLSR Model:"))
  summary(model)
  plot(model, ncomp = ncomp_onesigma, asp = 1, line = TRUE, main= paste(property,"validation"))
}


getSampleID <- function(refSet, predSet){
  if(predType == "predict"){
    sample_id <- data.frame(predSet$sample_id)
  }
  else{
    sample_id <- data.frame(refSet$sample_id)
  }
  names(sample_id) <- c("sample_id")
  return(sample_id)
}


getPredictions <- function(model, ncomp_onesigma, predType="fitted", predSet=NULL){
  if(predType=="valid"){
    sqrt_pred <- data.frame(model$validation$pred)[ncomp_onesigma]
  }
  if(predType == "fitted"){
    sqrt_pred <- data.frame(model$fitted.values)[ncomp_onesigma]
  }
  
  if(predType == "predict"){
    sqrt_pred <- c(predict(model, newdata = predSet$spc, ncomp=ncomp_onesigma))
  }
  predictions <- c(sqrt_pred^2)
  
  return(predictions)
}


getLabData <- function(model, predType="fitted", predSet=NULL, property=NULL){

  # Extract Lab Data
  if(predType=="predict"){
    lab_data <- c(predSet[,property])
  }else{
    lab_data <- data.frame(model$model$`sqrt(get(property))`^2)
  }
  
  return(lab_data)
}


saveModStats <- function(predictions, lab_data, property, ncomp_onesigma, modType=NA, predType=NA, validType=NA, predSet=NA, notes=FALSE){
  pred_obs <- data.frame(predictions, lab_data)
  names(pred_obs) <- c("pred","obs")
  
  # Regress predicted versus observed
  sum_perf <- summary(lm(pred_obs$obs ~ pred_obs$pred))
  R2 <- sum_perf$r.squared
  R2_adj <- sum_perf$adj.r.squared
  slope <- sum_perf$coefficients[2]
  RMSE <- round(sqrt(mean((pred_obs$pred - pred_obs$obs)^2)),2)
  bias <- round((sum(pred_obs$pred, na.rm=TRUE)- sum(pred_obs$obs, na.rm=TRUE))/length(pred_obs$obs),2)
  
  # Get {optional} Notes
  if(notes==TRUE){
    notes <- readline(prompt="Notes: ")
  }else{notes <- ""}

  # Assemble Row
  row <- paste(Sys.time(), property, modType, predType, validType, predSet, R2,R2_adj,slope,RMSE,bias,ncomp_onesigma, notes,"\r", sep=",")
  rownames <- c("Timestamp,Property,Model Type,Pred Type,Valid Type,Pred Data,R2,R2 Adj,Slope,RMSE,bias,ncomp, Notes\r")
  
  # Write Row 
  if(file.exists("./Models/model_performance.csv")==FALSE){
    write.table(rownames, file = "./Models/model_performance.csv", sep = ",", append = TRUE, quote = FALSE, row.names = FALSE, col.names=FALSE, eol="\r")
  }
  write.table(row, file = "./Models/model_performance.csv", sep = ",", append = TRUE, quote = FALSE, row.names = FALSE, col.names=FALSE, eol="\r")
  
}


saveModCoefs <- function(model, property){
  beta_coeff <- data.frame(model$coefficients)
  ncomp_onesigma <- selectNcomp(model, method = "onesigma")
  beta_coeff <- beta_coeff[ncomp.onesigma]
  write.csv(beta_coeff, paste("/Models/", property, "_plsr_coefs.csv", sep=""))
}