###---Performance Functions---###

# Get Model Results (Statistics & Predictions)

getModResults <- function(PROP, MODTYPE, MODNAME=NA, MODPATH=NA, PREDNAME=NA, PREDPATH=NA, SAVEPRED=TRUE, MODPERF=TRUE){
  
  # Load Model
  if(!is.na(MODPATH)){
    MODNAME <- load(MODPATH)
  }
  model <- get(MODNAME)
  
  # Load Prediction Set
  if(!is.na(PREDPATH)){
    PREDNAME <- load(PREDPATH)
  }
  predSet <- get(PREDNAME)
  
  # Extract Predictions
  if(MODTYPE=="MBL"){
    ncomp_onesigma <- NA
    pred_type <- NA
    predictions <- getPredMBL(model)
  }
  if(MODTYPE=="PLS"){
    # Find Optimal Number of Components
    ncomp_onesigma <- selectNcomp(model, method = "onesigma", plot=TRUE, main=PROP)
    # Get Predictions
    pred_type <- "predict"
    predictions <- getPredPLS(model, ncomp_onesigma, pred_type , predSet)
  }
  
  # Get Pred Versus Observations
  lab_data <- predSet[,PROP] # Lab Data
  predobs <- data.frame(predSet[,"sample_id"], predictions, lab_data)
  names(predobs) <- c("sample_id","pred","obs")
  
  # {Optional} Model Performance
  if(MODPERF==TRUE){
    modstats <- getModStats(PREDOBS= predobs, PROP=PROP, NCOMP=ncomp_onesigma, MODNAME=MODNAME, 
                            PREDTYPE= pred_type, PREDNAME=PREDNAME, SAVE=TRUE)
  }
  # {Optional} Save Predictions
  if(SAVEPRED==TRUE){
    savePredictions(predobs, PROP, MODTYPE, predSet, paste0(PREDNAME,"_predictions.csv"))
  }
  names(predobs) <- c("sample_id", paste0(PROP,".",MODTYPE), PROP)
  return(predobs)
}


# Get Model Statistics

getModStats <- function(PREDOBS, PROP=NA, NCOMP=NA, MODNAME=NA, PREDTYPE=NA, PREDNAME=NA, SAVE=FALSE){
  
  cat(paste("\n", PROP, "Summary:"))
  TIME <- as.character(Sys.time()[1])
  
  # Regress predicted versus observed
  PREDOBS <- na.omit(PREDOBS)
  reg_mod <- lm(PREDOBS$pred ~ PREDOBS$obs)
  sum_perf <- summary(reg_mod)
  
  # Get statistics
  R2 <- round(sum_perf$r.squared,4)
  R2_adj <- round(sum_perf$adj.r.squared,4)
  b0 <- round(sum_perf$coefficients[1], 2) # Y-Intercept
  b1 <- round(sum_perf$coefficients[2],2) # Slope
  RMSE <- round(sqrt(mean((PREDOBS$pred - PREDOBS$obs)^2)),2)
  bias <- round((sum(PREDOBS$pred, na.rm=TRUE)- sum(PREDOBS$obs, na.rm=TRUE))/length(PREDOBS$obs),2)
  std <- round(sd(PREDOBS$pred, na.rm=TRUE),2) # Standard Deviation
  rpd <- round(std / RMSE,2) # Residual Prediction Deviation
  
  # Assemble Row
  modStats <- data.frame(Timestamp=TIME, Property=PROP, Mod_Name=MODNAME, Pred_Type=PREDTYPE,
                         Pred_Data=PREDNAME, ncomp=NCOMP, R2=R2, R2_Adj=R2_adj, Y_Int=b0, Slope=b1,
                         RMSE=RMSE, bias=bias,STD=std, RPD=rpd)
  # Print Statistics
  print(t(modStats))
  
  # Write Row 
  if(SAVE==TRUE){saveModStats(modStats)}
  
  # Plot Pred Obs
  plotPred(PREDOBS$obs, PREDOBS$pred, modStats, paste(MODNAME,PREDNAME,"Predictions"), "")
  
  return(modStats)
  
} # End of getModStats


# Save Model Statistics

saveModStats <- function(MODSTATS){
  
  if(file.exists("./Predictions")==FALSE){dir.create("./Predictions")}
  
  modStats_file <- "Predictions/prediction_performance.csv"
  if(file.exists(modStats_file)==FALSE){
    write.csv(MODSTATS, file = modStats_file, row.names=FALSE)
  }else{
    save_table <- read.csv(modStats_file)
    save_table <- rbind(save_table,MODSTATS)
    write.csv(save_table, modStats_file, row.names=FALSE)
  }
  cat(paste("\nStatistics saved to", modStats_file))
}


# Plot Predictions v Observed (with equation & stats)

plotPred <- function(x,y, stats, name=NA, units=NA){
  max <- max(c(x,y))
  lims = c(0,(1.1*max))
  plot(y ~ x, 
       ylab = paste("Predicted", units), 
       xlab=paste("Observed", units), 
       xlim = lims,
       ylim=lims,main = name)
  
  
  reg_model <-lm(y~x)
  abline(reg_model)

  topstats <- bquote(R^2 == .(stats$R2) * "," ~~italic(bias)== .(stats$bias) * "," ~~ RMSE == .(stats$RMSE))
  text(min(x,y),max(x,y), topstats, pos = 4, col="blue")
  
  eqn <- bquote(y== .(stats$Slope) * "x"  * " + " * .(stats$Y_Int))
  text(min(x,y),max(x,y)-(max(x,y)-min(x,y))/10, eqn, pos = 4, col="blue")
}


# Make/Load Table to Save Predictions

getSavePredTable <- function(PREDSET, SAVENAME){
  
  if(file.exists("./Predictions")==FALSE){dir.create("./Predictions")}
  
  predSavePath <- paste0("./Predictions/", SAVENAME)
  if(file.exists(predSavePath) ){
    all_predictions <- read.csv(predSavePath)
  }else{
    all_predictions <- PREDSET[,-ncol(PREDSET)] # remove spectra, last column
  }
  
  return(all_predictions)
}


# Save Predictions

savePredictions <- function(PREDOBS, PROP, MODTYPE, PREDSET, SAVENAME){
  all_predictions <- getSavePredTable(PREDSET, SAVENAME) # Make/Load file to save predictions
  savename <- paste(PROP,MODTYPE,sep=".") # Ex: OC.PLSR column name
  
  if(!(savename %in% names(all_predictions))){
    
    all_predictions <- merge(all_predictions, PREDOBS[,1:2] , all.X=TRUE) # Merge with all_predictions
    ncolm <- ncol(all_predictions)
    names(all_predictions)[ncolm] <- savename
    savefile <- paste0("Predictions/", SAVENAME) # Set file savename
    write.csv(all_predictions, file=savefile, row.names=FALSE) # Save
    
    cat(paste("\nPredictions saved to", savefile)) # Print save location
    
  }else{
    cat("\nPrediction column already exists")
  }
  
  View(all_predictions)
}
  
  
