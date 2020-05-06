library(dplyr)
library(matrixStats)
library(bimixt)
library(car)
library(pls)
library(ranger)
library(resemble)
library(Cubist)
library(Rcpp)
library(doParallel)
library(hexbin)
library(RColorBrewer)
library(openxlsx)

setwd("/Volumes/Macintosh HD/Users/crivard/Desktop/PREDICT_vSelf")
source("functions_modelChoice.R")
source("functions_udev.R")

if(TRUE){
  
  #Initiate variables
  datname <- "Sample-Prep" #replaces modelname
  dat <- read.csv("predDat.csv")
  predDatName <- "Sample-Prep"
  property <- "CaCO3"
  mod_save_wkbk <- "*models.xlsx"
  pred_save_wkbk <- "*predictions.xlsx"
  
  ##-----------Model Data Processing---------------###
  
  if(TRUE){
    #check conditions...
    names(dat)[1:300]
    specCol1 <- 7 #first spectral column
    char_preNum <- 1 #cuts off spec.X
    
    #make data matrix
    n_col <- ncol(dat)
    datm <- dat[,1:specCol1-1] #matrix form of datatable
    datm$spc <- as.matrix(dat[, specCol1:n_col])
    
    #Get rid of CO2 sensitive bands
    col.names <- colnames(datm$spc)
    col.names <- as.numeric(substring(col.names,char_preNum+1))
    
    #adjust range of spectra (truncate)
    cutoff <- which(col.names <= 628)[1]
    datm$spc <- datm$spc[,-c(cutoff:length(col.names))]
    
    min.index <- which(col.names <= 2389)[1]
    max.index <- which(col.names <= 2268)[1]
    datm$spc <- datm$spc[,-c(min.index:max.index)] 
    
    #baseline transformation and divide data by soil properties
    datm$spc <- base_offset(datm$spc)
    datm <- datm[!is.na(datm[,property]),] #no NAs
    datm <- datm[which(datm[,property] > 0),] #no negative values
    
  }
  ## repeat this step for other soil properties -- Al, Ca, CO3, pH, Fe, Ca, BD, OCD, Clay.
  
  #remove outliers by selecting using a standard deviaton threshold
  #Functions in "functions_modelChoice.R"
  if(TRUE){
    #fit1 <- plsr(sqrt(get(property))~spc, ncomp=20, data = datm, valid="CV", segments = 50) #y, x, number components, data, cross validation, 
    fit1 <- plsr(sqrt(get(property))~spc, ncomp= 20, data = datm, valid="CV", segments = 50) #y, x, number components, data, cross validation, 
    pred <- c(predict(fit1, newdata = datm$spc,ncomp=20))^2
    sd.outlier <- optimum_sd_outlier(pred, datm[,property], seq(0.1,3, by =0.02))
    row.index <- outlier(pred, datm[,property], sd.outlier[1])
    if(length(row.index) > 0){
      datm <- datm[row.index,]
    }
  }
  
  #perform kennard stone to separate data into 80% calibration and 20% validation sets
  ks<- prospectr::kenStone(X = datm$spc, k = as.integer(0.8*nrow(datm)), metric = "mahal", pc = 10) ## repeat this step for other soil properties -- Al, Ca, CO3, pH, Fe, Ca, BD, OCD, Clay.
  calib <- datm[ks$model, ]
  valid <- datm[ks$test, ]
  
  #remove outliers of spectral data as well (those in valid that are not represented in calib)
  if(TRUE){
    fit2 <- plsr(sqrt(get(property))~spc, ncomp=20, data = calib, valid="CV", segments = 50)
    com.spc <- rbind(calib$spc, valid$spc)
    calval.scores <- fit2$scores #scores from the calibration set
    all.loads <- fit2$loadings #loadings for all data
    
    pred.scores <- predict(fit2, newdata= valid$spc, type ="scores")
    pred.loads <- fit2$loadings
    
    calval.scores <- rbind(calval.scores, pred.scores)
    
    test <- calval.scores %*% t(all.loads)
    
    all <- scale(com.spc, scale = FALSE) #Takes the log base e (ln)
    plot(test[1,], type = "l", ylim = c(-0.5,0.5))
    lines(all[1,], col="red")
    
    res <- (test - all)^2
    res <- sqrt(rowSums(res))
    
    for(i in 1:length(res)){
      sample.fratio <- (length(res)-1) * res^2/sum((res[-i])^2)
    }
    
    ok <- pf(sample.fratio, 1, length(res)) 
    rows <- which(ok>0.99)
    if(length(rows) > 0){
      print(rows)
      #datm <- datm[-rows,] #unsure if they should also be deleted here
      #valid <- valid[-rows,] #should be deleted from validation, but 
    }
    
  }
  
}#End of data pre-processing


##-----------PLSR Model---------------###
if (TRUE){
  ## 1. Creating PLSR model
  modelType <- "plsr"
  units <- "%"
  unit_adj <- 1
  
  fd="fd"
  
  if(fd=="fd"){
    datm$spc <- t(diff(t(datm$spc), differences =1)) #first derivative
    savename <- paste(property, modelType, fd, datname, sep=".")
  }else{
    savename <- paste(property, modelType, datname, sep=".")
  }
  
  validType ="LOO"
  plsr.model <- plsr(sqrt(get(property))~spc, ncomp=20, data = datm, valid=validType)
  save(plsr.model, file = paste(savename,".RData", sep=""))#Saving the model
  
  ###2. Applying PLSR model
  predDatName <- "datm"
  predDat <- get(predDatName)
  #savename is within fd/nofd conditionals
  ncomp.onesigma <- selectNcomp(plsr.model, method = "onesigma", plot = TRUE, ylim = c(0, 50))
  ncomp.onesigma
  
  if(fd=="fd" && predDatName!="datm"){
    predVals <- c(predict(plsr.model, newdata = t(diff(t(predDat$spc), differences =1)), ncomp=ncomp.onesigma))^2
    savename <- paste(property, modelType, fd, predDatName, paste("v",datname,sep=""), sep=".")
    print("took 1st deriv")
  }else{
    predVals <- c(predict(plsr.model, newdata = predDat$spc, ncomp=ncomp.onesigma))^2
    savename <- paste(property, modelType, predDatName, paste("v",datname,sep=""), sep=".")
    print("did not take 1st deriv")
  }
  
  #pred_obs table
  col.names <- colnames(predDat)
  propCol <- which(col.names == toString(property))[1]
  pred_obs <- data.frame(predDat[,1], predVals, (predDat[,propCol])*unit_adj)
  names(pred_obs) <- c("ID", "pred", "obs")
  
  #check results
  max <- max(pred_obs[,c("pred", "obs")])
  plot.plsr(pred_obs$obs, pred_obs$pred, property, c(0,(1.1*max)),units)
  #----
  
  ###3. Saving Results
  if(TRUE){
    #save model coefficients
    mod_tags <- c("OR")
    savename <- paste(property,modelType,paste(mod_tags,collapse="."),sep=".")
    mwb <- loadWorkbook(mod_save_wkbk)
    addWorksheet(mwb,savename)
    beta_coeff <- data.frame(plsr.model$coefficients)
    beta_coeff <- beta_coeff[ncomp.onesigma]
    names(beta_coeff) <- property
    writeData(mwb,savename,beta_coeff,startCol=1,startRow=1)
    saveWorkbook(mwb,mod_save_wkbk, overwrite=TRUE)
    
    #save model
    save(plsr.model,file=paste(savename,".RData",sep=""))
    
    #save predictions
    pwb <- loadWorkbook(pred_save_wkbk)
    addWorksheet(pwb,savename)
    writeData(pwb,savename,pred_obs,startCol=1,startRow=1)
    
    #save plot
    max <- max(pred_obs[,c("pred", "obs")])
    plot.plsr(pred_obs$obs, pred_obs$pred,property, c(0,(1.1*max)),units)
    insertPlot(pwb,savename,width=5,height=5,xy=NULL,startRow=2, startCol=5, fileType="png",units="in",dpi=300 )
    saveWorkbook(pwb,pred_save_wkbk, overwrite=TRUE)
    
    #save model performance & details
    sum_perf <- summary(lm(pred_obs$obs ~ pred_obs$pred))
    R2 <- sum_perf$r.squared
    R2_adj <- sum_perf$adj.r.squared
    slope <- sum_perf$coefficients[2]
    RMSE <- round(sqrt(mean((pred_obs$pred - pred_obs$obs)^2)),2)
    bias <- round((sum(pred_obs$pred, na.rm=TRUE)- sum(pred_obs$obs, na.rm=TRUE))/length(pred_obs$obs),2)
    notes <- "" #add special details here
    row <- paste(paste("\r",Sys.time(),sep=""), property, modelType, validType, datname, predDatName,notes,R2,R2_adj,slope,RMSE,bias,ncomp.onesigma,sep=",")
    rownames <- c("Timestamp","Property","Model Type","LOO/CV","Model Data","Pred Data","Notes","R2","R2 Adj","Slope","RMSE","bias","ncomp")
    write.table(row, file = "*models.csv", sep = ",", append = TRUE, quote = FALSE, row.names = FALSE, col.names=FALSE, eol="\r")
  }
  
}


##-----------MBL Model---------------###
if (TRUE){
  modelType <- "mbl"
  predDatName <- "valid"
  predDat <- get(predDatName)
  units <- ""
  unit_adj <- 1
  
  fd = ""
  if(fd=="fd"){
    calib$spc <- t(diff(t(calib$spc), differences =1)) #first derivative
    predDat$spc <- t(diff(t(predDat$spc), differences =1)) #first derivative
    print("took 1st derivative")
  }
  #savepath <- paste("./Models/USGS/Geochem/Models/", sep="")
  #savename <- paste(property, modelType, fd, predDatName, paste("v",datname,sep=""), sep=".")
  
  
  #1- Creating the model {option 1}
  if(FALSE){
    Xu <- predDat$spc
    Yu <- sqrt(predDat[,property]) 
    Yr <- sqrt(calib[,property])
    Xr <- calib$spc
    Xu <- Xu[!is.na(Yu),]
    Yu <- Yu[!is.na(Yu)]
    Xr <- Xr[!is.na(Yr),]
    Yr <- Yr[!is.na(Yr)]
    ctrl <- mblControl(sm = 'pc', pcSelection = list('opc', 50),
                       valMethod = 'loc_crossval',center=TRUE,scale=FALSE,allowParallel=FALSE)
    mbl.sqrt <- mbl(Yr = Yr, Xr = Xr, Yu = Yu, Xu = Xu,
                    mblCtrl = ctrl,
                    dissUsage = 'none',
                    #k = seq(40, 100, by = 20),
                    k = seq(10, 20, by = 2),
                    method = 'pls', pls.c = 6)
    #save(mbl.sqrt, file= paste(savepath, savename,".RData", sep=""))
    
    #2- Applying the model
    #savepath <- paste("./Models/USGS/Geochem/Predictions/", sep="")
    #savename <- paste(property, modelType, fd, predDatName, paste("v",datname,sep=""), sep=".")
    predVals <- c(mbl.sqrt$results$Nearest_neighbours_40$pred)^2
    #predVals <- c(mbl.sqrt$results$Nearest_neighbours_14$pred)^2
  }
  
  #{option 2-- weighted to closer neighbors}
  if(TRUE){
    
    Xu <- predDat$spc
    Yu <- sqrt(predDat[,property]) 
    Yr <- sqrt(calib[,property])
    Xr <- calib$spc
    Xu <- Xu[!is.na(Yu),]
    Yu <- Yu[!is.na(Yu)]
    Xr <- Xr[!is.na(Yr),]
    Yr <- Yr[!is.na(Yr)]
    
    dmetric = "pls"
    diss2test <- seq(0.3, 1, by=0.1)
    kminmax <- c(10, nrow(calib$spc))
    rmethod <- "wapls1"
    pls.f <- c(minpls=3, maxpls=20)
    
    ctrl <- mblControl(sm = dmetric, pcSelection = list("opc", 50), valMethod = "NNv", 
                       returnDiss = TRUE, scaled = FALSE, center = TRUE)
    
    mbl.sqrt <- mbl(Yr = Yr, Xr = Xr, Xu = Xu, mblCtrl = ctrl, dissUsage = "none", k.diss = diss2test, k.range = kminmax, 
                    pls.c = pls.f, method = rmethod)
    
    idx.best.ca <- which.min(mbl.sqrt$nnValStats$st.rmse)
    best.kdiss.ca <- mbl.sqrt$nnValStats$k.diss[idx.best.ca]
    
    ## Get the predicted values for the validation set
    predVals <- c(getPredictions(mbl.sqrt)[, idx.best.ca])^2
  }
  
  
  #Pred/Obs Table
  col.names <- colnames(predDat)
  propCol <- which(col.names == toString(property))[1]
  
  pred_obs <- data.frame(predDat[,1], predVals, (predDat[,propCol])*unit_adj)
  names(pred_obs) <- c("ID", "pred", "obs")
  
  #Plot
  max <- max(pred_obs[,c("pred", "obs")])
  plot.plsr(pred_obs$obs, pred_obs$pred,property, c(0,(1.1*max)),units)
  #--check results
  
  if(TRUE){
    
    #save model
    mod_tags = c("2","OR")
    savename = paste(property,modelType, paste(mod_tags,collapse="."), sep=".")
    save(mbl.sqrt, file= paste(savepath, savename,".RData", sep=""))
    
    
    #save predictions
    pwb <- loadWorkbook(pred_save_wkbk)
    addWorksheet(pwb,savename)
    writeData(pwb,savename,pred_obs,startCol=1,startRow=1)
    
    #save plot
    max <- max(pred_obs[,c("pred", "obs")])
    plot.plsr(pred_obs$obs, pred_obs$pred,property, c(0,(1.1*max)),units)
    insertPlot(pwb,savename,width=5,height=5,xy=NULL,startRow=2, startCol=5, fileType="png",units="in",dpi=300 )
    saveWorkbook(pwb,pred_save_wkbk, overwrite=TRUE)
    
    #save model performance & details
    sum_perf <- summary(lm(pred_obs$obs ~ pred_obs$pred))
    R2 <- sum_perf$r.squared
    R2_adj <- sum_perf$adj.r.squared
    slope <- sum_perf$coefficients[2]
    RMSE <- round(sqrt(mean((pred_obs$pred - pred_obs$obs)^2)),2)
    bias <- round((sum(pred_obs$pred, na.rm=TRUE)- sum(pred_obs$obs, na.rm=TRUE))/length(pred_obs$obs),2)
    notes <- "mbl option 2- both outlier removal strategies" #add special details here
    validType = "" # not applicable for mbl
    ncomp.onesigma = ""# not applicable for mbl
    row <- paste(paste("\r",Sys.time(),sep=""), property, modelType, validType, datname, predDatName,notes,R2,R2_adj,slope,RMSE,bias,ncomp.onesigma,sep=",")
    rownames <- c("Timestamp","Property","Model Type","LOO/CV","Model Data","Pred Data","Notes","R2","R2 Adj","Slope","RMSE","bias","ncomp")
    write.table(row, file = "./Models/USGS/Geochem/Models/*Models_USGS_GC.csv", sep = ",", append = TRUE, quote = FALSE, row.names = FALSE, col.names=FALSE, eol="\r")
  }
  
}

#dirname(rstudioapi::getSourceEditorContext()$path) #gets you to your current directory which is cool
