###-----MBL Functions-----###

# Run MBL Model
library(resemble)
runMBL <- function(PROP, REFNAME=NA, REFPATH=NA, PREDNAME=NA, PREDPATH=NA, SAVENAME=paste0("mbl.",PROP)){
  
  # Load Reference Set
  if(!is.na(REFPATH)){
    REFNAME <- load(REFPATH)
  }
  refSet <- get(REFNAME)
  
  # Load Prediction Set
  if(!is.na(PREDPATH)){
    PREDNAME <- load(PREDPATH)
  }
  predSet <- get(PREDNAME)
  
  # Define Input Datasets
  Xu <- predSet$spc              # Prediction Spectra 
  Yu <- sqrt(predSet[,PROP]) # Prediction Lab Data
  Yr <- sqrt(refSet[,PROP])  # Reference Spectra
  Xr <- refSet$spc               # Reference Lab Data
  
  # Get Rid of NAs
  Xu <- Xu[!is.na(Yu),]
  Xr <- Xr[!is.na(Yr),]
  Yu <- Yu[!is.na(Yu)]
  Yr <- Yr[!is.na(Yr)]
  
  # Set Control Parameters
  ctrl <- mblControl(sm = 'pc', pcSelection = list('opc', 50),
                     valMethod = 'NNv',center=TRUE,scale=FALSE,allowParallel=FALSE)
  
  # Run MBL Model
  mbl.sqrt <- mbl(Yr = Yr, Xr = Xr, Yu = Yu, Xu = Xu, mblCtrl = ctrl, dissUsage = 'none',
                  k = seq(40, 100, by = 20), method = 'pls', pls.c = 6)
  
  #mbl.sqrt <- mbl(Yr = Yr, Xr = Xr, Xu = Xu, mblCtrl = ctrl, dissUsage = "none", k.diss = seq(0.3, 1, by=0.1),
  #k.range = c(20, nrow(refSet)), pls.c = c(minpls=3, maxpls=20), method = "wapls1")
  
  
  # Save MBL Model
  if(SAVENAME != "none"){
    modelName <- paste("mbl", PROP, sep=".")
    assign(modelName, mbl.sqrt)
    savefile <- paste("./Models/mbl", PROP,"RData", sep=".")
    save(list= modelName, file = savefile)
    cat(paste("\nModel",modelName, "saved to", savefile))
  }
  
  return(mbl.sqrt)
  
}


# Finding Best Model Within

bestModMBL <- function(mbl.sqrt){
  
  valType <- mbl.sqrt$cntrlParam$valMethod
  
  if(valType=="NNv"){
    index_best_model <- which.min(mbl.sqrt$nnValStats$st.rmse)
  }
  if(valType=="loc_crossval"){
    index_best_model <- which.min(mbl.sqrt$localCrossValStats$st.rmse)
  }
  
  best_model_name <- names(mbl.sqrt$results)[index_best_model]
  
  return(best_model_name)
}


# Getting Predictions

getPredMBL <- function(mbl.sqrt, model_name=NULL){
  
  if(is.null(model_name)){
    model_name <- bestModMBL(mbl.sqrt)
  }
  sqrt_preds <- eval(parse( text=paste0("mbl.sqrt$results$", model_name,"$pred" )))
  
  predictions <- c(sqrt_preds)^2
  
  return(predictions)
}


# Getting Lab Data

getLabMBL <- function(mbl.sqrt){
  
  Yu <- mbl.sqrt$call$Yu
  
  if(!is.null(Yu)){
    sqrt_lab <- eval(parse( text=paste0("mbl.sqrt$results$",best_model_name,"$yu.obs" )))
  }else{
    sqrt_lab <- NULL
  }
  
  lab <- c(sqrt_lab)^2
  
  return(lab)
  
}