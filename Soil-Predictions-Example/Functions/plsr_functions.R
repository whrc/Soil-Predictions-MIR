###--- PLSR Functions ---###

# Make PLS Model
library(pls)

makePLSModel <- function(PROP, REFSET=NA, REFPATH=NA, VAL_TYPE="CV", SAVENAME=paste0("plsr.",PROP)){
  print(SAVENAME)
  # Load Data
  if(!is.na(REFPATH)){ 
    refSet.name <- load(REFPATH) # If REFPATH is given
    REFSET <- get(refSet.name) # load as variable refSet
  }
  
  # Create Model
  model <- plsr(sqrt(get(PROP))~spc, ncomp=20, data = REFSET , valid=VAL_TYPE) 
  
  # Save Model
  if(SAVENAME != "none"){
    if(file.exists("./Models")==FALSE){dir.create("./Models")} # Create Folder to Save Models
    assign(SAVENAME, model)
    save(list= SAVENAME, file = paste0("./Models/",SAVENAME,".RData")) #Ex: plsr.OC.RData
  }
  
  return(model)
}


# Get PLS Predictions

getPredPLS <- function(MODEL, NCOMP, PREDTYPE="fitted", PREDSET=NULL){
  if(PREDTYPE=="valid"){
    sqrt_pred <- unlist(data.frame(MODEL$validation$pred)[NCOMP])
  }
  if(PREDTYPE == "fitted"){
    sqrt_pred <- unlist(data.frame(MODEL$fitted.values)[NCOMP])
  }
  
  if(PREDTYPE == "predict"){
    sqrt_pred <- c(predict(MODEL, newdata = PREDSET$spc, ncomp=NCOMP))
  }
  predictions <- c(sqrt_pred^2)
  names(predictions) <- c(seq(1:length(predictions)))
  
  return(predictions)
}


#----------------Optional Functions------------------#

# Get Sample ID

getSampleID <- function(PROPERTY, PREDTYPE, PREDSET=NULL, REFSET=NULL){
  if(PREDTYPE == "predict"){
    sample_id <- data.frame(PREDSET$sample_id)
  }
  else{ # If predtype is "valid" or "fitted"
    sample_id <- data.frame(REFSET$sample_id)
  }
  names(sample_id) <- c("sample_id")
  return(sample_id)
}


# Get Lab Data

getLabPLS<- function(MODEL, PREDTYPE="fitted", PREDSET=NULL, PROP=NULL){

  # Extract Lab Data
  if(PREDTYPE=="predict"){
    lab_data <- c(PREDSET[,PROP])
  }else{
    Y_name <- names(MODEL$model)[1]
    print(Y_name)
    sqrt_lab <- eval(parse( text=paste0("MODEL$model$`", Y_name,"`")))
    #lab_data <- data.frame(MODEL$model$`sqrt(get(property))`^2)
    lab_data <- c(sqrt_lab)^2
  }
  
  return(lab_data)
}


# Get PLS Model Coefficients

getPLSCoefs <- function(MODEL, PROP, NCOMP){
  beta_coeff <- data.frame(MODEL$coefficients)[NCOMP]
  return(beta_coeff)
}

