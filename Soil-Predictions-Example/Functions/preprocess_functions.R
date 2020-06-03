###--- Preprocessing functions ---###

#----------------------------------------------#
# Get Spectral Library Functions #
#----------------------------------------------#

# Get Spectral Library
getSpecLib <- function(SPECPATH="/Data_Raw/SPECTRA", LABPATH="Data_Raw/LAB_DATA.csv", SAVENAME="none"){
  # Process Reference Set Spectra 
  spectra <- opus_to_dataset(SPECPATH)
  spectra$spc <- subset_spectral_range(spectra$spc)
  # (Where calibration transfer would occur)
  spectra$spc <- base_offset(spectra$spc)
  
  # Merge with Reference Set Lab Data
  lab <- data.frame(read.csv(LABPATH))
  speclib <- merge(lab, spectra, all.y=TRUE)
  
  # Save speclib after preprocessing
  if(SAVENAME!= "none"){
    assign(SAVENAME,speclib)
    speclib <- get(SAVENAME)
    if(file.exists("./Data_Processed")==FALSE){dir.create("./Data_Processed")}
    savepath <- paste0("Data_Processed/",SAVENAME,".RData")
    save(speclib, file=savepath)
    write.csv(speclib, savepath, row.names=FALSE)
    
  }
  return(speclib)
}


# Opus to Dataset
library(stringr) #used for str_sub
library(foreach)
source("Functions/simplerspc/gather-spc.R")
source("Functions/simplerspc/read-opus-universal.R")

opus_to_dataset <- function(SPECPATH ="/Data_Raw/SPECTRA", NWAVE=3017, SAVE=FALSE){
  #--- Converts OPUS files to RData ---#
  
  # List directory and files containing MIR in opus format
  dirs <- list.dirs(paste(getwd(),SPECPATH,sep=""), full.names=TRUE)
  all.files <- list.files(dirs, pattern= "*.0", recursive=TRUE,full.names=TRUE)
  all.files[1]
  
  # Read opus format files as list and gather spc
  spc_list <- read_opus_univ(fnames = all.files, extract = c("spc"))
  soilspec_tbl <- spc_list %>%
    gather_spc() # Gather list of spectra data into tibble data frame
  
  # Process to output spectra in desired format
  spc <- soilspec_tbl$spc
  spc.trun <- lapply(1:length(spc),function(x) spc[[x]][,1:NWAVE]) # Truncate at 3017, the minimal number of wavenumbers of spectra
  spc.df <- as.data.frame(matrix(unlist(spc.trun), nrow=length(spc.trun), byrow=T)) # Convert to dataframe
  
  # Assigns wavenumbers as column names
  colnames(spc.df) <- colnames(spc.trun[[1]]) 
  rownames(spc.df) <- as.character(seq(1,nrow(spc.df)))
  
  # Assigns sample_id from opus file names
  spc.df <- data.frame(sample_id = soilspec_tbl$sample_id, spc.df)
  spc.df$sample_id <- str_sub(spc.df$sample_id,1,str_length(spc.df$sample_id)-10)
  #spc.avg.df <- aggregate(.~sample_id, data = spc.df, FUN=mean, na.rm=TRUE)
  
  # Reformates dataframe with spectral matrix column
  spectra <- data.frame(spc.df[,1])
  spectra$spc <- as.matrix(spc.df[,2:ncol(spc.df)])
  colnames(spectra) <- c("sample_id", "spc")
  
  if(SAVE==TRUE){
    ##optional save the spectra here
    save(spectra, file="Data_Processed/ref-spectra_original.RData")
    write.csv(spectra, "Data_Processed/ref-spectra_original.csv", row.names=FALSE)
  }
  
  return(spectra)
  
}


# Subset Spectral Range

subset_spectral_range <- function(SPECTRA){
  #--- Subsets columns of the spectral range,
  #--- excluding CO2 sensitive region and truncating at 628.
  col.names <- colnames(SPECTRA)
  col.names <- as.numeric(substring(col.names,2))
  
  cutoff <- which(col.names <= 628)[1]
  SPECTRA <- SPECTRA[,-c(cutoff:length(col.names))] #truncate at >= 628
  
  min.index <- which(col.names <= 2389)[1]
  max.index <- which(col.names <= 2268)[1]
  SPECTRA <- SPECTRA[,-c(min.index:max.index)] #remove CO2 region
  
  return(SPECTRA)
}


# Base Offset

library(matrixStats)
base_offset <- function(x){
  row_mins <- rowMins(x)
  return(x-row_mins)
}


#----------------------------------------------#
# Refine Spectral Library Functions #
#----------------------------------------------#

# Refine Spectral Library
source("Functions/outlier_functions.R")

refineSpecLib <- function(SPECLIB, PROP=NA, OUTLIER=c("stdev"), LARGE=FALSE, SAVENAME=paste0("sub.", PROP)){
  
  # Remove rows with faulty lab data
  if(!is.na(PROP)){
    SPECLIB  <- noNA(SPECLIB , PROP) # Remove NAs
    SPECLIB  <- noNeg(SPECLIB , PROP) # Remove Negative
    if("stdev" %in% OUTLIER){
      SPECLIB  <- SPECLIB[-stdev_outliers(SPECLIB,PROP),] # Remove lab data outliers
    } 
  }
  
  # Remove spectral outliers
  if(!("fratio" %in% OUTLIER)){
    #SPECLIB  <- SPECLIB[-fratio_outliers(SPECLIB),] # Identified with fratio
  } 
  
  # Subset a large dataset to 15000
  if(LARGE==TRUE){
    SPECLIB$spc <- sub_large_set(SPECLIB) # Subset to 15000 samples
  }
  
  # Save the refined reference set for OC
  if(SAVENAME != "none"){
    if(file.exists("./Data_Processed")==FALSE){dir.create("./Data_Processed")}
    assign(SAVENAME, SPECLIB)
    save(list= SAVENAME, file=paste0("Data_Processed/", SAVENAME, ".RData"))
  }
  
  return(SPECLIB)
}


# Get Rid of NAs

noNA <- function(DATASET, column){
  return(DATASET[!is.na(DATASET[,column]),])
}


# Get Rid of Negative Values

noNeg <- function(DATASET, column){
  return(DATASET[which(DATASET[,column] > 0),])
}


# Subset Large Datasets

library(clhs)
sub_large_set <- function(SPECLIB, subcount=15000){
  #Conditional Latin Hypercube Sampling if the set exceeds 15000 samples
  spectra <- data.frame(SPECLIB$spc)
  subset <- clhs(spectra, size = subcount, progress = TRUE, iter = 500)
  SPECLIB <- SPECLIB[subset,] #double check
  
  return(SPECLIB)
}


# Split into Calibration and Validation Sets
library(prospectr)
calValSplit <- function(SPECLIB, PROP=NA, SAVEDIR="Data_Processed", FRAC=0.8){
  
  #perform kennard stone to separate data into 80% calibration and 20% validation sets
  ken_stone<- prospectr::kenStone(X = SPECLIB$spc, 
                                  k = as.integer(FRAC*nrow(SPECLIB)), 
                                  metric = "mahal", pc = 10)
  calib <- SPECLIB[ken_stone$model, ]
  valid <- SPECLIB[ken_stone$test, ]
  
  # {Optional} Save the Sets
  if(SAVEDIR != "none"){
    if(file.exists(SAVEDIR)==FALSE){dir.create(SAVEDIR)}
    if(!is.na(PROP)){
      cal_savename <- paste("calib", property, sep=".") # Ex: calib.OC
      val_savename <- paste("valid", property, sep=".") # Ex: valid.OC
    }else{
      cal_savename <- "calib.ALL"
      val_savename <- "valid.ALL"
    }
    assign(cal_savename, calib)
    save(list=cal_savename, file=paste0(SAVEDIR,"/", cal_savename,".RData"))
    assign(val_savename, valid)
    save(list=val_savename, file=paste0(SAVEDIR,"/", val_savename,".RData"))
  }
  
  return(list(calib=get(cal_savename), valid=get(val_savename)))
}









