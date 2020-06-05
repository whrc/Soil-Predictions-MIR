###--- Preprocessing functions ---###

#----------------------------------------------#
# Get Spectral Library Functions #
#----------------------------------------------#

# Get Spectral Library

getSpecLib <- function(SPECPATH="/Data_Raw/SPECTRA", 
                       LABPATH="Data_Raw/LAB_DATA.csv", SAVENAME="none"){
  # Extract OPUS Files
  spectra <- opus_to_dataset(SPECPATH)
  
  # Subset Spectral Range
  spectra$spc <- subset_spectral_range(spectra$spc)
  
  # Where calibration transfer would occur
  
  # Baseline Transformation
  spectra$spc <- base_offset(spectra$spc)
  
  # Merge with Lab Data
  lab <- data.frame(read.csv(LABPATH))
  speclib <- merge(lab, spectra, all.y=TRUE)
  
  # Optional Save after Processing
  if(SAVENAME!= "none"){
    assign(SAVENAME,speclib)
    if(file.exists("./Data_Processed")==FALSE){dir.create("./Data_Processed")}
    savepath <- paste0("./Data_Processed/",SAVENAME,".RData")
    save(list=SAVENAME, file=savepath)
    write.csv(get(SAVENAME), savepath, row.names=FALSE)
  }
  return(speclib)
}


# Opus to Dataset
library(stringr) #used for str_sub
library(foreach)
source("Functions/simplerspec/gather-spc.R")
source("Functions/simplerspec/read-opus-universal.R")

opus_to_dataset <- function(SPECPATH ="/Data_Raw/SPECTRA", NWAVE=3017, SAVENAME="none"){
  #--- Converts OPUS files to RData ---#
  
  #---List Files---#
  dirs <- list.dirs(paste(getwd(),SPECPATH,sep=""), full.names=TRUE)
  all.files <- list.files(dirs, pattern= "*.0", recursive=TRUE,full.names=TRUE)
  
  #---Extract Spectra---#
  spc_list <- read_opus_univ(fnames = all.files, extract = c("spc"))
  soilspec_tbl <- spc_list %>%
    gather_spc() # Gather list of spectra data into tibble data frame
  spc <- soilspec_tbl$spc
  
  #---Truncate Spectra---#
  spc.trun <- lapply(1:length(spc),function(x) spc[[x]][,1:NWAVE]) # Truncate at 3017 by default
  
  #---Process to Dataframe---#
  spc.df <- as.data.frame(matrix(unlist(spc.trun), nrow=length(spc.trun), byrow=T))
  colnames(spc.df) <- colnames(spc.trun[[1]]) 
  rownames(spc.df) <- as.character(seq(1,nrow(spc.df)))
  
  #---Assign sample_ids---#
  spc.df <- data.frame(sample_id = soilspec_tbl$sample_id, spc.df)
  spc.df$sample_id <- str_sub(spc.df$sample_id,1,str_length(spc.df$sample_id)-10)
  #spc.avg.df <- aggregate(.~sample_id, data = spc.df, FUN=mean, na.rm=TRUE)
  
  #---Reformat w/ Spectral Matrix Column---#
  spectra <- data.frame(spc.df[,1])
  spectra$spc <- as.matrix(spc.df[,2:ncol(spc.df)])
  colnames(spectra) <- c("sample_id", "spc")
  
  #---Optionally Saves---#
  if(SAVENAME != "none"){
    assign(SAVENAME, spectra)
    savefile <- paste0("Data_Processed/", SAVENAME, ".RData")
    save(list= SAVENAME, file= savefile)
    print(SAVENAME,"saved to", savefile)
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

refineSpecLib <- function(SPECLIB, PROP=NA, OUTLIER=c("stdev"), LARGE=FALSE, CALVAL=FALSE, SAVENAME="none"){
  
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
    SPECLIB  <- SPECLIB[-fratio_outliers(SPECLIB),] # Identified with fratio
  } 
  
  # Subset a large dataset to 15000
  if(LARGE==TRUE){
    SPECLIB$spc <- sub_large_set(SPECLIB) # Subset to 15000 samples
  }
  
  # Split calibration/validation sets
  if(CALVAL==TRUE){
    SPECLIB <- calValSplit(SPECLIB)
  }
  
  # Save the refined reference set for OC
  if(SAVENAME != "none"){
    if(file.exists("./Data_Processed")==FALSE){dir.create("./Data_Processed")}
    assign(SAVENAME, SPECLIB)
    savefile <- paste0("Data_Processed/", SAVENAME, ".RData")
    save(list= SAVENAME, file= savefile)
    print(paste(SAVENAME,"saved to", savefile))
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
                     
sub_large_set <- function(SPECLIB, SUBCOUNT=15000){
  #Conditional Latin Hypercube Sampling if the set exceeds 15000 samples
  spectra <- data.frame(SPECLIB$spc)
  subset <- clhs(spectra, size = SUBCOUNT, progress = TRUE, iter = 500)
  SPECLIB <- SPECLIB[subset,] #double check
  
  return(SPECLIB)
}


# Split into Calibration and Validation Sets
library(prospectr)
calValSplit <- function(SPECLIB, FRAC=0.8){
  
  #perform kennard stone to separate data into 80% calibration and 20% validation sets
  ken_stone<- prospectr::kenStone(X = SPECLIB$spc, 
                                  k = as.integer(FRAC*nrow(SPECLIB)), 
                                  metric = "mahal", pc = 10)
  
  subset <- data.frame(SPECLIB[,1, drop=F])
  subset$calib <- 0
  subset[ken_stone$model, "calib"] <- 1
  SPECLIB <- data.frame(subset,SPECLIB[,2:ncol(SPECLIB),drop=F])
  
  return(SPECLIB)
}









