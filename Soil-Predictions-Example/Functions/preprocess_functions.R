###--- Preprocessing functions ---###

# Packages: opus_to_dataset
library(stringr) #used for str_sub
library(foreach)
source("Functions/gather-spc.R")
source("Functions/read-opus-universal.R")

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
  
  # Assigns sample_id from opus file names
  spc.df <- data.frame(sample_id = soilspec_tbl$sample_id, spc.df)
  spc.df$sample_id <- str_sub(spc.df$sample_id,1,-11)
  
  # Reformates dataframe with spectral matrix column
  spectra <- data.frame(spc.df[,1])
  spectra$spc <- as.matrix(spc.df[,2:ncol(spc.df)])
  colnames(spectra) <- c("sample_id", "spc")
  
  if(SAVE==TRUE){
    ##optional save the spectra here
    save(spectra, file="Data_Processed/spectra_original.RData")
    write.csv(spectra, "Data_Processed/spectra_original.csv", row.names=FALSE)
  }
  
  return(spectra)
  
}

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

# Packages: base_offset
library(matrixStats)
base_offset <- function(x){
  test <- rowMins(x)
  return(x-test)
}

subset15000 <- function(spectra){
  #Conditional Latin Hypercube Sampling if the set exceeds 15000 samples
  subset <- clhs(spectra, size = 15000, progress = TRUE, iter = 500)
  spectra <- spectra[subset,] #double check
  return(spectra)
}

#----------------------------------------------#
        # Property Specific Functions #
#----------------------------------------------#

noNA <- function(dataset, column){
  return(dataset[!is.na(dataset[,column]),])
}

noNeg <- function(dataset, column){
  return(dataset[which(dataset[,column] > 0),])
}

source("Functions/outlier_functions.R")
noOut <- function(dataset, column){
  outliers <- stdev_outliers(dataset, column)
  return(dataset[-outliers,])
}

library(prospectr)
calValSplit <- function(dataset, column){
  #perform kennard stone to separate data into 80% calibration and 20% validation sets
  ken_stone<- prospectr::kenStone(X = dataset$spc, 
                                  k = as.integer(0.8*nrow(dataset)), 
                                  metric = "mahal", pc = 10)
  calib <- dataset[ken_stone$model, ]
  valid <- dataset[ken_stone$test, ]
  
  return(c(calib, valid))
}



