#Preprocessing functions

preprocess_spectra <- function(spectraPath="/Data_Raw/SPECTRA"){
  #--- Converts OPUS files to RData ---#
  
  # Load packages
  library(stringr) #used for str_sub
  library(foreach)
  source("Functions/gather-spc.R")
  source("Functions/read-opus-universal.R")
  
  # List directory and files containing MIR in opus format
  dirs <- list.dirs(paste(getwd(),spectraPath,sep=""), full.names=TRUE)
  all.files <- list.files(dirs, pattern= "*.0", recursive=TRUE,full.names=TRUE)
  all.files[1]
  
  # Read opus format files as list and gather spc
  spc_list <- read_opus_univ(fnames = all.files, extract = c("spc"))
  soilspec_tbl <- spc_list %>%
    # Gather list of spectra data into tibble data frame
    gather_spc()
  
  # Process to output spectra in desired format
  spc <- soilspec_tbl$spc
  spc.trun <- lapply(1:length(spc),function(x) spc[[x]][,1:3017]) #for neon we want to make sure to truncate here
  
  spc.df <- as.data.frame(matrix(unlist(spc.trun), nrow=length(spc.trun), byrow=T))
  colnames(spc.df) <- colnames(spc.trun[[1]])
  spc.df <- data.frame(sample_id = soilspec_tbl$sample_id, spc.df)
  spc.df$sample_id <- str_sub(spc.df$sample_id,1,9)
  
  spectra <- data.frame(spc.df[,1])
  spectra$spc <- as.matrix(spc.df[,2:ncol(spc.df)])
  colnames(spectra) <- c("sample_id", "spc")
  
  return(spectra)
  
}

