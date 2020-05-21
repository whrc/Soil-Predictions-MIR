### Data Preprocessing

#Convert OPUS files to RData
library(stringr) #used for str_sub
library(foreach)
source("Functions/gather-spc.R")
source("Functions/read-opus-universal.R")

#List directory and files containing MIR in opus format
spectraPath <- "/Data_Raw/SPECTRA"
dirs <- list.dirs(paste(getwd(),spectraPath,sep=""), full.names=TRUE)
all.files <- list.files(dirs, pattern= "*.0", recursive=TRUE,full.names=TRUE)
all.files[1]

#Read opus format files as list and gather spc
spc_list <- read_opus_univ(fnames = all.files, extract = c("spc"))
soilspec_tbl <- spc_list %>%
  #Gather list of spectra data into tibble data frame
  gather_spc()

#Process to output spectra in desired format
spc <- soilspec_tbl$spc
spc.trun <- lapply(1:length(spc),function(x) spc[[x]][,1:3017]) #for neon we want to make sure to truncate here

spc.df <- as.data.frame(matrix(unlist(spc.trun), nrow=length(spc.trun), byrow=T))
colnames(spc.df) <- colnames(spc.trun[[1]])
spc.df <- data.frame(sample_id = soilspec_tbl$sample_id, spc.df)
spc.df$sample_id <- str_sub(spc.df$sample_id,1,9)

spectra <- data.frame(spc.df[,1])
spectra$spc <- as.matrix(spc.df[,2:ncol(spc.df)])
colnames(spectra) <- c("sample_id", "spc")

##optional save the spectra here
save(spectra, file="Data_Processed/spectra_original.RData")
write.csv(spectra, "Data_Processed/spectra_original.csv", row.names=FALSE)

#Edit Spectral Columns
col.names <- colnames(spectra$spc)
col.names <- as.numeric(substring(col.names,2))

cutoff <- which(col.names <= 628)[1]
spectra$spc <- spectra$spc[,-c(cutoff:length(col.names))] #truncate at >= 628

min.index <- which(col.names <= 2389)[1]
max.index <- which(col.names <= 2268)[1]
spectra$spc <- spectra$spc[,-c(min.index:max.index)] #remove CO2 region

#Baseline Transformation
library(matrixStats)
base_offset <- function(x){
  test <- rowMins(x)
  return(x-test)
}
spectra$spc <- base_offset(spectra$spc)

#Optional save of spectra after processing
save(spectra, file="Data_Processed/spectra_processed.RData")
write.csv(spectra, "Data_Processed/spectra_processed.csv", row.names=FALSE)

#Merge with Lab Data
library(readr)
lab <- data.frame(read_csv("LAB_DATA.csv"))
all_data <- merge(lab, spectra, all.y=TRUE)

#Save all_data file after preprocessing
save(all_data, file="Data_Processed/spectra_lab_merge.RData")
write.csv(all_data, "Data_Processed/spectra_lab_merge.csv", row.names=FALSE)

#________________________________________________________#

#For EACH Soil Property:
properties <- c("OC","CLAY","SAND","SILT") #Column names of lab data

for(property in properties){
  
  #Get rid of rows with no validation data
  all_data <- all_data[!is.na(all_data[,property]),] #no NAs
  all_data <- all_data[which(all_data[,property] > 0),] #no negative values
  
  #outlier removal
  #remove outliers by selecting using a standard deviaton threshold
  library(pls)
  pls.fit <- plsr(sqrt(get(property))~spc, ncomp= 20, data = all_data, valid="CV", segments = 50) #y, x, number components, data, cross validation, 
  pred <- c(predict(pls.fit, newdata = all_data$spc,ncomp=20))^2
  
  source("Reference_Files/functions_modelChoice.R")
  sd.outlier <- optimum_sd_outlier(pred, all_data[,property], seq(0.1,3, by =0.02))
  row.index <- outlier(pred, all_data[,property], sd.outlier[1])
  if(length(row.index) > 0){
    all_data <- all_data[row.index,]
  }
  
  #Subset if dataset is very large
  if(FALSE){
    #{Optional} Conditional Latin Hypercube Sampling if the set exceeds 15000 samples
    spectra <- data.frame(all_data$spc)
    subset <- clhs(spectra, size = 15000, progress = TRUE, iter = 500)
    #another line for using the indices of the clhs to subset, double check
    all_data <- all_data[subset,]
  }
  
  #Save processed datasets for each property
  save(all_data, file=paste("Data_Processed/all",property,"RData", sep="."))
  
  #{Optional} Split calibration and validation sets
  if(FALSE){
    #perform kennard stone to separate data into 80% calibration and 20% validation sets
    ken_stone<- prospectr::kenStone(X = all_data$spc, 
                                    k = as.integer(0.8*nrow(all_data)), 
                                    metric = "mahal", pc = 10)
    calib <- all_data[ken_stone$model, ]
    valid <- all_data[ken_stone$test, ]
    
    #Save for property
    save(calib, file=paste("Data_Processed/calib",property,"RData", sep="."))
    save(valid, file=paste("Data_Processed/valid",property,"RData", sep="."))
  }
  
}
