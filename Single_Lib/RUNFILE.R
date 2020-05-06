### Full Script Data Preprocessing
setwd("~/Desktop/Example-Models/bookdown-MIR/Single_Lib")

#Converts opus files to RData
library(stringr) #used for str_sub
library(foreach)
source("reference_files/gather-spc.R")
source("reference_files/read-opus-universal.R")

#list directory and files containing MIR in opus format
spectraPath <- "/SPECTRA"
dirs <- list.dirs(paste(getwd(),spectraPath,sep=""), full.names=TRUE)
all.files <- list.files(dirs, pattern= "*.0", recursive=TRUE,full.names=TRUE)
all.files[1]

# read opus format files as list and gather spc
spc_list <- read_opus_univ(fnames = all.files, extract = c("spc"))
soilspec_tbl <- spc_list %>%
  # Gather list of spectra data into tibble data frame
  gather_spc()

#process to output spectra in desired format, average four replicates and save all spectra as a dataframe to build models
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
save(spectra, file="spectra_original.RData")
write.csv(spectra, "spectra_original.csv", row.names=FALSE)

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

#optional save of spectra after processing
save(spectra, file="spectra_processed.RData")
write.csv(spectra, "spectra_processed.csv", row.names=FALSE)

# Merge with Lab Data
lab <- data.frame(read_csv("LAB_DATA.csv"))
all_data <- merge(lab, spectra, all.y=TRUE)

#save all_data file after preprocessing
save(all_data, file="spectra_lab_merge.RData")
write.csv(all_data, "spectra_lab_merge.csv", row.names=FALSE)

#-----------------------------------------------#

#FOR EACH SOIL PROPERTY!
property <- "OC"

#Get rid of rows with no validation data
all_data <- all_data[!is.na(all_data[,property]),] #no NAs
all_data <- all_data[which(all_data[,property] > 0),] #no negative values

#outlier removal
#remove outliers by selecting using a standard deviaton threshold
#Functions in "functions_modelChoice.R"
library(pls)
pls.fit <- plsr(sqrt(get(property))~spc, ncomp= 20, data = all_data, valid="CV", segments = 50) #y, x, number components, data, cross validation, 
pred <- c(predict(pls.fit, newdata = all_data$spc,ncomp=20))^2

source("reference_files/functions_modelChoice.R")
sd.outlier <- optimum_sd_outlier(pred, all_data[,property], seq(0.1,3, by =0.02))
row.index <- outlier(pred, all_data[,property], sd.outlier[1])
if(length(row.index) > 0){
  all_data <- all_data[row.index,]
}

#Split calibration and validation sets
#perform kennard stone to separate data into 80% calibration and 20% validation sets
ken_stone<- prospectr::kenStone(X = all_data$spc, k = as.integer(0.8*nrow(all_data)), metric = "mahal", pc = 10) ## repeat this step for other soil properties -- Al, Ca, CO3, pH, Fe, Ca, BD, OCD, Clay.
calib <- all_data[ken_stone$model, ]
valid <- all_data[ken_stone$test, ]

#Save for property
save(calib, file=paste("calib",property,"RData", sep="."))
save(valid, file=paste("valid",property,"RData", sep="."))


#________________________________________________________#
properties <- c("CLAY","SAND","SILT")
#FOR EACH SOIL PROPERTY!
#property <- "OC"
for(property in properties){
  
  #Get rid of rows with no validation data
  all_data <- all_data[!is.na(all_data[,property]),] #no NAs
  all_data <- all_data[which(all_data[,property] > 0),] #no negative values
  
  #outlier removal
  #remove outliers by selecting using a standard deviaton threshold
  library(pls)
  pls.fit <- plsr(sqrt(get(property))~spc, ncomp= 20, data = all_data, valid="CV", segments = 50) #y, x, number components, data, cross validation, 
  pred <- c(predict(pls.fit, newdata = all_data$spc,ncomp=20))^2
  
  source("Single_Lib/reference_files/functions_modelChoice.R")
  sd.outlier <- optimum_sd_outlier(pred, all_data[,property], seq(0.1,3, by =0.02))
  row.index <- outlier(pred, all_data[,property], sd.outlier[1])
  if(length(row.index) > 0){
    all_data <- all_data[row.index,]
  }
  
  #Split calibration and validation sets
  #perform kennard stone to separate data into 80% calibration and 20% validation sets
  ken_stone<- prospectr::kenStone(X = all_data$spc, k = as.integer(0.8*nrow(all_data)), metric = "mahal", pc = 10) ## repeat this step for other soil properties -- Al, Ca, CO3, pH, Fe, Ca, BD, OCD, Clay.
  calib <- all_data[ken_stone$model, ]
  valid <- all_data[ken_stone$test, ]
  
  #Save for property
  save(calib, file=paste("Single_Lib/calib",property,"RData", sep="."))
  save(valid, file=paste("Single_Lib/valid",property,"RData", sep="."))
  
}





"https://stackoverflow.com/questions/750786/whats-the-best-way-to-use-r-scripts-on-the-command-line-terminal"


