#Converts opus files to RData
library(stringr) #used for str_sub
source(gather-spc.R)
source(read-opus-universal.R)

#list directory and files containing MIR in opus format
spectraPath <- "/spectra"
dirs <- list.dirs(paste(getwd(),spectraPath,sep=""), full.names=TRUE)
all.files <- list.files(dirs, pattern= "*.0", recursive=TRUE,full.names=TRUE)

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
spc.df$sample_id <- str_sub(spc.df$sample_id,5,9)

write.csv(spc.df, file="./Models/ASCC/opus-files/spec.csv", row.names=FALSE)

