#1-Resampling
library(readr)
library(prospectr)
RODL_Spectra <- read_csv("/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/spc/RODL_Spectra_nolabdata.csv", 
                         col_names = TRUE)
dataset <- RODL_Spectra
rm(RODL_Spectra)

nc <- ncol(dataset)
rod <- data.frame(dataset[,c(1,7)])
rod$spc <- as.matrix(dataset[,8:nc])
colnames(rod) <- c("ID","OC","spc")

wav <- as.numeric(substring(colnames(rod$spc),2))
new.wav <- seq(from = 4000, to = 628, by = -2)  #628 cut to
rod$spc <- resample2(rod$spc,wav,new.wav) # 

#2- Cal Transfer
library(ggfortify)
library(prospectr)
library(pls)
source("/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/functions_calTransfer.R")
load("/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/pls.moving.w2k.RData")

ncomps <- 1 # should be less than or equal to window size????
win.size <- 3
y.start <- ifelse(is.even(win.size)== TRUE, win.size/2, round2(win.size/2,0))
y.end <- ncol(rod$spc) - y.start +1
#set parameters for starting and ending col index
#these parameters only changes if window size changes
col.start.index <- y.start
col.end.index <-  ifelse(is.even(win.size)==TRUE, ncol(rod$spc)-y.start, ncol(rod$spc)-y.start+1)
tmp.trans.X1 <- predictPLSR(plsr.moving.w2k,rod$spc, 1, 3, y.start)
trans.X1 <- matrix(unlist(tmp.trans.X1), nrow=nrow(rod$spc))

rod$pds.spc <- trans.X1
rod$spc <- rod$spc[,-c(1,1687)]
colnames(rod$pds.spc) <- colnames(rod$spc)


plot(rod$spc[1,], type ="l", ylim =c(0,2.7))
lines(rod$pds.spc[1,], col="blue")

save(rod, file ="/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/spc/rod.spc.nolab.RData")
#write.csv(rod, file = "/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/rod.spc.csv")