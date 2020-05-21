#Indigo Extraction
load("./Data_Processed/spectra_lab_merge.RData")
sample_id = all_data$sample_id

lists.prop <- list.files("./Models", pattern = paste("*", property,".RData", sep=""), full.names =TRUE)
output.pred <- matrix(0,length(sample_id),length(lists.prop))
out.names <- strsplit(lists.anal4, "/")
out.names <- sapply(out.names, function(x) x[[11]])
out.names <- substr(out.names, 5, nchar(out.names)-6)
for(i in 1:length(lists.anal4)){
  pred<- get(load(lists.anal4[i]))
  opt.in <- which.min(pred$localCrossValStats$rmse)
  output.pred[,i] <- pred$results[[opt.in]]$pred^2
  
}
colnames(output.pred) <- out.names

#outlier flagging
output.pred <- data.frame(output.pred)
output.pred$outlier <- 0
oc <- read.csv("/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/output/indigo_output/fratio/bd.anal4.csv")
outs <- c(oc$x) 
outs <- outs[outs>10363]
outs <- outs-10363
output.pred$outlier[outs] <- 1
output.pred <- cbind(TERR_ID,output.pred)

write.csv(output.pred, file ="/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/output/indigo_output/all-predictions/pred.bd.anal4.csv")



#bd analyte 21
lists.anal21 <- list.files("/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/output/indigo_output/analyte21", pattern = "*.RData", full.names =TRUE)
output.pred <- matrix(0,661,length(lists.anal21))

out.names <- strsplit(lists.anal21, "/")
out.names <- sapply(out.names, function(x) x[[11]])
out.names <- substr(out.names, 5, nchar(out.names)-6)
for(i in 1:length(lists.anal21)){
  pred<- get(load(lists.anal21[i]))
  opt.in <- which.min(pred$localCrossValStats$rmse)
  output.pred[,i] <- pred$results[[opt.in]]$pred^2
  
}
colnames(output.pred) <- out.names

#outlier flagging
output.pred <- data.frame(output.pred)
output.pred$outlier <- 0
oc <- read.csv("/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/output/indigo_output/fratio/bd.anal21.csv")
outs <- c(oc$x) 
outs <- outs[outs>6890]
outs <- outs-6890
output.pred$outlier[outs] <- 1
output.pred <- cbind(TERR_ID,output.pred)

write.csv(output.pred, file ="/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/output/indigo_output/all-predictions/pred.bd.anal21.csv")



#combined bd
lists.anal21 <- list.files("/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/output/indigo_output/bd", pattern = "*.RData", full.names =TRUE)
output.pred <- matrix(0,661,length(lists.anal21))

out.names <- strsplit(lists.anal21, "/")
out.names <- sapply(out.names, function(x) x[[11]])
out.names <- substr(out.names, 5, nchar(out.names)-6)
for(i in 1:length(lists.anal21)){
  pred<- get(load(lists.anal21[i]))
  opt.in <- which.min(pred$localCrossValStats$rmse)
  output.pred[,i] <- pred$results[[opt.in]]$pred^2
  
}
colnames(output.pred) <- out.names

#outlier flagging
output.pred <- data.frame(output.pred)
output.pred$outlier <- 0
oc <- read.csv("/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/output/indigo_output/fratio/bd.csv")
outs <- c(oc$x) 
outs <- outs[outs>13667]
outs <- outs-13667
output.pred$outlier[outs] <- 1
output.pred <- cbind(TERR_ID,output.pred)

write.csv(output.pred, file ="/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/output/indigo_output/all-predictions/pred.bd.all.csv")


#CEC
lists.cec <- list.files("/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/output/indigo_output/cec", pattern = "*.RData", full.names =TRUE)
output.pred <- matrix(0,661,length(lists.cec))

out.names <- strsplit(lists.cec, "/")
out.names <- sapply(out.names, function(x) x[[11]])
out.names <- substr(out.names, 5, nchar(out.names)-6)
for(i in 1:length(lists.cec)){
  pred<- get(load(lists.cec[i]))
  opt.in <- which.min(pred$localCrossValStats$rmse)
  output.pred[,i] <- pred$results[[opt.in]]$pred^2
  
}
colnames(output.pred) <- out.names

#outlier flagging
output.pred <- data.frame(output.pred)
output.pred$outlier <- 0
oc <- read.csv("/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/output/indigo_output/fratio/cec.csv")
outs <- c(oc$x) 
outs <- outs[outs>15000]
outs <- outs-15000
output.pred$outlier[outs] <- 1
output.pred <- cbind(TERR_ID,output.pred)

write.csv(output.pred, file ="/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/output/indigo_output/all-predictions/pred.cec.csv")




#Clay
lists.clay <- list.files("/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/output/indigo_output/clay", pattern = "*.RData", full.names =TRUE)
output.pred <- matrix(0,661,length(lists.clay))

out.names <- strsplit(lists.clay, "/")
out.names <- sapply(out.names, function(x) x[[11]])
out.names <- substr(out.names, 5, nchar(out.names)-6)
for(i in 1:length(lists.clay)){
  pred<- get(load(lists.clay[i]))
  opt.in <- which.min(pred$localCrossValStats$rmse)
  output.pred[,i] <- pred$results[[opt.in]]$pred^2
  
}
colnames(output.pred) <- out.names

#outlier flagging
output.pred <- data.frame(output.pred)
output.pred$outlier <- 0
oc <- read.csv("/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/output/indigo_output/fratio/clay.csv")
outs <- c(oc$x) 
outs <- outs[outs>15000]
outs <- outs-15000
output.pred$outlier[outs] <- 1
output.pred <- cbind(TERR_ID,output.pred)

write.csv(output.pred, file ="/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/output/indigo_output/all-predictions/pred.clay.csv")




#organic carbon
lists.oc <- list.files("/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/output/indigo_output/oc", pattern = "*.RData", full.names =TRUE)
output.pred <- matrix(0,661,length(lists.oc))

out.names <- strsplit(lists.oc, "/")
out.names <- sapply(out.names, function(x) x[[11]])
out.names <- substr(out.names, 5, nchar(out.names)-6)
for(i in 1:length(lists.oc)){
  pred<- get(load(lists.oc[i]))
  opt.in <- which.min(pred$localCrossValStats$rmse)
  output.pred[,i] <- pred$results[[opt.in]]$pred^2
  
}
colnames(output.pred) <- out.names

#outlier flagging
output.pred <- data.frame(output.pred)
output.pred$outlier <- 0
oc <- read.csv("/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/output/indigo_output/fratio/oc.csv")
outs <- c(oc$x) 
outs <- outs[outs>15000]
outs <- outs-15000
output.pred$outlier[outs] <- 1
output.pred <- cbind(TERR_ID,output.pred)

write.csv(output.pred, file ="/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/output/indigo_output/all-predictions/pred.oc.csv")



#OCDensity
lists.ocden <- list.files("/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/output/indigo_output/ocden", pattern = "*.RData", full.names =TRUE)
output.pred <- matrix(0,661,length(lists.ocden))

out.names <- strsplit(lists.ocden, "/")
out.names <- sapply(out.names, function(x) x[[11]])
out.names <- substr(out.names, 5, nchar(out.names)-6)
for(i in 1:length(lists.ocden)){
  pred<- get(load(lists.ocden[i]))
  opt.in <- which.min(pred$localCrossValStats$rmse)
  output.pred[,i] <- pred$results[[opt.in]]$pred^2
  
}
colnames(output.pred) <- out.names

#outlier flagging
output.pred <- data.frame(output.pred)
output.pred$outlier <- 0
oc <- read.csv("/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/output/indigo_output/fratio/ocden.csv")
outs <- c(oc$x) 
outs <- outs[outs>12374]
outs <- outs-12374
output.pred$outlier[outs] <- 1
output.pred <- cbind(TERR_ID,output.pred)

write.csv(output.pred, file ="/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/output/indigo_output/all-predictions/pred.ocden.csv")



#Carbonates -- remember the number of samples -- we removed samples that have no detectable limits of co3
lists.co3 <- list.files("/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/output/indigo_output/co3", pattern = "*.RData", full.names =TRUE)
output.pred <- matrix(0,660,length(lists.co3))

out.names <- strsplit(lists.co3, "/")
out.names <- sapply(out.names, function(x) x[[11]])
out.names <- substr(out.names, 5, nchar(out.names)-6)
for(i in 1:length(lists.co3)){
  pred<- get(load(lists.co3[i]))
  opt.in <- which.min(pred$localCrossValStats$rmse)
  output.pred[,i] <- pred$results[[opt.in]]$pred^2
  
}
colnames(output.pred) <- out.names

#outlier flagging
output.pred <- data.frame(output.pred)
output.pred$outlier <- 0
oc <- read.csv("/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/output/indigo_output/fratio/caco3.trun.csv")
outs <- c(oc$x) 
outs <- outs[outs> 14003]
outs <- outs-14003
output.pred$outlier[outs] <- 1
output.pred <- cbind(TERR_ID,output.pred)

write.csv(output.pred, file ="/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/output/indigo_output/all-predictions/pred.co3.csv")
#15000, 14003 (trun), 869 val



#total Nitrogen
lists.totn <- list.files("/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/output/indigo_output/totN", pattern = "*.RData", full.names =TRUE)
output.pred <- matrix(0,661,length(lists.totn))

out.names <- strsplit(lists.totn, "/")
out.names <- sapply(out.names, function(x) x[[11]])
out.names <- substr(out.names, 5, nchar(out.names)-6)
for(i in 1:length(lists.totn)){
  pred<- get(load(lists.totn[i]))
  opt.in <- which.min(pred$localCrossValStats$rmse)
  output.pred[,i] <- pred$results[[opt.in]]$pred^2
  
}
colnames(output.pred) <- out.names

#outlier flagging
output.pred <- data.frame(output.pred)
output.pred$outlier <- 0
oc <- read.csv("/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/output/indigo_output/fratio/totN.csv")
outs <- c(oc$x) 
outs <- outs[outs>15000]
outs <- outs-15000
output.pred$outlier[outs] <- 1
output.pred <- cbind(TERR_ID,output.pred)

write.csv(output.pred, file ="/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/output/indigo_output/all-predictions/pred.totn.csv")


#total Sand
lists.sand <- list.files("/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/output/indigo_output/sand", pattern = "*.RData", full.names =TRUE)
output.pred <- matrix(0,661,length(lists.sand))

out.names <- strsplit(lists.sand, "/")
out.names <- sapply(out.names, function(x) x[[11]])
out.names <- substr(out.names, 5, nchar(out.names)-6)
for(i in 1:length(lists.sand)){
  pred<- get(load(lists.sand[i]))
  opt.in <- which.min(pred$localCrossValStats$rmse)
  output.pred[,i] <- pred$results[[opt.in]]$pred^2
  
}
colnames(output.pred) <- out.names

#outlier flagging
output.pred <- data.frame(output.pred)
output.pred$outlier <- 0
oc <- read.csv("/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/output/indigo_output/fratio/sand.csv")
outs <- c(oc$x) 
outs <- outs[outs>15000]
outs <- outs-15000
output.pred$outlier[outs] <- 1
output.pred <- cbind(TERR_ID,output.pred)

write.csv(output.pred, file ="/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/output/indigo_output/all-predictions/pred.sand.csv")




#total silt
lists.silt <- list.files("/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/output/indigo_output/silt", pattern = "*.RData", full.names =TRUE)
output.pred <- matrix(0,661,length(lists.silt))

out.names <- strsplit(lists.silt, "/")
out.names <- sapply(out.names, function(x) x[[11]])
out.names <- substr(out.names, 5, nchar(out.names)-6)
for(i in 1:length(lists.silt)){
  pred<- get(load(lists.silt[i]))
  opt.in <- which.min(pred$localCrossValStats$rmse)
  output.pred[,i] <- pred$results[[opt.in]]$pred^2
  
}
colnames(output.pred) <- out.names

#outlier flagging
output.pred <- data.frame(output.pred)
output.pred$outlier <- 0
oc <- read.csv("/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/output/indigo_output/fratio/silt.csv")
outs <- c(oc$x) 
outs <- outs[outs>15000]
outs <- outs-15000
output.pred$outlier[outs] <- 1
output.pred <- cbind(TERR_ID,output.pred)

write.csv(output.pred, file ="/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/output/indigo_output/all-predictions/pred.silt.csv")





## create a mean prediction with uncertainty estimates
lists.csv <- list.files("/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/output/indigo_output/all-predictions", pattern = "*.csv", full.names=TRUE, recursive=TRUE)
#remove co3 since samples does not match
#lists.csv <- lists.csv[9:11]
lists.csv <- lists.csv[-6]
library(Rmisc)

out.name <- strsplit(lists.csv,"/")
out.name <- sapply(out.name, function(x) x[[11]])
out.name <- substr(out.name, 1, nchar(out.name)-4)
y <- NULL
i <- 4
CI <- 0.95
#for(i in 1:length(out.name)){
dat <- read.csv(lists.csv[i])
dat.trun <- dat[,-1]
dat.trun <- dat.trun[-1]
dat.trun <- na.omit(dat.trun)
confint <- t(apply(dat.trun, 1, CI))
colnames(confint) <- c(paste0("upper.", out.name[i]), paste0("mean.", out.name[i]), paste0("lower.", out.name[i]))
y <- cbind(y, confint)
#}


#add all sample information for validation sets here
ok.all <- cbind(terra.spc[,1:2], y)
write.csv(ok.all, file ="/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/output/indigo_output/all-predictions/pred.soilproperties.csv")
write.csv(ok.all, file ="/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/output/indigo_output/all-predictions/pred.ceconly.csv")
#write.csv(ok.all, file ="/mnt/data2/disk1/sdangal/KSSL/Terraton/finalout/pred.sandsilttotN.csv")


#get co3 separately
co3 <- read.csv("/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/output/indigo_output/all-predictions/pred.co3.csv")
co3.trun <- co3[,-1]
confint <- t(apply(co3.trun,1,CI))
colnames(confint) <- c("upper.co3", "mean.co3", "lower.co3")
#merge this with co3 cal
load("/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/spc/indigo.val.caco3.RData")
ok.co3 <- cbind(val.caco3[,1:2], confint)
write.csv(ok.co3, file ="/mnt/data2/disk1/soi