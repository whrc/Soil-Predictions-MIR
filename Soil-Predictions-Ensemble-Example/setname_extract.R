#Rodale Extract
setwd("/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/ensemble_example/")

# Get the ID column
load("./spc/rod.pds.spc.RData")
N <-nrow(rod)
WHRC_ID <- rod[,1]

# Make csv of all organic carbon predictions (column per mbl model)
lists.oc <- list.files("./output/oc", pattern = "*.RData", full.names =TRUE)
output.pred <- matrix(0,N,length(lists.oc))

out.names <- strsplit(lists.oc, "/")
out.names <- sapply(out.names, function(x) x[[4]])
out.names <- substr(out.names, 5, nchar(out.names)-6)
for(i in 1:length(lists.oc)){
  pred<- get(load(lists.oc[i]))
  opt.in <- which.min(pred$localCrossValStats$rmse)
  output.pred[,i] <- pred$results[[opt.in]]$pred^2
}
colnames(output.pred) <- out.names

# Make a column flagging outliers (1 or 0)
output.pred <- data.frame(output.pred)
output.pred$outlier <- 0  # Make a column for tagging outliers
oc <- read.csv("./fratio/oc.csv")
outs <- c(oc$x) 
outs <- outs[outs>15000]
outs <- outs-15000
output.pred$outlier[outs] <- 1
output.pred <- cbind(WHRC_ID,output.pred)
write.csv(output.pred, file ="./output/pred.oc.csv", row.names=FALSE)

#.... {repeat for additional soil properties}

#----------Final Comprehensive File- All Soil Properties------------#

## create a mean prediction with uncertainty estimates
library(Rmisc)
lists.csv <- list.files("./output", pattern = "*.csv", full.names=TRUE, recursive=TRUE)
out.name <- strsplit(lists.csv,"/")
out.name <- sapply(out.name, function(x) x[[3]]) #3 is number of directories
out.name <- substr(out.name, 1, nchar(out.name)-4)
y <- NULL
CI <- 0.95
for(i in 1:length(out.name)){
  dat <- read.csv(lists.csv[i])
  dat.trun <- dat[,-1]
  confint <- t(apply(dat.trun, 1, CI))
  colnames(confint) <- c(paste0("upper.", out.name[i]), paste0("mean.", out.name[i]), paste0("lower.", out.name[i]))
  y <- cbind(y, confint)
}

#add all sample information for validation sets here
ok.all <- cbind(WHRC_ID, y)
write.csv(ok.all, file ="./output/all-predictions.csv", row.names=FALSE)

