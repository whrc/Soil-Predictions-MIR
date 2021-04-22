#detect spectral outliers associated with each soil property
setwd("/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/ensemble_example/")
source("functions_udev.R")

lists <- list.files("./spc", pattern = "cal.*.RData", full.names=TRUE)
load("./spc/rod.pds.spc.RData")
out.names <- strsplit(lists, "/")
out.names <- sapply(out.names, function(x) x[[3]]) # 3 is number of directories
out.names <- substr(out.names, 5, nchar(out.names)-6)

for(i in 1:length(lists)){
  test<- get(load(lists[i]))$kssl.spc
  test1 <- rbind(test, rod$pds.spc)
  test.pca <- prcomp(test1)
  all.scores <- test.pca$x
  all.loads <- test.pca$rotation
  
  ##### save pca plot #### for developing maps
  assign(paste0("pca.",out.names[i]), test.pca)
  save(list=paste0("pca.",out.names[i]), 
       file = paste0("./fratio/pca.",out.names[i],".RData"))
  ################################################################
  
  
  spc.recons <- all.scores %*% t(all.loads)
  
  obs.spc <- rbind(test, rod$pds.spc)
  spc.orig.all <- scale(obs.spc, scale =FALSE) #only center the data because prcomp centers the data before pca analysis
  
  spc.res <- (spc.recons - spc.orig.all)^2
  spc.res <- sqrt(rowSums(spc.res))
  
  for(j in 1:length(spc.res)){
    sample.fratio <- (length(spc.res)-1) * spc.res^2/sum((spc.res[-j])^2)
  }
  
  #create the probability distriution of the residuals
  pf.resid <- pf(sample.fratio, 1, length(spc.res)) 
  rows <- which(pf.resid>0.99)
  write.csv(rows, file = paste0("./fratio/",out.names[i], ".csv"))
  cat("Finished", out.names[i], "\n")
  
}

