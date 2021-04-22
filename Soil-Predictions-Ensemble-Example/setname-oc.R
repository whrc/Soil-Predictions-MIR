library(matrixStats)
library(resemble)
setwd("/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/ensemble_example/")
source("functions_soilp.R")


#------Ensemble Combinations------#
## define parameters and method for local model -- these parameters are same across all combination of models
## static parameters
diss2test <- seq(0.5, 3, by = 0.5)
kminmax <- c(100, 300)
pls.f <- c(minpls = 7, maxpls = 12)

k.neigh <- seq(50,200,by = 25)
pls.f1 <- 12

#define methods for regression model
## Regression method: Weighted average pls (wapls)
rmethod <- c("wapls1","pls")

#define paramters that change for different combination of models
# the dynamics are implemented as different approaches to test the similarity
sm2test <- c("euclid", "cosine", "cor", "pc", "pls")

#define selection of predictor variables
diss2use <- c("none", "predictors")

para.grid <- expand.grid(sm2test, diss2use, rmethod)
out.names <- sapply(1:nrow(para.grid), function(x) {paste0("mbl.", para.grid[x,1], ".", para.grid[x,2],".", para.grid[x,3])})
para.grid # <--- To view the combinations the script will run


#--------Inputs----------#
load("./spc/cal.oc.RData")
all.set.cal <- cal.oc
rm(cal.oc)
load("./spc/rod.pds.spc.RData")
all.set.val <- rod
rm(rod)

#remove CO2 sensitive region calibration
col.names <- as.numeric(colnames(all.set.cal$kssl.spc))
min.index <- which(col.names <= 2389)[1]
max.index <- which(col.names <= 2268)[1]
all.set.cal$kssl.spc <- all.set.cal$kssl.spc[,-c(min.index:max.index)] 
all.set.val$pds.spc <- all.set.val$pds.spc[,-c(min.index:max.index)]

#base offset
all.set.cal$kssl.spc <- base_offset(all.set.cal$kssl.spc) #The calibration set (KSSL) spectra
all.set.val$pds.spc <- base_offset(all.set.val$pds.spc) #The transformed spectra to make predictions from

#inputs matrix -- these matrix are also same across all combination of models
#start building a local model
Xu <- all.set.val$pds.spc
#Yu <- sqrt(all.set.val$OC) #If we have lab data
Yr <- sqrt(all.set.cal$OC)
Xr <- all.set.cal$kssl.spc

#Xu <- Xu[!is.na(Yu),]
Xr <- Xr[!is.na(Yr),]
#Yu <- Yu[!is.na(Yu)]
Yr <- Yr[!is.na(Yr)]


#--------Running MBL Models----------#
#parallelize the models
#library(parallel)
#library(doSNOW)
#cores <- detectCores() - 1 
#cl <-makeCluster(cores) 
#print(cl, "\n")
#registerDoSNOW(cl)


#pred.all <- foreach (i = 1:nrow(para.grid), .packages=c("resemble"), .combine ='cbind')%dopar%{
  #loop through parameters
for(i in 1:nrow(para.grid)){
#for(i in 1:2){
  sm2test <- paste0(para.grid[i,1])
  diss2use <- paste0(para.grid[i,2])
  rmethod <- paste0(para.grid[i,3])
  
  ctrl <- mbl_control(validation_type = 'local_cv',allow_parallel=FALSE) #Resemble 2.0
  #ctrl <- mblControl(sm = sm2test, pcSelection = list('opc', 40),
                     #valMethod = 'loc_crossval',center=TRUE,scale=FALSE,allowParallel=FALSE) 
                      
  if(rmethod == "wapls1"){
    mbl.out <- mbl(Yr = Yr, Xr = Xr, Yu = NULL, Xu = Xu,
                   mblCtrl = ctrl, method=local_fit_wapls(min_pls_c = 7, max_pls_c = 12),
                   diss_usage = diss2use, k_diss = diss2test, k_range = kminmax, diss_method = sm2test, center=TRUE, scale=FALSE)
    
    #mbl.out <- mbl(Yr = Yr, Xr = Xr, Yu = NULL, Xu = Xu,
                   #mblCtrl = ctrl,
                   #dissUsage = diss2use, k.diss = diss2test, k.range = kminmax, 
                   #pls.c = pls.f, method = rmethod)
  }else {
    mbl.out <- mbl(Yr = Yr, Xr = Xr, Yu = NULL, Xu = Xu,
                   mblCtrl = ctrl, method=local_fit_pls(pls_c = pls.f1),
                   diss_usage = diss2use, diss_method = sm2test, k = k.neigh, center=TRUE, scale=FALSE)
    
    #mbl.out <- mbl(Yr = Yr, Xr = Xr, Yu = NULL, Xu = Xu,
                   #mblCtrl = ctrl,
                   #dissUsage = diss2use, k = k.neigh, 
                   #pls.c = pls.f1, method = rmethod)
  }
  
  #assign(paste0("mbl.",sm2test[i]), mbl.out)
  assign(paste0("mbl.",rmethod,".",sm2test,".",diss2use), mbl.out)
  
  save(list=paste0("mbl.",rmethod,".",sm2test,".",diss2use), 
       file = paste0("./output/oc/mbl.",rmethod,".",sm2test,".",diss2use,".RData"))
  opt.ind <- which.min(mbl.out$localCrossValStats$rmse)
  pred.oc <- c(mbl.out$results[[opt.ind]]$pred^2)
}

#colnames(pred.all) <- out.names
#pred.all <- cbind(pred.all, obs.oc =Yu^2)
#write.csv(pred.all, file = "/mnt/data2/disk1/soilcarbon/crivard/predEnsemble/pred.mbl.whrc.pds.csv")
#stopCluster(cl)

