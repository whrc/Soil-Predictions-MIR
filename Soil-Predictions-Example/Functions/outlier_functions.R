###--- Outlier Functions ---###

#getOutliers <- function(dataset, column, method){
  #outliers<- NULL
  #if(method=="stdev"){
    #outliers <- stdev_outlier(dataset, column)
  #}
  #return(outliers)
#}

# Packages: stdev_outlier
library(pls)
stdev_outliers <- function(dataset, column){
  #remove outliers by selecting using a standard deviaton threshold
  pls.fit <- plsr(sqrt(get(column))~spc, ncomp= 20, data = dataset, valid="CV", segments = 50) #y, x, number components, data, cross validation, 
  pred <- c(predict(pls.fit, newdata = dataset$spc,ncomp=20))^2
  
  sd.outlier <- optimum_sd_outlier(pred, dataset[,column], seq(0.1,3, by =0.02))
  outlier.indices <- outlier_indices(pred, dataset[,column], sd.outlier[1])
  
  print(paste(column,"outliers"))
  print(outlier.indices)
  
  return(outlier.indices)
}

optimum_sd_outlier <- function(x,y, temp.sd,.....){
  mod <- lm(y~x)
  reg <- fitted(mod)
  stdev <- sd(reg)
  vec <- vector(mode="double", length= length(temp.sd))
  len.outl <- vector(mode = "numeric", length=length(temp.sd))
  for(i in 1:length(temp.sd)) 
  {
    lwr <- reg - temp.sd[i] * stdev
    upr <- reg + temp.sd[i] * stdev
    tmp.index <- which(y < upr & y > lwr)
    len.outl[i] <- length(y) - length(tmp.index)
    vec[i] <- len.outl[i]/length(y) * 100
  }
  sd.index <- which(vec <= 1)[1]
  sd.value <- temp.sd[sd.index]
  #cat(sd.value, len.outl[sd.index],  "\n")
  return(c(sd.value, len.outl[sd.index]))
}

outlier_indices <- function(x,y,sd,.....)
{
  mod <- lm(y~x)
  reg <- fitted(mod)
  stdev <- sd(reg)
  lwr <- reg - sd*stdev
  upr <- reg + sd*stdev
  outlier.rows<- which(y > upr | y < lwr)
  #newd <- datum[datum$sqrt.oc.10.comps<upr & datum$sqrt.oc.10.comps>lwr,]
  return(outlier.rows)
}

