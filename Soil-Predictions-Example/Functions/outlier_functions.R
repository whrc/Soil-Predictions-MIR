###--- Outlier Functions ---###

#----------------------------------------------#
# Standard Deviation Outlier Functions #
#----------------------------------------------#

# Standard Deviation/ Lab Data Outliers
library(pls)

stdev_outliers <- function(SPECLIB, PROP, SHOW=TRUE, PLOT=TRUE){
  
  # Create a PLS model with the data
  pls.fit <- plsr(sqrt(get(PROP))~spc, ncomp= 20, data = SPECLIB, valid="CV", segments = 50) #y, x, number components, data, cross validation, 
  pred <- c(predict(pls.fit, newdata = SPECLIB$spc,ncomp=20))^2
  
  # Identify outliers using a standard deviation threshold
  sd.outlier <- optimum_sd_outlier(pred, SPECLIB[,PROP], seq(0.1,3, by =0.02))
  outliers <- outlier_indices(pred, SPECLIB[,PROP], sd.outlier[1])
  
  # Display outlier identification
  if(SHOW==TRUE){
    
    # Show outlier prediction versus observed values
    predobs <- data.frame(pred, SPECLIB[,PROP])
    names(predobs) <- c("pred", "obs")
    print(predobs[outliers,])
    
    # Plot with Outliers
    plot(pred ~ obs, data=predobs[-outliers,], main="Outliers")
    points(pred ~ obs, data=predobs[outliers,], col="red", pch=24, bg="red")
    
  }
  
  return(outliers)
}


# Find Standard Deviation that Gets Outer 1% of samples

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


# Get Indices of Outer 1% of Samples

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


#----------------------------------------------#
# FRatio Outlier Functions #
#----------------------------------------------#

