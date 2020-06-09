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
  sd.threshold <- optimum_sd_outlier(pred, SPECLIB[,PROP], seq(0.1,3, by =0.02))
  outliers <- outlier_indices(pred, SPECLIB[,PROP], sd.threshold[1])

  # Display outlier identification
  if(SHOW==TRUE){

    # Show outlier prediction versus observed values
    predobs <- data.frame(SPECLIB[,"sample_id"], round(pred,2), SPECLIB[,PROP])
    names(predobs) <- c("sample_id", "pred", "obs")
    
    cat("\nLab Outliers:\n")
    print(predobs[outliers,])

    # Plot with Outliers
    plot(pred ~ obs, data=predobs[-outliers,], main="Lab Outliers")
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

library(plot3D) # For 3D plot
fratio_outliers <- function(SPECLIB, P=0.99, SHOW=TRUE, PLOT=TRUE){
  
  # Get Principle Component Analysis
  pca <- prcomp(SPECLIB$spc)
  
  # Get Scores
  scores <- pca$x
  
  # Get Loadings
  loadings <- pca$rotation
  
  # Get Predicted Spectra
  pred_spc <- scores %*% t(loadings) 
  
  # Scale Spectra
  spc <- scale(SPECLIB$spc,scale=FALSE)
  
  # Get Residuals
  res <- (pred_spc - spc)^2
  res <- sqrt(rowSums(res))
  
  # Get Fratio
  for(i in 1:length(res)){
    sample.fratio <- (length(res)-1) * res^2/sum((res[-i])^2)
  }
  
  # Get Samples that Exceed the Threshold, P
  ok <- pf(sample.fratio, 1, length(res)) 
  outliers <- which(ok>P)
  
  # Show Results
  if(length(outliers)>0){
    if(PLOT==TRUE){
      # 3D Plot
      x <- pc1 <- scores[,1]
      y <- pc2 <- scores[,2]
      z <- pc3 <- scores[,3]
      scatter3D(x,y,z, col="black", cex = 0.5, main="Spectral Outliers")
      points3D(x[outliers], y[outliers], z[outliers], col="red", pch=16, add=TRUE)
    }
    
    if(SHOW==TRUE){
      # Print Outliers
      cat("\nSpectral Outliers:\n")
      print(data.frame(sample_id=SPECLIB[outliers,1], PF=round(ok[outliers],6)))
    }
    
    return(outliers)
    
  }else{
    print("No Spectral Outliers")
  }
  
}


