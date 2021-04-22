

##functions for prediction of soil properties
##baseoffset -- requires matrixStats package to be uploaded
base_offset <- function(x){
  test <- rowMins(x)
  return(x-test)
  #apply(x,1,function(y) na.omit(y)-min(y, na.rm=TRUE))
}


#process to screenout outliers from the model 
#this function provides optimum sd value to retain the threshold number of samples (in %)
#in this case theshold level is set to 1 (vec <1)
optimum_sd_outlier <- function(x,y, temp.sd,.....)
{
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

#data sets with fitted line +/- 1sd are only retained
outlier <- function(x,y,sd,.....)
{
  mod <- lm(y~x)
  reg <- fitted(mod)
  stdev <- sd(reg)
  lwr <- reg - sd*stdev
  upr <- reg + sd*stdev
  tmp.index <- which(y < upr & y > lwr)
  #newd <- datum[datum$sqrt.oc.10.comps<upr & datum$sqrt.oc.10.comps>lwr,]
  return(tmp.index)
}


fit.local <- function(x,y, name, axis.min, axis.max,x.lim, y.lim){
  lm.model <- lm(x~y)
  print(summary(lm.model))
  #axis.max <- max(x,y, na.rm=TRUE)
  #axis.min <- min(x,y, na.rm=TRUE)
  bias <- round((sum(y, na.rm=TRUE)- sum(x, na.rm=TRUE))/length(x),2)
  plot(x~y, ylab = "Measured", xlab="Modeled", col= "black", xlim = c(axis.min, axis.max), ylim =c(axis.min, axis.max),main = paste(name))
  rmse <- round(sqrt(mean(resid(lm.model)^2)), 2)
  coefs <- coef(lm.model)
  b0 <- round(coefs[1], 2)
  b1 <- round(coefs[2],2)
  r2 <- round(summary(lm.model)$r.squared, 2)
  std <- round(sd(x, na.rm=TRUE),2)
  rpd <- round(std / rmse,2)
  eqn <- bquote(italic(y) == .(b1)*italic(x) + .(b0) * "," ~~ 
                  R^2 == .(r2) * "," ~~ RMSE == .(rmse) * ","~~ RPD== .(rpd))
  abline(lm.model)
  abline(0,1)
  text(x.lim, y.lim, eqn, pos = 4, col="blue")
  cat(paste("y = ", b1, "x", "+", b0, ", R2 = ", r2, ", RMSE = ", rmse, ", sd =", std), "\n")
  cat("N=", length(x), "bias =", bias)
}
