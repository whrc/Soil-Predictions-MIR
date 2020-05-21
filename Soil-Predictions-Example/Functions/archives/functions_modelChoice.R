
##baseoffset
base_offset <- function(x){
  test <- rowMins(x)
  return(x-test)
  #apply(x,1,function(y) na.omit(y)-min(y, na.rm=TRUE))
}

plot.plsr <- function(x,y, name, xaxes, units){
  x <- as.vector(x)
  y <- as.vector(y)
  lm.model <- lm(x~y)
  RMSE = round(sqrt(mean((y - x)^2)),2)
  print(summary(lm.model))
  axis.max <- max(x,y, na.rm=TRUE)
  axis.min <- min(x,y, na.rm=TRUE)
  bias <- round((sum(y, na.rm=TRUE)- sum(x, na.rm=TRUE))/length(x),2)
  plot(x~y, ylab = paste("Lab Estimated", units), xlab=paste("Predicted", units), xlim = xaxes,ylim=xaxes,main = paste(name))
  coefs <- coef(lm.model)
  b0 <- round(coefs[1], 2)
  b1 <- round(coefs[2],2)
  r2 <- round(summary(lm.model)$r.squared, 2)
  std <- round(sd(x, na.rm=TRUE),2)
  rpd <- round(std / RMSE,2)
  #eqn <- bquote(italic(y) == .(b1)*italic(x) + .(b0) * "," ~~ 
  #                R^2 == .(r2) * "," ~~ RMSE == .(RMSE) * ","~~ RPD== .(rpd))
  eqn <- bquote(italic(bias)== .(bias) * ","~~R^2 == .(r2) * "," ~~ RMSE == .(RMSE))
  abline(lm.model)
  text(min(x,y),max(x,y)-(max(x,y)-min(x,y))/10, eqn, pos = 4, col="blue")
  cat(min(x, na.rm=TRUE), min(y, na.rm=TRUE), max(x, na.rm=TRUE), max(y, na.rm=TRUE), "\n")
  cat(paste("y = ", b1, "x", "+", b0, ", R2 = ", r2, ", RMSE = ", RMSE, ", RPD = ", rpd , ", sd =", std), "\n")
  cat("N=", length(x), "bias =", bias)
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

##functions for boxcox
##functions for boxcox transformation
for.trans <- function(x){
  lam <- powerTransform(x[x>1e-30], family="bcPower")$lambda[[1]]
  test1 <- boxcox(x+1, lambda = lam )
  return(list("var" = test1, "lam"=lam))
}

back.trans <- function(x,lam){   #x is a transformed vector, while is an untransformed vector
  x1 <- boxcox.inv(x, lambda = lam) - 1
  return(x1)
}

##log transformation
for.log <- function(x){
  trans.tmp <- log10(x+1)
  return(trans.tmp)
}

back.log <- function(x){
  back.tmp <- (10^x) - 1 
  return(back.tmp)
}

getsummary <- function(x,y){
  x <- as.vector(x)
  y <- as.vector(y)
  lm.model <- lm(x~y)
  rss <- c(crossprod(lm.model$residuals))
  mse <- rss / length(lm.model$residuals)
  RMSE <- round(sqrt(mse),2)
  bias1 <- round((sum(y, na.rm=TRUE)- sum(x, na.rm=TRUE))/length(x),2)
  bias2 <- sum(y-x, na.rm=TRUE)/length(x)
  r2 <- round(summary(lm.model)$r.squared, 2)
  coefs <- coef(lm.model)
  b0 <- round(coefs[1], 2)
  b1 <- round(coefs[2],2)
  cat("N=", length(x), ", bias1 =", bias1, ", bias2 = ", bias2,", slope =", b1,", offset =", b0, ", R2 = ", r2, ", RMSE = ", RMSE, "\n")
  
}

pfun.lm <- function(x,y, ...){
  panel.hexbinplot(x,y, ...)  
  panel.lmline(x,y,lty=1,lw=2,col="black")
  lm.model <- lm(y~x)
  #lm.model <- lm(x~y)
  rss <- c(crossprod(lm.model$residuals))
  mse <- rss / length(lm.model$residuals)
  rmse <- round(sqrt(mse),2)
  bias <- round(sum(y-x)/length(x),2) #mean of obs - pred
  r2 <- round(summary(lm.model)$r.squared, 2)
  eqn <- bquote(italic(bias)== .(bias) * ","~~R^2 == .(r2) * "," ~~ RMSE == .(rmse))
  panel.text(min(x),max(y),eqn, pos =4, col="blue")
  #panel.text(0,0.25,eqn, pos =4, col="blue")
}