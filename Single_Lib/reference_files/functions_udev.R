
#1 Xresidual in validation objects
getResXValSamp <- function(x.val.mat,x.cal.mat,x.val.scores,x.cal.loadings){
  nobj <- dim(x.val.mat)[1]
  ncomp <- dim(x.val.scores)[2]
  npred <- dim(x.cal.loadings)[1]
  res.val <- matrix(0, nrow=nobj, ncol=ncomp)
  Xmeans.cal <- colMeans(x.cal.mat) #compare predictors average from calib object with predictors from each object
  # to get predictor redisuals.
  X.center <- x.val.mat - matrix(rep(Xmeans.cal, each = nobj), nrow=nobj)
  for(i in 1:ncomp){
    x.fac.load.wts <- x.val.scores[,1:i, drop=FALSE] %*% t(x.cal.loadings[,1:i,drop=FALSE])
    #res.val[,i] <- rowMeans((-x.fac.load.wts + X.center)^2)
    res.val[,i] <- rowSums((-x.fac.load.wts + X.center)^2)/(npred-i)
    
  }
  return(res.val)
}



#2 Total X residual in Validation sets
#extracted from getXvalres function by using a colMeans function 
#use residuals based on xcal rather than cross-validation to make it more deterministic

getTotResXCal <- function(x.cal.mat,x.cal.scores,x.cal.loadings){
  nobj <- dim(x.cal.mat)[1]
  ncomp <- dim(x.cal.scores)[2]
  npred <- dim(x.cal.loadings)[1]
  res.val <- matrix(0, nrow=nobj, ncol=ncomp)
  #Xmeans.cal <- colMeans(x.cal.mat) #compare predictors average from calib object with predictors from each object
  # to get predictor redisuals.
  #X.center <- x.cal.mat - matrix(rep(Xmeans.cal, each = nobj), nrow=nobj)
  X.center <- scale(x.cal.mat, scale=FALSE)
  for(i in 1:ncomp){
    x.fac.load.wts <- x.cal.scores[,1:i, drop=FALSE] %*% t(x.cal.loadings[,1:i,drop=FALSE])
    #res.val[,i] <- rowMeans((-x.fac.load.wts + X.center)^2)
    res.val[,i] <- rowSums((-x.fac.load.wts + X.center)^2)/(npred-i)
    
  }
  tot.res <- colMeans(res.val)
  return(tot.res)
}


#3
#note that pred.mat should be m * n matrix with m the number of object and n the number of components
getResYValVar <- function(val.resp, pred.mat){
  ncomp <- dim(pred.mat)[2]
  nobj <- dim(pred.mat)[1]
  res.val <- matrix(0, nrow=nobj, ncol=ncomp)
  #Y.center <- scale(val.resp, scale=FALSE)
  res.val <- sapply(1:ncomp, function(x){mean((pred.mat[,x] - val.resp)^2)})
  return(res.val)
} 

#4. Leverage corresponding to best PLS -- This function gives the leverage for each component. The total leverage is
# the sum of individual component leverage for each prediction samples
#************ Gives a cumulative leverage based on the number of PC used ***************#
getLeverage <- function(scores.calib, scores.valid){
  ta.calib <- diag(crossprod(scores.calib))
  ta.calib1 <- matrix(rep(ta.calib, each = nrow(scores.valid)), nrow=nrow(scores.valid))
  ncal <- dim(scores.calib)[1]
  Hi <- scores.valid^2 / ta.calib1
  ncomp <- dim(scores.valid)[2]
  nobj <- dim(scores.valid)[1]
  Hi.pr <- matrix(0, nrow=nobj, ncol=ncomp)
  for(i in 1:ncomp){
    if(i == 1){
      Hi.pr[,1] <- Hi[,1]
    }
    else {
      Hi.pr[,i] <- rowSums(Hi[,1:i]) 
    }
  }
  Hi.pr <- Hi.pr + (1/ncal)
  return(Hi.pr)
}



#5 : Compute prediction error ydev
getYdev <- function(ResYValVar, ResXValSamp, ResXValTot, Hi.pr, ncalobj){
  nobj <- dim(ResXValSamp)[1]
  ncomp <- dim(ResXValSamp)[2]
  ydev <- matrix(0, nrow= nobj, ncol=ncomp)
  for( i in 1:ncomp){  
    ydev[,i] <- sqrt(ResYValVar[i] * (ResXValSamp[,i]/ResXValTot[i] + Hi[,i] + 1/ncalobj) * (1- (i+1)/ncalobj))
  }
  return(ydev)
}

#6: Function to plot observed vs predicted
fit.local <- function(x,y, name, x.lim, y.lim){
  lm.model <- lm(x~y)
  print(summary(lm.model))
  axis.max <- max(x,y, na.rm=TRUE)
  axis.min <- min(x,y, na.rm=TRUE)
  bias <- round((sum(y, na.rm=TRUE)- sum(x, na.rm=TRUE))/length(x),2)
  plot(x~y, ylab = "Measured", xlab="Modeled", pch=16, col= "blue", xlim = c(axis.min, axis.max),main = paste(name))
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
  text(x.lim, y.lim, eqn, pos = 3, col="blue")
  cat(paste("y = ", b1, "x", "+", b0, ", R2 = ", r2, ", RMSE = ", rmse, ", RPD = ", rpd , ", sd =", std), "\n")
  cat("N=", length(x), "bias =", bias)
  #return(lm.model)
}
