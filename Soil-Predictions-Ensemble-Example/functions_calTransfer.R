
#***** functions ***************#
round2 = function(x, n) {
  posneg = sign(x)
  z = abs(x)*10^n
  z = z + 0.5
  z = trunc(z)
  z = z/10^n
  z*posneg
}

mean_center <- function(x){
  apply(x,2,function(y) na.omit(y)-mean(y, na.rm=TRUE))
}

#base_offset <- function(x){
#  apply(x,1,function(y) na.omit(y)-min(y, na.rm=TRUE))
#}

##baseoffset
base_offset <- function(x){
  test <- rowMins(x)
  return(x-test)
  #apply(x,1,function(y) na.omit(y)-min(y, na.rm=TRUE))
}


#define a function
movingWindowPLSR <- function(X,Y,ncomp,width,ystart, ...) {
  nr <- ncol(as.matrix(X))
  if (width <= 0 || width >= nr)
    stop(paste("width must be in the range 1,...,", nr, sep=""))
  nreg = nr - width + 1
  base = 0:(width - 1)
  sumrys <- lapply(1:nreg,
                   function(st) {
                     plsr(Y[,st+ystart-1]~as.matrix(X[,base+st]), ncomp)
                   })
  sumrys
}



predictPLSR <- function(plsr.moving,X,ncomp,width,ystart, ...) {
  nr <- ncol(as.matrix(X))
  if (width <= 0 || width >= nr)
    stop(paste("width must be in the range 1,...,", nr, sep=""))
  nreg = nr - width + 1
  print(nreg)
  base = 0:(width - 1)
  sumrys <- lapply(1:nreg,
                   function(st) {
                     predict(plsr.moving[[st+ystart-1]], newdata = as.matrix(X[,base+st]), ncomp)
                   })
  
  sumrys
}



is.even <- function(x) x%%2 == 0



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



#define new functions for transformations
#function to transform the data
for.trans <- function(x){
  x1 <- sqrt(x)
  lam <- powerTransform(x1[x1>1e-30], family="bcPower")$lambda[[1]]
  test1 <- boxcox(x1+ 1 - min(x1, na.rm=TRUE), lambda = lam )
  return(test1)
}

back.trans <- function(x,y){   #x is a transformed vector, while is an untransformed vector
  y1 <- sqrt(y)
  lam <- powerTransform(y1[y1>1e-30], family="bcPower")$lambda[[1]]
  x1 <- (boxcox.inv(x, lambda = lam) -1  + min(y1, na.rm=TRUE))
  x2 <- x1^2
  return(x2)
}






#define a function
predictPLSR <- function(obj.name,X,ncomp,width,ystart, ...) {
  nr <- ncol(as.matrix(X))
  width = as.integer(width)[1]
  if (width <= 0 || width >= nr)
    stop(paste("width must be in the range 1,...,", nr, sep=""))
  nreg = nr - width + 1
  print(nreg)
  base = 0:(width - 1)
  sumrys <- lapply(1:nreg,
                   function(st) {
                     predict(obj.name[[st]], newdata = X[,base+st], ncomp)
                   })
  
  sumrys
}

predictPLSR_onewin <- function(plsr.moving,X,ncomp,width,ystart, ...) {
  nr <- ncol(as.matrix(X))
  if (width <= 0 || width >= nr)
    stop(paste("width must be in the range 1,...,", nr, sep=""))
  nreg = nr - width + 1
  print(nreg)
  base = 0:(width - 1)
  sumrys <- lapply(1:nreg,
                   function(st) {
                     predict(plsr.moving[[st]], newdata = as.matrix(X[,st]), ncomp)
                   })
  
  sumrys
}



DWindowPLSR <- function(X,Y,ncomp,width1, width2,ystart, ...) {
  nr <- ncol(as.matrix(X))
  width1 = as.integer(width1)[1]
  width2 = as.integer(width2)[1]
  if (width1 <= 0 || width2 <= 0 || width1 >= nr || width2 >= nr)
    stop(paste("width must be in the range 1,...,", nr, sep=""))
  nreg = nr - width1 + 1
  for(i in seq(from=width2, to = ncol(X), by = width2)){
    base1 = 0:(width2 - 1)
    base2 = 0:(width1-1)
    sumrys <- lapply(1:width2,
                     function(st) {
                       plsr(Y3a[,i+st-1]~as.matrix(X3[,base2+st+i-width2]), ncomp)
                     })
    win.start <- ifelse(is.even(win2)== TRUE, win2/2, round2(win2/2,0))
    for (j in win.start:width2){
      if (i<=(width2+11)){ print(win.start+i-1)
        print(base1+win.start+j-1)
      }
    }
  }
}


splitWithOverlap <- function(vec, seg.length, overlap) {
  starts = seq(1, length(vec), by=seg.length-overlap)
  ends   = starts + seg.length -1
  ends[ends > length(vec)] = length(vec)
  #lapply(1:length(starts), function(i) vec[starts[i]:ends[i]])
  lapply(1:(length(starts)-seg.length+overlap+1), function(i) vec[starts[i]:ends[i]])
  #print(starts)
  #print(ends)
  # print(length(starts)-seg.length+overlap+1)
  #print(length(starts))
}

getsummary <- function(x,y){
  x <- as.vector(x)
  y <- as.vector(y)
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]
  lm.model <- lm(y~x)
  rss <- c(crossprod(lm.model$residuals))
  mse <- rss / length(lm.model$residuals)
  RMSE <- round(sqrt(mse),2)
  rmse.obs <- round(sqrt(mean((x-y)^2)),2)
  rmse1 <- sqrt(mse)
  std <- sd(x)
  rpd <- round(std/rmse.obs,2)
  bias1 <- round((sum(y, na.rm=TRUE)- sum(x, na.rm=TRUE))/length(x),2)
  bias2 <- round(sum(y-x, na.rm=TRUE)/length(x),2)
  r2 <- round(summary(lm.model)$r.squared, 2)
  coefs <- coef(lm.model)
  b0 <- round(coefs[1], 2)
  b1 <- round(coefs[2],2)
  cat("N=", length(x), ", bias1 =", bias1, ", bias2 = ", bias2,", slope =", b1,", offset =", b0, ", R2 = ", r2, ", RMSE = ", RMSE, ", rmse1 = ", rmse.obs, ", RPD", rpd,"\n")
  
}
