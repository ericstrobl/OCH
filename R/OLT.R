OLT <- function(tR,xR,yR,tO,xO,yO,xT){
  
  xR = as.matrix(xR)
  xT = as.matrix(xT)
  yR = as.matrix(yR)
  xO = as.matrix(xO)
  yO = as.matrix(yO)
  
  dO = ncol(xO)
  nO = nrow(xO)
  nR = nrow(xR)
  nT = nrow(xT)
  
  xn = rbind(xO,xR,xT) #####
  xn = normalizeAB2(xn,0.01,1)
  xOn = xn[1:nO,,drop=FALSE]
  xRn = xn[(nO+1):(nO+nR),,drop=FALSE]
  xTn = xn[(nO+nR+1):(nO+nR+nT),,drop=FALSE]
  
  kernel = INK_spline2n(xOn,rbind(xOn,xRn,xTn))
  
  ridges = c(1E-1,1E-2,1E-3,1E-4,1E-5,1E-6,1E-7,1E-8); ##
  
  err = Inf
  
  gRf = matrix(0,nR,2)
  gTf = matrix(0,nT,2)
  for (t in 0:1){
    iX = which(tO==t)
    nX = length(iX)
    err = Inf
    for (r in 1:length(ridges)){
      iK = spdinv(kernel[iX,iX] + nX*ridges[r]*diag(nX))
      alpha = t(yO[iX]) %*% iK
      gO = t(alpha %*% kernel[iX,iX])
      
      pL = LOO_KRR_fast(yO[iX], yO[iX] - gO, iK %*% kernel[iX,iX])
      
      errn = mean( (yO[iX] - pL)^2 )
      
      if (errn < err){
        err = errn
        gRf[,t+1] = t(alpha %*% kernel[iX,(nO+1):(nO+nR)])
        gTf[,t+1] = t(alpha %*% kernel[iX,(nO+nR+1):(nO+nR+nT)])
      }
    }
  }
  gRf = gRf[,2]-gRf[,1]
  gTf = gTf[,2]-gTf[,1]
  
  

  if (sd(gRf)<1E-20){
    return(gTf)
  } else{
    mod = lm.fit(cbind(gRf,1),yR)
    return(cbind(gTf,1) %*% mod$coefficients)
  }
  
}
