OBS1 <- function(tO,xO,yO,xT){
  
  xT = as.matrix(xT)
  xO = as.matrix(xO)
  yO = as.matrix(yO)
  
  dO = ncol(xO)
  nO = nrow(xO)
  nT = nrow(xT)
  
  # create RBF kernels
  xn = rbind(xO,xT) #####
  xn = normalizeAB2(xn,0.01,1)
  xOn = xn[1:nO,,drop=FALSE]
  xTn = xn[(nO+1):(nO+nT),,drop=FALSE]
  
  kernel = INK_spline2n(xOn,rbind(xOn,xTn))
  
  ridges = c(1E-1,1E-2,1E-3,1E-4,1E-5,1E-6,1E-7,1E-8); ##
  
  err = Inf
  
  gTf = matrix(0,nT,2)
  # print(xn[501,])
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
        gTf[,t+1] = t(alpha %*% kernel[iX,(nO+1):(nO+nT)])
      }
    }
  }
  gTf = gTf[,2]-gTf[,1]
  
  return(gTf)
  
}
