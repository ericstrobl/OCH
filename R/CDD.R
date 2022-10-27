CDD <- function(tO,xO,yO,mO,
                xT){
  
  #iI = included individuals
  #xT = X test
  #mO = time step
  
  xO = as.matrix(xO)
  yO = as.matrix(yO)
  xT = as.matrix(xT)
  
  dO = ncol(xO)
  nO = nrow(xO)
  nT = nrow(xT)
  
  # create RBF kernels
  xn = rbind(xO,xT) #####
  xn = normalizeAB2(xn,0,1)
  xOn = xn[1:nO,,drop=FALSE]
  xTn = xn[(nO+1):(nO+nT),,drop=FALSE]

  ridges = c(1E-1,1E-2,1E-3,1E-4,1E-5,1E-6,1E-7,1E-8); ##
  
  ## CATT
  
  mm = length(unique(mO))
  gTf = matrix(1,nT,2*mm)
  cnt=0
  for(m in 0:(mm-1)){
    for (t in 0:1){
      cnt = cnt + 1
      iX = which(mO==m & tO==t)
      nX = length(iX)
      
      kernel = INK_spline2n(xOn[iX,,drop=FALSE],rbind(xOn[iX,,drop=FALSE],xTn))
      
      err = Inf
      for (r in 1:length(ridges)){
        
        iK = spdinv(kernel[1:nX,1:nX] + nX*ridges[r]*diag(nX)) 
        alpha = t(yO[iX]) %*% iK
        gO = t(alpha %*% kernel[1:nX,1:nX])
        
        pL = LOO_KRR_fast(yO[iX], yO[iX] - gO, iK %*% kernel[1:nX,1:nX])
        
        errn = mean( (yO[iX] - pL)^2 )
        
        if (errn < err){
          err = errn
          gTf[,cnt] = t(alpha %*% kernel[1:nX,(nX+1):(nX+nT)]) #final prediction model
        }
      }
    }
  }
  
  pF = (gTf[,4]-gTf[,3]) - (gTf[,2]-gTf[,1])
  
  
  return(pF)
  
}
