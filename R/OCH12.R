OCH12 <- function(tR,xR,yR,
                   tO,xO,yO,mO,
                   xT){
  
  #iI = included individuals
  #xT = X test
  #mO = time step
  
  xR = as.matrix(xR)
  xO = as.matrix(xO)
  xT = as.matrix(xT)
  
  dO = ncol(xO)
  nO = nrow(xO)
  nR = nrow(xR)
  nT = nrow(xT)
 
  xn = rbind(xR,xO,xT) #####
  xn = normalizeAB2(xn,0,1)
  xRn = xn[1:nR,,drop=FALSE]
  xOn = xn[(nR+1):(nR+nO),,drop=FALSE]
  xTn = xn[(nO+nR+1):(nO+nR+nT),,drop=FALSE]
  
  kernel = INK_spline2n(xRn,xRn)
  # kernel_eig = eigen(kernel)
  # iV = which(kernel_eig$values>1E-10)
  
  ridges = c(1E-1,1E-2,1E-3,1E-4,1E-5,1E-6,1E-7,1E-8); ##
  
  ## CATE with RCT data
  gRft = matrix(0,nR,2)
  for (t in 0:1){
    iX = which(tR==t)
    nX = length(iX)
    err = Inf
    for (r in 1:length(ridges)){
      iK = spdinv(kernel[iX,iX] + nX*ridges[r]*diag(nX))
      # iK = t(kernel_eig$vectors[iV,]) %*% diag(1/kernel_eig$values[iV]) %*% t(kernel_eig$vectors[,iV])
      alpha = t(yR[iX]) %*% iK
      gR = t(alpha %*% kernel[iX,iX])
      
      pL = LOO_KRR_fast(yR[iX], yR[iX] - gR, iK %*% kernel[iX,iX])
      
      errn = mean( (yR[iX] - pL)^2 )
      
      if (errn < err){
        err = errn
        gRft[,t+1] = t(alpha %*% kernel[iX,])
      }
    }
  }
  gRft = gRft[,2]-gRft[,1]
  
  # plot(xR,gRft)
  
  ## Observational components  
  gOf = matrix(1,nR,3)
  gTf = matrix(1,nT,3)
  cnt=0
  for(m in 0:1){
    
    if (m == 0){
      
      cnt = cnt + 1
      iX = which(mO==m)
      nX = length(iX)
      
      kernel = INK_spline2n(xOn[iX,,drop=FALSE],rbind(xOn[iX,,drop=FALSE],xRn,xTn))
      
      err = Inf
      for (r in 1:length(ridges)){
        
        iK = spdinv(kernel[1:nX,1:nX] + nX*ridges[r]*diag(nX))
        alpha = t(yO[iX]) %*% iK
        gO = t(alpha %*% kernel[1:nX,1:nX])
        
        pL = LOO_KRR_fast(yO[iX], yO[iX] - gO, iK %*% kernel[1:nX,1:nX])
        
        errn = mean( (yO[iX] - pL)^2 )
        
        if (errn < err){
          err = errn
          gOf[,cnt] = t(alpha %*% kernel[1:nX,(nX+1):(nX+nR)]) #final prediction model
          gTf[,cnt] = t(alpha %*% kernel[1:nX,(nX+nR+1):(nX+nR+nT)]) #final prediction model
        }
        
      }
      
      # plot(xR,gOf[,1])
      
    } else {
      for (t in 0:1){
        cnt = cnt + 1
        iX = which(mO==m & tO==t)
        nX = length(iX)
        
        kernel = INK_spline2n(xOn[iX,,drop=FALSE],rbind(xOn[iX,,drop=FALSE],xRn,xTn))
        
        err = Inf
        for (r in 1:length(ridges)){
          
          iK = spdinv(kernel[1:nX,1:nX] + nX*ridges[r]*diag(nX))
          alpha = t(yO[iX]) %*% iK
          gO = t(alpha %*% kernel[1:nX,1:nX])
          
          pL = LOO_KRR_fast(yO[iX], yO[iX] - gO, iK %*% kernel[1:nX,1:nX])
          
          errn = mean( (yO[iX] - pL)^2 )
          
          if (errn < err){
            err = errn
            gOf[,cnt] = t(alpha %*% kernel[1:nX,(nX+1):(nX+nR)]) #final prediction model
            gTf[,cnt] = t(alpha %*% kernel[1:nX,(nX+nR+1):(nX+nR+nT)]) #final prediction model
          }
          
        }
      }
    }
  }
  
  # plot(xR,gOf[,3])
  
  gOff = matrix(1,nR,2)
  gTff = matrix(1,nT,2)
  
  gOff[,1] = gOf[,3] - gOf[,1]
  gOff[,2] = -gOf[,2] + gOf[,1]

  gTff[,1] = gTf[,3] - gTf[,1]
  gTff[,2] = -gTf[,2] + gTf[,1]
  
  # gOff[,1] = gOf[,3]
  # gOff[,2] = -gOf[,2]
  # 
  # gTff[,1] = gTf[,3]
  # gTff[,2] = -gTf[,2]
  
  # plot(xR,gOff[,1]*0.3449493 + gOff[,2]*0.7701808)
  
  Ii = diag(2)

  Dmat = t(gOff) %*% gOff
  Amat = t(rbind(Ii,-Ii))
  bvec = c(0,0,-1,-1)
  dvec = t(gRft) %*% gOff
  
  beta = as.matrix(solve.QP(Dmat, dvec, Amat, bvec)$solution)
  
  # print(beta)
  
  # plot(xR,gOff[,1]*beta[1] + gOff[,2]*beta[2])
  
  pF= gTff %*% beta
  
  #final prediction model
  # print(gTf)
  # 
  # 
  # 
  # print(beta + betas)
  # points(rowSums(xT),pF,col="blue")
  
  
  return (pF)
  
}
