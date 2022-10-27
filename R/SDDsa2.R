SDDsa2 <- function(tR,xR,yR,
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
  # 
  # y = c(yR, yO)
  # ym = mean(y)
  # ys = sd(y)
  # y = (y - ym) / ys
  # yR = as.matrix(y[1:nR])
  # yO = as.matrix(y[(nR+1):(nR +nO)])
  # 
  # create RBF kernels
  # xn = normalize2(rbind(xR,xO,xT))
  xn = rbind(xR,xO,xT) #####
  xn = normalizeAB2(xn,0,1)
  xRn = xn[1:nR,,drop=FALSE]
  xOn = xn[(nR+1):(nR+nO),,drop=FALSE]
  xTn = xn[(nO+nR+1):(nO+nR+nT),,drop=FALSE]
  
  kernel = INK_spline2n(xRn,xRn)
  
  ridges = c(1E-1,1E-2,1E-3,1E-4,1E-5,1E-6,1E-7,1E-8); ##
  
  ## CATE
  gRft = matrix(0,nR,2)
  for (t in 0:1){
    iX = which(tR==t)
    nX = length(iX)
    err = Inf
    for (r in 1:length(ridges)){
      iK = spdinv(kernel[iX,iX] + nX*ridges[r]*diag(nX))
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
  
  ## CATT
  
  mm = length(unique(mO))
  gOf = matrix(1,nR,2*mm)
  gTf = matrix(1,nT,2*mm)
  cnt=0
  for(m in 0:(mm-1)){
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

  # xTs = rowSums(xT)
  # points(xTs,gTf[,1],col="red")
  # points(xTs,gTf[,2],col="red")
  # points(xTs,gTf[,3],col="red")
  # points(xTs,gTf[,4],col="red")
  
  #linear regression model
  cnt = 0
  for(m in 0:(mm-1)){
    for (t in 0:1){
      cnt = cnt + 1
      if (t==1 & m != (mm-1)){
        gOf[,cnt] = -gOf[,cnt]
        gTf[,cnt] = -gTf[,cnt]
      } else if (t==0 & m==(mm-1)){
        gOf[,cnt] = -gOf[,cnt]
        gTf[,cnt] = -gTf[,cnt]
      }
    }
  }
  
  gOt = gOf[,2*mm]
  gOf = gOf[,-(2*mm),drop=FALSE]
  
  gTt = gTf[,2*mm]
  gTf = gTf[,-(2*mm),drop=FALSE]
  
  # gOf[,2] = -gOf[,2]
  # gOf[,3] = -gOf[,3]
  # gTf[,2] = -gTf[,2]
  # gTf[,3] = -gTf[,3]
  # 
  pF = rep(0,nT)
  err = Inf
  Ii = diag(mm*2-1)*1
  
  # betas = lm.fit(gOf,gRft)$coefficients
  # print(betas)
  # pF = gTf %*% betas
  
  betas = matrix(1,mm*2-1,1)
  gRft2 = gRft - gOt - gOf %*% betas
  # gR = yR

  nF = 5
  fld = rep(1:5,length.out=nR)
  # ridges = c(1,1E-1,1E-2,1E-3,1E-4,1E-5,1E-6); ##
  for (r in 1:length(ridges)){
    for (f in 1:nF){
      iX = !(fld==f)
      nX = sum(iX)

      cOf = t(gOf[iX,]) %*% gOf[iX,]
      Dmat = cOf + nR*ridges[r]*Ii
      Amat = t(rbind(Ii,-Ii))
      bvec = c(rep(-1,mm*2-1),rep(-2,mm*2-1))
      dvec = t(gOf[iX,]) %*% gRft2[iX]

      beta = as.matrix(solve.QP(Dmat, dvec, Amat, bvec)$solution)

      gR[fld==f] = gOf[fld==f,] %*% beta
    }

    errn = mean( (gRft2 - gR)^2 )

    if (errn < err){
      err = errn
      rf = r
    }

  }

  cOf = t(gOf) %*% gOf
  Dmat = cOf + nR*ridges[rf]*Ii
  Amat = t(rbind(Ii,-Ii))
  bvec = c(rep(-1,mm*2-1),rep(-2,mm*2-1))
  dvec = t(gOf) %*% gRft2
  beta = as.matrix(solve.QP(Dmat, dvec, Amat, bvec)$solution)
  pF= gTf %*% (beta + betas) + gTt #final prediction model
  # print(gTf)
  # 
  # 
  # 
  # print(beta + betas)
  # points(rowSums(xT),pF,col="blue")
  
  
  return (pF)
  
}
