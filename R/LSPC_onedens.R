LSPC_onedens <- function(tR,xR,yR,
                        xT,cy){
  
  xR = as.matrix(xR)
  xT = as.matrix(xT)
  
  nR = nrow(xR)
  nT = nrow(xT)
  
  ### normalize predictors to [0,1] and compute INK spline kernel
  
  xn = rbind(xR,xT)
  xn = normalizeAB2(xn,0,1)
  xRn = xn[1:nR,,drop=FALSE]
  xTn = xn[(nR+1):(nR+nT),,drop=FALSE]
  
  K = INK_spline2n(xRn,xn)
  
  ### normalize outcome
  
  yRm = matrix(0,nR,length(cy))
  for (i in 1:length(yR)){
    yRm[i,which(yR[i]==cy)]=1
  }

  
  ### compute indices of RCT, observational and test data
  
  iR = 1:nR
  iR0 = which(tR==0)
  iR1 = which(tR==1)
  
  iT = (nR+1):(nR+nT)
  
  ### estimate p(Y_1(1)|X)
  
  p1 = LSPC_DCH(K[iR1,iR1],yR[iR1],yRm[iR1,],K[iR1,iT],cy)$dens
  
  ### estimate p(Y_1(0)|X)
  
  p0 = LSPC_DCH(K[iR0,iR0],yR[iR0],yRm[iR0,],K[iR0,iT],cy)$dens
  
  
  return (list(y = cy, dens1 = p1, dens0 = p0))
}
