DMD_onedens <- function(tR,xR,yR,
                xT,lb=-Inf,ub=Inf,discrete=FALSE){
  
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
  
  y_all = yR
  my = mean(y_all);
  sy = sd(y_all);
  y_all = (y_all - my)/sy
  lb = (lb - my)/sy
  ub = (ub - my)/sy
  if (discrete == FALSE){
    cy = seq(min(y_all),max(y_all),length.out=500);
  } else{
    cy = sort(unique(y_all))
  }
  
  ### compute indices of RCT, observational and test data
  
  iR = 1:nR
  iR0 = which(tR==0)
  iR1 = which(tR==1)
  
  iT = (nR+1):(nR+nT)
  
  ### estimate p(Y_1(1)|X)
  
  p1 = DDR_DCH(K[iR1,iR1],y_all[iR1],K[iR1,iT],cy,lb,ub,discrete)$dens
  
  ### estimate p(Y_1(0)|X)
  
  p0 = DDR_DCH(K[iR0,iR0],y_all[iR0],K[iR0,iT],cy,lb,ub,discrete)$dens
  
  cy = sy*cy + my;
  if (discrete){
    dens1 = p1
    dens0 = p0
  } else{
    dens1 = p1 / sy
    dens0 = p0 / sy
  }
  
  return (list(y = cy, dens1 = dens1, dens0 = dens0))
}
