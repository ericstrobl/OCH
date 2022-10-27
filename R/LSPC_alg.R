LSPC_alg <- function(tR,xR,yR,
                tO,xO,yO,mO,
                xT,cy){
  
  xR = as.matrix(xR)
  xO = as.matrix(xO)
  xT = as.matrix(xT)
  
  dO = ncol(xO)
  nO = nrow(xO)
  nR = nrow(xR)
  nT = nrow(xT)
  
  ### normalize predictors to [0,1] and compute INK spline kernel
  
  xn = rbind(xR,xO,xT)
  xn = normalizeAB2(xn,0,1)
  xRn = xn[1:nR,,drop=FALSE]
  xOn = xn[(nR+1):(nR+nO),,drop=FALSE]
  xTn = xn[(nO+nR+1):(nO+nR+nT),,drop=FALSE]
  
  K = INK_spline2n(xOn,xn)
  
  ### normalize outcome
  
  y_all = c(yR,yO)
  yOm = matrix(0,nO,length(cy))
  for (i in 1:length(yO)){
    yOm[i,which(yO[i]==cy)]=1
  }
  
  
  ### compute indices of RCT, observational and test data
  
  iR = 1:nR
  iR0 = which(tR==0)
  iR1 = which(tR==1)
  
  iO = (nR+1):(nR+nO)
  iT = (nR+nO+1):(nR+nO+nT)
  
  ### estimate p(Y_0|X)
  
  iX0 = which(mO==0)
  iX0a = iO[iX0]
  
  p0 = LSPC_DCH(K[iX0,iX0a],y_all[iX0a],yOm[iX0,],K[iX0,c(iR,iT)],cy)$dens
  
  ### estimate p(Y_1(0)|X,T=0)
  
  iX10 = which(mO==1 & tO==0)
  iX10a = iO[iX10]
  
  p10 = LSPC_DCH(K[iX10,iX10a],y_all[iX10a],yOm[iX10,],K[iX10,c(iR,iT)],cy)$dens
  
  ### estimate p(Y_1(1)|X,T=1)
  
  iX11 = which(mO==1 & tO==1)
  iX11a = iO[iX11]
  
  p11 = LSPC_DCH(K[iX11,iX11a],y_all[iX11a],yOm[iX11,],K[iX11,c(iR,iT)],cy)$dens
  
  # ### verify DDR
  # lines(sy*cy + my,p11[(nR+1),]/sy,col="green")
  # lines(sy*cy + my,p10[(nR+1),]/sy,col="red")
  # lines(sy*cy + my,p0[(nR+1),]/sy,col="orange")
  
  
  # # ### toy data to verify optimization
  # 
  # nn = length(c(iR,iT))
  # p11 = matrix(0,nn,length(cy))
  # p10 = p11; p0 = p11;
  # for (i in 1:nn){
  #   p11[i,] = b1*dnorm(cy,yOt[i,4],sqrt(0.1)) + (1-b1)*(0.5*dnorm(cy,yOt[i,2],sqrt(0.1))+0.5*dnorm(cy,yOt[i,1],sqrt(0.1)))
  #   p10[i,] = dnorm(cy,-1,0.1)
  #   p0[i,] = dnorm(cy,0,0.1)
  # }
  # y_all[iR1] = sample(c(rnorm(length(iR1),1,0.1),rnorm(length(iR1),0,0.1)),length(iR1))
  # y_all[iR0] = sample(c(rnorm(length(iR0),-1,0.1),rnorm(length(iR0),0,0.1)),length(iR0))
  # y_all[iR1] = rnorm(length(iR1),0,0.1)
  # y_all[iR0] = rnorm(length(iR0),0,0.1)
  
  
  #### estimate mu1
  

  p0_sq = mean(rowSums(p0[iR,]^2)) ## AUC for p(Y_0|X)^2 
  p11_sq = mean(rowSums(p11[iR,]^2))
  p11_p0 = mean(rowSums(p11[iR,]*p0[iR,]))
 
  Dmat1 = matrix(0,1,1)
  Dmat1[1,1] = p11_sq + p0_sq - 2*p11_p0
  
  nn_cy1 = get.knnx(cy,y_all[iR1],k=1)$nn.index
  
  dvec1 = mean(p11[cbind(iR1,nn_cy1)])-mean(p0[cbind(iR1,nn_cy1)]) + p0_sq - p11_p0
  
  Amat = matrix(0,1,2)
  Amat[1,1] = 1; Amat[1,2]=-1;
  bvec = c(0,-1)
  
  mu1 = as.matrix(solve.QP(Dmat1, dvec1, Amat, bvec)$solution)[1]
  
  dens1 = p11[(nR+1):(nR+nT),]*mu1 + p0[(nR+1):(nR+nT),]*(1-mu1) ## test data only
  
  #### estimate mu0
  
  p10_sq = mean(rowSums(p10[iR,]^2))
  p10_p0 = mean(rowSums(p10[iR,]*p0[iR,]))    
  
  Dmat0 = matrix(0,1,1)
  Dmat0[1,1] = p10_sq + p0_sq - 2*p10_p0
  
  nn_cy0 = get.knnx(cy,y_all[iR0],k=1)$nn.index
  
  dvec0 = mean(p10[cbind(iR0,nn_cy0)])-mean(p0[cbind(iR0,nn_cy0)])+ p0_sq - p10_p0
  
  mu0 = as.matrix(solve.QP(Dmat0, dvec0, Amat, bvec)$solution)[1]
  
  dens0 = p10[(nR+1):(nR+nT),]*mu0 + p0[(nR+1):(nR+nT),]*(1-mu0) ## test data only
  
  dens1 = dens1
  dens0 = dens0
  
  return (list(y = cy, dens1 = dens1, dens0 = dens0, mu1 = mu1, mu0 = mu0))
}
