get_RCT_OBS_data <- function(nR,nO,d,prop){
  
  
  xR = matrix(runif(nR*d)*2-1,nR,d)
  xR[,1] = runifI2(nR,-1+prop*2,1) ###
  tR = c(rep(0,nR/2),rep(1,nR/2))
  
  xO = matrix(runif(nO*d)*2-1,nO,d)
  tO = c(rep(0,nO/2),rep(1,nO/2))
  mO = c(rep(0,nO/4),rep(1,nO/4),rep(0,nO/4),rep(1,nO/4))
  
  ff = sample(1:4,4,replace=TRUE)
  yOt = matrix(0,nO,4)
  yRt = matrix(0,nR,4)
  
  xRs = rowSums(xR)
  xOs = rowSums(xO)
  for (f in 1:length(ff)){
    if (ff[f]==1){
      yRt[,f] = xRs
      yOt[,f] = xOs
    } else if (ff[f]==2){
      yRt[,f] = xRs*pnorm(xRs,0,1)
      yOt[,f] = xOs*pnorm(xOs,0,1)
    } else if (ff[f]==3){
      yRt[,f] = exp(-xRs^2)
      yOt[,f] = exp(-xOs^2)
    } else if (ff[f]==4){
      yRt[,f] = tanh(xRs)
      yOt[,f] = tanh(xOs)
    }
  }
  
  yO = rnorm(nO)
  for (f in 1:length(ff)){
    if (f==1){
      iX = which(tO == 0  & mO == 0)
      yO[iX] = yOt[iX,1]
    } else if (f==2){
      iX = which(tO == 1  & mO == 0)
      yO[iX] = yOt[iX,2] ###
    } else if (f==3){
      iX = which(tO == 0 & mO == 1)
      yO[iX] = yOt[iX,3]
    } else if (f==4){
      iX = which(tO == 1 & mO == 1)
      yO[iX] = yOt[iX,4]
    }
    
  }
  
  
  ## mixing coefficients
  b1 = runifI2(1,0,1)
  b0 = runifI2(1,0,1)
  
  # b1 = 0.5
  # b0 = 0.5
  
  yO = yO + sqrt(0.1)*rnorm(nO)
  
  
  ##
  yR = rep(0,nR)
  
  iX1 = which(tR==1)
  iX0 = which(tR==0)
  
  i1 = rbinom(length(iX1),1,b1)
  i0 = rbinom(length(iX0),1,b0)
  
  yR[iX1[i1==1]] = yRt[iX1[i1==1],4] + sqrt(0.1)*rnorm(length(iX1[i1==1])) # T=1, effective
  ii = rbinom(length(iX1[i1==0]),1,0.5)
  yR[(iX1[i1==0])[ii==1]] = yRt[(iX1[i1==0])[ii==1],2] + sqrt(0.1)*rnorm(length((iX1[i1==0])[ii==1]))# T=1, ineffective
  yR[(iX1[i1==0])[ii==0]] = yRt[(iX1[i1==0])[ii==0],1] + sqrt(0.1)*rnorm(length((iX1[i1==0])[ii==0]))# T=1, ineffective
  
  
  yR[iX0[i0==1]] = yRt[iX0[i0==1],3] + sqrt(0.1)*rnorm(length(iX0[i0==1])) # T=0, effective
  ii = rbinom(length(iX0[i0==0]),1,0.5)
  yR[(iX0[i0==0])[ii==1]] = yRt[(iX0[i0==0])[ii==1],2] + sqrt(0.1)*rnorm(length((iX0[i0==0])[ii==1]))# T=1, ineffective
  yR[(iX0[i0==0])[ii==0]] = yRt[(iX0[i0==0])[ii==0],1] + sqrt(0.1)*rnorm(length((iX0[i0==0])[ii==0]))# T=1, ineffective
  
  
  yOR = yOt[,4]*b1 - yOt[,3]*b2 + (yOt[,2]*0.5 + yOt[,1]*0.5)*(b2-b1)
  
  
  return(list(xR=xR,yR=yR,tR=tR,xO=xO,yO=yO,tO=tO,mO=mO))
}
