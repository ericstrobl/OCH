library(proxy)
library(Rfast)
library(quadprog)

props = c(0, 0.25, 0.5, 0.75, 0.9, 0.95) # percent excluded
ds = c(1,3,6,10)
nR = 100
nO = 1000

MSE = array(0,c(10,1000,length(ds),length(props)))
Time = array(0,c(10,1000,length(ds),length(props)))
for (t in 1:500){
  print(t)
  for (d in 1:length(ds)){
    for (p in 1:length(props)){
      xR = matrix(runif(nR*ds[d])*2-1,nR,ds[d])
      xR[,1] = runifI2(nR,-1+props[p]*2,1) ###
      tR = c(rep(0,nR/2),rep(1,nR/2))
      
      xO = matrix(runif(nO*ds[d])*2-1,nO,ds[d])
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
      
      delta1 = yOt[,4] - (yOt[,2]*0.5 + yOt[,1]*0.5) 
      delta0 = yOt[,3] - (yOt[,2]*0.5 + yOt[,1]*0.5) 
      
      iX1 = which(tR==1)
      iX0 = which(tR==0)
      
      i1 = rbinom(length(iX1),1,b1)
      i0 = rbinom(length(iX0),1,b0)
      
      pick = sample(c(0,1),1)
      if (pick == 0){
        yR[iX1[i1==1]] = (yRt[iX1[i1==1],4] + delta1[iX1[i1==1]]) + sqrt(0.1)*rnorm(length(iX1[i1==1])) # T=1, M=1 plus
        yR[iX1[i1==0]] = yRt[iX1[i1==0],4] + sqrt(0.1)*rnorm(length(iX1[i1==0])) # T=1, M=1
        
        yR[iX0[i0==1]] = (yRt[iX0[i0==1],3] + delta0[iX0[i0==1]]) + sqrt(0.1)*rnorm(length(iX0[i0==1])) # T=0, M=1 plus
        yR[iX0[i0==0]] = yRt[iX0[i0==0],3] + sqrt(0.1)*rnorm(length(iX0[i0==0])) # T=0, M=1
        
        yOR = ((yOt[,4] + delta1)*b1 + yOt[,4]*(1-b1)) - ((yOt[,3] + delta0)*b0 + yOt[,3]*(1-b0))
      } else{
        ii = rbinom(length(iX1[i1==1]),1,0.5)
        yR[(iX1[i1==1])[ii==1]] = yRt[(iX1[i1==1])[ii==1],2] + sqrt(0.1)*rnorm(length((iX1[i1==1])[ii==1]))# T=1, M=1, mix = 1
        yR[(iX1[i1==1])[ii==0]] = yRt[(iX1[i1==1])[ii==0],1] + sqrt(0.1)*rnorm(length((iX1[i1==1])[ii==0]))# T=1, M=1, mix = 1
        ii = rbinom(length(iX1[i1==0]),1,0.5)
        yR[(iX1[i1==0])[ii==1]] = yRt[(iX1[i1==0])[ii==1],2] - delta1[(iX1[i1==0])[ii==1]]/2 + sqrt(0.1)*rnorm(length((iX1[i1==0])[ii==1]))# T=1, M=1 minus, mix = 0
        yR[(iX1[i1==0])[ii==0]] = yRt[(iX1[i1==0])[ii==0],1] - delta1[(iX1[i1==0])[ii==0]]/2 + sqrt(0.1)*rnorm(length((iX1[i1==0])[ii==0]))# T=1, M=1 minus, mix = 0
        
        ii = rbinom(length(iX0[i0==1]),1,0.5)
        yR[(iX0[i0==1])[ii==1]] = yRt[(iX0[i0==1])[ii==1],2] + sqrt(0.1)*rnorm(length((iX0[i0==1])[ii==1]))# T=0, M=1, mix = 1
        yR[(iX0[i0==1])[ii==0]] = yRt[(iX0[i0==1])[ii==0],1] + sqrt(0.1)*rnorm(length((iX0[i0==1])[ii==0]))# T=0, M=1, mix = 1
        ii = rbinom(length(iX0[i0==0]),1,0.5)
        yR[(iX0[i0==0])[ii==1]] = yRt[(iX0[i0==0])[ii==1],2] - delta0[(iX0[i0==0])[ii==1]]/2 + sqrt(0.1)*rnorm(length((iX0[i0==0])[ii==1]))# T=1, M=1 minus, mix = 0
        yR[(iX0[i0==0])[ii==0]] = yRt[(iX0[i0==0])[ii==0],1] - delta0[(iX0[i0==0])[ii==0]]/2 + sqrt(0.1)*rnorm(length((iX0[i0==0])[ii==0]))# T=1, M=1 minus, mix = 0
        
        yOR = ((yOt[,2]*0.5+yOt[,1]*0.5 - delta1)*b1 + (yOt[,2]*0.5+yOt[,1]*0.5)*(1-b1)) - ((yOt[,2]*0.5+yOt[,1]*0.5 - delta0)*b0 + (yOt[,2]*0.5+yOt[,1]*0.5)*(1-b0))
      }
      
      ###
      ptm = proc.time()
      
      MSE[1,t,d,p] = mean( (yOR - OCH12(tR,xR,yR,tO,xO,yO,mO,xT=xO) )^2)
      Time[1,t,d,p] = (proc.time() - ptm)[3]
      
      ptm = proc.time()
      MSE[2,t,d,p] = mean( (yOR - CDD(tO,xO,yO,mO,xO))^2)
      Time[2,t,d,p] = (proc.time() - ptm)[3]
      
      m1 = which(mO==1)
      
      ptm = proc.time()
      MSE[3,t,d,p] = mean( (yOR - TwoStep(tR,xR,yR,tO[m1],xO[m1,],yO[m1],xO))^2)
      Time[3,t,d,p] = (proc.time() - ptm)[3]
      
      ptm = proc.time()
      MSE[4,t,d,p] = mean( (yOR - OLT(tR,xR,yR,tO[m1],xO[m1,],yO[m1],xO))^2)
      Time[4,t,d,p] = (proc.time() - ptm)[3]
      
      ptm = proc.time()
      MSE[5,t,d,p] = mean( (yOR - OBS1(tO[m1],xO[m1,],yO[m1],xO))^2)
      Time[5,t,d,p] = (proc.time() - ptm)[3]
      
      ptm = proc.time()
      MSE[6,t,d,p] = mean( (yOR - OBS1(tR,xR,yR,xO))^2)
      Time[6,t,d,p] = (proc.time() - ptm)[3]
      
      ptm = proc.time()
      MSE[7,t,d,p] = mean( (yOR - SDDsa2(tR,xR,yR,
                                         tO,xO,yO,mO,
                                         xT=xO) )^2)
      Time[7,t,d,p] = (proc.time() - ptm)[3]
      
      
      mOt = mO[m1]
      ptm = proc.time()
      MSE[8,t,d,p] = mean( (yOR - OCH12(tR,xR,yR,
                                      tO[m1],xO[m1,],yO[m1],mOt,
                                      xT=xO))^2)
      Time[8,t,d,p] = (proc.time() - ptm)[3]
      
      
      ptm = proc.time()
      MSE[9,t,d,p] = mean( (yOR - OCH12_LM(tR,xR,yR,
                                         tO,xO,yO,mO,
                                         xT=xO))^2)
      Time[9,t,d,p] = (proc.time() - ptm)[3]
      
      ptm = proc.time()
      MSE[10,t,d,p] = mean( (yOR - OCH12_LM(tR,xR,yR,
                                          tO[m1],xO[m1,],yO[m1],mOt,
                                          xT=xO))^2)
      Time[10,t,d,p] = (proc.time() - ptm)[3]
    }
  }
  save(MSE,Time, file="ExpSynth_Viol_mix.RData")
}

for (pp in 1:6){
  print(median(MSE[6,1:500,,pp]))
  # print(skewness(MSE[6,1:500,,pp]))
  # print(median.test(MSE[1,1:500,,1],MSE[8,1:500,,1]))
}

library(e1071)  
for (pp in 1:4){
  # print(median(MSE[1,1:500,pp,]))
  # print(skewness(MSE[10,1:500,,pp]))
  print(median.test(MSE[1,1:500,pp,],MSE[8,1:500,pp,]))
}
