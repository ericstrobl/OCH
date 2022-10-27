library(proxy)
library(Rfast)
library(quadprog)
library(MASS)
library(DDR)
library(matrixcalc)

props = c(0, 0.25, 0.5, 0.75, 0.9, 0.95) # percent excluded
ds = c(1,3,6,10)
nR = 500
nO = 1000

MISE = array(0,c(10,1000,length(ds),length(props)))
Time = array(0,c(10,1000,length(ds),length(props)))
for (t in 89:1000){ #907
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
      
      ##
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
      

      
      ####
      numCores <- detectCores()-1; registerDoParallel(numCores) # set up parallel computing
      ptm = proc.time()
      out = DMD(tR,xR,yR,tO,xO,yO,mO,xT=xO)
      MISE[1,t,d,p] = compute_MISE(out$ddens,out$y[2]-out$y[1],out$y,yOt,b1,b0)
      Time[1,t,d,p] = (proc.time() - ptm)[3]
      stopImplicitCluster()
      # # 
      # # ptm = proc.time()
      # # MSE[2,t,d,p] = mean( (yOR - CDD(tO,xO,yO,mO,xO))^2)
      # # Time[2,t,d,p] = (proc.time() - ptm)[3]
      # 
      # m1 = which(mO==1)
      # 
      # numCores <- detectCores()-1; registerDoParallel(numCores)
      # ptm = proc.time()
      # out = DDR_ddens1(tO[m1],xO[m1,],yO[m1],xO)
      # MISE[2,t,d,p] = compute_MISE(out$ddens,out$y[2]-out$y[1],out$y,yOt,b1,b0)
      # Time[2,t,d,p] = (proc.time() - ptm)[3]
      # stopImplicitCluster()
      # 
      numCores <- detectCores()-1; registerDoParallel(numCores)
      ptm = proc.time()
      # out = DDR_ddens1(tR,xR,yR,xO)
      out1 = DDR(xR[tR==1,,drop=FALSE],yR[tR==1],xO)
      out0 = DDR(xR[tR==0,,drop=FALSE],yR[tR==0],xO)
      out = DDR_ddens(tR[tR==1],xR[tR==1],yR[tR==1],xO)
      
      i=1
      cy = out$y
      true_dens1 = b1*dnorm(cy,yOt[i,4],sqrt(0.1)) + (1-b1)*(0.5*dnorm(cy,yOt[i,2],sqrt(0.1))+0.5*dnorm(cy,yOt[i,1],sqrt(0.1)))
      true_dens0 = b0*dnorm(cy,yOt[i,3],sqrt(0.1)) + (1-b0)*(0.5*dnorm(cy,yOt[i,2],sqrt(0.1))+0.5*dnorm(cy,yOt[i,1],sqrt(0.1)))
      plot(cy,true_dens1)
      # lines(cy,out1$dens[1,],col="green")
      lines(cy,true_dens0)
      lines(cy,out$dens1[1,],col="yellow")
      lines(cy,out$dens0[1,],col="purple")
      
      # MISE[3,t,d,p] = compute_MISE(out$ddens,out$y[2]-out$y[1],out$y,yOt,b1,b0)
      # Time[3,t,d,p] = (proc.time() - ptm)[3]
      stopImplicitCluster()
      
      numCores <- detectCores()-1; registerDoParallel(numCores)
      ptm = proc.time()
      out = DCH_ddens2_LM(tR,xR,yR,tO,xO,yO,mO,xT=xO)
      MISE[4,t,d,p] = compute_MISE(out$ddens,out$y[2]-out$y[1],out$y,yOt,b1,b0)
      Time[4,t,d,p] = (proc.time() - ptm)[3]
      stopImplicitCluster()
    }
  }
  save(MSE,Time, file="ExpSynth_DCH_ddens2.RData")
}


for (pp in 1:6){
  print(median(MISE[2,1:100,,pp]))
  # print(skewness(MISE[4,1:100,,pp]))
  # print(median.test(MISE[1,1:110,,pp],MISE[4,1:110,,pp]))
}

ms = c()
for (t in 1:4){
  bt = sample(1:500,500,replace=TRUE)
  # ms = c(ms, median(MSE[8,bt,1:4,1]))
  
}