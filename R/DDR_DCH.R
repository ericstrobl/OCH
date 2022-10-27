DDR_DCH <- function(K_tr,y_tr,K_te,cy,lb=-Inf,ub=Inf,discrete=FALSE){

  n = nrow(K_tr)
  
  K = cbind(K_tr,K_te);
  m = ncol(K)
  
  nn_cy = get.knnx(cy,y_tr,k=1)$nn.index #for fast computation of term 2
  
  cym = matrix(cy); cym = repmat(t(cy),n,1);
  
  #compute RBF kernels
  lambdas = c(1E-1,1E-2,1E-3,1E-4,1E-5,1E-6);
  # lambdas = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
  Hs = lapply(1,function(x) vector("list", length = length(lambdas)));
  for (l in seq_len(length(lambdas))){
    K_tr1 = K[1:n,1:n]

    Hs[[1]][[l]] = spdinv( K_tr + m*diag(nrow(K_tr))*lambdas[l]) %*%
      K[1:n,1:m];
  }
  
  # perform DDR in parallel (Step 1)
  out=CV_KRR_dotX_SQ3_fast(t(K),y_tr,Hs,cym,nn_cy,lb,ub,discrete)
  h_star = out$h_star
  
  # sharpen estimate (Step 3)
  densf = sharpen_density(out$densf,cym[1,2]-cym[1,1],nn_cy,discrete)
  dens = densf[(n+1):m,]
  
  return(list(y = cy, dens = dens))
  
}
