LSPC_DCH <- function(K_tr,y_tr,y_trm,K_te,cy){
  
  n = nrow(K_tr)
  
  K = cbind(K_tr,K_te);
  m = ncol(K)
  
  nn_cy = get.knnx(cy,y_tr,k=1)$nn.index #for fast computation of term 2
  
  cym = matrix(cy); cym = repmat(t(cy),n,1);
  
  # invert kernels
  lambdas = c(1E-1,1E-2,1E-3,1E-4,1E-5,1E-6);
  Hs = vector("list", length = length(lambdas))
  for (l in seq_len(length(lambdas))){
    K_tr1 = K[1:n,1:n]
    
    Hs[[l]] = spdinv( K_tr + m*diag(nrow(K_tr))*lambdas[l]) %*%
      K[1:n,1:m];
  }
  
  # perform LSPC (Step 1)
  out=CV_LSPC_dotX_SQ3_fast(t(K),y_trm,Hs,cym,nn_cy)
  
  # sharpen estimate (Step 3)
  densf = sharpen_density(out$densf,cym[1,2]-cym[1,1],nn_cy,discrete=TRUE)
  dens = densf[(n+1):m,]
  
  return(list(y = cy, dens = dens))
  
}