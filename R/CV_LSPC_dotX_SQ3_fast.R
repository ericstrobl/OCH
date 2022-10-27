CV_LSPC_dotX_SQ3_fast <- function(dotX,y_trm,Hs,cym,nn_cy){
  
  ## estimate conditional expectations
  n = nrow(y_trm); m = nrow(dotX)
  err = Inf; #maximum error allowed
  for (l in seq_len(length(Hs))){
      
      KRRs = t(t(y_trm) %*% Hs[[l]])
      res = LOO_LSPC(y_trm, y_trm-KRRs[1:n,],Hs[[l]])$pre
      dens = rbind(res,KRRs[(n+1):m,])
      
      dens = pmax(dens,0) + 1E-10;
      dens = dens/rowSums(dens)
      term1 = 0.5*mean(rowSums(dens[1:n,]^2))
      
      term2 = mean(dens[cbind(1:n,nn_cy)])
      errn = term1 - term2
      
      if (errn<err){
        err = errn
        densf = dens
      }
  }
  
  return(list(densf=densf,err=err))
  
}