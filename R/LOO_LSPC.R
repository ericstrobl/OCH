LOO_LSPC <- function(tr_y, tr_e, H){
  dH = diag(H)
  err = tr_e/(1-dH)
  pre = tr_y - err
  # # 
  prep = pmax(pre,0) + 1E-10
  pre = prep/rowSums(prep)

  err = tr_y - pre
  
  list( pre=pre, err=mean(err^2), resid=err )
}