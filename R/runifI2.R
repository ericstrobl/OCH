runifI2 <- function(nsamps,a,b){
  
  ss = (runif(nsamps)*(b-a) + a)
  # ss=1
  return(ss)
  
}