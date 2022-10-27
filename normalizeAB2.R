normalizeAB2 <- function(x,a,b){
  
  x = as.matrix(x)
  for (i in 1:ncol(x)){
    x[,i] = (b-a)*(x[,i]-min(x[,i]))/(max(x[,i]) - min(x[,i]))  + a
  }
  
  x
}