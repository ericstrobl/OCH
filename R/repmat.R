repmat = function(X,m,n){
  ##R equivalent of repmat (matlab)
  X=matrix2(X);
  mx = dim(X)[1]
  nx = dim(X)[2]
  matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
}