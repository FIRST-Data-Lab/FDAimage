FDAimage.est.pc <- function(B,X,Y){
  n=nrow(Y)
  nx=ncol(X)
  npix=ncol(Y)
  B=as.matrix(B)
  J=ncol(B)

  WW=kronecker(crossprod(X),crossprod(B))
  rhs=as.vector(crossprod(B,t(Y))%*%X)
  theta=solve(WW,rhs)
  theta.mtx=matrix(theta,J,nx)
  beta=B%*%theta.mtx
  Yhat=tcrossprod(X,beta)

  ssei=apply((Y-Yhat)^2,1,sum)
  sse=sum(ssei)
  lambdac=NA;

  list(beta=beta,gamma=theta.mtx,sse=sse)
}
