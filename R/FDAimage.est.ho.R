FDAimage.est.ho <- function(B,Q2,K,X,Y,lambda){
  n=nrow(Y)
  nx=ncol(X)
  npix=ncol(Y)
  BQ2=B%*%Q2
  BQ2=as.matrix(BQ2)
  K=as.matrix(K)
  J=ncol(BQ2)

  WW=kronecker(crossprod(X),crossprod(BQ2))
  rhs=rowSums(sapply(1:n, function(iter) as.matrix(kronecker(X[iter,],crossprod(BQ2,Y[iter,])))))
  P=as.matrix(crossprod(Q2,K)%*%Q2)

  lambda=as.matrix(lambda)
  nlam=nrow(lambda)
  gamma_all=list(nlam); beta_all=list(nlam); sse_all=c();
  for(il in 1:nlam){
    Lam=diag(rep(lambda[il],nx))
    Dlam=as.matrix(kronecker(Lam,P))
    lhs=WW+Dlam
    theta=solve(lhs,rhs)
    theta.mtx=matrix(theta,J,nx)
    gamma=Q2%*%theta.mtx
    gamma_all[[il]]=gamma
    beta=BQ2%*%theta.mtx
    beta_all[[il]]=beta
    Yhat=tcrossprod(X,beta)
    ssei=apply((Y-Yhat)^2,1,sum)
    sse=sum(ssei)
    sse_all=c(sse_all,sse)
  }
  list(beta=beta_all,gamma=gamma_all,sse=sse)
}
