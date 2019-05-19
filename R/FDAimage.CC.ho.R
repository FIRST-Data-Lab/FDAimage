FDAimage.CC.ho <- function(B,Q2,K,X,Y,lambda,alpha0=0.05,nboot=100,ksi=0.01){
  n <- nrow(Y)
  nx <- ncol(X)
  npix <- ncol(Y)
  BQ2 <- B%*%Q2
  BQ2 <- as.matrix(BQ2)
  K <- as.matrix(K)
  J <- ncol(BQ2)

  ######################################################################
  # Step 0: Preparation
  # Step 0.1: For Estimation
  WW <- kronecker(crossprod(X),crossprod(BQ2))
  P <- as.matrix(crossprod(Q2,K)%*%Q2)

  # Step 0.2: For Band
  VV <- crossprod(BQ2)
  VV.inv <- chol2inv(chol(VV))
  P.band <- BQ2%*%VV.inv%*%t(BQ2)
  tmp <- chol2inv(chol(crossprod(X)/n))
  XX.inv <- matrix(diag(tmp),ncol=1)

  # Step 0.3: For Alpha Adjusting
  a.grid <- c(10^(seq(-10,-4,by=2)),(2:9)*1e-4,seq(0.001,alpha0,0.004))
  Zalpha <- matrix(qnorm((1-a.grid/2)),nrow=1)
  Zalpha <- split(Zalpha,rep(1:ncol(Zalpha),each=nrow(Zalpha)))

  ######################################################################
  # Step 1: Estimation
  rhs <- as.vector(crossprod(BQ2,t(Y))%*%X)
  flag <- (rankMatrix(WW)<(J*nx))
  if(!flag){
    Ainv <- chol(WW,pivot=TRUE)
  }
  if(flag){
    warning("The solution may not be correct.")
    VV <- WW+diag(rep(1e-12),nx*J)
    Ainv <- chol(VV,pivot=TRUE)
  }
  A <- solve(t(Ainv))
  D <- as.matrix(kronecker(diag(rep(1,nx)),P))
  ADA <- A%*%D%*%t(A)
  eigs <- eigen(ADA)
  C <- eigs$values

  lambda <- as.matrix(lambda)
  nlam <- nrow(lambda)
  sse_all <- c(); df_all <- c(); gcv_all <- c();
  for(il in 1:nlam){
    Lam <- diag(rep(lambda[il],nx))
    Dlam <- as.matrix(kronecker(Lam,P))
    lhs <- WW+Dlam
    theta <- solve(lhs,rhs)
    theta.mtx <- matrix(theta,J,nx)
    beta <- BQ2%*%theta.mtx
    Yhat <- tcrossprod(X,beta)

    ssei <- apply((Y-Yhat)^2,1,sum)
    sse <- sum(ssei)
    sse_all <- c(sse_all,sse)

    df <- sum(1/(1+C*lambda[il]))
    df_all <- c(df_all,df)
    gcv <- n*npix*sse/(n*npix-df)^2
    gcv_all <- c(gcv_all,gcv)
  }
  j <- which.min(gcv_all)
  lamc <- lambda[j]
  df <- df_all[j]
  sse <- sse_all[j]
  gcv <- gcv_all[j]
  Lam <- diag(rep(lamc,nx))
  Dlam <- as.matrix(kronecker(Lam,P))
  lhs <- WW+Dlam
  WD.inv <- chol2inv(chol(lhs)) # for future use
  theta <- WD.inv%*%rhs
  theta.mtx <- matrix(theta,J,nx)
  gamma <- as.matrix(Q2%*%theta.mtx)
  beta <- as.matrix(BQ2%*%theta.mtx)
  Yhat <- tcrossprod(X,beta)
  R <- Y-Yhat

  ######################################################################
  # Step 2: Get Eta.hat
  eta <- t(tcrossprod(P.band,R))
  eps <- R-eta
  # sig2 <- apply(eps^2,2,mean)

  # Step 2.1: calculate Geta
  M <- lapply(1:n, function(iter) WD.inv%*%kronecker(X[iter,],VV)%*%VV.inv)

  eta.theta <- t(R%*%BQ2)
  A <- lapply(1:n, function(iter) tcrossprod(M[[iter]]%*%eta.theta))
  A <- Reduce("+",A)
  A <- as.matrix(A)

  Geta <- c()
  for(iter in 1:nx){
    a.begin <- (iter-1)*J+1
    a.end <- iter*J
    A.temp <- A[a.begin:a.end,a.begin:a.end]
    Geta <- cbind(Geta,apply(BQ2,1,function(x) x%*%A.temp%*%x))
  }
  Sigma <- Geta/n

  ######################################################################
  # Step 3: Bootstrap
  cover_all <- list()
  for(ib in 1:nboot){
    # delta1 <- rnorm(n)
    # delta2 <- matrix(rnorm(n*npix),nrow=n,ncol=npix)
    delta1 <- sample(c(-1,1),n,replace=TRUE)
    delta2 <- matrix(sample(c(-1,1),n*npix,replace=TRUE),nrow=n,ncol=npix)
    Yb <- Yhat+diag(delta1)%*%eta+delta2*eps

    rhsb <- as.vector(crossprod(BQ2,t(Yb))%*%X)
    thetab <- as.matrix(WD.inv%*%rhsb)
    thetab.mtx <- matrix(thetab,J,nx)
    betab <- as.matrix(BQ2%*%thetab.mtx)
    Ybhat <- tcrossprod(X,betab)
    Rb <- Yb-Ybhat

    eta.theta <- t(Rb%*%BQ2)
    A <- lapply(1:n, function(iter) tcrossprod(M[[iter]]%*%eta.theta))
    A <- Reduce("+",A)
    A <- as.matrix(A)

    Geta <- c()
    for(iter in 1:nx){
      a.begin <- (iter-1)*J+1
      a.end <- iter*J
      A.temp <- A[a.begin:a.end,a.begin:a.end]
      Geta <- cbind(Geta,apply(BQ2,1,function(x) x%*%A.temp%*%x))
    }
    Sigmab <- Geta/n

    coverb <- lapply(Zalpha,cover_boot,beta.hat=beta,
                  betab.hat=betab,Sigmab=Sigmab)
    if(ib==1){
      cover_all <- coverb
    }else{
      cover_all <- mapply("+",cover_all,coverb,SIMPLIFY=FALSE)
    }
  }

  ######################################################################
  # Step 4: Adjust Alpha Level
  tau_all <- list()
  tau_hat <- c()
  for(ii in 1:nx){
    tmp1 <- lapply(1:length(a.grid),function(x) cover_all[[x]][,ii])
    tau_all[[ii]] <- do.call(cbind,tmp1)/nboot
    tmp2 <- tau_all[[ii]]-1+alpha0
    tmp2[tmp2<0] <- 1e6
    tmp3 <- apply(tmp2,1,function(x) max(which(x==min(x))))
    tau_hat <- cbind(tau_hat,tmp3)
  }
  a.adjust <- a.grid[apply(tau_hat,2,quantile,prob=ksi)]
  Za <- matrix(qnorm((1-a.adjust/2)),nrow=1)
  Za <- matrix(rep(Za,each=npix),nrow=npix)

  cc.l <- beta-sqrt(Sigma)*Za
  cc.u <- beta+sqrt(Sigma)*Za
  list(cc.l=cc.l,cc.u.=cc.u,beta.hat=beta,alpha.adj=a.adjust)
}
