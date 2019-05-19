FDAimage.CC.pc <- function(B,X,Y,tria,alpha0=0.05,nboot=100,ksi=0.01){
  n <- nrow(Y)
  nx <- ncol(X)
  npix <- ncol(Y)
  B <- as.matrix(B)
  J <- ncol(B)
  ######################################################################
  # Step 0: Preparation
  # Step 0.1: For Estimation
  WW <- kronecker(crossprod(X),crossprod(B))
  WW.inv <- chol2inv(chol(WW))

  # Step 0.2: For Band
  VV <- crossprod(B)
  VV.inv <- chol2inv(chol(VV))
  P.band <- B%*%VV.inv%*%t(B)
  tmp <- chol2inv(chol(crossprod(X)/n))
  XX.inv <- matrix(diag(tmp),ncol=1)

  # Step 0.3: For Alpha Adjusting
  a.grid <- c(10^(seq(-10,-4,by=2)),(2:9)*1e-4,seq(0.001,alpha0,0.004))
  Zalpha <- matrix(qnorm((1-a.grid/2)),nrow=1)
  Zalpha <- split(Zalpha,rep(1:ncol(Zalpha),each=nrow(Zalpha)))

  ######################################################################
  # Step 1: Estimation
  rhs <- as.vector(crossprod(B,t(Y))%*%X)
  theta <- as.matrix(WW.inv%*%rhs)
  theta.mtx <- matrix(theta,J,nx)
  beta <- as.matrix(B%*%theta.mtx)
  Yhat <- tcrossprod(X,beta)
  R <- Y-Yhat

  ######################################################################
  # Step 2: Get Eta.hat
  eta <- t(tcrossprod(P.band,R))
  eps <- R-eta
  sig2 <- apply(eps^2,2,mean)
  Geta <- apply((eta^2),2,mean)+sig2/npix/tria
  Sigma <- t(kronecker(XX.inv,t(Geta)))/n

  ######################################################################
  # Step 3: Bootstrap
  cover_all <- list()
  for(ib in 1:nboot){
    delta1 <- sample(c(-1,1),n,replace=TRUE)
    delta2 <- matrix(sample(c(-1,1),n*npix,replace=TRUE),nrow=n,ncol=npix)
    Yb <- Yhat+diag(delta1)%*%eta+delta2*eps

    rhsb <- as.vector(crossprod(B,t(Yb))%*%X)
    thetab <- as.matrix(WW.inv%*%rhsb)
    thetab.mtx <- matrix(thetab,J,nx)
    betab <- as.matrix(B%*%thetab.mtx)
    Ybhat <- tcrossprod(X,betab)
    Rb <- Yb-Ybhat

    etab <- t(tcrossprod(P.band,Rb))
    epsb <- Rb-etab
    sigb2 <- apply(epsb^2,2,mean)
    Getab <- apply(etab^2,2,mean)+sigb2/npix/tria
    Sigmab <- t(kronecker(XX.inv,t(Getab)))/n

    coverb <- lapply(Zalpha,cover_boot,beta.hat=beta,
                  betab.hat=betab,Sigmab=Sigmab)
    if(ib==1){
      cover_all <- coverb
    }else{
      cover_all <- mapply("+",cover_all,coverb,SIMPLIFY=FALSE)
    }
  }

  ######################################################################
  # Step 4: Sdjust Alpha Level
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
