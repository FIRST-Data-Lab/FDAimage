#' Generating simulation examples
#'
#' This function is used to generate the simulation examples for image-on-scalar regression.
#'
#' @importFrom MASS mvrnorm
#' @param n The number of images.
#' \cr
#' @param Z The cooridinates of dimension \code{npix} by two. Each row is the coordinates of a pixel/voxel.
#' \cr
#' @param ind.inside The indices of the pixels/voxels which are inside the boundary.
#' \cr
#' @param sigma The level of white noise.
#' \cr
#' @param rho The correlation used to generate the explanatory variales -- default is 0.5.
#' \cr
#' @param iter The seed used to generate the simulation examples -- default is 2019.
#' \cr
#' @param lambda1&lambda2 The eigenvalues used to generate the simulation examples -- default are 0.03 and 0.006, respectively.
#' \cr
#' @return
#' \item{Y}{The response images.}
#' \item{X}{The covariates.}
#' \item{Z}{The coordinates of the pixels/voxels.}
#' \item{beta.true}{The true coefficient functions.}
#' \cr
#' @details This R package is the implementation program for manuscript entitled ``Multivariate Spline Estimation and Inference for Varying Coeffiient Models with Imaging Data" by Shan Yu, Guannan Wang, Li Wang and Lijian Yang.
#' @export
#'
data.FDAimage <- function(n,Z,ind.inside,sigma,rho=0.5,iter=2019,lambda1=0.03,
                        lambda2=0.006){
  xx <- Z[,1]
  n1 <- length(unique(xx))
  yy <- Z[,2]
  n2 <- length(unique(yy))
  npix <- n1*n2;

  uu <- seq(0,1,length.out=n1)
  vv <- seq(0,1,length.out=n2)

  # generate true function
  beta0 <- 5*((xx-0.5)^2+(yy-0.5)^2)
  beta1 <- -xx^3+yy^3
  beta2 <- exp(-15*((xx-0.5)^2+(yy-0.5)^2))
  beta.true <- cbind(beta0,beta1,beta2)

  ind.outside <- setdiff(1:(n1*n2),ind.inside)
  beta.true[ind.outside,] <- NA

  set.seed(iter)
  mu <- rep(0,2)
  Sigma <- diag(rep(1,2))
  for(i in 1:2){
    for(j in i:2){
      Sigma[i,j] <- rho^(j-i)
      Sigma[j,i] <- rho^(j-i)
    }
  }
  X <- mvrnorm(n,mu,Sigma)
  X <- apply(X,2,scale)
  X[X>3] <- 3
  X[X<(-3)] <- -3
  X <- cbind(1,X)

  eps1 <- rnorm(n,mean=0,sd=sqrt(lambda1))
  eps2 <- rnorm(n,mean=0,sd=sqrt(lambda2))
  phi1 <- 1.738*sin(2*pi*xx/n1)
  phi2 <- 1.939*cos(2*pi*yy/n2)
  Y <- c()
  signal <- c()
  noise <- c()
  eta.true <- c()
  for(i in 1:n){
    epsi <- rnorm(npix,mean=0,sd=sigma)
    signali <- colSums(X[i,]*t(beta.true))
    etai <- eps1[i]*phi1+eps2[i]*phi2
    noisei <- etai+epsi
    Yi <- signali+noisei
    signal <- rbind(signal,signali)
    noise <- rbind(noise,noisei)
    eta.true <- rbind(eta.true,etai)
    Y <- rbind(Y,Yi)
  }
  Y[,ind.outside] <- NA
  signal[,ind.outside] <- NA
  noise[ind.outside] <- NA

  Geta.true <- t(eta.true)%*%eta.true/n

  list(Y=Y,X=X,Z=Z,beta.true=beta.true,signal=signal,noise=noise,eta.true=eta.true)
}
