#' Fitting multivariate varying coefficient models
#'
#' This function is used to fit the multivariate varying coefficient models.
#'
#' @importFrom Matrix Matrix
#' @importFrom BPST basis
#' @param Y The response images of dimension \code{n} by \code{npix}, where \code{n} is the number of observed images and \code{npix} is the number of pixels/voxels in each image. Each row is an observation of the images.
#' \cr
#' @param X The design matrix of dimension \code{n} by \code{p}, with an intercept. Each row is an observation vector.
#' \cr
#' @param Z The cooridinates of dimension \code{npix} by two. Each row is the coordinates of a pixel/voxel.
#' \cr
#' @param V The \code{N} by two matrix of vertices of a triangulation, where \code{N} is the number of vertices. Each row is the coordinates for a vertex.
#' \cr
#' @param Tr The triangulation matrix of dimention \code{nT} by three, where \code{nT} is the number of triangles in the triangulation. Each row is the indices of vertices in \code{V}.
#' \cr
#' @param d The degree of piecewise polynomials -- default is 5, and -1 represents piecewise constant.
#' \cr
#' @param r The smoothness parameter -- default is 1, and 0 \eqn{\le} \code{r} \eqn{<} \code{d}.
#' \cr
#' @param lambda The vector of the candidates of penalty parameter -- default is grid points of 10 to the power of a sequence from -6 to 6 by 0.5.
#' \cr
#' @return An object with S3 class "FDAimage", including
#' \item{beta}{The estimated coefficient functions.}
#' \item{Yhat}{The estimated images.}
#' \item{sse}{Sum of squared errors.}
#' \item{df}{Effective degree of freedom.}
#' \item{gcv}{Generalized cross-validation (GCV).}
#' \item{lamc}{Selected tuning parameter for bivariate penalized spline based on GCV.}
#'
#' @details This R package is the implementation program for manuscript entitled ``Multivariate Spline Estimation and Inference for Varying Coeffiient Models with Imaging Data" by Shan Yu, Guannan Wang, Li Wang and Lijian Yang.
#'
#' @examples
#' # Libraries and source files needed
#' library(devtools)
#' install_github("funstatpackages/BPST")
#' library(BPST)
#' install_github("funstatpackages/Triangulation")
#' library(Triangulation)
#' # Location information
#' n1=95; n2=79;
#' u1=seq(0,1,length.out=n1)
#' v1=seq(0,1,length.out=n2)
#' uu=rep(u1,each=n2)
#' vv=rep(v1,times=n1)
#' uu.mtx=matrix(uu,n2,n1)
#' vv.mtx=matrix(vv,n2,n1)
#' Z=as.matrix(cbind(uu,vv))
#' # Parameters
#' d=5; r=1; rho=0.5; nfold=10; alpha0=0.05;
#' # Triangulation
#' data(V1); data(Tr1); # triangulation_1 in the paper;
#' # data(brain_boundary); # brain imaging boundary of slide 48;
#' V=V1; Tr=Tr1;
#' ind=inVT(V,Tr,Z[,1],Z[,2])
#' ind.inside=ind$ind.inside; ind=ind$ind;
#' n=50; sigma=0.5; lambda1=0.03; lambda2=0.006;
#' dat=data.FDAimage(n,Z,ind.inside,sigma,rho,2019,lambda1,lambda2)
#' Y=dat$Y; X=dat$X; Z=dat$Z;
#' Y=Y[,ind.inside] # remove the points which are outside the boundary;
#' lambda=10^(seq(-6,6,by=1))
#' cv=cv.FDAimage(Y,X,Z,V,Tr,d,r,lambda,nfold)
#' lamc=cv$lamc
#' mfit0=fit.FDAimage(Y,X,Z,V,Tr,d,r,lamc)
#'
#' @export
#'
fit.FDAimage <- function(Y,X,Z,V,Tr,d=5,r=1,lambda=0.1){
  if(!is.matrix(Y)){
    warning("The response variable, Y, should be a matrix with each row represents an image.")
    Y <- as.matrix(Y)
  }
  if(!is.matrix(X)){
    warning("The explanatory variable, X, should be a matrix.")
    X <- as.matrix(X)
  }
  np <- ncol(X)
  if(!is.matrix(Z)){
    warning("The coordinates of each pixel/voxel, Z, should be a matrix.")
    Z <- as.matrix(Z)
  }
  lambda <- as.matrix(lambda)
  if(nrow(lambda)>1){
    warning("The tuning parameter, lambda, should be a scalar. Instead, the default 0.1 is used.")
    lambda <- as.matrix(0.1)
  }
  this.call <- match.call()
  Ball <- basis(V,Tr,d,r,Z)
  K <- Ball$K
  Q2 <- Ball$Q2
  B <- Ball$B
  ind.inside <- Ball$Ind.inside
  npix <- length(ind.inside)
  tria.all <- Ball$tria.all

  if(d>1){
    mfit <- FDAimage.est.ho(B,Q2,K,X,Y,lambda)
  }
  if(d== -1){
    mfit <- FDAimage.est.pc(B,X,Y)
  }
  mfit$beta <- matrix(unlist(mfit$beta),nrow=npix,ncol=np);
  mfit$gamma <- matrix(unlist(mfit$gamma),ncol=np);
  mfit$V <- V; mfit$Tr <- Tr; mfit$d <- d; mfit$r <- r;
  mfit$B <- B; mfit$Q2 <- Q2; mfit$K <- K;
  mfit$ind.inside <- ind.inside; mfit$tria.all <- tria.all;
  mfit$lamc <- lambda;
  mfit$X <- X; mfit$Y <- Y; mfit$Z <- Z;
  mfit$call <- this.call;
  class(mfit) <- "FDAimage"
  return(mfit)
}
