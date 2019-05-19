#' Cross-validation for Imga-on-Scalar Regression
#'
#' This function implements k-fold cross-validation for image-on-scalar regression, and returns the mean squared prediction error.
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
#' @param nfolds The number of folds -- default is 10. Although \code{nfold} can be as large as the sample size (leave-one-out CV), it is not recommended for large datasets. Smallest value allowable for \code{nfolds} is 3.
#' \cr
#' @param iter The seed used for cross-validation.
#' \cr
#' @return
#' \item{mspe}{The mean squared prediction error based on k-fold cross validation (CV).}
#' \item{lamc}{The tuning parameter selected by k-fold cross validation (CV).}
#' \cr
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
#' cv.FDAimage(Y,X,Z,V,Tr,d,r,lambda,nfold,iter)
#'
#' @export
#'
cv.FDAimage <- function(Y,X,Z,V,Tr,d=5,r=1,lambda=10^seq(-6,6,by=0.5),nfold=10,iter=2019){
  if(nfold<3){
    warning("The number of folds in CV is too small. Instead, the default 10-fold CV is used.")
    nfold <- 10
  }
  if(!is.matrix(Y)){
    warning("The response variable, Y, should be a matrix with each row represents an image.")
    Y <- as.matrix(Y)
  }
  if(!is.matrix(X)){
    warning("The explanatory variable, X, should be a matrix.")
    X <- as.matrix(X)
  }
  if(!is.matrix(Z)){
    warning("The coordinates of each pixel/voxel, Z, should be a matrix.")
    Z <- as.matrix(Z)
  }
  Ball <- basis(V,Tr,d,r,Z)
  K <- Ball$K
  Q2 <- Ball$Q2
  B <- Ball$B
  ind.inside <- Ball$Ind.inside
  npix <- length(ind.inside)
  tria.all <- Ball$tria.all

  n <- nrow(Y)
  npix <- ncol(Y)
  sfold <- round(n/nfold)
  set.seed(iter)
  Test <- sample(1:n)
  cv.error <- c()
  for(ii in 1:nfold){
    if(ii<nfold){
      Test.set <- sort(Test[((ii-1)*sfold+1):(ii*sfold)])
    }
    if(ii ==nfold){
      Test.set <- sort(Test[((ii-1)*sfold+1):n])
    }
    Train.set <- setdiff(1:n,Test.set)
    X.test <- X[Test.set,]
    X.train <- X[Train.set,]
    Y.test <- Y[Test.set,]
    Y.train <- Y[Train.set,]

    if(d>1){
      mfit.ii <- FDAimage.est.ho(B,Q2,K,X.train,Y.train,lambda)
      nl <- length(lambda)
      pred.error <- lapply(1:nl,function(il) sum((Y.test-tcrossprod(X.test,mfit.ii$beta[[il]]))^2))
      pred.error <- unlist(pred.error)
      cv.error <- cbind(cv.error,pred.error)
    }
    if(d== -1){
      mfit.ii <- FDAimage.est.pc(B,X.train,Y.train)
      Ypred.ii <- tcrossprod(X.test,mfit.ii$beta)
      pred.error <- sum((Y.test-Ypred.ii)^2)
      cv.error <- c(cv.error,pred.error)
    }
  }
  if(d>1){
    j <- which.min(apply(cv.error,1,sum))
    lamc <- lambda[j]
    cv.error <- cv.error[j,]
  }
  list(mspe=sum(cv.error)/n/npix,lamc=lamc)
}
