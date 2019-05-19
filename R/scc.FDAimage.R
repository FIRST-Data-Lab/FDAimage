#' Simultaneous confidence corridors for Image-on-Scalar Regression
#'
#' This function is used to construct the simultaneous confidence corridors (SCC) for the coefficient functions in image-on-scalar regression.
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
#' @param alpha0 The nominal level of the SCC -- default is 0.05.
#' \cr
#' @param nboot The number of bootstraps -- default is 100.
#' @return A list of the lower SCC and upper SCC for the coefficient functions.
#'
#' @details This R package is the implementation program for manuscript entitled "Multivariate Spline Estimation and Inference for Image-on-Scalar Regression" by Shan Yu, Guannan Wang, Li Wang and Lijian Yang.
#'
#' @export
#'
scc.FDAimage <- function(Y,X,Z,V,Tr,d=5,r=1,lambda=10^seq(-6,6,by=0.5),
                        alpha0=0.05,nboot=100,ksi=0.01){
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

  if(d>1){
    CC <- FDAimage.CC.ho(B,Q2,K,X,Y,lambda,alpha0,nboot,ksi)
  }
  if(d== -1){
    CC <- FDAimage.CC.pc(B,X,Y,tria.all,alpha0,nboot,ksi)
  }
  return(CC)
}
