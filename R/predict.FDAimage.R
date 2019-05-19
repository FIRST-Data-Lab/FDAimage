#' Make predictions from a fitted multivariate spline image-on-scalar regression model.
#'
#' This function is used to make predictions of the multivariate spline image-on-scalar regression models.
#'
#' @importFrom Matrix Matrix
#' @importFrom BPST basis
#'
#' @param mfit Fitted ``FDAimage" object.
#' \cr
#' @param Xpred The design matrix for prediction.
#' \cr
#' @param Zpred The cooridinates for prediction.
#' \cr
#' @return A matrix of predicted images is returned.
#'
#' @details This R package is the implementation program for manuscript entitled "Multivariate Spline Estimation and Inference for Image-on-Scalar Regression" by Shan Yu, Guannan Wang, Li Wang and Lijian Yang.
#'
#' @export
#'
predict.FDAimage <- function(mfit,Xpred,Zpred){
  if(all.equal(Zpred,Z)){
    Ypred <- tcrossprod(Xpred,mfit$beta)
  }else{
    V <- mfit$V; Tr <- mfit$Tr; d <- mfit$d; r <- mfit$r;
    Ball <- basis(V,Tr,d,r,Zpred,FALSE,FALSE)
    Bpred <- Ball$B
    ind.inside.pred <- Ball$Ind.inside
    beta.pred <- Bpred%*%mfit$gamma
    Ypred <- tcrossprod(Xpred,beta.pred)
  }
  return(Ypred)
}
