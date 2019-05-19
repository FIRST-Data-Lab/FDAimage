cover_boot <- function(Zalpha,beta.hat,betab.hat,Sigmab){
  cc.l <- betab.hat-Zalpha*sqrt(Sigmab)
  cc.u <- betab.hat+Zalpha*sqrt(Sigmab)
  cover <- (beta.hat>=cc.l & beta.hat<=cc.u)
  return(cover)
}
