
#' @export
Calculate_ratio = function( params=c("K","M"), Cov_gjj=Estimate_database$Cov_gjj, Mean_gj=Estimate_database$ParHat$beta_gj,
  ParentChild_gz=Estimate_database$ParentChild_gz ){

  # Rotation matrix
  RotateM = function( angle ) matrix( c(cos(angle),sin(angle),-sin(angle),cos(angle)), 2,2 )
  # Rotate covariance matrix via eigen-decomposition
  RotateCov = function(angle, cov){
    Eigen = eigen(cov)
    newcov = RotateM(angle)%*%Eigen$vectors %*% diag(Eigen$values) %*% solve(RotateM(angle)%*%Eigen$vectors)
    return( newcov )
  }

  # Calculate ratio
  Ratio_gz = array(NA, dim=c(nrow(ParentChild_gz),2), dimnames=list(NULL,c("median","logSD")) )
  for( gI in 1:nrow(Mean_gj)){
    Cov_hat = Cov_gjj[gI,params,params]
    #Angle_hat = atan(eigen(Cov_hat)$vector[2,1] / eigen(Cov_hat)$vector[1,1])
    Cov = RotateCov(angle=-3/4*pi, Cov_hat)
    MargVar = sum(abs(eigen(Cov)$vectors[1,] * eigen(Cov)$values))
    Median = Mean_gj[gI,params]
    Ratio_gz[gI,] = c( exp(Median[2]-Median[1]), sqrt(MargVar) )
  }
  return( Ratio_gz )
}
