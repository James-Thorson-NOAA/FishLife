
#' Calculate ratio of two modeled traits
#'
#' \code{Calculate_ratio} projects the covariance among traits onto a single difference (interpreted as a ratio for log-scaled traits)
#'
#' @param params Parameter names to use in calculation
#' @param Cov_gjj Array of estimated covariance for each taxonomic group \code{g} and trait \code{j}
#' @param Mean_gj Matrix of Empirical Bayes predictions of traits for each taxonomic group
#' @param ParentChild_gz Matrix representing taxonomic tree for analyzed data

#' @return Matrix of exponentiated-differences (ratios) and log-standard deviations for projection
#' \describe{
#'   \item{median}{median difference for taxonomic group}
#'   \item{N_obsfactors}{log-standard deviation of difference for taxonomic group}
#' }

#' @export
Calculate_ratio = function( params=c("K","M"), Cov_gjj=FishLife::database$Cov_gjj, Mean_gj=FishLife::database$ParHat$beta_gj,
  ParentChild_gz=FishLife::database$ParentChild_gz ){

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
