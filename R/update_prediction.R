
#' Update predictions from a model
#'
#' Updates the predictive distribution for each trait analytically as
#'   precision-weighted average of the existing prediction and new sampled values. Presumably
#'   the function will be run such that samples are weighted by the estimated measurement-covariance.
#'   The measurement-covariance is edited in each sample to account for missing values,
#'   which are inputted as NA values;  this approximates the analytic calculation
#'   without requiring integrating across missing samples.
#'
#' @param predmean_j vector of predictions for a given taxon
#' @param predcov_jj matrix of estimated covariance for predictions
#' @param obscov_jj matrix of estimated sampling imprecision for new samples
#' @param Ynew_ij matrix of new samples (NAs are ignored and have essentially no effect)
#'
#' @examples
#' \dontrun{
#'  # New values
#'  linf = c(61.1, 50.2, 41.8, 35.4, 52.8, 46.3, 53.4, 54.1, 43.9, 45.5, 50.7, 55.1)
#'  vbk = c(0.28, 0.45, 0.57, 0.56, 0.1426, 0.3669, 0.182, 0.3764, 0.34, 0.31, 0.454, 0.246)
#'
#'  # Extract estimates
#'  predict_GP = Plot_taxa( Search_species(Genus="Macquaria",Species="ambigua", add_ancestors=FALSE)$match_taxonomy, mfrow=c(3,2) )
#'
#'  # Format new
#'  Ynew_ij = matrix(NA, nrow=length(linf), ncol=length(predict_GP[[1]]$Mean_pred) )
#'  colnames(Ynew_ij) = names(predict_GP[[1]]$Mean_pred)
#'  Ynew_ij[,"Loo"] = log(linf)
#'  Ynew_ij[,"K"] = log(vbk)
#'
#'  # Update
#'  which_cols = which( names(predict_GP[[1]]$Mean_pred) %in% colnames(FishLife::FishBase_and_RAM$obsCov_jj) )
#'  Update = update_prediction( predmean_j = predict_GP[[1]]$Mean_pred[which_cols],
#'                     predcov_jj = predict_GP[[1]]$Cov_pred[which_cols,which_cols],
#'                     obscov_jj = FishLife::FishBase_and_RAM$obsCov_jj,
#'                     Ynew_ij = Ynew_ij[,which_cols] )
#'
#'  # Check
#'  cbind( "Orig"=Update$updatemean_j, "New"=colMeans(Ynew_ij[,which_cols],na.rm=TRUE), "Updated"=predict_GP[[1]]$Mean_pred[which_cols] )
#' }
#'
#' @export
update_prediction <-
function( predmean_j,
          predcov_jj,
          obscov_jj,
          Ynew_ij ){

  # Local function
  replace_missing = function(b, v, largevar=max(c(1000,diag(v)*1000)) ){
    R = cov2cor(v)
    SD = sqrt(diag(v))
    SD = ifelse( is.na(b), sqrt(largevar), SD )
    return( diag(SD) %*% R %*% diag(SD) )
  }

  #
  if(is.vector(Ynew_ij)){
    Ynew_ij = matrix(Ynew_ij,nrow=1)
  }

  # Calculate running precision
  Q_jj = solve(predcov_jj)
  b_j = Q_jj %*% predmean_j
  for( i in 1:nrow(Ynew_ij) ){
    Q_jj = Q_jj + solve(replace_missing( b=Ynew_ij[i,], v=obscov_jj ))
    b_j = b_j + solve(replace_missing( b=Ynew_ij[i,], v=obscov_jj )) %*% ifelse( is.na(Ynew_ij[i,]), predmean_j, Ynew_ij[i,] )
  }
  V_jj = solve(Q_jj)
  b_j = V_jj %*% b_j

  # return stuff
  return(list(updatemean_j=b_j, updatecov_jj=V_jj))
}
