
#' Resize parameter list
#'
#' Takes a list of estimated parameters, and changes the size to use as starting value for new fit
#'
#' @export
Resize_estimates = function( ParList, Data, N_factors, N_obsfactors ){

  # Local functions
  rmatrix = function( nrow, ncol, mean=0, sd=1 ) matrix( rnorm(nrow*ncol,mean=mean,sd=sd), nrow=nrow, ncol=ncol )
  generate = function( L_val, n_rows, n_old, n_new ){
    L_z = cov_matrix( L_val=L_val, n_rows=n_rows, n_cols=n_old, output="L_rc" )
    L_z = cbind( L_z, rmatrix(nrow=n_rows, ncol=n_rows-abs(n_old), sd=0.001) )
    L_z[upper.tri(L_z)] = NA
    L_z = L_z[,seq(1,abs(n_new),length=abs(n_new))]
    L_z = as.vector(t(L_z))
    L_z = L_z[which(!is.na(L_z))]
    if( n_new <= 0 ){
      if( n_old <= 0 ) L_z = c( L_z, L_val[c(length(L_val)-n_rows:1+1)] )
      if( n_old > 0 ) L_z = c( L_z, rep(log(0.001),n_rows) )
    }
    return(L_z)
  }

  # Combine
  ParList_new = ParList
  ParList_new[["L_z"]] = generate( L_val=ParList$L_z, n_rows=ncol(ParList$beta_gj), n_old=Data$Options_vec['n_factors'], n_new=N_factors )
  ParList_new[["obsL_z"]] = generate( L_val=ParList$obsL_z, n_rows=ncol(ParList$beta_gj), n_old=Data$Options_vec['n_obsfactors'], n_new=N_obsfactors )
  # L_val=ParList$obsL_z; n_rows=ncol(ParList$beta_gj); n_old=Data$Options_vec['n_obsfactors']; n_new=N_obsfactors

  # Other changes
  if( "L_logmult_col" %in% names(ParList) ){
    ParList_new[["L_logmult_col"]] = c( ParList[["L_logmult_col"]], rep(0,ncol(ParList$beta_gj)-length(ParList[["L_logmult_col"]])) )
    ParList_new[["L_logmult_col"]] = ParList_new[["L_logmult_col"]][1:abs(N_factors)]
  }
  if( "obsL_logmult_col" %in% names(ParList) ){
    ParList_new[["obsL_logmult_col"]] = c( ParList[["obsL_logmult_col"]], rep(0,ncol(ParList$beta_gj)-length(ParList[["obsL_logmult_col"]])) )
    ParList_new[["obsL_logmult_col"]] = ParList_new[["obsL_logmult_col"]][1:abs(N_obsfactors)]
  }

  # insibile return
  return( ParList_new )
}
