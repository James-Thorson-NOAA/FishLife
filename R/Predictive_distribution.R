
#' @export
Predictive_distribution = function( mean_vec, process_cov, obs_cov, include_obscov=FALSE, check_bounds=TRUE, check_names=FALSE, lowerbound_MLSPS=0, rho_option=0, include_r=FALSE ){

  var_names = names(mean_vec)
  n_j = length(mean_vec)

  if( all(c("ln_var","rho","ln_MASPS") %in% var_names) ){
    # Define objects
    var_names = union( var_names, c("rho", "ln_margsd", "h", "logitbound_h", "ln_Fmsy_over_M", "ln_Fmsy") )
    if( include_r==TRUE ){
      var_names = union( var_names, c("ln_r", "r", "ln_G", "G") )
    }
  }
  if( check_names==TRUE ) return( var_names )

  if( all(c("ln_var","rho","ln_MASPS") %in% var_names) ){
    n_v = length(var_names)
    pred_cov = array(NA, dim=c(n_v,n_v), dimnames=list(var_names,var_names) )
    pred_mean = array(NA, dim=n_v, dimnames=list(var_names) )
    # Sample
    Samp_zv = mvtnorm::rmvnorm( n=1000, sigma=process_cov + obs_cov*include_obscov, mean=mean_vec )
    MLSPS = lowerbound_MLSPS + exp(Samp_zv[,'ln_MASPS']) / (1-exp(-exp(Samp_zv[,'M'])))
    if( rho_option %in% c(1,2) ){
      Samp_zv = cbind( Samp_zv, "rho"=(2*plogis(Samp_zv[,'logit_rho']))-1 )
    }
    Samp_zv = cbind( Samp_zv, "ln_margsd"=log(exp(Samp_zv[,'ln_var'])/( 1-ifelse(Samp_zv[,'rho']>0.99,NA,Samp_zv[,'rho'])^2 ))/2 )
    Samp_zv = cbind( Samp_zv, "h"=MLSPS/(4+MLSPS) )
    Samp_zv = cbind( Samp_zv, "logitbound_h"=qlogis((ifelse(Samp_zv[,'h']<0.2,NA,Samp_zv[,'h'])-0.2)*5/4) )
    Samp_zv = cbind( Samp_zv, "ln_Fmsy_over_M"=log(sqrt((4*Samp_zv[,'h'])/(1-Samp_zv[,'h'])) - 1) )  #  From Mangel et al. 2013 Eq. 13
    Samp_zv = cbind( Samp_zv, "ln_Fmsy"=Samp_zv[,'ln_Fmsy_over_M'] + Samp_zv[,'M'] )

    if( include_r==TRUE ){
      # Linf=exp(vec['Loo']); K=exp(vec['K']); t0=-0.1; W_a=0.001; W_b=3.04; tm=exp(vec['tm']); dm=tm/4; minage=0; maxage=ceiling(min(100,2*exp(vec['tmax']))); h=vec['h']; M=exp(vec['M'])
      Mat_zv = t(apply( Samp_zv, MARGIN=1, FUN=get_r ))
      Samp_zv = cbind( Samp_zv, "ln_r"=log(Mat_zv[,'intrinsic_growth_rate']), 'r'=Mat_zv[,'intrinsic_growth_rate'], 'ln_G'=log(Mat_zv[,'generation_time']), 'G'=Mat_zv[,'generation_time'] )
    }

    # Check for problems
    if( check_bounds==TRUE ){
      Stop = FALSE
      if( mean(Samp_zv[,'rho']>0.99)>0.5 ){ message("Median sampled rho is above upper bound of 0.99"); Stop = TRUE }
      if( mean(Samp_zv[,'h']<0.2)>0.9 ){ message("90% of sampled h is below upper bound of 0.20"); Stop = TRUE }
      if( Stop==TRUE ) return( Samp_zv )
    }
    # Ensure that var_names and Samp_zv have same order
    Samp_zv = Samp_zv[,var_names]
    # Save results
    pred_mean[1:n_j] = mean_vec
    pred_mean[n_j+1:(n_v-n_j)] = colMeans( Samp_zv[,n_j+1:(n_v-n_j)], na.rm=TRUE )
    pred_cov[,] = cov( Samp_zv, use="pairwise.complete" )
    pred_cov[1:n_j,1:n_j] = process_cov + obs_cov*include_obscov
  }else{
    pred_cov = process_cov + obs_cov*include_obscov
    pred_mean = mean_vec
  }

  Return = list("pred_cov"=pred_cov, "pred_mean"=pred_mean)
  return( Return )
}
