
#' @export
Predictive_distribution <-
function( mean_vec,
          process_cov,
          obs_cov,
          Samp_rj,
          group_j,
          include_obscov = FALSE,
          check_bounds = TRUE,
          check_names = FALSE,
          lowerbound_MLSPS = 0,
          rho_option = 0,
          include_r = FALSE ){

  # Pull out inputs
  var_orig = var_names = names(mean_vec)
  n_j = length(mean_vec)
  groupR_j = group_j + 1     # Switch from CPP to R indexing

  # Informative errors
  if( include_obscov!=FALSE ){
    stop("Input `include_obscov` is only implemented for FALSE")
  }
  if( any(table(group_j)>1) & all(c("ln_var","rho","ln_MASPS") %in% var_names) ){
    stop("`Predictive_distribution` not yet built to work with both factors and stock-recruit values")
  }

  ################
  # Part 1:   Extract variable names
  ################

  # Expand to include factors
  if( any(table(group_j)>1) ){
    # Define objects
    var_new = vector()
    groupfull_j = vector()
    for( varI in 1:length(unique(group_j)) ){
      if( sum(groupR_j==varI) > 1 ){
        var_new = c( var_new, "base", var_names[which(groupR_j==varI)] )
        groupfull_j = c( groupfull_j, groupR_j[which(groupR_j==varI)][1], groupR_j[which(groupR_j==varI)] )
      }else{
        var_new = c( var_new, var_names[which(groupR_j==varI)] )
        groupfull_j = c( groupfull_j, groupR_j[which(groupR_j==varI)] )
      }
    }
    var_names = var_new
  }
  # Expand to include stock-recruit stuff
  if( all(c("ln_var","rho","ln_MASPS") %in% var_names) ){
    # Define objects
    var_names = union( var_names, c("rho", "ln_margsd", "h", "logitbound_h", "ln_Fmsy_over_M", "ln_Fmsy") )
    if( include_r==TRUE ){
      var_names = union( var_names, c("ln_r", "r", "ln_G", "G") )
    }
  }
  if( check_names==TRUE ) return( var_names )

  ################
  # Part 2:   Expand values
  ################

  # Expand to include factors
  n_v = length(var_names)
  if( any(table(group_j)>1) ){
    pred_cov = array(NA, dim=c(n_v,n_v), dimnames=list(var_names,var_names) )
    pred_mean = array(NA, dim=n_v, dimnames=list(var_names) )
    # Experiment with projection matrix
    Proj = data.frame( "varname" = var_names )
    Proj = cbind(Proj, "colnum"=match(var_names, var_orig) )
    Proj = cbind(Proj, "group"=groupfull_j )
    Proj = cbind(Proj, "is_factor"=ifelse( Proj$group %in% which(table(Proj$group)>1), TRUE, FALSE ) )

    # Set up projected values beta_gw
    Samp_rv = array( NA, dim=c(nrow(Samp_rj),nrow(Proj)), dimnames=list(NULL,Proj[,'varname']) )

    # Add values
    Samp_rv[] = Samp_rj[ , Proj[,'colnum'] ]

    # Transform factor values
    Samp_rv[] = ifelse( is.na(Samp_rv), 0, Samp_rv)
    Samp_rv[,which(Proj$is_factor)] = exp(Samp_rv[,which(Proj$is_factor)])

    # Add missing default values
    for( groupI in 1:max(Proj$group) ){
      if( all(Proj$is_factor[which(Proj$group==groupI)])==TRUE ){
        tmp_gw = Samp_rv[,which(Proj$group==groupI),drop=FALSE]
        tmp_gw = tmp_gw / outer(rowSums(tmp_gw),rep(1,ncol(tmp_gw)))
        Samp_rv[,which(Proj$group==groupI)] = tmp_gw
      }
    }
    # Save results
    pred_mean = colMeans( Samp_rv, na.rm=TRUE )
    pred_cov[,] = cov( Samp_rv, use="pairwise.complete" )
  }else
  # Expand to include stock-recruit stuff
  if( all(c("ln_var","rho","ln_MASPS") %in% var_names) ){
    pred_cov = array(NA, dim=c(n_v,n_v), dimnames=list(var_names,var_names) )
    pred_mean = array(NA, dim=n_v, dimnames=list(var_names) )
    # Sample
    if( missing(Samp_rj) ){
      Samp_rj = mvtnorm::rmvnorm( n=1000, sigma=process_cov + obs_cov*include_obscov, mean=mean_vec )
    }
    MLSPS_r = lowerbound_MLSPS + exp(Samp_rj[,'ln_MASPS']) / (1-exp(-exp(Samp_rj[,'M'])))

    # Add columns
    Samp_rv = Samp_rj
    if( rho_option %in% c(1,2) ){
      Samp_rv = cbind( Samp_rv, "rho"=(2*plogis(Samp_rv[,'logit_rho']))-1 )
    }
    Samp_rv = cbind( Samp_rv, "ln_margsd"=log(exp(Samp_rv[,'ln_var'])/( 1-ifelse(Samp_rv[,'rho']>0.99,NA,Samp_rv[,'rho'])^2 ))/2 )
    Samp_rv = cbind( Samp_rv, "h"=MLSPS_r/(4+MLSPS_r) )
    Samp_rv = cbind( Samp_rv, "logitbound_h"=qlogis((ifelse(Samp_rv[,'h']<0.2,NA,Samp_rv[,'h'])-0.2)*5/4) )
    Samp_rv = cbind( Samp_rv, "ln_Fmsy_over_M"=log(sqrt((4*Samp_rv[,'h'])/(1-Samp_rv[,'h'])) - 1) )  #  From Mangel et al. 2013 Eq. 13
    Samp_rv = cbind( Samp_rv, "ln_Fmsy"=Samp_rv[,'ln_Fmsy_over_M'] + Samp_rv[,'M'] )

    if( include_r==TRUE ){
      Mat_zv = t(apply( Samp_rv, MARGIN=1, FUN=get_r ))
      Samp_rv = cbind( Samp_rv, "ln_r"=log(Mat_zv[,'intrinsic_growth_rate']), 'r'=Mat_zv[,'intrinsic_growth_rate'], 'ln_G'=log(Mat_zv[,'generation_time']), 'G'=Mat_zv[,'generation_time'] )
    }

    # Check for problems
    if( check_bounds==TRUE ){
      Stop = FALSE
      if( mean(Samp_rv[,'rho']>0.99)>0.5 ){ message("Median sampled rho is above upper bound of 0.99"); Stop = TRUE }
      if( mean(Samp_rv[,'h']<0.2)>0.9 ){ message("90% of sampled h is below upper bound of 0.20"); Stop = TRUE }
      if( Stop==TRUE ) return( Samp_rv )
    }
    # Ensure that var_names and Samp_rv have same order
    Samp_rv = Samp_rv[,var_names]
    # Save results
    pred_mean = colMeans( Samp_rv, na.rm=TRUE )
    pred_cov[,] = cov( Samp_rv, use="pairwise.complete" )
    pred_mean[1:n_j] = mean_vec
    if( !missing(process_cov) & !missing(obs_cov) ){
      pred_cov[1:n_j,1:n_j] = process_cov + obs_cov*include_obscov
    }
  }else{
    if(missing(process_cov)){
      process_cov = cov( Samp_rj, use="pairwise.complete" )
    }
    if(missing(obs_cov)){
      obs_cov = matrix( 0, nrow=ncol(Samp_rj), ncol=ncol(Samp_rj) )
    }
    pred_cov = process_cov + obs_cov*include_obscov
    pred_mean = mean_vec
  }

  Return = list("pred_cov"=pred_cov, "pred_mean"=pred_mean)
  return( Return )
}
