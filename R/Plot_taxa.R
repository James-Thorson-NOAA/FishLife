
#' Plots trait predictions for taxa
#'
#' Plots bivariate predictive intervals for specified parameters and taxa
#'
#' @param Taxa One or more taxonomic names (matching rows from \code{ParentChild_gz[,'ChildName']})
#' @param params matrix (2 columns) of parameters to plot, where each row gives parameter names for a given plot
#' @param drop_pred Boolean whether to drop "_predictive" from legend
#' @param mfrow numeric vector giving number of rows and columns for panel figure
#' @param legendnum which panel to use for legend
#' @inheritParams Plot_ellipse
#' @inheritParams Calculate_ratio
#' @inheritParams Fit_model
#' @inheritParams Plot_trait

#' @export
Plot_taxa = function( Taxa, prob=0.95, params=matrix(c("K","M","Winfinity","Loo","tmax","tm","Lm","Temperature"),ncol=2,byrow=TRUE),
  Cov_gjj=FishLife::database$Cov_gjj, Mean_gj=FishLife::database$ParHat$beta_gj, ParentChild_gz=FishLife::database$ParentChild_gz,
  Y_ij=FishLife::database$Y_ij, xlim=log(c(0.01,2)), ylim=xlim, ticks=c(0,5), logticks=c(1,2,5), partial_match=FALSE, drop_pred=TRUE,
  mfrow=c(nrow(params),1), legendnum=2, verbose=FALSE, plot_lines=FALSE, lcol=rainbow(length(Taxa)), lty=rep("solid",length(Taxa)),
  xaxt="s", yaxt="s", ... ){

  # Loop through parameter-pairs
  par( mfrow=mfrow, mar=c(3,3,0,0), mgp=c(1.75,0.25,0), tck=-0.02, oma=c(1,1,1,1) )
  Pred_taxa = NULL
  for( rowI in 1:nrow(params)){
    # Loop through specified taxa
    for( uniqueI in 1:length(unique(Taxa)) ){
      Pred_taxa[[uniqueI]] = Plot_trait( Taxon=Taxa[uniqueI], params=params[rowI,], Cov_gjj=Cov_gjj, Mean_gj=Mean_gj, ParentChild_gz=ParentChild_gz, Y_ij=Y_ij, add=ifelse(uniqueI==1,FALSE,TRUE), xlim=quantile(Mean_gj[,params[rowI,1]],na.rm=TRUE,c(0,1)), ylim=quantile(Mean_gj[,params[rowI,2]],na.rm=TRUE,prob=c(0,1)), partial_match=partial_match, main="", lcol=lcol[uniqueI], ticks=ticks, logticks=logticks, plot_lines=plot_lines, verbose=verbose, prob=prob, lty=lty[uniqueI], xaxt=xaxt, yaxt=yaxt, ... )
      for( aI in 1:2 ){
        Text = switch( params[rowI,aI], "Loo"="Asymptotic length (L_inf)", "K"="Relative growth rate (K)", "Winfinity"="Asymptotic mass (W_inf)", "tmax"="Maximum age (A_max)", "tm"="Age at maturity (A_mat)", "M"="Mortality rate (M)", "Lm"="Length at maturity (L_mat)", "Temperature"="Average temperature", "ln_var"="Conditional recruitment variance", "rho"="Recruitment autocorrelation (rho)", "ln_MASPS"="Maximum annual spawners per spawner", "ln_margsd"="SD of recruitment (Sigma_R)", "h"="Steepness (h)", "logitbound_h"="Steepness (h)", "ln_Fmsy_over_M"="Ratio of F_msy and M", "ln_Fmsy"="Fishing mortality rate at MSY", "ln_r"="Intrinsic growth rate (r)", "ln_G"="Generation time", "r"="Intrinsic growth rate (r)", "G"="Generation time", params[rowI,aI] )
        mtext( side=aI, text=Text, line=1.5 )
      }
    }
    # Add legend
    if( rowI%in%legendnum ){
      Legend = Taxa
      if( drop_pred==TRUE ) Legend = gsub(x=Taxa, pattern="_predictive", replacement="")
      legend("topleft", legend=Legend, fill=lcol, bty="n")
    }
  }

  # insibile return
  return( invisible(Pred_taxa) )
}
