
#' Plots trait predictions
#'
#' Plots predictive interval for a given taxon and trait-pair
#'
#' @param Taxon A taxon (matching row from \code{ParentChild_gz[,'ChildName']})
#' @param params character vector (length of two) of parameters to plot
#' @param SpeciesMatch Species for which to plot data
#' @param partial_match Should \code{Taxon} be an partial match or not (exact match)
#' @param verbose Boolean whether to print name matches or not
#' @param g_i Vector that associates every observation with a level of the taxonomic tree
#' @inheritParams Plot_ellipse
#' @inheritParams Fit_model
#' @inheritParams Calculate_ratio

#' @export
Plot_trait = function( Taxon, params=c('K','M'), Cov_gjj=FishLife::database$Cov_gjj, Mean_gj=FishLife::database$ParHat$beta_gj,
  ParentChild_gz=FishLife::database$ParentChild_gz, Y_ij=FishLife::database$Y_ij, g_i=FishLife::database$g_i,
  SpeciesMatch=NULL, prob=0.95, add=FALSE, xlim=log(c(0.01,2)), ylim=xlim, partial_match=TRUE, main="", xlab="", ylab="",
  lcol="black", plot_lines=FALSE, verbose=FALSE, ticks=c(0,5), logticks=c(1,2,5), obsCov_jj=FishLife::database$obsCov_jj,
  include_obscov=FALSE, lty="solid", xaxt="s", yaxt="s", ... ){

  # Match taxon
  if(partial_match==TRUE) Which = grep(Taxon, ParentChild_gz[,'ChildName'])
  if(partial_match==FALSE) Which = which(Taxon == ParentChild_gz[,'ChildName'])
  if( length(Which)!=1 ) stop( paste0("'Taxon' ",Taxon," input matches more or less than one element") )
  if(verbose==TRUE) print( ParentChild_gz[Which,] )

  # Plot ellipse
  #Plot_ellipse( Cov=Cov_gjj[Which,params,params], Mean=Mean_gj[Which,params], add=add, whichlog=paste(c("x","y")[which(!params%in%c("Temperature","h"))],collapse=""), xlim=xlim, ylim=ylim, main=main, xlab=xlab, lcol=lcol, plot_lines=plot_lines, ticks=ticks, logticks=logticks, prob=prob, ... )
  Plot_ellipse( Cov=Cov_gjj[Which,params,params], Mean=Mean_gj[Which,params], add=add, axis_scale=sapply(params,FUN=switch,"Temperature"="natural","h"="natural","rho"="natural","G"="natural","r"="natural","logitbound_h"="logit_0.2_1.0","log"), xlim=xlim, ylim=ylim, main=main, xlab=xlab, lcol=lcol, plot_lines=plot_lines, ticks=ticks, logticks=logticks, prob=prob, lty=lty, xaxt=xaxt, yaxt=yaxt, ... )

  # Plot observations
  if( !is.null(SpeciesMatch) ){
    Which = grep(SpeciesMatch, ParentChild_gz[,'ChildName'])
    Which = which( g_i %in% Which )
    points( x=Y_ij[Which,params[1]], y=Y_ij[Which,params[2]] )
  }

  # insibile return
  return( invisible(list("Cov_pred"=Cov_gjj[Which,,], "Mean_pred"=Mean_gj[Which,])) )
}
