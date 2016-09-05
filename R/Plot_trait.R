
#' @export
Plot_trait = function( Taxon, params=c('K','M'), Cov_gjj=Estimate_database$Cov_gjj, Mean_gj=Estimate_database$ParHat$beta_gj,
  ParentChild_gz=Estimate_database$ParentChild_gz, SpeciesMatch=NULL, prob=0.95, add=FALSE,
  xlim=log(c(0.01,2)), ylim=xlim, partial_match=TRUE, main="", xlab="", ylab="",
  lcol="black", plot_lines=FALSE, verbose=FALSE, ticks=c(0,5), logticks=c(1,2,5), ... ){

  # Match taxon
  if(partial_match==TRUE) Which = grep(Taxon, ParentChild_gz[,'ChildName'])
  if(partial_match==FALSE) Which = which(Taxon == ParentChild_gz[,'ChildName'])
  if( length(Which)!=1 ) stop( paste0("'Taxon' ",Taxon," input matches more or less than one element") )
  if(verbose==TRUE) print( ParentChild_gz[Which,] )

  # Plot ellipse
  Plot_ellipse( Cov=Cov_gjj[Which,params,params], Mean=Mean_gj[Which,params], add=add, whichlog=paste(c("x","y")["Temperature"!=params],collapse=""), xlim=xlim, ylim=ylim, main=main, xlab=xlab, lcol=lcol, plot_lines=plot_lines, ticks=ticks, logticks=logticks, ... )

  # Plot observations
  if( !is.null(SpeciesMatch) ){
    Which = grep(SpeciesMatch, ParentChild_gz[,'ChildName'])
    Which = which( g_i %in% Which )
    points( x=Y_ij[Which,'K'], y=Y_ij[Which,'M'] )
  }

  # insibile return
  return( invisible(list("Cov_pred"=Cov_gjj[Which,,])) )
}
