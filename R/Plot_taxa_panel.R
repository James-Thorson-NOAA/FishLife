
#' @export
Plot_taxa_panel = function( lowerTaxa, upperTaxa=lowerTaxa, prob=0.95, params=names(Estimate_database$Y_ij),
  Mean_gj=Estimate_database$ParHat$beta_gj, xlim=log(c(0.01,2)), ylim=xlim, ticks=c(1,2,5), partial_match=FALSE, verbose=FALSE, ... ){

  # Loop through all pairs of parameters
  par( mfrow=c(length(params),length(params)), mar=c(2,2,0,0), mgp=c(1.75,0.25,0), tck=-0.02, oma=c(2,3,0,0))
  for( rowI in 1:length(params)){
  for( colI in 1:length(params)){
    if(partial_match==FALSE) lowerWhich = match(lowerTaxa, ParentChild_gz[,'ChildName'])
    if(partial_match==FALSE) upperWhich = match(upperTaxa, ParentChild_gz[,'ChildName'])
    if(partial_match==TRUE) lowerWhich = grep(lowerTaxa, ParentChild_gz[,'ChildName'])
    if(partial_match==TRUE) upperWhich = grep(upperTaxa, ParentChild_gz[,'ChildName'])
    lowerY = Mean_gj[lowerWhich,params[rowI]]
    lowerX = Mean_gj[lowerWhich,params[colI]]
    upperY = Mean_gj[upperWhich,params[rowI]]
    upperX = Mean_gj[upperWhich,params[colI]]
    if( rowI==colI ){
      # If on diagonal, plot histogram
      if( params[rowI]!="Temperature" ) lowerX = exp(lowerX)
      hist( lowerX, breaks=25, col=rgb(0,0,0,0.2), yaxs="i", xlab="", ylab="", main="" )
    }else{
      # If on off-diagonal, plot ellipse for pairwise covariance
      Params = c(params[colI], params[rowI]) # first: X, second: Y
      Taxa = list(lowerTaxa, upperTaxa)[[ifelse( rowI>colI, 1, 2 )]]
      for( uniqueI in 1:length(unique(Taxa)) ){
        Plot_trait( Taxon=Taxa[uniqueI], params=Params, add=ifelse(uniqueI==1,FALSE,TRUE), xlim=NULL, lcol="black", lty="solid", lwd=1, partial_match=partial_match, main="", verbose=verbose, ... )
      }
    }
    if( colI==1 ) mtext(side=2, text=params[rowI], line=2, cex=1.3 )
    if( rowI==length(params) ) mtext(side=1, text=params[colI], line=2, cex=1.3 )
  }}
}
