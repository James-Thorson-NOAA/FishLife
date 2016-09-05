
#' @export
Plot_taxa = function( Taxa, prob=0.95, params=matrix(c("K","M","Winfinity","Loo","tmax","tm","Lm","Temperature"),ncol=2,byrow=TRUE),
  Cov_gjj=Estimate_database$Cov_gjj, Mean_gj=Estimate_database$ParHat$beta_gj, ParentChild_gz=Estimate_database$ParentChild_gz,
  Y_ij=Estimate_database$Y_ij, xlim=log(c(0.01,2)), ylim=xlim, ticks=c(0,5), logticks=c(1,2,5), partial_match=FALSE, drop_pred=TRUE,
  mfrow=c(nrow(params),1), legendnum=2, verbose=FALSE, plot_lines=FALSE, ... ){

  # Loop through parameter-pairs
  par( mfrow=mfrow, mar=c(3,3,0,0), mgp=c(1.75,0.25,0), tck=-0.02, oma=c(0,0,0,0))
  for( rowI in 1:nrow(params)){
    # Loop through specified taxa
    for( uniqueI in 1:length(unique(Taxa)) ){
      Plot_trait( Taxon=Taxa[uniqueI], params=params[rowI,], Cov_gjj=Cov_gjj, Mean_gj=Mean_gj, ParentChild_gz=ParentChild_gz, add=ifelse(uniqueI==1,FALSE,TRUE), xlim=range(Y_ij[,params[rowI,1]],na.rm=TRUE), ylim=range(Y_ij[,params[rowI,2]],na.rm=TRUE), partial_match=partial_match, main="", lcol=rainbow(length(Taxa))[uniqueI], ticks=ticks, logticks=logticks, plot_lines=plot_lines, verbose=verbose, ... )
      mtext( side=1:2, text=params[rowI,], line=1.5 )
    }
    # Add legend
    if( rowI==legendnum ){
      Legend = Taxa
      if( drop_pred==TRUE ) Legend = gsub(x=Taxa, pattern="_predictive", replacement="")
      legend("topleft", legend=Legend, fill=rainbow(length(Taxa)), bty="n")
    }
  }
}
