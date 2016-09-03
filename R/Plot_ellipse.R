
#' @export
Plot_ellipse = function( Taxon, SpeciesMatch=NULL, prob=0.95, params=c('K','M'), add=FALSE, xlim=log(c(0.01,2)), ylim=xlim, ticks=c(1,2,5), partial_match=TRUE, main="", xlab="", ylab="", lcol="black", verbose=TRUE, plot_lines=FALSE, ... ){
  # Match taxon
  if(partial_match==TRUE) Which = grep(Taxon, ParentChild_gz[,'ChildName'])
  if(partial_match==FALSE) Which = which(Taxon == ParentChild_gz[,'ChildName'])
  if( length(Which)!=1 ) stop( paste0("'Taxon' ",Taxon," input matches more or less than one element") )
  if(verbose==TRUE) print( ParentChild_gz[Which,] )

  # Calculate eigen-decomposition
  # http://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/
  Eigen = eigen( Cov_gjj[Which,params,params] )
  Major_radius = sqrt(qchisq(prob,df=2) * Eigen$values[1])
  Minor_radius = sqrt(qchisq(prob,df=2) * Eigen$values[2])

  # Calculate bounds
  if( is.null(xlim) ){
    xlim = ParHat$beta_gj[Which,params][1] + c(-1.2,1.2)*max(c(Major_radius,Minor_radius)*abs(Eigen$vectors[1,]))
    ylim = ParHat$beta_gj[Which,params][2] + c(-1.2,1.2)*max(c(Major_radius,Minor_radius)*abs(Eigen$vectors[2,]))
  }

  # Plot log-scale ellipse
  if(add==FALSE){
    f = list( function(a){a}, log )
    plot( 1, type="n", xlim=xlim, ylim=ylim, xaxt="n", yaxt="n", main=main, xlab=xlab, ylab=ylab)    #
    axis( side=1, at=f[[ifelse(params[1]=="Temperature",1,2)]](as.vector(outer(ticks,10^(-10:10)))), labels=as.vector(outer(ticks,10^(-10:10))) )
    axis( side=2, at=f[[ifelse(params[2]=="Temperature",1,2)]](as.vector(outer(ticks,10^(-10:10)))), labels=as.vector(outer(ticks,10^(-10:10))) )
  }
  shape::plotellipse( rx=Major_radius, ry=Minor_radius, mid=ParHat$beta_gj[Which,params], angle=atan(Eigen$vector[2,1]/Eigen$vector[1,1])/(2*pi)*360, lcol=lcol, ...)
  points( x=ParHat$beta_gj[Which,params[1]], y=ParHat$beta_gj[Which,params[2]], col=lcol, pch=20 )

  # Plot observations
  if( !is.null(SpeciesMatch) ){
    Which = grep(SpeciesMatch, ParentChild_gz[,'ChildName'])
    Which = which( g_i %in% Which )
    points( x=Y_ij[Which,'K'], y=Y_ij[Which,'M'] )
  }

  # Plot lines for OLS and MA regression
  calc_coef = function( Mean_z, Cov_zz, Type="Y|X" ){
    # http://www.unc.edu/courses/2007spring/biol/145/001/docs/lectures/Nov5.html
    Coef = matrix( NA, nrow=3, ncol=3, dimnames=list(c("Y|X","major_axis","X|Y"),c("plot_intercept","slope","param_intercept")) )
    Coef["Y|X",'slope'] = cov2cor(Cov_gjj[Which,params,params])[1,2] * sqrt(Cov_zz[2,2]/Cov_zz[1,1])
    Coef["major_axis",'slope'] = Eigen$vectors[2,1] / Eigen$vectors[1,1]
    Coef["X|Y",'slope'] = 1/cov2cor(Cov_gjj[Which,params,params])[1,2] * sqrt(Cov_zz[2,2]/Cov_zz[1,1])
    Coef[,'plot_intercept'] = Mean_z[2] - Coef[,'slope']*Mean_z[1]
    Coef[,'param_intercept'] = exp(Mean_z[2] - Mean_z[1])
    return( Coef )
  }
  if( plot_lines==TRUE ){
    Coef = calc_coef(Mean_z=ParHat$beta_gj[Which,params], Cov_zz=Cov_gjj[Which,params,params])
    abline( coef=Coef["Y|X",c("plot_intercept","slope")], col="red" )
    abline( coef=Coef["major_axis",c("plot_intercept","slope")], col="black" )
    abline( coef=Coef["X|Y",c("plot_intercept","slope")], col="blue" )
    if( !("Temperature" %in% params) ){
      legend( "topleft", fill=c("red","black","blue"), legend=formatC(Coef[,'slope'],format="f",digits=2), bty="n", title="Slope" )
      legend( "bottomright", fill=c("black"), legend=formatC(Coef["major_axis",'param_intercept'],format="f",digits=2), bty="n", title="Intercept" )
    }
  }

  # insibile return
  return( invisible(list("Cov_pred"=Cov_gjj[Which,,])) )
}
