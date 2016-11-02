#' Plot an ellipse
#'
#' Plots an ellipse representing the bivariate predictive inteval
#'
#' @param Cov 2x2 matrix representing covariance
#' @param Mean numeric vector (length of 2) representing Empirical Bayes prediction of median
#' @param add Boolean whether to add ellipse to existing plot
#' @param probability for defining predictive interval
#' @param xlim x-limits if \code{add=FALSE}
#' @param ylim y-limits if \code{add=FALSE}
#' @param logticks ticks to use if plotting on a log-scale
#' @param ticks ticks to use if plotting on a natural-scale
#' @param whichlog which axes are log-scale
#' @param main legend for each panel
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param lcol line color for ellipse
#' @param plot_lines whether to plot lines representing "major axis" and "OLS" regression parameters
#'
#' @return integer of row numbers of \code{ParentChild_gz} matching \code{genus_species}
#'
#' @export
Plot_ellipse = function( Cov, Mean=rep(0,2), add=FALSE, prob=0.95, xlim=log(c(0.01,2)), ylim=xlim, logticks=c(1,2,5), ticks=c(0,5),
  whichlog="xy", main="", xlab="", ylab="", lcol="black", plot_lines=FALSE, ... ){

  # Calculate eigen-decomposition
  # http://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/
  Eigen = eigen( Cov )
  Major_radius = sqrt(qchisq(prob,df=2) * Eigen$values[1])
  Minor_radius = sqrt(qchisq(prob,df=2) * Eigen$values[2])
  Angle = atan(Eigen$vector[2,1]/Eigen$vector[1,1])/(2*pi)*360

  # Plot log-scale ellipse
  if(add==FALSE){
  # Calculate bounds
    if( is.null(xlim) ){
      xlim = Mean[1] + c(-1.2,1.2)*max(c(Major_radius,Minor_radius)*abs(Eigen$vectors[1,]))
      ylim = Mean[2] + c(-1.2,1.2)*max(c(Major_radius,Minor_radius)*abs(Eigen$vectors[2,]))
    }
    plot( 1, type="n", xlim=xlim, ylim=ylim, xaxt="n", yaxt="n", main=main, xlab=xlab, ylab=ylab)    #
    Logticks = as.vector(outer(logticks,10^(-10:10),FUN="*"))
    Ticks = as.vector(outer(ticks,10*(-10:10),FUN="+"))
    if(length(grep("x",whichlog))==1){ axis(side=1, at=log(Logticks), labels=Logticks) }else{ axis(side=1, at=Ticks, labels=Ticks) }
    if(length(grep("y",whichlog))==1){ axis(side=2, at=log(Logticks), labels=Logticks) }else{ axis(side=2, at=Ticks, labels=Ticks) }
  }
  shape::plotellipse( rx=Major_radius, ry=Minor_radius, mid=Mean, angle=Angle, lcol=lcol, ...)

  # Plot lines for OLS and MA regression
  calc_coef = function( Mean, Cov, Type="Y|X" ){
    # http://www.unc.edu/courses/2007spring/biol/145/001/docs/lectures/Nov5.html
    Coef = matrix( NA, nrow=3, ncol=3, dimnames=list(c("Y|X","major_axis","X|Y"),c("plot_intercept","slope","param_intercept")) )
    Coef["Y|X",'slope'] = cov2cor(Cov)[1,2] * sqrt(Cov[2,2]/Cov[1,1])
    Coef["major_axis",'slope'] = Eigen$vectors[2,1] / Eigen$vectors[1,1]
    Coef["X|Y",'slope'] = 1/cov2cor(Cov)[1,2] * sqrt(Cov[2,2]/Cov[1,1])
    Coef[,'plot_intercept'] = Mean[2] - Coef[,'slope']*Mean[1]
    Coef[,'param_intercept'] = exp(Mean[2] - Mean[1])
    return( Coef )
  }
  if( plot_lines==TRUE ){
    Coef = calc_coef(Mean=Mean, Cov=Cov)
    abline( coef=Coef["Y|X",c("plot_intercept","slope")], col="red" )
    abline( coef=Coef["major_axis",c("plot_intercept","slope")], col="black" )
    abline( coef=Coef["X|Y",c("plot_intercept","slope")], col="blue" )
    if( whichlog=="xy" ){
      legend( "topleft", fill=c("red","black","blue"), legend=formatC(Coef[,'slope'],format="f",digits=2), bty="n", title="Slope" )
      legend( "bottomright", fill=c("black"), legend=formatC(Coef["major_axis",'param_intercept'],format="f",digits=2), bty="n", title="Intercept" )
    }
  }
  return( invisible(c("Major_radius"=Major_radius,"Minor_radius"=Minor_radius,"Angle"=Angle)) )
}
