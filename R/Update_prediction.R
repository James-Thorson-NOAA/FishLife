#' Update predictions with new data
#'
#' Combines the taxonomic predictions from FishBase with new user-supplied data
#'
#' This function allows the user to update predictions using new data from one or more studies.  A study for a single parameter is entered as a new
#' row of \code{Ynew_ij}, where all elements are \code{NA} except for the parameter of interest.
#'
#' @param Ynew_ij a data frame of new trait values, perhaps log-scaled, with rows for records and tagged-columns for traits
#' @param obsCov_jj a matrix of the observation covariance for new data
#' @param partial_match Boolean, whether to require a partial match to a taxon-name \code{partial_match=TRUE} , or an exact match \code{partial_match=FALSE}
#' @param Taxon Boolean, whether to require a partial match to a taxon-name \code{partial_match=TRUE} , or an exact match \code{partial_match=FALSE}
#' @inheritParams Plot_ellipse
#' @inheritParams Fit_model
#' @inheritParams Search_species
#' @inheritParams Calculate_ratio
#'
#' @return Tagged list containing old and new distribution for the \code{Taxon}
#' \describe{
#'   \item{predMean_j}{Predictive mean using data from FishBase}
#'   \item{predCov_jj}{Predictive covariance using data from FishBase}
#'   \item{updateMean_j}{Updated mean after assimilating new information in \code{Ynew_ij}}
#'   \item{updateCov_jj}{Updated covariance after assimilating new information in \code{Ynew_ij}}
#' }
#'
#' @export
Update_prediction = function( Taxon, Ynew_ij, partial_match=TRUE, verbose=FALSE, ParentChild_gz=FishLife::database$ParentChild_gz,
  Cov_gjj=FishLife::database$Cov_gjj, Mean_gj=FishLife::database$ParHat$beta_gj, obsCov_jj=FishLife::database$obsCov_jj,
  Version="Update_v1_0_0", TmbDir=system.file("executables",package="FishLife"), RunDir=tempfile(pattern="run_",tmpdir=tempdir(),fileext="/")){

  # Match taxon
  if(partial_match==TRUE) Which = grep(Taxon, ParentChild_gz[,'ChildName'])
  if(partial_match==FALSE) Which = which(Taxon == ParentChild_gz[,'ChildName'])
  if( length(Which)!=1 ) stop( paste0("'Taxon' ",Taxon," input matches more or less than one element") )
  if(verbose==TRUE) print( ParentChild_gz[Which,] )

  # Extract stuff
  predMean_j = Mean_gj[Which,]
  predCov_jj = Cov_gjj[Which,,]

  #####################
  # TMB inputs
  #####################

  # Figure out missingness
  Missing_az = matrix( nrow=0, ncol=2)
  for(jI in 1:ncol(Ynew_ij)){
    whichRow = which(is.na(Ynew_ij[,jI]))
    if( length(whichRow)>0 ){
      Missing_az = rbind(Missing_az, cbind(whichRow,jI))
    }
  }
  Y_a = rnorm(nrow(Missing_az))
  if( nrow(Missing_az)==0 ){
    Missing_az = rbind(Missing_az, cbind(1,1))
    Y_a = as.vector( Ynew_ij )
  }

  # Data
  if(Version%in%"Update_v1_0_0") Data = list("Ynew_ij"=as.matrix(Ynew_ij), "Missing_az"=Missing_az-1, "predMean_j"=predMean_j, "predCov_jj"=predCov_jj, "obsCov_jj"=obsCov_jj)

  # Parameters
  if(Version%in%"Update_v1_0_0") Params = list( "beta_j"=rep(0,length(predMean_j)), "Y_a"=Y_a )

  # Random
  if(Version%in%"Update_v1_0_0") Random = c("Y_a")

  # Map
  Map = list()
  if( all(!is.na(Ynew_ij)) ) Map[["Y_a"]] = factor(rep(NA,length(Params$Y_a)))

  #####################
  # Build and run
  #####################

  Test = require(TMB, quietly=TRUE)
  if( Test==FALSE ){
    stop("Please install package `TMB` from CRAN")
  }

  # Compile TMB software
  dir.create( RunDir )
  #dyn.unload( paste0(RunDir,"/",TMB::dynlib(Version)) )          #
  file.copy( from=paste0(TmbDir,"/",Version,".cpp"), to=paste0(RunDir,"/",Version,".cpp"), overwrite=FALSE)
  setwd( RunDir )
  compile( paste0(Version,".cpp") )

  # Build
  dyn.load( paste0(RunDir,"/",TMB::dynlib(Version)) )          #
  Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map, DLL=Version, silent=TRUE )
  Report = Obj$report()

  # Optimize                         #  , startpar=opt$par[-grep("alpha",names(opt$par))]
  Opt = TMBhelper::Optimize( obj=Obj, savedir=RunDir, getJointPrecision=TRUE, control=list(eval.max=10000,iter.max=10000,trace=FALSE) )
  Report = Obj$report()

  # SE
  ParHat = Obj$env$parList( Opt$par )
  ParHat_SE = as.list( Opt$SD, what="Std" )

  # Updates
  updateMean_j = Report$beta_j
  updateCov_jj = Opt$SD$cov.fixed
  names(updateMean_j) = rownames(updateCov_jj) = colnames(updateCov_jj) = names(predMean_j)

  # Return match
  Return = list( "predMean_j"=predMean_j, "predCov_jj"=predCov_jj, "updateMean_j"=updateMean_j, "updateCov_jj"=updateCov_jj )
  return( Return )
}


