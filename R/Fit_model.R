
#' Fit and predict fish traits
#'
#' \code{Fit_model} estimates parameters and predicts values from a multivariate random-walk model for fish traits
#'
#' @param N_factors Number of factors in decomposotion of covariance for random-walk along evolutionary tree (same format as \code{N_obsfactors})
#' @param N_obsfactors Number of factors in decomposition of observation covariance (0 means a diagonal but unequal covariance; negative is the sum of a factor decomposition and a diagonal-but-unequal covariance)
#' @param Use_REML, OPTIONAL boolean whether to use maximum marginal likelihood or restricted maximum likelihood (termed "REML")
#' @param Y_ij a data frame of trait values (perhaps log-scaled) with rows for records, and tagged-columns for traits
#' @param Z_ik a data frame of taxonomic classification for each row of \code{Y_ij}
#' @param Version TMB version number
#' @param TmbDir Directory containing pre-compiled TMB script
#' @param RunDir Directory to use when compiling and running TMB script (different to avoid problems with read-write restrictions)

#' @return Tagged list containing objects from FishTraits run
#' \describe{
#'   \item{Obj}{The built TMB object}
#'   \item{Opt}{Output from optimization}
#'   \item{Report}{tagged list of report-file from TMB}
#'   \item{ParHat_SE}{Estimated/predicted standard errors for fixed/random effects}
#'   \item{Estimate_database}{Tagged list containing database of parameter estimates in format for other functions}
#' }

#' @export
Fit_model = function( N_factors, N_obsfactors, Use_REML=TRUE, Y_ij=Estimate_database$Y_ij, Z_ik=Estimate_database$Z_ik,
  Version="Taxon_v1_1_0", TmbDir=system.file("executables",package="FishTraits"), RunDir=getwd() ){

  #####################
  # Pre-process data
  #####################

  # Figure out missingness
  Missing_az = NULL
  for(jI in 1:ncol(Y_ij)){
    Which = which(is.na(Y_ij[,jI]))
    if( length(Which)>0 ){
      Missing_az = rbind(Missing_az, cbind(Which,jI))
    }
  }

  # Figure out network structure for taxonomy
  ParentChild_gz = NULL
  # 1st column: child taxon name
  # 2nd column: parent taxon name
  # 3rd column: parent row-number in ParentChild_gz
  # 4th column: Taxon level
  # Loop through
  for( colI in 1:ncol(Z_ik)){
    Taxa_Names = apply( Z_ik[,1:colI,drop=FALSE], MARGIN=1, FUN=paste, collapse="_")
    Unique_Taxa = unique(Taxa_Names)
    for( uniqueI in 1:length(Unique_Taxa) ){
      Which = which( Taxa_Names == Unique_Taxa[uniqueI] )
      if( colI==1 ){
        ParentChild_gz = rbind( ParentChild_gz, c(Unique_Taxa[uniqueI], NA, NA, colI) )
      }else{
        if( length(unique(Z_ik[Which,colI-1]))>1 ) stop("Taxa has multiple parents")
        ChildName = Unique_Taxa[uniqueI]
        ParentName = paste(rev(rev(strsplit(ChildName,"_")[[1]])[-1]),collapse="_")
        ParentChild_gz = rbind( ParentChild_gz, c(ChildName, ParentName, match(ParentName,ParentChild_gz[,1]), colI) )
      }
    }
  }
  # Loop through again to add predictive elements
  for( colI in 1:(ncol(Z_ik)-1)){
    Taxa_Names = apply( Z_ik[,1:colI,drop=FALSE], MARGIN=1, FUN=paste, collapse="_")
    Unique_Taxa = unique(Taxa_Names)
    for( uniqueI in 1:length(Unique_Taxa) ){
      ParentName = Unique_Taxa[uniqueI]
      for( predI in 1:(ncol(Z_ik)-colI) ){
        ChildName = paste0(ParentName,"_predictive")
        ParentChild_gz = rbind( ParentChild_gz, c(ChildName, ParentName, match(ParentName,ParentChild_gz[,1]), colI+predI) )
        ParentName = ChildName
      }
    }
  }
  # Add top predictive
  ParentChild_gz = rbind( ParentChild_gz, c("predictive", NA, NA, 1) )
  for( colI in 2:ncol(Z_ik)) ParentChild_gz = rbind( ParentChild_gz, c(paste(rep("predictive",colI),collapse="_"), paste(rep("predictive",colI-1),collapse="_"), match(paste(rep("predictive",colI-1),collapse="_"),ParentChild_gz[,1]), colI) )
  # Relabel
  ParentChild_gz = data.frame( ParentChild_gz )
  colnames(ParentChild_gz) = c("ChildName", "ParentName", "ParentRowNumber", "ChildTaxon")
  ParentChild_gz[,'ParentRowNumber'] = as.numeric(as.character(ParentChild_gz[,'ParentRowNumber']))
  ParentChild_gz[,'ChildTaxon'] = as.numeric(as.character(ParentChild_gz[,'ChildTaxon']))

  # Identify location for every observation
  Taxa_Names = apply( Z_ik, MARGIN=1, FUN=paste, collapse="_")
  g_i = match( Taxa_Names, ParentChild_gz[,'ChildName'] )
  n_j = ncol(Y_ij)
  n_k = ncol(Z_ik)
  n_g = nrow(ParentChild_gz)

  #####################
  # TMB inputs
  #####################

  # Data
  if(Version%in%"Taxon_v1_0_0") Data = list("Options_vec"=c("n_obsfactors"=N_obsfactors), "Z_ik"=as.matrix(Z_ik)-1, "Y_ij"=as.matrix(Y_ij), "Missing_az"=Missing_az-1)
  if(Version%in%"Taxon_v1_1_0") Data = list("Options_vec"=c("n_obsfactors"=N_obsfactors,"n_factors"=N_factors), "Y_ij"=as.matrix(Y_ij), "Missing_az"=Missing_az-1, "PC_gz"=as.matrix(ParentChild_gz[,c('ParentRowNumber','ChildTaxon')])-1, "g_i"=g_i-1)

  # Parameters
  rmatrix = function( nrow, ncol, mean=0, sd=1 ) matrix( rnorm(nrow*ncol,mean=mean,sd=sd), nrow=nrow, ncol=ncol )
  rloadings = function( n_row, n_col ){
    param = vector(length=0)
    if( n_col!=0 ) param = c(param, rnorm(sum(n_row:(n_row-abs(n_col)+1))))
    if( n_col<=0 ) param = c(param, rnorm(n_row))
    #if( n_row == -1*n_col ) stop("Illogical inputs")
    return(param)
  }
  if(Version%in%"Taxon_v1_0_0") Params = list( "alpha_j"=rep(0,n_j), "obsL_z"=rloadings(n_row=n_j, n_col=Data$Options_vec['n_obsfactors']), "Y_a"=rnorm(nrow(Data$Missing_az)) )
  if(Version%in%"Taxon_v1_1_0") Params = list( "alpha_j"=rep(0,n_j), "L_z"=rloadings(n_row=n_j, n_col=Data$Options_vec['n_factors']), "obsL_z"=rloadings(n_row=n_j, n_col=Data$Options_vec['n_obsfactors']), "beta_gj"=rmatrix(nrow=n_g,ncol=n_j), "Y_a"=rnorm(nrow(Data$Missing_az)) )

  # Random
  if(Version%in%"Taxon_v1_0_0") Random = c("Y_a")
  if(Version%in%"Taxon_v1_1_0") Random = c("Y_a", "beta_gj")
  if(Use_REML==TRUE) Random = c(Random, "alpha_j")

  #####################
  # Build and run
  #####################

  # Compile TMB software
  #dyn.unload( paste0(RunDir,"/",dynlib(TMB:::getUserDLL())) ) # random=Random,
  file.copy( from=paste0(TmbDir,"/",Version,".cpp"), to=paste0(RunDir,"/",Version,".cpp"), overwrite=FALSE)
  setwd( RunDir )
  compile( paste0(Version,".cpp") )

  # Build
  dyn.load( paste0(RunDir,"/",TMB::dynlib(Version)) )          #
  Obj = MakeADFun( data=Data, parameters=Params, random=Random, DLL=Version )
  Report = Obj$report()

  # Optimize                         #  , startpar=opt$par[-grep("alpha",names(opt$par))]
  Opt = TMBhelper::Optimize( obj=Obj, savedir=RunDir, getJointPrecision=TRUE ) # jointPrecision is used below, and is too big to invert whole
  Report = Obj$report()

  # SE
  ParHat = Obj$env$parList()
  ParHat_SE = as.list( Opt$SD, what="Std" )
  colnames(ParHat$beta_gj) = colnames(ParHat_SE$beta_gj) = colnames(Y_ij)
  dyn.load( paste0(RunDir,"/",TMB::dynlib(Version)) )          #

  ####################
  # Interpret results
  ####################
  # Approximate joint precision
  Prec_zz = Opt$SD$jointPrecision[ grep("beta_gj",names(unlist(Params))), ]
  Prec_zz = Prec_zz[ , grep("beta_gj",names(unlist(Params))) ]
  PartialCorr_gjj = Corr_gjj = Cov_gjj = Prec_gjj = array(NA, dim=c(n_g,n_j,n_j), dimnames=list(ParentChild_gz[,'ChildName'],colnames(Y_ij),colnames(Y_ij)) )
  for( gI in 1:n_g ){
    Indices = as.vector( outer(seq(1,n_g*8,by=n_g)-1, Find_ancestors(child_num=gI, ParentChild_gz=ParentChild_gz), FUN="+") )
    Full_Precision = matrix(Prec_zz[Indices,Indices],length(Indices),length(Indices))
    Prec_gjj[gI,,] = Full_Precision[1:n_j,1:n_j]
    PartialCorr_gjj[gI,,] = -1*cov2cor( Prec_gjj[gI,,] )
    Cov_gjj[gI,,] = solve( Full_Precision )[1:n_j,1:n_j]
    Corr_gjj[gI,,] = cov2cor( Cov_gjj[gI,,] )
  }

  # Return stuff
  Estimate_database = list("N_factors"=N_factors, "N_obsfactors"=N_obsfactors, "Use_REML"=Use_REML, "Cov_gjj"=Cov_gjj, "ParentChild_gz"=ParentChild_gz, "ParHat"=ParHat, "g_i"=g_i, "Y_ij"=Y_ij, "Z_ik"=Z_ik)
  Return = list("Obj"=Obj, "Opt"=Opt, "Report"=Report, "ParHat_SE"=ParHat_SE, "Estimate_database"=Estimate_database)
  return( Return )
}
