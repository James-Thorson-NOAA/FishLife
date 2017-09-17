#' Fit and predict fish traits
#'
#' \code{Fit_model} estimates parameters and predicts values from a multivariate random-walk model for fish traits
#'
#' @param N_factors Number of factors in decomposition of covariance for random-walk along evolutionary tree (0 means a diagonal but unequal covariance; negative is the sum of a factor decomposition and a diagonal-but-unequal covariance)
#' @param N_obsfactors Number of factors in decomposotion of observation covariance (same format as \code{N_obsfactors})
#' @param Use_REML, OPTIONAL boolean whether to use maximum marginal likelihood or restricted maximum likelihood (termed "REML")
#' @param Y_ij a data frame of trait values (perhaps log-scaled) with rows for records, and tagged-columns for traits
#' @param Z_ik a data frame of taxonomic classification for each row of \code{Y_ij}
#' @param Version TMB version number
#' @param Process_cov Whether process-error covariance is equal or differs multiplicatively for different taxonomic levels (Options:  "Equal" or "Unequal")
#' @param TmbDir Directory containing pre-compiled TMB script
#' @param RunDir Directory to use when compiling and running TMB script (different to avoid problems with read-write restrictions)
#' @param Params optional list of parameter estimates to use as starting values (Default \code{Params="Generate"} starts from random values)
#' @param verbose Boolean whether to print diagnostics to terminal
#' @param ... other paramers passed to \code{TMBhelper::Optimize}
#'
#' @return Tagged list containing objects from FishLife run (first 9 slots constitute list 'Estimate_database' for archiving results)
#' \describe{
#'   \item{N_factors}{Number of factors used for evolution in life-history model}
#'   \item{N_obsfactors}{Number of factors used for measurent-error in life-history model}
#'   \item{Use_REML}{Boolean, whether REML was used for model}
#'   \item{Cov_gjj}{Covariance among traits for every taxon in tree}
#'   \item{ParentChild_gz}{Record of taxonomic tree}
#'   \item{ParHat}{Parameter estimates and predictions}
#'   \item{g_i}{Associates every observation with a level of the taxonomic tree}
#'   \item{Y_ij}{Raw data}
#'   \item{Z_ik}{Taxonomy for each datum}
#'   \item{Obj}{The built TMB object}
#'   \item{Opt}{Output from optimization}
#'   \item{Report}{tagged list of report-file from TMB}
#'   \item{ParHat_SE}{Estimated/predicted standard errors for fixed/random effects}
#' }
#'
#' @export
Fit_model = function( N_factors, N_obsfactors, Use_REML=TRUE, Y_ij=FishLife::database$Y_ij, Z_ik=FishLife::database$Z_ik,
  Version="Taxon_v1_2_0", Process_cov="Equal", TmbDir=system.file("executables",package="FishLife"),
  RunDir=tempfile(pattern="run_",tmpdir=tempdir(),fileext="/"), Params="Generate", verbose=FALSE, ... ){

  #####################
  # Check for potential problems
  #####################

  # Code uses "_" so throw error if included in names
  if( any(apply(Z_ik, MARGIN=2, FUN=function(vec){length(grep("_",vec))})>0) ){
    stop("Please do not use character '_' in names included in 'Z_ik'")
  }

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
  if(Version%in%c("Taxon_v1_2_0","Taxon_v1_1_0")) Data = list("Options_vec"=c("n_obsfactors"=N_obsfactors,"n_factors"=N_factors), "Y_ij"=as.matrix(Y_ij), "Missing_az"=Missing_az-1, "PC_gz"=as.matrix(ParentChild_gz[,c('ParentRowNumber','ChildTaxon')])-1, "g_i"=g_i-1)

  # Parameters
  if( Params[1]=="Generate" ){
    rmatrix = function( nrow, ncol, mean=0, sd=1 ) matrix( rnorm(nrow*ncol,mean=mean,sd=sd), nrow=nrow, ncol=ncol )
    rloadings = function( n_row, n_col ){
      param = vector(length=0)
      if( is.na(n_col) ){
        param = c(param, rep(0.0001,n_row))
      }else{
        if( n_col!=0 ) param = c(param, rnorm(sum(n_row:(n_row-abs(n_col)+1))))
        if( n_col<=0 ) param = c(param, rnorm(n_row))
        #if( n_row == -1*n_col ) stop("Illogical inputs")
      }
      return(param)
    }
    if(Version%in%"Taxon_v1_0_0") Params = list( "alpha_j"=rep(0,n_j), "obsL_z"=rloadings(n_row=n_j, n_col=Data$Options_vec['n_obsfactors']), "Y_a"=rnorm(nrow(Data$Missing_az)) )
    if(Version%in%"Taxon_v1_1_0") Params = list( "alpha_j"=rep(0,n_j), "L_z"=rloadings(n_row=n_j, n_col=Data$Options_vec['n_factors']), "obsL_z"=rloadings(n_row=n_j, n_col=Data$Options_vec['n_obsfactors']), "beta_gj"=rmatrix(nrow=n_g,ncol=n_j), "Y_a"=rnorm(nrow(Data$Missing_az)) )
    if(Version%in%c("Taxon_v1_2_0")) Params = list( "alpha_j"=rep(0,n_j), "L_z"=rloadings(n_row=n_j, n_col=Data$Options_vec['n_factors']), "obsL_z"=rloadings(n_row=n_j, n_col=Data$Options_vec['n_obsfactors']), "cov_logmult_z"=rep(0,max(Data$PC_gz[,'ChildTaxon'])+1), "beta_gj"=rmatrix(nrow=n_g,ncol=n_j), "Y_a"=rnorm(nrow(Data$Missing_az)) )
  }

  # Random
  if(Version%in%"Taxon_v1_0_0") Random = c("Y_a")
  if(Version%in%c("Taxon_v1_2_0","Taxon_v1_1_0")) Random = c("Y_a", "beta_gj")
  if(Use_REML==TRUE) Random = c(Random, "alpha_j")

  # Map
  Map = list()
  if(Version%in%c("Taxon_v1_2_0")){
    if(Process_cov=="Equal") Map[["cov_logmult_z"]] = factor( rep(NA,length(Params[["cov_logmult_z"]])) )
    if(Process_cov=="Unequal") Map[["cov_logmult_z"]] = factor( c(NA,1:(length(Params[["cov_logmult_z"]])-1)) )
  }
  if(is.na(Data$Options_vec['n_obsfactors'])){
    Map[["obsL_z"]] = factor( rep(NA,length(Params[["obsL_z"]])) )
  }

  # Change
  if( is.na(Data$Options_vec[1]) ){
    Data$Options_vec[1] = 0
  }

  #####################
  # Build and run
  #####################

  Test = require(TMB, quietly=TRUE)
  if( Test==FALSE ){
    stop("Please install package `TMB` from CRAN")
  }

  # Compile TMB software
  #dyn.unload( paste0(RunDir,"/",dynlib(TMB:::getUserDLL())) ) # random=Random,
  dir.create( RunDir )
  file.copy( from=paste0(TmbDir,"/",Version,".cpp"), to=paste0(RunDir,"/",Version,".cpp"), overwrite=FALSE)
  origwd = getwd()
  on.exit(setwd(origwd),add=TRUE)
  setwd( RunDir )
  compile( paste0(Version,".cpp") )

  # Build
  dyn.load( paste0(RunDir,"/",TMB::dynlib(Version)) )          #
  Obj = MakeADFun( data=Data, parameters=Params, random=Random, DLL=Version, map=Map )
  Report = Obj$report()

  # Print to screen
  if( verbose==TRUE ){
    print( "Numer of fixed effects:")
    print( table(names(Obj$par)) )
    print( "Numer of random effects:")
    print( table(names(Obj$env$last.par[Obj$env$random])) )
  }

  # Optimize                         #  , startpar=opt$par[-grep("alpha",names(opt$par))]
  Opt = TMBhelper::Optimize( obj=Obj, savedir=RunDir, getJointPrecision=TRUE, ... ) # jointPrecision is used below, and is too big to invert whole
  Report = Obj$report()

  # SE
  ParHat = Obj$env$parList()
  ParHat_SE = as.list( Opt$SD, what="Std" )
  colnames(ParHat$beta_gj) = colnames(ParHat_SE$beta_gj) = colnames(Y_ij)
  dyn.unload( paste0(RunDir,"/",TMB::dynlib(Version)) )          #

  ####################
  # Interpret results
  ####################

  ### Approximate joint precision
  # Extract predictive covariance for species-specific traits (necessary to do separately for rows and columns)
  Prec_zz = Opt$SD$jointPrecision[ , grep("beta_gj",colnames(Opt$SD$jointPrecision)) ]
    Prec_zz = Prec_zz[ grep("beta_gj",rownames(Opt$SD$jointPrecision)), ]
  # Extract and invert
  PartialCorr_gjj = Corr_gjj = Cov_gjj = Prec_gjj = array(NA, dim=c(n_g,n_j,n_j), dimnames=list(ParentChild_gz[,'ChildName'],colnames(Y_ij),colnames(Y_ij)) )
  for( gI in 1:n_g ){
    # Extract precision for species and ancestors
    Indices = as.vector( outer(seq(1,n_g*n_j,by=n_g)-1, Find_ancestors(child_num=gI, ParentChild_gz=ParentChild_gz), FUN="+") )
    Full_Precision = matrix(Prec_zz[Indices,Indices],length(Indices),length(Indices))
    # Record
    Prec_gjj[gI,,] = Full_Precision[1:n_j,1:n_j]
    PartialCorr_gjj[gI,,] = -1*cov2cor( Prec_gjj[gI,,] )
    # Invert approximate cov and corr
    Cov_gjj[gI,,] = solve( Full_Precision )[1:n_j,1:n_j]
    Corr_gjj[gI,,] = cov2cor( Cov_gjj[gI,,] )
  }

  # Return stuff
  Return = list("N_factors"=N_factors, "N_obsfactors"=N_obsfactors, "Use_REML"=Use_REML, "Cov_gjj"=Cov_gjj, "ParentChild_gz"=ParentChild_gz, "ParHat"=ParHat, "g_i"=g_i, "Y_ij"=Y_ij, "Z_ik"=Z_ik, "Obj"=Obj, "Opt"=Opt, "Report"=Report, "ParHat_SE"=ParHat_SE, "obsCov_jj"=Report$obsCov_jj)
  dimnames(Return$obsCov_jj) = list(colnames(Y_ij),colnames(Y_ij))
  if(Version %in% c("Taxon_v1_1_0","Taxon_v1_0_0")) Return = c(Return, list("obsCov_jj"=Report$obsCov_jj, "Cov_jj"=Report$Cov_jj))
  if(Version %in% c("Taxon_v1_2_0")){
    Return = c(Return, list("obsCov_jj"=Report$obsCov_jj, "Cov_jjz"=Report$Cov_jj %o% exp(ParHat$cov_logmult_z)))
    dimnames(Return$Cov_jjz) = list(colnames(Y_ij),colnames(Y_ij),colnames(Z_ik))
  }

  return( Return )
}
