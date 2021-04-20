#' Fit and predict fish traits
#'
#' \code{Fit_model} estimates parameters and predicts values from a multivariate random-walk model for fish traits
#'
#' @param N_factors Number of factors in decomposition of covariance for random-walk along evolutionary tree (0 means a diagonal but unequal covariance; negative is the sum of a factor decomposition and a diagonal-but-unequal covariance)
#' @param N_obsfactors Number of factors in decomposotion of observation covariance (same format as \code{N_obsfactors})
#' @param Database, Whether to use results for both adult and stock-recruit parameters, \code{Database="FishBase_and_RAM"}, or just adult parameters, \code{Database="FishBase"}
#' @param Use_REML, OPTIONAL boolean whether to use maximum marginal likelihood or restricted maximum likelihood (termed "REML")
#' @param Y_ij a data frame of trait values (perhaps log-scaled) with rows for records, and tagged-columns for traits
#' @param Z_ik a data frame of taxonomic classification for each row of \code{Y_ij}
#' @param SR_obs Stock-recruit records from RAM Legacy stock-recruit database
#' @param StockData Auxiliary information for every stock with stock-recruit information
#' @param Version TMB version number
#' @param Process_cov Whether process-error covariance is equal or differs multiplicatively for different taxonomic levels (Options:  "Equal" or "Unequal")
#' @param TmbDir Directory containing pre-compiled TMB script
#' @param RunDir Directory to use when compiling and running TMB script (different to avoid problems with read-write restrictions)
#' @param Params optional list of parameter estimates to use as starting values (Default \code{Params="Generate"} starts from random values)
#' @param verbose Boolean whether to print diagnostics to terminal
#' @param run_model Boolean indicating whether to run the model or just return the built TMB object
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
Fit_model = function( N_factors,
          N_obsfactors,
          Use_REML = TRUE,
          Database = FishLife::FishBase_and_RAM,
          Y_ij = Database$Y_ij,
          Z_ik = Database$Z_ik,
          SR_obs = Database$SR_obs,
          StockData = Database$StockData,
          Version = "Taxon_v2_14_0",
          Process_cov = "Unequal",
          TmbDir = system.file("executables",package="FishLife"),
          RunDir = tempfile(pattern = "run_",tmpdir=tempdir(),fileext="/"),
          verbose = FALSE,
          debug_mode = FALSE,
          j_SR = ncol(Y_ij)-3:1,
          zerocovTF_j = rep(FALSE,ncol(Y_ij)),
          additional_variance = c(0,0),
          invertTF = FALSE,
          SD_b_stock = 10,
          b_type = 0,
          Cov_design = NULL,
          Cov_RAM = c("M"=FALSE),
          Turn_off_taxonomy = FALSE,
          Pen_lowvar_lnRhat = 1,
          lowerbound_MLSPS = 1,
          Use_RAM_Mvalue_TF = TRUE,
          rho_space = "natural",
          include_r = TRUE,
          PredTF_stock = NULL,
          extract_covariance = TRUE,
          run_model = TRUE,
          Params = "Generate",
          Random = "Generate",
          Map = "Generate",
          ... ){

  #### Local function
  # Sample from GMRF using sparse precision
  rmvnorm_prec <- function(mu, prec, n.sims, seed=1) {
    set.seed( seed )
    z <- matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
    L <- Matrix::Cholesky(prec, super=TRUE)
    z <- Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
    z <- Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
    z <- as.matrix(z)
    mu + z
  }

  #####################
  # Check for potential problems
  #####################
  start_time = Sys.time()

  # Code uses "_" so throw error if included in names
  if( any(apply(Z_ik, MARGIN=2, FUN=function(vec){length(grep("_",vec))})>0) ){
    stop("Please do not use character '_' in names included in 'Z_ik'")
  }

  # Process stock-recruit inputs
  if( is.null(SR_obs) | is.null(StockData) ){
    ### TO ADD LATER
    Nstock = 0
    Nobs = 0
    SR_obs = cbind("R_obs"=1, "SSB_obs"=1, "AR_Index"=1, "StockNum"=1)
    include_r = FALSE
    PredTF_stock = vector()
    StockData = cbind( "SPRF0"=vector(), "M"=vector(), "SSBmax"=vector(), "Rmax"=vector(), "Stock_to_i"=vector() )
  }else{
    if( !all( c("R_obs","SSB_obs","AR_Index","StockNum") %in% colnames(SR_obs)) ) stop("Check `SR_obs` input")
    if( !all( c("M","SPRF0","Stock_to_i") %in% colnames(StockData)) ) stop("Check `StockData` input")
    Nobs = nrow(SR_obs)
    Nstock = nrow(StockData)
    if( !("SSBmax" %in% colnames(StockData)) ){
      StockData = cbind(StockData, "SSBmax"=tapply(SR_obs[,'SSB_obs'], INDEX=SR_obs[,'StockNum'], FUN=max, na.rm=TRUE) )
    }
    if( !("Rmax" %in% colnames(StockData)) ){
      StockData = cbind(StockData, "Rmax"=tapply(SR_obs[,'R_obs'], INDEX=SR_obs[,'StockNum'], FUN=max, na.rm=TRUE) )
    }
    if( any(is.na(StockData)) ) stop("Check `StockData`")
    # Rename rho if necessary
    if( rho_space%in%c("logit","logit_with_jacobian") & ("rho" %in% colnames(Y_ij)) ){
      colnames(Y_ij) = sapply( 1:ncol(Y_ij), FUN=function(num){switch(colnames(Y_ij)[num], "rho"="logit_rho", colnames(Y_ij)[num])} )
    }
    if( is.null(PredTF_stock) ){
      PredTF_stock = rep(FALSE, Nstock)
    }
  }

  if( is.null(Cov_design) ){
    Cov_pz = matrix( c(0,0), ncol=2, nrow=1 )
    n_p = 0
  }else{
    Cov_pz = Cov_design
    n_p = nrow(Cov_pz)
    if( !is.matrix(Cov_pz) | ncol(Cov_pz)!=2 ) stop("Check `Cov_pz` type and dimensions")
    if( any(Cov_pz[,1]>=Cov_pz[,2]) ) stop("Check entries in `Cov_pz`")
    if( any(Cov_pz<0 | Cov_pz>(ncol(Y_ij)-1)) ) stop("Entries in `Cov_pz` exceed bounds")
  }

  # Check for errors
  if(nrow(Y_ij)!=nrow(Z_ik)) stop("Check number of rows in `Y_ij` and `Z_ik`")

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
  n_k = ncol(Z_ik)
  n_j = ncol(Y_ij)
  n_g = nrow(ParentChild_gz)

  # Context-specific inputs
  if("M" %in% colnames(Y_ij)){
    j_logM = which(colnames(Y_ij)=="M")-1
  }else{
    j_logM = 0
  }

  #####################
  # TMB inputs
  #####################

  # Data
  if(Version%in%"Taxon_v1_0_0") Data = list("Options_vec"=c("n_obsfactors"=N_obsfactors), "Z_ik"=as.matrix(Z_ik)-1, "Y_ij"=as.matrix(Y_ij), "Missing_az"=Missing_az-1)
  if(Version%in%c("Taxon_v1_2_0","Taxon_v1_1_0")) Data = list("Options_vec"=c("n_obsfactors"=N_obsfactors,"n_factors"=N_factors), "Y_ij"=as.matrix(Y_ij), "Missing_az"=Missing_az-1, "PC_gz"=as.matrix(ParentChild_gz[,c('ParentRowNumber','ChildTaxon')])-1, "g_i"=g_i-1)
  if(Version%in%c("Taxon_v2_2_0","Taxon_v2_1_0","Taxon_v2_0_0")){
    Data = list("Options_vec"=c("n_obsfactors"=N_obsfactors,"n_factors"=N_factors,"invertTF"=invertTF), "Y_ij"=as.matrix(Y_ij), "Missing_az"=Missing_az-1, "PC_gz"=as.matrix(ParentChild_gz[,c('ParentRowNumber','ChildTaxon')])-1, "g_i"=g_i-1)
    Data = c(Data, list("Nobs"=Nobs, "Nstock"=Nstock, "Obs2Stock"=SR_obs[,'StockNum']-1, "AR_Index"=SR_obs[,'AR_Index'], "ln_R_obs"=log(SR_obs[,'R_obs']), "SSB_obs"=SR_obs[,'SSB_obs'], "SPRF0_stock"=StockData[,'SPRF0'], "M_stock"=StockData[,'M'], "j_SR"=j_SR, "i_stock"=StockData[,'Stock_to_i']-1) )
  }
  if(Version%in%c("Taxon_v2_5_0","Taxon_v2_4_0","Taxon_v2_3_0")){
    Data = list("Options_vec"=c("n_obsfactors"=N_obsfactors,"n_factors"=N_factors,"invertTF"=invertTF), "Options"=c("minvar_obsfactors"=additional_variance[1],"minvar_factors"=additional_variance[2],"SD_b_stock"=SD_b_stock), "Y_ij"=as.matrix(Y_ij), "Missing_az"=Missing_az-1, "PC_gz"=as.matrix(ParentChild_gz[,c('ParentRowNumber','ChildTaxon')])-1, "g_i"=g_i-1)
    Data = c(Data, list("Nobs"=Nobs, "Nstock"=Nstock, "Obs2Stock"=SR_obs[,'StockNum']-1, "AR_Index"=SR_obs[,'AR_Index'], "ln_R_obs"=log(SR_obs[,'R_obs']), "SSB_obs"=SR_obs[,'SSB_obs'], "SPRF0_stock"=StockData[,'SPRF0'], "M_stock"=StockData[,'M'], "j_SR"=j_SR, "i_stock"=StockData[,'Stock_to_i']-1) )
  }
  if(Version%in%c("Taxon_v2_6_0")){
    Data = list("Options_vec"=c("n_obsfactors"=N_obsfactors,"n_factors"=N_factors,"invertTF"=invertTF,"b_type"=b_type), "Options"=c("minvar_obsfactors"=additional_variance[1],"minvar_factors"=additional_variance[2],"SD_b_stock"=SD_b_stock), "Y_ij"=as.matrix(Y_ij), "Missing_az"=Missing_az-1, "PC_gz"=as.matrix(ParentChild_gz[,c('ParentRowNumber','ChildTaxon')])-1, "g_i"=g_i-1)
    Data = c(Data, list("Nobs"=Nobs, "Nstock"=Nstock, "Obs2Stock"=SR_obs[,'StockNum']-1, "AR_Index"=SR_obs[,'AR_Index'], "ln_R_obs"=log(SR_obs[,'R_obs']), "SSB_obs"=SR_obs[,'SSB_obs'], "SPRF0_stock"=StockData[,'SPRF0'], "M_stock"=StockData[,'M'], "SSBmax_stock"=StockData[,'SSBmax'], "Rmax_stock"=StockData[,'Rmax'], "j_SR"=j_SR, "i_stock"=StockData[,'Stock_to_i']-1) )
  }
  if(Version%in%c("Taxon_v2_8_0","Taxon_v2_7_0")){
    Data = list("Options_vec"=c("n_obsfactors"=N_obsfactors,"n_factors"=N_factors,"invertTF"=invertTF,"b_type"=b_type), "Options"=c("minvar_obsfactors"=additional_variance[1],"minvar_factors"=additional_variance[2],"SD_b_stock"=SD_b_stock), "Cov_pz"=Cov_pz, "Y_ij"=as.matrix(Y_ij), "Missing_az"=Missing_az-1, "PC_gz"=as.matrix(ParentChild_gz[,c('ParentRowNumber','ChildTaxon')])-1, "g_i"=g_i-1)
    Data = c(Data, list("Nobs"=Nobs, "Nstock"=Nstock, "Obs2Stock"=SR_obs[,'StockNum']-1, "AR_Index"=SR_obs[,'AR_Index'], "ln_R_obs"=log(SR_obs[,'R_obs']), "SSB_obs"=SR_obs[,'SSB_obs'], "SPRF0_stock"=StockData[,'SPRF0'], "M_stock"=StockData[,'M'], "SSBmax_stock"=StockData[,'SSBmax'], "Rmax_stock"=StockData[,'Rmax'], "j_SR"=j_SR, "i_stock"=StockData[,'Stock_to_i']-1) )
  }
  if(Version%in%c("Taxon_v2_9_0")){
    Data = list("Options_vec"=c("n_obsfactors"=N_obsfactors,"n_factors"=N_factors,"invertTF"=invertTF,"b_type"=b_type,"Turn_off_taxonomy"=Turn_off_taxonomy), "Options"=c("minvar_obsfactors"=additional_variance[1],"minvar_factors"=additional_variance[2],"SD_b_stock"=SD_b_stock), "Cov_pz"=Cov_pz, "Y_ij"=as.matrix(Y_ij), "Missing_az"=Missing_az-1, "PC_gz"=as.matrix(ParentChild_gz[,c('ParentRowNumber','ChildTaxon')])-1, "g_i"=g_i-1)
    Data = c(Data, list("Nobs"=Nobs, "Nstock"=Nstock, "Obs2Stock"=SR_obs[,'StockNum']-1, "AR_Index"=SR_obs[,'AR_Index'], "ln_R_obs"=log(SR_obs[,'R_obs']), "SSB_obs"=SR_obs[,'SSB_obs'], "SPRF0_stock"=StockData[,'SPRF0'], "M_stock"=StockData[,'M'], "SSBmax_stock"=StockData[,'SSBmax'], "Rmax_stock"=StockData[,'Rmax'], "j_SR"=j_SR, "i_stock"=StockData[,'Stock_to_i']-1) )
  }
  if(Version%in%c("Taxon_v2_10_0")){
    Data = list("Options_vec"=c("n_obsfactors"=N_obsfactors,"n_factors"=N_factors,"invertTF"=invertTF,"b_type"=b_type,"Turn_off_taxonomy"=Turn_off_taxonomy), "Options"=c("minvar_obsfactors"=additional_variance[1],"minvar_factors"=additional_variance[2],"SD_b_stock"=SD_b_stock,"Pen_lowvar_lnRhat"=Pen_lowvar_lnRhat), "Cov_pz"=Cov_pz, "Y_ij"=as.matrix(Y_ij), "Missing_az"=Missing_az-1, "PC_gz"=as.matrix(ParentChild_gz[,c('ParentRowNumber','ChildTaxon')])-1, "g_i"=g_i-1)
    Data = c(Data, list("Nobs"=Nobs, "Nstock"=Nstock, "Obs2Stock"=SR_obs[,'StockNum']-1, "AR_Index"=SR_obs[,'AR_Index'], "ln_R_obs"=log(SR_obs[,'R_obs']), "SSB_obs"=SR_obs[,'SSB_obs'], "SPRF0_stock"=StockData[,'SPRF0'], "M_stock"=StockData[,'M'], "SSBmax_stock"=StockData[,'SSBmax'], "Rmax_stock"=StockData[,'Rmax'], "j_SR"=j_SR, "i_stock"=StockData[,'Stock_to_i']-1) )
  }
  if(Version%in%c("Taxon_v2_11_0")){
    Data = list("Options_vec"=c("n_obsfactors"=N_obsfactors,"n_factors"=N_factors,"invertTF"=invertTF,"b_type"=b_type,"Turn_off_taxonomy"=Turn_off_taxonomy), "Options"=c("minvar_obsfactors"=additional_variance[1],"minvar_factors"=additional_variance[2],"SD_b_stock"=SD_b_stock,"Pen_lowvar_lnRhat"=Pen_lowvar_lnRhat, "lowerbound_MLSPS"=lowerbound_MLSPS), "Cov_pz"=Cov_pz, "Y_ij"=as.matrix(Y_ij), "Missing_az"=Missing_az-1, "PC_gz"=as.matrix(ParentChild_gz[,c('ParentRowNumber','ChildTaxon')])-1, "g_i"=g_i-1)
    Data = c(Data, list("Nobs"=Nobs, "Nstock"=Nstock, "Obs2Stock"=SR_obs[,'StockNum']-1, "AR_Index"=SR_obs[,'AR_Index'], "ln_R_obs"=log(SR_obs[,'R_obs']), "SSB_obs"=SR_obs[,'SSB_obs'], "SPRF0_stock"=StockData[,'SPRF0'], "M_stock"=StockData[,'M'], "SSBmax_stock"=StockData[,'SSBmax'], "Rmax_stock"=StockData[,'Rmax'], "j_SR"=j_SR, "i_stock"=StockData[,'Stock_to_i']-1) )
  }
  if(Version%in%c("Taxon_v2_12_0")){
    Data = list("Options_vec"=c("n_obsfactors"=N_obsfactors,"n_factors"=N_factors,"invertTF"=invertTF,"b_type"=b_type,"Turn_off_taxonomy"=Turn_off_taxonomy,"Use_RAM_Mvalue_TF"=Use_RAM_Mvalue_TF), "Options"=c("minvar_obsfactors"=additional_variance[1],"minvar_factors"=additional_variance[2],"SD_b_stock"=SD_b_stock,"Pen_lowvar_lnRhat"=Pen_lowvar_lnRhat, "lowerbound_MLSPS"=lowerbound_MLSPS), "Cov_pz"=Cov_pz, "Y_ij"=as.matrix(Y_ij), "Missing_az"=Missing_az-1, "PC_gz"=as.matrix(ParentChild_gz[,c('ParentRowNumber','ChildTaxon')])-1, "g_i"=g_i-1)
    Data = c(Data, list("Nobs"=Nobs, "Nstock"=Nstock, "Obs2Stock"=SR_obs[,'StockNum']-1, "AR_Index"=SR_obs[,'AR_Index'], "ln_R_obs"=log(SR_obs[,'R_obs']), "SSB_obs"=SR_obs[,'SSB_obs'], "SPRF0_stock"=StockData[,'SPRF0'], "M_stock"=StockData[,'M'], "SSBmax_stock"=StockData[,'SSBmax'], "Rmax_stock"=StockData[,'Rmax'], "j_SR"=j_SR, "j_logM"=j_logM, "i_stock"=StockData[,'Stock_to_i']-1) )
  }
  if(Version%in%c("Taxon_v2_13_0")){
    Data = list("Options_vec"=c("n_obsfactors"=N_obsfactors,"n_factors"=N_factors,"invertTF"=invertTF,"b_type"=b_type,"Turn_off_taxonomy"=Turn_off_taxonomy,"Use_RAM_Mvalue_TF"=Use_RAM_Mvalue_TF,"rho_option"=switch(rho_space,"natural"=0,"logit"=1,"logit_with_jacobian"=2)), "Options"=c("minvar_obsfactors"=additional_variance[1],"minvar_factors"=additional_variance[2],"SD_b_stock"=SD_b_stock,"Pen_lowvar_lnRhat"=Pen_lowvar_lnRhat, "lowerbound_MLSPS"=lowerbound_MLSPS), "Cov_pz"=Cov_pz, "Y_ij"=as.matrix(Y_ij), "Missing_az"=Missing_az-1, "PC_gz"=as.matrix(ParentChild_gz[,c('ParentRowNumber','ChildTaxon')])-1, "g_i"=g_i-1)
    Data = c(Data, list("Nobs"=Nobs, "Nstock"=Nstock, "Obs2Stock"=SR_obs[,'StockNum']-1, "AR_Index"=SR_obs[,'AR_Index'], "ln_R_obs"=log(SR_obs[,'R_obs']), "SSB_obs"=SR_obs[,'SSB_obs'], "SPRF0_stock"=StockData[,'SPRF0'], "M_stock"=StockData[,'M'], "SSBmax_stock"=StockData[,'SSBmax'], "Rmax_stock"=StockData[,'Rmax'], "j_SR"=j_SR, "j_logM"=j_logM, "i_stock"=StockData[,'Stock_to_i']-1) )
  }
  if(Version%in%c("Taxon_v2_14_0")){
    Data = list("Options_vec"=c("n_obsfactors"=N_obsfactors,"n_factors"=N_factors,"invertTF"=invertTF,"b_type"=b_type,"Turn_off_taxonomy"=Turn_off_taxonomy,"Use_RAM_Mvalue_TF"=Use_RAM_Mvalue_TF,"rho_option"=switch(rho_space,"natural"=0,"logit"=1,"logit_with_jacobian"=2)), "Options"=c("minvar_obsfactors"=additional_variance[1],"minvar_factors"=additional_variance[2],"SD_b_stock"=SD_b_stock,"Pen_lowvar_lnRhat"=Pen_lowvar_lnRhat, "lowerbound_MLSPS"=lowerbound_MLSPS), "Cov_pz"=Cov_pz, "Y_ij"=as.matrix(Y_ij), "Missing_az"=Missing_az-1, "PC_gz"=as.matrix(ParentChild_gz[,c('ParentRowNumber','ChildTaxon')])-1, "g_i"=g_i-1)
    Data = c(Data, list("Nobs"=Nobs, "Nstock"=Nstock, "Obs2Stock"=SR_obs[,'StockNum']-1, "AR_Index"=SR_obs[,'AR_Index'], "ln_R_obs"=log(SR_obs[,'R_obs']), "SSB_obs"=SR_obs[,'SSB_obs'], "PredTF_stock"=PredTF_stock, "SPRF0_stock"=StockData[,'SPRF0'], "M_stock"=StockData[,'M'], "SSBmax_stock"=StockData[,'SSBmax'], "Rmax_stock"=StockData[,'Rmax'], "j_SR"=j_SR, "j_logM"=j_logM, "i_stock"=StockData[,'Stock_to_i']-1) )
  }

  # Fix potential issues
  if( "j_logM"%in%names(Data) && length(Data$j_logM)==0 ){
    Data$j_logM = 0
  }

  # Parameters
  if( Params[1]=="Generate" ){
    # Local functions
    rmatrix = function( nrow, ncol, mean=0, sd=1 ) matrix( rnorm(nrow*ncol,mean=mean,sd=sd), nrow=nrow, ncol=ncol )
    rloadings = function( n_row, n_col, sd=1, mean=0 ){
      param = vector(length=0)
      if( is.na(n_col) ){
        param = c(param, rep(0.0001,n_row))
      }else{
        if( n_col!=0 ) param = c(param, rnorm(sum(n_row:(n_row-abs(n_col)+1)),mean=mean,sd=sd))
        if( n_col<=0 ) param = c(param, rnorm(n_row,mean=mean,sd=sd))
        #if( n_row == -1*n_col ) stop("Illogical inputs")
      }
      return(param)
    }

    ### Get informative starting values for RAM variables
    # Latent variables -- rho ~= 0.5
    Y_a = rnorm( nrow(Data$Missing_az) ,sd=0.1)
    #Y_a = ifelse( colnames(Y_ij)[Data$Missing_az[,2]+1] == "rho", plogis(Y_a), Y_a )
    alpha_j = sapply( colnames(Y_ij), FUN=switch, "ln_MASPS"=-1, "ln_b"=2, 0 )

    # Generate parameter list
    if(Version%in%"Taxon_v1_0_0") Params = list( "alpha_j"=alpha_j, "obsL_z"=rloadings(n_row=n_j, n_col=Data$Options_vec['n_obsfactors']), "Y_a"=Y_a )
    if(Version%in%"Taxon_v1_1_0") Params = list( "alpha_j"=alpha_j, "L_z"=rloadings(n_row=n_j, n_col=Data$Options_vec['n_factors']), "obsL_z"=rloadings(n_row=n_j, n_col=Data$Options_vec['n_obsfactors']), "beta_gj"=rmatrix(nrow=n_g,ncol=n_j), "Y_a"=Y_a )
    if(Version%in%c("Taxon_v1_2_0")) Params = list( "alpha_j"=alpha_j, "L_z"=rloadings(n_row=n_j, n_col=Data$Options_vec['n_factors']), "obsL_z"=rloadings(n_row=n_j, n_col=Data$Options_vec['n_obsfactors']), "cov_logmult_z"=rep(0,max(Data$PC_gz[,'ChildTaxon'])+1), "beta_gj"=rmatrix(nrow=n_g,ncol=n_j), "Y_a"=Y_a )
    if(Version%in%c("Taxon_v2_0_0")){
      Params = list( "alpha_j"=alpha_j, "L_z"=rloadings(n_row=n_j, n_col=Data$Options_vec['n_factors']), "obsL_z"=rloadings(n_row=n_j, n_col=Data$Options_vec['n_obsfactors']), "cov_logmult_z"=rep(0,max(Data$PC_gz[,'ChildTaxon'])+1), "beta_gj"=rmatrix(nrow=n_g,ncol=n_j), "Y_a"=Y_a )
      Params = c( Params, list("ln_b_stock"=rep(0,ifelse(Nstock>0,Nstock,1))) )
    }
    if(Version%in%c("Taxon_v2_5_0","Taxon_v2_4_0","Taxon_v2_3_0","Taxon_v2_2_0","Taxon_v2_1_0")){
      Params = list( "alpha_j"=alpha_j, "L_z"=rloadings(n_row=n_j,n_col=Data$Options_vec['n_factors'],mean=0.1,sd=0.1), "obsL_z"=rloadings(n_row=n_j,n_col=Data$Options_vec['n_obsfactors'],mean=0.1,sd=0.1), "L_logmult_col"=rep(0,ifelse(N_factors==0,1,abs(N_factors))), "obsL_logmult_col"=rep(0,ifelse(N_obsfactors==0,1,abs(N_obsfactors))), "cov_logmult_z"=rep(0,max(Data$PC_gz[,'ChildTaxon'])+1), "beta_gj"=rmatrix(nrow=n_g,ncol=n_j), "Y_a"=Y_a )
      Params = c( Params, list("ln_b_stock"=rep(0,ifelse(Nstock>0,Nstock,1))) )
    }
    if(Version%in%c("Taxon_v2_6_0")){
      Params = list( "alpha_j"=alpha_j, "L_z"=rloadings(n_row=n_j,n_col=Data$Options_vec['n_factors'],mean=0.1,sd=0.1), "obsL_z"=rloadings(n_row=n_j,n_col=Data$Options_vec['n_obsfactors'],mean=0.1,sd=0.1), "L_logmult_col"=rep(0,ifelse(N_factors==0,1,abs(N_factors))), "obsL_logmult_col"=rep(0,ifelse(N_obsfactors==0,1,abs(N_obsfactors))), "cov_logmult_z"=rep(0,max(Data$PC_gz[,'ChildTaxon'])+1), "beta_gj"=rmatrix(nrow=n_g,ncol=n_j), "Y_a"=Y_a )
      Params = c( Params, list("bparam_stock"=rep(0,ifelse(Nstock>0,Nstock,1))) )
    }
    if(Version%in%c("Taxon_v2_7_0")){
      Params = list( "alpha_j"=alpha_j, "L_z"=rloadings(n_row=n_j,n_col=Data$Options_vec['n_factors'],mean=0.1,sd=0.1), "obsL_z"=rloadings(n_row=n_j,n_col=Data$Options_vec['n_obsfactors'],mean=0.1,sd=0.1), "L_logmult_col"=rep(0,ifelse(N_factors==0,1,abs(N_factors))), "obsL_logmult_col"=rep(0,ifelse(N_obsfactors==0,1,abs(N_obsfactors))), "cov_logmult_z"=rep(0,max(Data$PC_gz[,'ChildTaxon'])+1), "beta_gj"=rmatrix(nrow=n_g,ncol=n_j), "Y_a"=Y_a )
      Params = c( Params, list("bparam_stock"=rep(0,ifelse(Nstock>0,Nstock,1)), "gamma_p"=rep(0,ifelse(n_p==0,1,n_p))) )
    }
    if(Version%in%c("Taxon_v2_14_0","Taxon_v2_13_0","Taxon_v2_12_0","Taxon_v2_11_0","Taxon_v2_10_0","Taxon_v2_9_0","Taxon_v2_8_0")){
      Params = list( "alpha_j"=alpha_j, "L_z"=rloadings(n_row=n_j,n_col=Data$Options_vec['n_factors'],mean=1,sd=0.1), "obsL_z"=rloadings(n_row=n_j,n_col=Data$Options_vec['n_obsfactors'],mean=1,sd=0.1), "L_logmult_col"=rep(0,ifelse(N_factors==0,1,abs(N_factors))), "obsL_logmult_col"=rep(0,ifelse(N_obsfactors==0,1,abs(N_obsfactors))), "cov_logmult_z"=rep(0,max(Data$PC_gz[,'ChildTaxon'])+1), "beta_gj"=rmatrix(nrow=n_g,ncol=n_j), "Y_a"=Y_a )
      Params = c( Params, list("bparam_stock"=rep(0,ifelse(Nstock>0,Nstock,1)), "gamma_p"=rep(0,ifelse(n_p==0,1,n_p)), "theta_q"=0) )
    }
  }

  # Random
  if( Random[1]=="Generate" ){
    if(Version%in%"Taxon_v1_0_0") Random = c("Y_a")
    if(Version%in%c("Taxon_v2_14_0","Taxon_v2_13_0","Taxon_v2_12_0","Taxon_v2_11_0","Taxon_v2_10_0","Taxon_v2_9_0","Taxon_v2_8_0","Taxon_v2_7_0","Taxon_v2_6_0","Taxon_v2_5_0","Taxon_v2_4_0","Taxon_v2_3_0","Taxon_v2_2_0","Taxon_v2_1_0","Taxon_v2_0_0","Taxon_v1_2_0","Taxon_v1_1_0")) Random = c("Y_a", "beta_gj")
    if(Use_REML==TRUE){
      Random = c(Random, "alpha_j")
      if("ln_b_stock" %in% names(Params)) Random = union(Random, "ln_b_stock")
      if("bparam_stock" %in% names(Params)) Random = union(Random, "bparam_stock")
      if("gamma_p" %in% names(Params)) Random = union(Random, "gamma_p")
      if("theta_q" %in% names(Params)) Random = union(Random, "theta_q")
    }
  }

  # Map
  if( Map[1] == "Generate" ){
    Map = list()
    if( "cov_logmult_z" %in% names(Params) ){
      if(Process_cov=="Equal") Map[["cov_logmult_z"]] = factor( rep(NA,length(Params[["cov_logmult_z"]])) )
      if(Process_cov=="Unequal") Map[["cov_logmult_z"]] = factor( c(NA,1:(length(Params[["cov_logmult_z"]])-1)) )
    }
    if(is.na(Data$Options_vec['n_obsfactors'])){
      Map[["obsL_z"]] = factor( rep(NA,length(Params[["obsL_z"]])) )
    }
    if( Nstock==0 ){
      if("ln_b_stock" %in% names(Params)) Map[["ln_b_stock"]] = factor(rep(NA,length(Params[["ln_b_stock"]])))
      if("bparam_stock" %in% names(Params)) Map[["bparam_stock"]] = factor(rep(NA,length(Params[["bparam_stock"]])))
    }
    if( N_factors==0 | TRUE ){
      if( "L_logmult_col" %in% names(Params) ) Map[["L_logmult_col"]] = factor(rep(NA,length(Params[["L_logmult_col"]])))
    }
    if( N_obsfactors==0 | TRUE ){
      if( "obsL_logmult_col" %in% names(Params) ) Map[["obsL_logmult_col"]] = factor(rep(NA,length(Params[["obsL_logmult_col"]])))
    }
    if( length(j_SR)==4 ){
      if("ln_b_stock" %in% names(Params)) Map[["ln_b_stock"]] = factor(rep(NA,length(Params[["ln_b_stock"]])))
      if("bparam_stock" %in% names(Params)) Map[["bparam_stock"]] = factor(rep(NA,length(Params[["bparam_stock"]])))
    }
    if( n_p==0 ){
      if( "gamma_p" %in% names(Params) ) Map[["gamma_p"]] = factor(NA)
    }
    if( Cov_RAM[["M"]]==FALSE ){
      if( "theta_q" %in% names(Params) ) Map[["theta_q"]] = factor(NA)
    }
    # Map off species where PredTF_stock == TRUE
    if( Nstock > 0 ){
      Map[["bparam_stock"]] = factor( ifelse(PredTF_stock==TRUE, NA, 1:Data$Nstock) )
    }

    # Turn off taxonomic hierarchy, such that all taxa have expected LH params = alpha_j
    if( Turn_off_taxonomy==TRUE ){
      Map[["beta_gj"]] = factor( array(NA, dim=dim(Params[["beta_gj"]])) )
      Map[["L_z"]] = factor( rep(NA,length(Params[["L_z"]])) )
    }

    # Change
    if( is.na(Data$Options_vec[1]) ){
      Data$Options_vec[1] = 0
    }

    # Zero out covariance for some variables
    if( N_factors!=0 ){
      Mat = diag(ncol(Y_ij))[1:abs(N_factors),,drop=FALSE]
      Which = upper.tri(Mat,diag=TRUE)
      rownum = col(Mat)[Which]
      Which = which( rownum %in% which(zerocovTF_j) )
      Params[["L_z"]][Which] = 0
      Map[["L_z"]] = 1:length(Params[["L_z"]])
      Map[["L_z"]][Which] = NA
      Map[["L_z"]] = factor(Map[["L_z"]])
    }
    if( N_obsfactors!=0 ){
      Mat = diag(ncol(Y_ij))[1:abs(N_obsfactors),,drop=FALSE]
      Which = upper.tri(Mat,diag=TRUE)
      rownum = col(Mat)[Which]
      Which = which( rownum %in% which(zerocovTF_j) )
      Params[["obsL_z"]][Which] = 0
      Map[["obsL_z"]] = 1:length(Params[["obsL_z"]])
      Map[["obsL_z"]][Which] = NA
      Map[["obsL_z"]] = factor(Map[["obsL_z"]])
    }
  }

  # Check for problems
  if( rho_space%in%c("logit","logit_with_jacobian") & !("rho_option"%in%names(Data$Options_vec)) ){
    stop("Using non-default value for `rho_space` doesn't make sense with the version you are using")
  }

  # Debugging option
  if( debug_mode==TRUE ){
    on.exit( assign("Params", Params, envir=.GlobalEnv), add=TRUE )
    on.exit( assign("Data", Data, envir=.GlobalEnv), add=TRUE )
    on.exit( assign("Map", Map, envir=.GlobalEnv), add=TRUE )
    on.exit( assign("Random", Random, envir=.GlobalEnv), add=TRUE )
  }

  #####################
  # Build and run
  #####################

  Test = require(TMB, quietly=TRUE)
  if( Test==FALSE ){
    stop("Please install package `TMB` from CRAN")
  }

  # Compile TMB software
  #dyn.unload( paste0(RunDir,"/",dynlib(Version)) ) # random=Random,
  dir.create( RunDir )
  file.copy( from=paste0(TmbDir,"/",Version,".cpp"), to=paste0(RunDir,"/",Version,".cpp"), overwrite=FALSE)
  origwd = getwd()
  on.exit(setwd(origwd),add=TRUE)
  setwd( RunDir )
  compile( paste0(Version,".cpp") )

  # Build
  dyn.load( paste0(RunDir,"/",TMB::dynlib(Version)) )          #
  Obj = MakeADFun( data=Data, parameters=Params, DLL=Version, map=Map, random=Random )
  Report = Obj$report()

  # Print number of parameters
  ThorsonUtilities::list_parameters( Obj )

  #
  if( run_model==FALSE ){
    Return = list("Data"=Data, "Params"=Params, "Random"=Random, "Map"=Map, "Obj"=Obj )
    return(Return)
  }

  # Debugging option
  if( debug_mode==TRUE ){
    #on.exit( assign("ParHat",Obj$env$parList(),envir=.GlobalEnv), add=TRUE )
    on.exit( assign("Obj",Obj,envir=.GlobalEnv), add=TRUE )
  }

  # Print to screen
  if( verbose==TRUE ){
    cat( "Number of fixed effects:")
    print( table(names(Obj$par)) )
    cat( "Number of random effects:")
    print( table(names(Obj$env$last.par[Obj$env$random])) )
  }

  # Optimize                         #  , startpar=opt$par[-grep("alpha",names(opt$par))]
  # JointPrecision is used below, and is too big to invert whole;  must have getReportCovariance=TRUE to get JointPrecision
  Opt = TMBhelper::fit_tmb( obj=Obj, savedir=RunDir, getJointPrecision=TRUE, getReportCovariance=TRUE, ... )

  #
  Report = Obj$report()
  colnames(Report$Ycomplete_ij) = colnames(Report$beta_gj) = colnames(Y_ij)

  # Debugging option
  if( debug_mode==TRUE ){
    #on.exit( assign("ParHat",Obj$env$parList(),envir=.GlobalEnv), add=TRUE )
    on.exit( assign("Opt",Opt,envir=.GlobalEnv), add=TRUE )
  }

  # return if necessary
  if( "h" %in% names(Opt)){
    return( list("Opt"=Opt, "Obj"=Obj) )
  }

  # Debugging
  if( FALSE ){
    plot( x=Report$ln_R_obs_hat, y=Data$ln_R_obs )
    plot( x=Report$mu_obs, y=Data$ln_R_obs )
    Table = cbind( Report$SD_stock, Report$ro_stock, Report$Ycomplete_ij[Data$i_stock+1,Data$j_SR+1] )
  }

  # SE
  ParHat = Obj$env$parList()
  ParHat_SE = as.list( Opt$SD, what="Std" )
  colnames(ParHat$beta_gj) = colnames(ParHat_SE$beta_gj) = colnames(Report$Ycomplete_ij) = colnames(Y_ij)
  #dyn.unload( paste0(RunDir,"/",TMB::dynlib(Version)) )          #

  # Return stuff
  Return = list("N_factors"=N_factors, "N_obsfactors"=N_obsfactors, "Use_REML"=Use_REML,
    "ParentChild_gz"=ParentChild_gz, "ParHat"=ParHat, "g_i"=g_i,
    "Y_ij"=Y_ij, "Z_ik"=Z_ik, "SR_obs"=SR_obs, "StockData"=StockData,
    "Obj"=Obj, "Opt"=Opt, "Report"=Report, "ParHat_SE"=ParHat_SE, "obsCov_jj"=Report$obsCov_jj)
  dimnames(Return$obsCov_jj) = list(colnames(Y_ij),colnames(Y_ij))
  if("Cov_jj" %in% names(Report)){
    Return = c(Return, list("Cov_jj"=Report$Cov_jj))
  }
  if("cov_logmult_z" %in% names(ParHat)){
    Return = c(Return, list("Cov_jjz"=Report$Cov_jj %o% exp(ParHat$cov_logmult_z)))
    dimnames(Return$Cov_jjz) = list(colnames(Y_ij),colnames(Y_ij),colnames(Z_ik))
  }

  ####################
  # Interpret results
  ####################

  if( extract_covariance==TRUE ){
    ### Approximate joint precision
    message( "Extracting predictive variance of life-history parameters for each taxon...")
    # Extract predictive covariance for species-specific traits (necessary to do separately for rows and columns)
    Prec_zz = Opt$SD$jointPrecision[ , grep("beta_gj",colnames(Opt$SD$jointPrecision)) ]
      Prec_zz = Prec_zz[ grep("beta_gj",rownames(Opt$SD$jointPrecision)), ]
    # Predict
    u_zr = rmvnorm_prec( mu=Obj$env$last.par.best, prec=Opt$SD$jointPrecision, n.sims=1000 )
    # Extract and invert
    VarNames = Predictive_distribution( mean_vec=Y_ij[1,], process_cov=NULL, obs_cov=NULL, check_names=TRUE, include_r=include_r )
    n_v = length(VarNames)
    Prec_gjj = PartialCorr_gjj = array(NA, dim=c(n_g,n_j,n_j), dimnames=list(ParentChild_gz[,'ChildName'],colnames(Y_ij),colnames(Y_ij)) )
    Prec_gvv = PartialCorr_gvv = Corr_gvv = Cov_gvv = array(NA, dim=c(n_g,n_v,n_v), dimnames=list(ParentChild_gz[,'ChildName'],VarNames,VarNames) )
    beta_gv = array(NA, dim=c(n_g,n_v), dimnames=list(ParentChild_gz[,'ChildName'],VarNames) )
    for( gI in 1:n_g ){
      # Method with full precision
      Indices = gI + ( seq(1,n_g*n_j,by=n_g)-1 )
      Samp_rj = t(u_zr[grep("beta_gj",colnames(Opt$SD$jointPrecision))[Indices],])
      colnames(Samp_rj) = colnames(Y_ij)
      Pred = Predictive_distribution( mean_vec=Report$beta_gj[gI,], Samp_rj=Samp_rj, include_obscov=FALSE, check_bounds=FALSE, include_r=include_r, lowerbound_MLSPS=lowerbound_MLSPS, rho_option=switch(rho_space,"natural"=0,"logit"=1,"logit_with_jacobian"=2) )
      beta_gv[gI,] = Pred$pred_mean
      Cov_gvv[gI,,] = Pred$pred_cov
      Corr_gvv[gI,,] = cov2cor( Cov_gvv[gI,,] )
      if( (gI%%1000) == 0 ) message( "Finished processing predictive variance for ", gI, " of ",n_g," taxa" )
    }

    # Return stuff
    if("Cov_jj" %in% names(Report)){
      Return = c(Return, list("Cov_gvv"=Cov_gvv, "beta_gv"=beta_gv))
    }
  }

  Return[["total_runtime"]] = Sys.time() - start_time
  return( Return )
}
