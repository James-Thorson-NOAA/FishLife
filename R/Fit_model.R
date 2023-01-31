#' Fit and predict fish traits
#'
#' \code{Fit_model} estimates parameters and predicts values from a multivariate random-walk model for fish traits
#'
#' @inheritParams sem::sem
#'
#' @param text structural equation model structure, passed to either \code{\link[sem]{specifyModel}}
#'        and then parsed to control the set of path coefficients and variance-covariance parameters
#' @param tree phylogenetic structure, using class \code{\link[ape]{as.phylo}}.  If \code{tree=NULL}
#'        then argument \code{Z_ik} is instead converted to a tree using 
#'        \code{ape::as.phylo( ~Class/Order/Family/Genus/GenusSpecies, ... )} to convert it to a tree
#' @param min_replicate_measurements specified threshold for the number of measurements for a given continuous traits,
#'        where any continuous trait having fewer replicated measurements for at least one taxon will have the 
#'        measurement variance fixed at an arbitrarily low value, such that estimated traits are forced
#'        to approach the unreplicated measurements
#' @param N_factors Number of factors in decomposition of covariance for random-walk along evolutionary tree (0 means a diagonal but unequal covariance; negative is the sum of a factor decomposition and a diagonal-but-unequal covariance)
#' @param N_obsfactors Number of factors in decomposotion of observation covariance (same format as \code{N_obsfactors})
#' @param Database Whether to use results for both adult and stock-recruit parameters, \code{Database="FishBase_and_RAM"}, or just adult parameters, \code{Database="FishBase"}
#' @param Use_REML OPTIONAL boolean whether to use maximum marginal likelihood or restricted maximum likelihood (termed "REML")
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
#' @param ... other paramers passed to \code{\link[TMBhelper]{fit_tmb}}
#'
#' @return Tagged list containing objects from FishLife run (first 9 slots constitute list 'Estimate_database' for archiving results)
#' \describe{
#'   \item{tree}{The phylogenetic tree used for analysis, whether inputted or generated from taxonomy based on \code{Z_ik}}
#'   \item{SEM_model}{The phylogenetic tree used for analysis, whether inputted or generated from taxonomy based on \code{Z_ik}}
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
#' @examples
#' \dontrun{
#' # Load data set
#' library(FishLife)
#'
#' # Simulate data
#' tree = ape::rtree(n=100)
#' xfit = x = 1 + 0.3 * phylolm::rTrait(n = 1, phy=tree)
#' yfit = y = 2 + 1*x + 0.3 * phylolm::rTrait(n = 1, phy=tree)
#' drop_x = sample( 1:ape::Ntip(tree), replace=FALSE, size=round(ape::Ntip(tree)*0.3) )
#' xfit[drop_x] = NA
#' drop_y = sample( 1:ape::Ntip(tree), replace=FALSE, size=round(ape::Ntip(tree)*0.3) )
#' yfit[drop_y] = NA
#'
#' # Fit model
#' Fit = Fit_model(
#'   text = "x -> y, p",
#'   Database = NULL,
#'   Use_REML = FALSE,
#'   Y_ij = data.frame( x=xfit, y=yfit ),
#'   tree = tree,
#'   min_replicate_measurements = Inf )
#'
#' # S3-defaults
#' print(Fit)
#' coef(Fit)
#'
#' # Convet and plot using sem
#' mysem = as(Fit,"sem")
#' sem::pathDiagram( model = mysem,
#'                   style = "traditional",
#'                   edge.labels = "values" )
#' myplot = semPlot::semPlotModel( as(Fit,"sem") )
#' semPlot::semPaths( myplot,
#'                    nodeLabels = myplot@Vars$name )
#'
#' # Convert and plot using phylobase / phylosignal
#' library(phylobase)
#' plot( as(Fit,"phylo4d") )
#' }
#'
#' @export
Fit_model <-
function( text = NULL,
          Database = FishLife::FishBase_and_RAM,
          tree = Database$tree,
          Y_ij = Database$Y_ij,
          Z_ik = Database$Z_ik,
          Use_REML = TRUE,
          N_factors = 0,
          N_obsfactors = 0,
          min_replicate_measurements = 0,
          SR_obs = Database$SR_obs,
          StockData = Database$StockData,
          group_j = 1:ncol(Y_ij)-1,
          Version = "Taxon_v3_0_0",
          TmbDir = system.file("executables",package="FishLife"),
          RunDir = tempfile(pattern = "run_",tmpdir=tempdir(),fileext="/"),
          verbose = FALSE,
          debug_mode = FALSE,
          j_SR = ncol(Y_ij)-3:1,
          additional_variance = c(0,0),
          SD_b_stock = 10,
          b_type = 0,
          Turn_off_taxonomy = FALSE,
          Pen_lowvar_lnRhat = 1,
          lowerbound_MLSPS = 1,
          Use_RAM_Mvalue_TF = TRUE,
          rho_space = "natural",
          n_sims = 1000,
          n_batches = NULL,
          include_r = TRUE,
          PredTF_stock = NULL,
          extract_covariance = FALSE,
          run_model = TRUE,
          multinomial_for_factors = FALSE,
          Params = "Generate",
          Random = "Generate",
          Map = "Generate",
          add_predictive = FALSE,
          ... ){

  #save.image( paste0(RunDir,"Save.RData") )
  #stop()

  #### Local function
  # Sample from GMRF using sparse precision
  rmvnorm_prec <- function(mu, prec, n_sims, seed=1, varnames="", n_batches=NULL) {
    set.seed( seed )
    which_vars = grep(varnames, names(mu))
    if(is.null(n_batches)) n_batches = ceiling( length(which_vars) * n_sims * 8 / 4e9 )
    if( (length(which_vars) * n_sims * 8 / n_batches) > 4e9 ) stop("Decrease `n_sims` to avoid R memory limits")
    batch_index = cut(1:n_sims, n_batches)
    out = mu[which_vars] %o% rep(1,n_sims)
    for( batch in levels(batch_index) ){
      which_index = which(batch_index == batch )
      z <- matrix(rnorm(length(mu) * length(which_index)), ncol=length(which_index))
      L <- Matrix::Cholesky(prec, super=TRUE)
      z <- Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
      z <- Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
      out[,which_index] <- out[,which_index] + as.matrix(z)[which_vars,]
    }
    return(out)
  }

  # Function that converts SEM model to a RAM, see `?sem` for more context
  build_ram = function( model, vars ){
    vars = sapply( vars, FUN=function(char){gsub("-", "", gsub(" ", "", char))} )
    n.paths = nrow(model)
    par.names = model[, 2]
    startvalues = model[,3]

    # EXCERPT FROM `getAnywhere("sem.semmod")`
    heads = from = to = rep(0, n.paths)
    for (p in 1:n.paths) {
      path = sem:::parse.path(model[p, 1])
      heads[p] = abs(path$direction)
      to[p] = path$second
      from[p] = path$first
      if (path$direction == -1) {
        to[p] = path$first
        from[p] = path$second
      }
    }
    missing_vars = setdiff( c(from,to), vars )
    if( length(missing_vars) > 0 ) stop( "Check `build_ram`:", paste0(missing_vars,sep=", ") )

    ram = data.frame(matrix(0, nrow=p, ncol=5))
    pars = na.omit(unique(par.names))
    ram[, 1] = heads
    ram[, 2] = apply(outer(vars, to, "=="), 2, which)
    ram[, 3] = apply(outer(vars, from, "=="), 2, which)
    par.nos = apply(outer(pars, par.names, "=="), 2, which)
    if(length(par.nos) > 0){
      ram[, 4] = unlist(lapply(par.nos, function(x) if (length(x)==0){0}else{x}))
    }
    ram[, 5] = startvalues
    colnames(ram) = c("heads", "to", "from", "parameter", "start")
    return(ram)
  }

  #####################
  # Check for potential problems
  #####################
  start_time = Sys.time()

  # Force setting
  Process_cov = "Equal"

  # Build RAM
  if(!is.null(text)){
    SEM_model = sem::specifyModel( text=text, exog.variances=TRUE, endog.variances=TRUE, covs=colnames(Y_ij) )
    RAM = build_ram( SEM_model, colnames(Y_ij) )
  }else{
    RAM = array( NA, dim=c(0,5), dimnames=list(NULL,c("heads","to","from","parameter","start")) )
  }

  # Process tree if supplied
  if( !is.null(Z_ik) & !is.null(tree) ){
    stop("Please only supply either `Z_ik` or `tree`")
  }
  if( is.null(Z_ik) & !is.null(tree) ){
    message("Using supplied tree")
    if(!all(rownames(Y_ij) %in% tree$tip.label)) stop("Some `rownames(Y_ij)` not in `tree$tip.label`")
    if(isTRUE(add_predictive)) stop("`add_predictive` is not designed to work when specifying a tree")
  }
  if( !is.null(Z_ik) & is.null(tree) ){
    # Code uses "_" so throw error if included in names
    if( any(apply(Z_ik, MARGIN=2, FUN=function(vec){length(grep("_",vec))})>0) ){
      stop("Please do not use character '_' in names included in 'Z_ik'")
    }
    message("Converting taxonomy `Z_ik` to tree of class `ape::phylo`")

    # Step 1 -- Build tree
    Data = Z_ik
    Data$GenusSpecies = paste0( Data$Genus, " ", Data$Species )
    for(colI in 1:ncol(Data)) Data[,colI] = factor(Data[,colI])
    Data = unique(Data)
    Formula = ~ Class/Order/Family/Genus/GenusSpecies
    tree = ape::as.phylo(Formula, data = Data, collapse=FALSE) # collapse=FALSE preserves taxa-hierarchy to get edge-length right
    tree$edge.length = rep(1,nrow(tree$edge))
    if( isFALSE(add_predictive) ){
      tree = ape::collapse.singles(tree)
    }

    # Step 2
    if(nrow(Y_ij)!=nrow(Z_ik)) stop("Check number of rows in `Y_ij` and `Z_ik`")
    #return(list(Y_ij=Y_ij, Data=Data))
    rownames(Y_ij) = paste0( Z_ik$Genus, " ", Z_ik$Species )

    message("Done with conversion")
  }
  if( is.null(Z_ik) & is.null(tree) ){
    stop("Must supply either `Z_ik` or `tree`")
  }

  # Step 3 -- Add "predictive" tip for each node
  if( isTRUE(add_predictive) ){
    message("Adding tips for predictive distribution")
    tree_full = tree
    depth = ape::node.depth.edgelength(tree_full)
    names(depth) = c( tree$tip.label, tree$node.label )
    tips_to_add = 5 - depth
    for( node in 1:length(tips_to_add) ){
      if( tips_to_add[node] > 0 ){
        Match = match( names(tips_to_add)[node], tree_full$node.label )
        Where = ape::Ntip(tree_full) + Match
        tip = ape::read.tree(text=paste0(paste0(rep("(",tips_to_add[node]),collapse=""),"predictive",paste0(rep(")",tips_to_add[node]-1),collapse=""),",b);"))
        tip = ape::drop.tip( phy=tip, tip="b", collapse.singles=FALSE, trim.internal=ifelse(tips_to_add[node]==1,TRUE,FALSE) )
        tip$edge.length = rep(1, nrow(tip$edge))
        tip$node.label = rep("predictive", tip$Nnode)
        tree_full = ape::bind.tree( tree_full, tip, where=Where )
      }
    }
    tree = tree_full
  }
  
  # Test naming in FishLife::FishBase$ParentChild_gz[,1]
  if( FALSE ){
    x = as.character(FishLife::FishBase$ParentChild_gz[,1])
    # x = Names_full
    Len = sapply( x, FUN=function(char){length(strsplit(char,"_")[[1]])} )
    Len = sapply( x, FUN=function(char){length(strsplit(char," ")[[1]])} )
  }

  # Step 4 -- Fix names
  if( isTRUE(add_predictive) ){
    message("Re-formatting node and tip names to match `Search_species` format")
    Names = Names_full = c( tree$tip.label, tree$node.label )
    Names[seq_len(ape::Ntip(tree))] = sapply( Names[seq_len(ape::Ntip(tree))], FUN=function(x){
      y = strsplit(x," ")[[1]]
      if(length(y)>=2) x = paste0(y[2:length(y)],collapse=" ")
      return(x)
    })
    Path = ape::nodepath(tree)
    for( nI in seq_len(ape::Ntip(tree)) ){
      Names_full[nI] = paste0( Names[Path[[nI]]][-1], collapse="_" )
    }
    root_index = ape::Ntip(tree) + 1
    for( nI in ape::Ntip(tree)+seq_len(ape::Nnode(tree)) ){
      Path = ape::nodepath(tree, from=nI, to=root_index )
      Names_full[nI] = paste0( Names[rev(Path)][-1], collapse="_" )
      #if( FALSE ){
      #  Which = which( tree$edge[,2] == nI )
      #  if(length(Which)==1){
      #    Names_full[nI] = paste( setdiff(Names[tree$edge[Which,1]],""), Names_full[nI], sep="_" )
      #  }else if(length(Which)>1){
      #    stop()
      #  }
      #}
    }
    # remote root
    Names_full = Names_full[-root_index]
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

  Cov_pz = matrix( c(0,0), ncol=2, nrow=1 )
  n_p = 0

  #####################
  # Pre-process data
  #####################

  # Figure out missingness
  Missing_az = matrix( nrow=0, ncol=2 )
  for(jI in 1:ncol(Y_ij)){
    Which = which(is.na(Y_ij[,jI]))
    if( length(Which)>0 ){
      Missing_az = rbind(Missing_az, cbind(Which,jI))
    }
  }

  # Figure out network structure for taxonomy
  #ParentChild_gz = NULL
  ## 1st column: child taxon name
  ## 2nd column: parent taxon name
  ## 3rd column: parent row-number in ParentChild_gz
  ## 4th column: Taxon level
  ## Loop through
  #for( colI in 1:ncol(Z_ik)){
  #  Taxa_Names = apply( Z_ik[,1:colI,drop=FALSE], MARGIN=1, FUN=paste, collapse="_")
  #  Unique_Taxa = unique(Taxa_Names)
  #  for( uniqueI in 1:length(Unique_Taxa) ){
  #    Which = which( Taxa_Names == Unique_Taxa[uniqueI] )
  #    if( colI==1 ){
  #      ParentChild_gz = rbind( ParentChild_gz, c(Unique_Taxa[uniqueI], NA, NA, colI) )
  #    }else{
  #      if( length(unique(Z_ik[Which,colI-1]))>1 ) stop("Taxa has multiple parents")
  #      ChildName = Unique_Taxa[uniqueI]
  #      ParentName = paste(rev(rev(strsplit(ChildName,"_")[[1]])[-1]),collapse="_")
  #      ParentChild_gz = rbind( ParentChild_gz, c(ChildName, ParentName, match(ParentName,ParentChild_gz[,1]), colI) )
  #    }
  #  }
  #}
  ## Loop through again to add predictive elements
  #for( colI in 1:(ncol(Z_ik)-1)){
  #  Taxa_Names = apply( Z_ik[,1:colI,drop=FALSE], MARGIN=1, FUN=paste, collapse="_")
  #  Unique_Taxa = unique(Taxa_Names)
  #  for( uniqueI in 1:length(Unique_Taxa) ){
  #    ParentName = Unique_Taxa[uniqueI]
  #    for( predI in 1:(ncol(Z_ik)-colI) ){
  #      ChildName = paste0(ParentName,"_predictive")
  #      ParentChild_gz = rbind( ParentChild_gz, c(ChildName, ParentName, match(ParentName,ParentChild_gz[,1]), colI+predI) )
  #      ParentName = ChildName
  #    }
  #  }
  #}
  ## Add top predictive
  #ParentChild_gz = rbind( ParentChild_gz, c("predictive", NA, NA, 1) )
  #for( colI in 2:ncol(Z_ik)) ParentChild_gz = rbind( ParentChild_gz, c(paste(rep("predictive",colI),collapse="_"), paste(rep("predictive",colI-1),collapse="_"), match(paste(rep("predictive",colI-1),collapse="_"),ParentChild_gz[,1]), colI) )
  ## Relabel
  #ParentChild_gz = data.frame( ParentChild_gz )
  #colnames(ParentChild_gz) = c("ChildName", "ParentName", "ParentRowNumber", "ChildTaxon")
  #ParentChild_gz[,'ParentRowNumber'] = as.numeric(as.character(ParentChild_gz[,'ParentRowNumber']))
  #ParentChild_gz[,'ChildTaxon'] = as.numeric(as.character(ParentChild_gz[,'ChildTaxon']))

  #root_index = Ntip(tree) + 1
  #ParentChild_gz = data.frame( ChildName=Names_full[tree$edge[,2]], ParentName=Names_full[tree$edge[,1]],
  #              ParentRowNumber=tree$edge[,1], ChildTaxon=node.depth.edgelength(tree)[tree$edge[,2]] )
  #ParentChild_gz[which(ParentChild_gz[,'ParentRowNumber']==root_index),"ParentRowNumber"] = NA
  #ParentChild_gz[which(ParentChild_gz[,'ParentRowNumber']==root_index),"ParentName"] = NA

  root_node = ape::Ntip(tree) + 1
  Edge = tree$edge
    Edge[,1] = ifelse( Edge[,1]==root_node, NA, Edge[,1] )
    Edge[,1:2] = ifelse( Edge>root_node, Edge-1, Edge )
  #Depth = ape::node.depth.edgelength(tree)
  #  Depth = Depth[-root_node]
  Length = tree$edge.length
  Order = order(Edge[,2])
  Depth = ifelse( is.na(Edge[,1]), 1, 2 )
  if( any(Length==0) ){
    min_nonzero_length = min( Length[Length!=0] )
    warning("Fixing ", sum(Length==0), " edges with zero length")
    Length = ifelse( Length==0, min(0.001,0.1*min_nonzero_length), Length)
  }
  ParentChild_gz = data.frame( ChildName = Names_full[Edge[Order,2]], 
                               ParentName = Names_full[Edge[Order,1]],
                               ParentRowNumber = Edge[Order,1], 
                               ChildTaxon = Depth[Order],     # Force ChildTaxon = 2
                               EdgeLength = Length[Order] )
  PC_gz = as.matrix(ParentChild_gz[,c('ParentRowNumber','ChildTaxon')]) - 1

  if( !all(1:nrow(Edge) %in% unique(as.vector(Edge))) ){
    stop("Check ParentChild_gz:  Some node is not included in graph")
  }
  #if( max(PC_gz[,1],na.rm=TRUE) >= nrow(betainput_gj) ){
  #  stop("Check PC_gz")
  #}

  # Identify location for every observation
  # ParentChild_gz is ordered with rownumber 1-10 for observations 1-10, etc
  #g_i = match( rownames(Y_ij), tree$tip.label[Order] )
  g_i = match( rownames(Y_ij), tree$tip.label )

  n_j = ncol(Y_ij)
  n_g = nrow(ParentChild_gz)

  # Context-specific inputs
  if( "M" %in% colnames(Y_ij) ){
    j_logM = which(colnames(Y_ij)=="M")-1
  }else{
    j_logM = 0
  }

  #####################
  # TMB inputs
  #####################

  # Data
  if(Version%in%"Taxon_v1_0_0") Data = list("Options_vec"=c("n_obsfactors"=N_obsfactors), "Z_ik"=as.matrix(Z_ik)-1, "Y_ij"=as.matrix(Y_ij), "Missing_az"=Missing_az-1)
  if(Version%in%c("Taxon_v1_2_0","Taxon_v1_1_0")) Data = list("Options_vec"=c("n_obsfactors"=N_obsfactors,"n_factors"=N_factors), "Y_ij"=as.matrix(Y_ij), "Missing_az"=Missing_az-1, "PC_gz"=PC_gz, "g_i"=g_i-1)
  if(Version%in%c("Taxon_v2_2_0","Taxon_v2_1_0","Taxon_v2_0_0")){
    Data = list("Options_vec"=c("n_obsfactors"=N_obsfactors,"n_factors"=N_factors,"invertTF"=FALSE), "Y_ij"=as.matrix(Y_ij), "Missing_az"=Missing_az-1, "PC_gz"=PC_gz, "g_i"=g_i-1)
    Data = c(Data, list("Nobs"=Nobs, "Nstock"=Nstock, "Obs2Stock"=SR_obs[,'StockNum']-1, "AR_Index"=SR_obs[,'AR_Index'], "ln_R_obs"=log(SR_obs[,'R_obs']), "SSB_obs"=SR_obs[,'SSB_obs'], "SPRF0_stock"=StockData[,'SPRF0'], "M_stock"=StockData[,'M'], "j_SR"=j_SR, "i_stock"=StockData[,'Stock_to_i']-1) )
  }
  if(Version%in%c("Taxon_v2_5_0","Taxon_v2_4_0","Taxon_v2_3_0")){
    Data = list("Options_vec"=c("n_obsfactors"=N_obsfactors,"n_factors"=N_factors,"invertTF"=FALSE), "Options"=c("minvar_obsfactors"=additional_variance[1],"minvar_factors"=additional_variance[2],"SD_b_stock"=SD_b_stock), "Y_ij"=as.matrix(Y_ij), "Missing_az"=Missing_az-1, "PC_gz"=PC_gz, "g_i"=g_i-1)
    Data = c(Data, list("Nobs"=Nobs, "Nstock"=Nstock, "Obs2Stock"=SR_obs[,'StockNum']-1, "AR_Index"=SR_obs[,'AR_Index'], "ln_R_obs"=log(SR_obs[,'R_obs']), "SSB_obs"=SR_obs[,'SSB_obs'], "SPRF0_stock"=StockData[,'SPRF0'], "M_stock"=StockData[,'M'], "j_SR"=j_SR, "i_stock"=StockData[,'Stock_to_i']-1) )
  }
  if(Version%in%c("Taxon_v2_6_0")){
    Data = list("Options_vec"=c("n_obsfactors"=N_obsfactors,"n_factors"=N_factors,"invertTF"=FALSE,"b_type"=b_type), "Options"=c("minvar_obsfactors"=additional_variance[1],"minvar_factors"=additional_variance[2],"SD_b_stock"=SD_b_stock), "Y_ij"=as.matrix(Y_ij), "Missing_az"=Missing_az-1, "PC_gz"=PC_gz, "g_i"=g_i-1)
    Data = c(Data, list("Nobs"=Nobs, "Nstock"=Nstock, "Obs2Stock"=SR_obs[,'StockNum']-1, "AR_Index"=SR_obs[,'AR_Index'], "ln_R_obs"=log(SR_obs[,'R_obs']), "SSB_obs"=SR_obs[,'SSB_obs'], "SPRF0_stock"=StockData[,'SPRF0'], "M_stock"=StockData[,'M'], "SSBmax_stock"=StockData[,'SSBmax'], "Rmax_stock"=StockData[,'Rmax'], "j_SR"=j_SR, "i_stock"=StockData[,'Stock_to_i']-1) )
  }
  if(Version%in%c("Taxon_v2_8_0","Taxon_v2_7_0")){
    Data = list("Options_vec"=c("n_obsfactors"=N_obsfactors,"n_factors"=N_factors,"invertTF"=FALSE,"b_type"=b_type), "Options"=c("minvar_obsfactors"=additional_variance[1],"minvar_factors"=additional_variance[2],"SD_b_stock"=SD_b_stock), "Cov_pz"=Cov_pz, "Y_ij"=as.matrix(Y_ij), "Missing_az"=Missing_az-1, "PC_gz"=PC_gz, "g_i"=g_i-1)
    Data = c(Data, list("Nobs"=Nobs, "Nstock"=Nstock, "Obs2Stock"=SR_obs[,'StockNum']-1, "AR_Index"=SR_obs[,'AR_Index'], "ln_R_obs"=log(SR_obs[,'R_obs']), "SSB_obs"=SR_obs[,'SSB_obs'], "SPRF0_stock"=StockData[,'SPRF0'], "M_stock"=StockData[,'M'], "SSBmax_stock"=StockData[,'SSBmax'], "Rmax_stock"=StockData[,'Rmax'], "j_SR"=j_SR, "i_stock"=StockData[,'Stock_to_i']-1) )
  }
  if(Version%in%c("Taxon_v2_9_0")){
    Data = list("Options_vec"=c("n_obsfactors"=N_obsfactors,"n_factors"=N_factors,"invertTF"=FALSE,"b_type"=b_type,"Turn_off_taxonomy"=Turn_off_taxonomy), "Options"=c("minvar_obsfactors"=additional_variance[1],"minvar_factors"=additional_variance[2],"SD_b_stock"=SD_b_stock), "Cov_pz"=Cov_pz, "Y_ij"=as.matrix(Y_ij), "Missing_az"=Missing_az-1, "PC_gz"=PC_gz, "g_i"=g_i-1)
    Data = c(Data, list("Nobs"=Nobs, "Nstock"=Nstock, "Obs2Stock"=SR_obs[,'StockNum']-1, "AR_Index"=SR_obs[,'AR_Index'], "ln_R_obs"=log(SR_obs[,'R_obs']), "SSB_obs"=SR_obs[,'SSB_obs'], "SPRF0_stock"=StockData[,'SPRF0'], "M_stock"=StockData[,'M'], "SSBmax_stock"=StockData[,'SSBmax'], "Rmax_stock"=StockData[,'Rmax'], "j_SR"=j_SR, "i_stock"=StockData[,'Stock_to_i']-1) )
  }
  if(Version%in%c("Taxon_v2_10_0")){
    Data = list("Options_vec"=c("n_obsfactors"=N_obsfactors,"n_factors"=N_factors,"invertTF"=FALSE,"b_type"=b_type,"Turn_off_taxonomy"=Turn_off_taxonomy), "Options"=c("minvar_obsfactors"=additional_variance[1],"minvar_factors"=additional_variance[2],"SD_b_stock"=SD_b_stock,"Pen_lowvar_lnRhat"=Pen_lowvar_lnRhat), "Cov_pz"=Cov_pz, "Y_ij"=as.matrix(Y_ij), "Missing_az"=Missing_az-1, "PC_gz"=PC_gz, "g_i"=g_i-1)
    Data = c(Data, list("Nobs"=Nobs, "Nstock"=Nstock, "Obs2Stock"=SR_obs[,'StockNum']-1, "AR_Index"=SR_obs[,'AR_Index'], "ln_R_obs"=log(SR_obs[,'R_obs']), "SSB_obs"=SR_obs[,'SSB_obs'], "SPRF0_stock"=StockData[,'SPRF0'], "M_stock"=StockData[,'M'], "SSBmax_stock"=StockData[,'SSBmax'], "Rmax_stock"=StockData[,'Rmax'], "j_SR"=j_SR, "i_stock"=StockData[,'Stock_to_i']-1) )
  }
  if(Version%in%c("Taxon_v2_11_0")){
    Data = list("Options_vec"=c("n_obsfactors"=N_obsfactors,"n_factors"=N_factors,"invertTF"=FALSE,"b_type"=b_type,"Turn_off_taxonomy"=Turn_off_taxonomy), "Options"=c("minvar_obsfactors"=additional_variance[1],"minvar_factors"=additional_variance[2],"SD_b_stock"=SD_b_stock,"Pen_lowvar_lnRhat"=Pen_lowvar_lnRhat, "lowerbound_MLSPS"=lowerbound_MLSPS), "Cov_pz"=Cov_pz, "Y_ij"=as.matrix(Y_ij), "Missing_az"=Missing_az-1, "PC_gz"=PC_gz, "g_i"=g_i-1)
    Data = c(Data, list("Nobs"=Nobs, "Nstock"=Nstock, "Obs2Stock"=SR_obs[,'StockNum']-1, "AR_Index"=SR_obs[,'AR_Index'], "ln_R_obs"=log(SR_obs[,'R_obs']), "SSB_obs"=SR_obs[,'SSB_obs'], "SPRF0_stock"=StockData[,'SPRF0'], "M_stock"=StockData[,'M'], "SSBmax_stock"=StockData[,'SSBmax'], "Rmax_stock"=StockData[,'Rmax'], "j_SR"=j_SR, "i_stock"=StockData[,'Stock_to_i']-1) )
  }
  if(Version%in%c("Taxon_v2_12_0")){
    Data = list("Options_vec"=c("n_obsfactors"=N_obsfactors,"n_factors"=N_factors,"invertTF"=FALSE,"b_type"=b_type,"Turn_off_taxonomy"=Turn_off_taxonomy,"Use_RAM_Mvalue_TF"=Use_RAM_Mvalue_TF), "Options"=c("minvar_obsfactors"=additional_variance[1],"minvar_factors"=additional_variance[2],"SD_b_stock"=SD_b_stock,"Pen_lowvar_lnRhat"=Pen_lowvar_lnRhat, "lowerbound_MLSPS"=lowerbound_MLSPS), "Cov_pz"=Cov_pz, "Y_ij"=as.matrix(Y_ij), "Missing_az"=Missing_az-1, "PC_gz"=PC_gz, "g_i"=g_i-1)
    Data = c(Data, list("Nobs"=Nobs, "Nstock"=Nstock, "Obs2Stock"=SR_obs[,'StockNum']-1, "AR_Index"=SR_obs[,'AR_Index'], "ln_R_obs"=log(SR_obs[,'R_obs']), "SSB_obs"=SR_obs[,'SSB_obs'], "SPRF0_stock"=StockData[,'SPRF0'], "M_stock"=StockData[,'M'], "SSBmax_stock"=StockData[,'SSBmax'], "Rmax_stock"=StockData[,'Rmax'], "j_SR"=j_SR, "j_logM"=j_logM, "i_stock"=StockData[,'Stock_to_i']-1) )
  }
  if(Version%in%c("Taxon_v2_13_0")){
    Data = list("Options_vec"=c("n_obsfactors"=N_obsfactors,"n_factors"=N_factors,"invertTF"=FALSE,"b_type"=b_type,"Turn_off_taxonomy"=Turn_off_taxonomy,"Use_RAM_Mvalue_TF"=Use_RAM_Mvalue_TF,"rho_option"=switch(rho_space,"natural"=0,"logit"=1,"logit_with_jacobian"=2)), "Options"=c("minvar_obsfactors"=additional_variance[1],"minvar_factors"=additional_variance[2],"SD_b_stock"=SD_b_stock,"Pen_lowvar_lnRhat"=Pen_lowvar_lnRhat, "lowerbound_MLSPS"=lowerbound_MLSPS), "Cov_pz"=Cov_pz, "Y_ij"=as.matrix(Y_ij), "Missing_az"=Missing_az-1, "PC_gz"=PC_gz, "g_i"=g_i-1)
    Data = c(Data, list("Nobs"=Nobs, "Nstock"=Nstock, "Obs2Stock"=SR_obs[,'StockNum']-1, "AR_Index"=SR_obs[,'AR_Index'], "ln_R_obs"=log(SR_obs[,'R_obs']), "SSB_obs"=SR_obs[,'SSB_obs'], "SPRF0_stock"=StockData[,'SPRF0'], "M_stock"=StockData[,'M'], "SSBmax_stock"=StockData[,'SSBmax'], "Rmax_stock"=StockData[,'Rmax'], "j_SR"=j_SR, "j_logM"=j_logM, "i_stock"=StockData[,'Stock_to_i']-1) )
  }
  if(Version%in%c("Taxon_v2_14_0")){
    Data = list("Options_vec"=c("n_obsfactors"=N_obsfactors,"n_factors"=N_factors,"invertTF"=FALSE,"b_type"=b_type,"Turn_off_taxonomy"=Turn_off_taxonomy,"Use_RAM_Mvalue_TF"=Use_RAM_Mvalue_TF,"rho_option"=switch(rho_space,"natural"=0,"logit"=1,"logit_with_jacobian"=2)), "Options"=c("minvar_obsfactors"=additional_variance[1],"minvar_factors"=additional_variance[2],"SD_b_stock"=SD_b_stock,"Pen_lowvar_lnRhat"=Pen_lowvar_lnRhat, "lowerbound_MLSPS"=lowerbound_MLSPS), "Cov_pz"=Cov_pz, "Y_ij"=as.matrix(Y_ij), "Missing_az"=Missing_az-1, "PC_gz"=PC_gz, "g_i"=g_i-1)
    Data = c(Data, list("Nobs"=Nobs, "Nstock"=Nstock, "Obs2Stock"=SR_obs[,'StockNum']-1, "AR_Index"=SR_obs[,'AR_Index'], "ln_R_obs"=log(SR_obs[,'R_obs']), "SSB_obs"=SR_obs[,'SSB_obs'], "PredTF_stock"=PredTF_stock, "SPRF0_stock"=StockData[,'SPRF0'], "M_stock"=StockData[,'M'], "SSBmax_stock"=StockData[,'SSBmax'], "Rmax_stock"=StockData[,'Rmax'], "j_SR"=j_SR, "j_logM"=j_logM, "i_stock"=StockData[,'Stock_to_i']-1) )
  }
  if(Version%in%c("Taxon_v2_15_0")){
    Data = list("Options_vec"=c("n_obsfactors"=N_obsfactors,"n_factors"=N_factors,"invertTF"=FALSE,"b_type"=b_type,"Turn_off_taxonomy"=Turn_off_taxonomy,"Use_RAM_Mvalue_TF"=Use_RAM_Mvalue_TF,"rho_option"=switch(rho_space,"natural"=0,"logit"=1,"logit_with_jacobian"=2)), "Options"=c("minvar_obsfactors"=additional_variance[1],"minvar_factors"=additional_variance[2],"SD_b_stock"=SD_b_stock,"Pen_lowvar_lnRhat"=Pen_lowvar_lnRhat, "lowerbound_MLSPS"=lowerbound_MLSPS), "Cov_pz"=Cov_pz, "Y_ij"=as.matrix(Y_ij), "Missing_az"=Missing_az-1, "PC_gz"=PC_gz, "g_i"=g_i-1, "group_j"=group_j)
    Data = c(Data, list("Nobs"=Nobs, "Nstock"=Nstock, "Obs2Stock"=SR_obs[,'StockNum']-1, "AR_Index"=SR_obs[,'AR_Index'], "ln_R_obs"=log(SR_obs[,'R_obs']), "SSB_obs"=SR_obs[,'SSB_obs'], "PredTF_stock"=PredTF_stock, "SPRF0_stock"=StockData[,'SPRF0'], "M_stock"=StockData[,'M'], "SSBmax_stock"=StockData[,'SSBmax'], "Rmax_stock"=StockData[,'Rmax'], "j_SR"=j_SR, "j_logM"=j_logM, "i_stock"=StockData[,'Stock_to_i']-1) )
  }
  if(Version%in%c("Taxon_v3_0_0")){
    Data = list("Options_vec"=c("n_obsfactors"=N_obsfactors,"n_factors"=N_factors,"invertTF"=FALSE,"b_type"=b_type,"Turn_off_taxonomy"=Turn_off_taxonomy,"Use_RAM_Mvalue_TF"=Use_RAM_Mvalue_TF,"rho_option"=switch(rho_space,"natural"=0,"logit"=1,"logit_with_jacobian"=2),"multinomial_for_factors"=multinomial_for_factors), "Options"=c("minvar_obsfactors"=additional_variance[1],"minvar_factors"=additional_variance[2],"SD_b_stock"=SD_b_stock,"Pen_lowvar_lnRhat"=Pen_lowvar_lnRhat, "lowerbound_MLSPS"=lowerbound_MLSPS), "RAM"=cbind(RAM[,1],RAM[,2],RAM[,3],RAM[,4]), "Y_ij"=as.matrix(Y_ij), "Missing_az"=Missing_az-1, "PC_gz"=PC_gz, "distance_g"=ParentChild_gz[,'EdgeLength'], "g_i"=g_i-1, "group_j"=group_j)
    Data = c(Data, list("Nobs"=Nobs, "Nstock"=Nstock, "Obs2Stock"=SR_obs[,'StockNum']-1, "AR_Index"=SR_obs[,'AR_Index'], "ln_R_obs"=log(SR_obs[,'R_obs']), "SSB_obs"=SR_obs[,'SSB_obs'], "PredTF_stock"=PredTF_stock, "SPRF0_stock"=StockData[,'SPRF0'], "M_stock"=StockData[,'M'], "SSBmax_stock"=StockData[,'SSBmax'], "Rmax_stock"=StockData[,'Rmax'], "j_SR"=j_SR, "j_logM"=j_logM, "i_stock"=StockData[,'Stock_to_i']-1) )
  }

  # Fix potential issues
  if( "j_logM"%in%names(Data) && length(Data$j_logM)==0 ){
    Data$j_logM = 0
    message("test")
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

    #
    if(!is.null(text)){
      L_z = rep( 0.1, max(RAM[,4]) )
    }else{
      L_z = rloadings( n_row=n_j, n_col=Data$Options_vec['n_factors'], mean=1, sd=0.1 )
    }

    ### Get informative starting values for RAM variables
    # Latent variables -- rho ~= 0.5
    if( N_obsfactors == 0 ){
      Y_a = vector()
    }else{
      Y_a = rnorm( nrow(Data$Missing_az), sd=0.1)
    }
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
    if(Version%in%c("Taxon_v2_15_0")){
      Params = list( "alpha_j"=alpha_j, "L_z"=rloadings(n_row=n_j,n_col=Data$Options_vec['n_factors'],mean=1,sd=0.1), "obsL_z"=rloadings(n_row=n_j,n_col=Data$Options_vec['n_obsfactors'],mean=1,sd=0.1), "L_logmult_col"=rep(0,ifelse(N_factors==0,1,abs(N_factors))), "obsL_logmult_col"=rep(0,ifelse(N_obsfactors==0,1,abs(N_obsfactors))), "cov_logmult_z"=rep(0,max(Data$PC_gz[,'ChildTaxon'])+1), "betainput_gj"=rmatrix(nrow=n_g,ncol=n_j), "Y_a"=Y_a )
      Params = c( Params, list("bparam_stock"=rep(0,ifelse(Nstock>0,Nstock,1)), "gamma_p"=rep(0,ifelse(n_p==0,1,n_p)), "theta_q"=0) )
    }
    if(Version%in%c("Taxon_v3_0_0")){
      Params = list( "alpha_j"=alpha_j, "L_z"=L_z, "obsL_z"=rloadings(n_row=n_j,n_col=Data$Options_vec['n_obsfactors'],mean=1,sd=0.1), "betainput_gj"=rmatrix(nrow=n_g,ncol=n_j), "Y_a"=Y_a )
      Params = c( Params, list("bparam_stock"=rep(0,ifelse(Nstock>0,Nstock,0))) )
    }
  }

  # Random
  if( Random[1]=="Generate" ){
    if(Version%in%"Taxon_v1_0_0"){
      Random = c("Y_a")
    }
    if(Version%in%c("Taxon_v2_14_0","Taxon_v2_13_0","Taxon_v2_12_0","Taxon_v2_11_0","Taxon_v2_10_0","Taxon_v2_9_0","Taxon_v2_8_0","Taxon_v2_7_0","Taxon_v2_6_0","Taxon_v2_5_0","Taxon_v2_4_0","Taxon_v2_3_0","Taxon_v2_2_0","Taxon_v2_1_0","Taxon_v2_0_0","Taxon_v1_2_0","Taxon_v1_1_0")){
      Random = c("Y_a", "beta_gj")
    }
    if(Version%in%c("Taxon_v3_0_0","Taxon_v2_15_0")){
      Random = c("Y_a", "betainput_gj")
    }
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
    if( "L_logmult_col" %in% names(Params) ) Map[["L_logmult_col"]] = factor(rep(NA,length(Params[["L_logmult_col"]])))
    if( "obsL_logmult_col" %in% names(Params) ) Map[["obsL_logmult_col"]] = factor(rep(NA,length(Params[["obsL_logmult_col"]])))
    if( length(j_SR)==4 ){
      if("ln_b_stock" %in% names(Params)) Map[["ln_b_stock"]] = factor(rep(NA,length(Params[["ln_b_stock"]])))
      if("bparam_stock" %in% names(Params)) Map[["bparam_stock"]] = factor(rep(NA,length(Params[["bparam_stock"]])))
    }
    if( n_p==0 ){
      if( "gamma_p" %in% names(Params) ) Map[["gamma_p"]] = factor(NA)
    }
    if( "theta_q" %in% names(Params) ) Map[["theta_q"]] = factor(NA)
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
  }

  if( (min_replicate_measurements>0) & (N_obsfactors!=0) ){
    stop("`min_replicate_measurements` must be zero except when `N_obsfactors=0`")
  }else{
    #Genus_species = paste0( Z_ik[,'Genus'], "_", Z_ik[,'Species'] )
    Genus_species = rownames(Y_ij)
    MaxNum_j = apply( Y_ij, MARGIN=2, FUN=function(x){max(tapply(x,INDEX=Genus_species,FUN=function(y){sum(!is.na(y))}))} )
    Map$obsL_z = seq_len(length(Params$obsL_z))
    Map$obsL_z[MaxNum_j <= min_replicate_measurements] = NA
    Map$obsL_z = factor(Map$obsL_z)
    Params$obsL_z[MaxNum_j <= min_replicate_measurements] = 0.01
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
  # SEE https://github.com/kaskr/adcomp/issues/321 for flags argument
  if( packageVersion("TMB") >= numeric_version("1.8.0") ){
    compile( paste0(Version,".cpp"), flags="-Wno-ignored-attributes -O2 -mfpmath=sse -msse2 -mstackrealign", framework="TMBad" )
  }else{
    compile( paste0(Version,".cpp"), flags="-Wno-ignored-attributes -O2 -mfpmath=sse -msse2 -mstackrealign" )
  }

  # Build
  dyn.load( paste0(RunDir,"/",TMB::dynlib(Version)) )          #
  #return( list(Params=Params, Data=Data, Map=Map, Random=Random) )
  Obj = MakeADFun( data=Data, parameters=Params, DLL=Version, map=Map, random=Random )
  Obj$env$beSilent()

  # Print number of parameters
  ThorsonUtilities::list_parameters( Obj )

  #
  if( run_model==FALSE ){
    Return = list("Data"=Data, "Params"=Params, "Random"=Random, "Map"=Map, "Obj"=Obj, "SEM_model"=SEM_model, "RAM"=RAM, "text"=text, "tree"=tree )
    return(Return)
  }

  # Debugging option
  if( debug_mode==TRUE ){
    #on.exit( assign("ParHat",Obj$env$parList(),envir=.GlobalEnv), add=TRUE )
    on.exit( assign("Obj",Obj,envir=.GlobalEnv), add=TRUE )
  }

  # Print to screen
  if( verbose==TRUE ){
    #cat( "Number of fixed effects:")
    #print( table(names(Obj$par)) )
    #cat( "Number of random effects:")
    #print( table(names(Obj$env$last.par[Obj$env$random])) )
  }

  # Optimize                         #  , startpar=opt$par[-grep("alpha",names(opt$par))]
  # JointPrecision is used below, and is too big to invert whole;  must have getReportCovariance=TRUE to get JointPrecision
  Opt = TMBhelper::fit_tmb( obj = Obj,
                            savedir = RunDir,
                            getJointPrecision = TRUE,
                            getReportCovariance = TRUE,
                            control = list(eval.max=10000, iter.max=10000, trace=1),
                            ... )

  # Debugging option
  if( debug_mode==TRUE ){
    #on.exit( assign("ParHat",Obj$env$parList(),envir=.GlobalEnv), add=TRUE )
    on.exit( assign("Opt",Opt,envir=.GlobalEnv), add=TRUE )
  }

  #
  Report = Obj$report()
  colnames(Report$Ycomplete_ij) = colnames(Report$beta_gj) = colnames(Y_ij)
  if("betainput_gj" %in% names(Report)) colnames(Report$betainput_gj) = colnames(Y_ij)

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
  #ParHat_SE = as.list( Opt$SD, what="Std" )
  colnames(Report$Ycomplete_ij) = colnames(Y_ij)
  if("beta_gj" %in% names(ParHat)) colnames(ParHat$beta_gj) = colnames(Y_ij)
  if("betainput_gj" %in% names(ParHat)) colnames(ParHat$betainput_gj) = colnames(Y_ij)
  #colnames(ParHat_SE$beta_gj) = colnames(Y_ij)
  #dyn.unload( paste0(RunDir,"/",TMB::dynlib(Version)) )          #

  # Return stuff
  Return = list( "N_factors" = N_factors,
                 "N_obsfactors" = N_obsfactors,
                 "Use_REML" = Use_REML,
                 "ParentChild_gz" = ParentChild_gz,
                 "ParHat" = ParHat,
                 "g_i" = g_i,
                 "Y_ij" = Y_ij,
                 "Z_ik" = Z_ik,
                 "SR_obs" = SR_obs,
                 "StockData" = StockData,
                 "Obj" = Obj,
                 "Opt" = Opt,
                 "Report" = Report,
                 #"ParHat_SE" = ParHat_SE,
                 "obsCov_jj" = Report$obsCov_jj,
                 "tree" = tree )
  dimnames(Return$obsCov_jj) = list(colnames(Y_ij),colnames(Y_ij))
  if("Cov_jj" %in% names(Report)){
    Return = c(Return, list("Cov_jj"=Report$Cov_jj))
  }
  #if("cov_logmult_z" %in% names(ParHat)){
  #  Return = c(Return, list("Cov_jjz"=Report$Cov_jj %o% exp(ParHat$cov_logmult_z)))
  #  dimnames(Return$Cov_jjz) = list(colnames(Y_ij),colnames(Y_ij),colnames(Z_ik))
  #}
  if( !is.null(text) ){
    Return = c(Return, list("SEM_model"=SEM_model, "RAM"=RAM, "text"=text))
  }

  ####################
  # Interpret results
  ####################

  # Get variable names
  if(Version%in%c("Taxon_v3_0_0","Taxon_v2_15_0")){
    varname = "betainput_gj"
  }else{
    varname = "beta_gj"
  }
  # Extract and invert
  VarNames = Predictive_distribution( mean_vec = Y_ij[1,],
                                      process_cov = NULL,
                                      group_j = group_j,
                                      obs_cov = NULL,
                                      check_names = TRUE,
                                      include_r = include_r )
  n_v = length(VarNames)

  if( extract_covariance==TRUE ){
    ### Approximate joint precision
    message( "Extracting predictive variance of life-history parameters for each taxon: ", Sys.time())
    # Extract predictive covariance for species-specific traits (necessary to do separately for rows and columns)
    #Prec_zz = Opt$SD$jointPrecision[ , grep("beta_gj",colnames(Opt$SD$jointPrecision)) ]
    #  Prec_zz = Prec_zz[ grep("beta_gj",rownames(Opt$SD$jointPrecision)), ]
    # Predict
    # Don't form full u_zr to avoid memory load
    #u_zr = rmvnorm_prec( mu=Obj$env$last.par.best, prec=Opt$SD$jointPrecision, n_sims=n_sims )
    u_zr = rmvnorm_prec( mu = Obj$env$last.par.best,
                         prec = Opt$SD$jointPrecision,
                         n_sims = n_sims,
                         varnames = varname,
                         n_batches = n_batches )
    Prec_gjj = PartialCorr_gjj = array(NA, dim=c(n_g,n_j,n_j), dimnames=list(ParentChild_gz[,'ChildName'],colnames(Y_ij),colnames(Y_ij)) )
    Prec_gvv = PartialCorr_gvv = Corr_gvv = Cov_gvv = array(NA, dim=c(n_g,n_v,n_v), dimnames=list(ParentChild_gz[,'ChildName'],VarNames,VarNames) )
    beta_gv = array(NA, dim=c(n_g,n_v), dimnames=list(ParentChild_gz[,'ChildName'],VarNames) )
    for( gI in 1:n_g ){
      # Method with full precision
      Indices = gI + ( seq(1,n_g*n_j,by=n_g)-1 )
      Samp_rj = t(u_zr[Indices,])
      #Samp_rj = matrix(NA, nrow=n_sims, ncol=n_j)
      #for( rI in 1:n_sims ){
      #  u_z = rmvnorm_prec( mu=Obj$env$last.par.best, prec=Opt$SD$jointPrecision, n_sims=1, seed=1+rI )[,1]
      #  Samp_rj[rI,] = u_z[ grep("beta_gj",colnames(Opt$SD$jointPrecision))[Indices] ]
      #}
      colnames(Samp_rj) = colnames(Y_ij)
      Pred = Predictive_distribution( mean_vec = Report[[varname]][gI,],
                                      Samp_rj = Samp_rj,
                                      group_j = group_j,
                                      include_obscov = FALSE,
                                      check_bounds = FALSE,
                                      include_r = include_r,
                                      lowerbound_MLSPS = lowerbound_MLSPS,
                                      rho_option = switch(rho_space, "natural"=0, "logit"=1, "logit_with_jacobian"=2) )
      beta_gv[gI,] = Pred$pred_mean
      Cov_gvv[gI,,] = Pred$pred_cov
      Corr_gvv[gI,,] = cov2cor( Cov_gvv[gI,,] )
      if( (gI%%1000) == 0 ) message( "Finished processing predictive variance for ", gI, " of ",n_g," taxa: ", Sys.time() )
    }
    if("Cov_jj" %in% names(Report)){
      Return = c(Return, list("Cov_gvv"=Cov_gvv, "beta_gv"=beta_gv))
    }
  }else{
    beta_gv = array(NA, dim=c(n_g,n_v), dimnames=list(ParentChild_gz[,'ChildName'],VarNames) )
    for( gI in 1:n_g ){
      Pred = Predictive_distribution( mean_vec = Report[[varname]][gI,],
                                      Samp_rj = matrix(Report[[varname]][gI,],nrow=1),
                                      group_j = group_j,
                                      include_obscov = FALSE,
                                      check_bounds = FALSE,
                                      include_r = include_r,
                                      lowerbound_MLSPS = lowerbound_MLSPS,
                                      rho_option = switch(rho_space, "natural"=0, "logit"=1, "logit_with_jacobian"=2) )
      beta_gv[gI,] = Pred$pred_mean
      if( (gI%%1000) == 0 ) message( "Finished processing plug-in mean estimator for ", gI, " of ",nrow(beta_gv)," taxa: ", Sys.time() )
    }
    Return = c(Return, list("beta_gv"=beta_gv))
  }

  Return[["total_runtime"]] = Sys.time() - start_time
  class(Return) = "FishLife"
  return( Return )
}

#' Print parameter estimates and standard errors.
#'
#' @title Print parameter estimates
#' @param x Output from \code{\link{Fit_model}}
#' @param ... Not used
#' @return NULL
#' @method print FishLife
#' @export
print.FishLife <- function(x, ...)
{
  cat("Fit_model(.) result\n")
  if( "Opt" %in% names(x) ){
    print( x$Opt )
  }else{
    cat("`parameter_estimates` not available in `Fit_model`\n")
  }
  invisible(x$Opt)
}

#' Extract path coefficients.
#'
#' @title Extract path coefficients
#' @param x Output from \code{\link{Fit_model}}
#' @param ... Not used
#' @return NULL
#' @method coef FishLife
#' @export
coef.FishLife = function( x ){
  beta_z = x$Opt$par[names(x$Opt$par)=="L_z"]
  RAM = x$Obj$env$data$RAM
  if(nrow(RAM) != nrow(x$SEM_model)) stop("Check assumptions")
  # Report variances
  for( i in which(RAM[,1]==2) ){
    if( RAM[i,3] == RAM[i,2] ){
      beta_z[i] = beta_z[i]^2
    }
  }
  SEM_params = beta_z[ifelse(RAM[,4]==0, NA, RAM[,4])]
  SEM_params = ifelse( is.na(SEM_params), as.numeric(x$SEM_model[,3]), SEM_params )
  return( data.frame(Path=x$SEM_model[,1], Parameter=x$SEM_model[,2], Estimate=SEM_params ) )
}
