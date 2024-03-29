% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Fit_model.R
\name{Fit_model}
\alias{Fit_model}
\title{Fit and predict fish traits}
\usage{
Fit_model(
  text = NULL,
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
  group_j = 1:ncol(Y_ij) - 1,
  Version = "Taxon_v3_0_0",
  TmbDir = system.file("executables", package = "FishLife"),
  RunDir = tempfile(pattern = "run_", tmpdir = tempdir(), fileext = "/"),
  verbose = FALSE,
  debug_mode = FALSE,
  j_SR = ncol(Y_ij) - 3:1,
  additional_variance = c(0, 0),
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
  ...
)
}
\arguments{
\item{text}{structural equation model structure, passed to either \code{\link[sem]{specifyModel}}
and then parsed to control the set of path coefficients and variance-covariance parameters}

\item{Database}{Whether to use results for both adult and stock-recruit parameters, \code{Database="FishBase_and_RAM"}, or just adult parameters, \code{Database="FishBase"}}

\item{tree}{phylogenetic structure, using class \code{\link[ape]{as.phylo}}.  If \code{tree=NULL}
then argument \code{Z_ik} is instead converted to a tree using 
\code{ape::as.phylo( ~Class/Order/Family/Genus/GenusSpecies, ... )} to convert it to a tree}

\item{Y_ij}{a data frame of trait values (perhaps log-scaled) with rows for records, and tagged-columns for traits}

\item{Z_ik}{a data frame of taxonomic classification for each row of \code{Y_ij}}

\item{Use_REML}{OPTIONAL boolean whether to use maximum marginal likelihood or restricted maximum likelihood (termed "REML")}

\item{N_factors}{Number of factors in decomposition of covariance for random-walk along evolutionary tree (0 means a diagonal but unequal covariance; negative is the sum of a factor decomposition and a diagonal-but-unequal covariance)}

\item{N_obsfactors}{Number of factors in decomposotion of observation covariance (same format as \code{N_obsfactors})}

\item{min_replicate_measurements}{specified threshold for the number of measurements for a given continuous traits,
where any continuous trait having fewer replicated measurements for at least one taxon will have the 
measurement variance fixed at an arbitrarily low value, such that estimated traits are forced
to approach the unreplicated measurements}

\item{SR_obs}{Stock-recruit records from RAM Legacy stock-recruit database}

\item{StockData}{Auxiliary information for every stock with stock-recruit information}

\item{Version}{TMB version number}

\item{TmbDir}{Directory containing pre-compiled TMB script}

\item{RunDir}{Directory to use when compiling and running TMB script (different to avoid problems with read-write restrictions)}

\item{verbose}{Boolean whether to print diagnostics to terminal}

\item{run_model}{Boolean indicating whether to run the model or just return the built TMB object}

\item{Params}{optional list of parameter estimates to use as starting values (Default \code{Params="Generate"} starts from random values)}

\item{...}{other paramers passed to \code{\link[TMBhelper]{fit_tmb}}}

\item{Process_cov}{Whether process-error covariance is equal or differs multiplicatively for different taxonomic levels (Options:  "Equal" or "Unequal")}
}
\value{
Tagged list containing objects from FishLife run (first 9 slots constitute list 'Estimate_database' for archiving results)
\describe{
  \item{tree}{The phylogenetic tree used for analysis, whether inputted or generated from taxonomy based on \code{Z_ik}}
  \item{SEM_model}{The phylogenetic tree used for analysis, whether inputted or generated from taxonomy based on \code{Z_ik}}
  \item{N_factors}{Number of factors used for evolution in life-history model}
  \item{N_obsfactors}{Number of factors used for measurent-error in life-history model}
  \item{Use_REML}{Boolean, whether REML was used for model}
  \item{Cov_gjj}{Covariance among traits for every taxon in tree}
  \item{ParentChild_gz}{Record of taxonomic tree}
  \item{ParHat}{Parameter estimates and predictions}
  \item{g_i}{Associates every observation with a level of the taxonomic tree}
  \item{Y_ij}{Raw data}
  \item{Z_ik}{Taxonomy for each datum}
  \item{Obj}{The built TMB object}
  \item{Opt}{Output from optimization}
  \item{Report}{tagged list of report-file from TMB}
  \item{ParHat_SE}{Estimated/predicted standard errors for fixed/random effects}
}
}
\description{
\code{Fit_model} estimates parameters and predicts values from a multivariate random-walk model for fish traits
}
\examples{
\dontrun{
# Load data set
library(FishLife)

# Simulate data
tree = ape::rtree(n=100)
xfit = x = 1 + 0.3 * phylolm::rTrait(n = 1, phy=tree)
yfit = y = 2 + 1*x + 0.3 * phylolm::rTrait(n = 1, phy=tree)
drop_x = sample( 1:ape::Ntip(tree), replace=FALSE, size=round(ape::Ntip(tree)*0.3) )
xfit[drop_x] = NA
drop_y = sample( 1:ape::Ntip(tree), replace=FALSE, size=round(ape::Ntip(tree)*0.3) )
yfit[drop_y] = NA

# Fit model
Fit = Fit_model(
  text = "x -> y, p",
  Database = NULL,
  Use_REML = FALSE,
  Y_ij = data.frame( x=xfit, y=yfit ),
  tree = tree,
  min_replicate_measurements = Inf )

# S3-defaults
print(Fit)
coef(Fit)

# Convet and plot using sem
mysem = as(Fit,"sem")
sem::pathDiagram( model = mysem,
                  style = "traditional",
                  edge.labels = "values" )
myplot = semPlot::semPlotModel( as(Fit,"sem") )
semPlot::semPaths( myplot,
                   nodeLabels = myplot@Vars$name )

# Convert and plot using phylobase / phylosignal
library(phylobase)
plot( as(Fit,"phylo4d") )
}

}
