#' Database of size, growth, maturity, and mortality parameters
#'
#' Output from `Fit_model` applied to database scraped from www.FishBase.org
#' using `rfishbase`
#'
#' @format A tagged list containing data and predictions
#' \describe{
#'   \item{N_factors}{Number of factors used for evolution in life-history model}
#'   \item{N_obsfactors}{Number of factors used for measurent-error in life-history model}
#'   \item{beta_gv}{Predictive mean (in transformed space) among traits for every taxon in tree}
#'   \item{Cov_gvv}{Covariance among traits for every taxon in tree}
#'   \item{Use_REML}{Boolean, whether REML was used for model}
#'   \item{ParentChild_gz}{Record of taxonomic tree}
#'   \item{ParHat}{Parameter estimates and predictions}
#'   \item{g_i}{Associates every observation with a level of the taxonomic tree}
#'   \item{Y_ij}{Life-history parameters from FishBase}
#'   \item{Z_ik}{Taxonomy for each datum}
#'   ...
#' }
"FishBase"

#' Database of stock-recruit, population-dynamics, size, growth, maturity, and mortality parameters
#'
#' Output from `Fit_model` applied to database scraped from www.FishBase.org
#' using `rfishbase` as well as RAM Legacy stock-recruit database
#'
#' @format A tagged list containing data and predictions
#' \describe{
#'   \item{N_factors}{Number of factors used for evolution in life-history model}
#'   \item{N_obsfactors}{Number of factors used for measurent-error in life-history model}
#'   \item{beta_gv}{Predictive mean (in transformed space) among traits for every taxon in tree}
#'   \item{Cov_gvv}{Covariance among traits for every taxon in tree}
#'   \item{Use_REML}{Boolean, whether REML was used for model}
#'   \item{ParentChild_gz}{Record of taxonomic tree}
#'   \item{ParHat}{Parameter estimates and predictions}
#'   \item{g_i}{Associates every observation with a level of the taxonomic tree}
#'   \item{Y_ij}{Life-history parameters from FishBase}
#'   \item{Z_ik}{Taxonomy for each datum}
#'   \item{SR_obs}{Stock-recruit records from RAM Legacy stock-recruit database}
#'   \item{StockData}{Auxiliary information for every stock with stock-recruit information}
#'   ...
#' }
"FishBase_and_RAM"

#' Database of FishBase (population-dynamics, size, growth, maturity, mortality parameters),
#'          habitat, trophic, reproductive, and morphometric traits
#'
#' Output from `Fit_model` applied to database scraped from www.FishBase.org
#' using `rfishbase` as well as many other variables
#'
#' @format A tagged list containing data and predictions
#' \describe{
#'   \item{text}{text file specifying the SEM structure}
#'   \item{N_obsfactors}{Number of factors used for measurent-error in life-history model}
#'   \item{beta_gv}{Predictive mean (in transformed space) among traits for every taxon in tree}
#'   \item{Cov_gvv}{Covariance among traits for every taxon in tree}
#'   \item{Use_REML}{Boolean, whether REML was used for model}
#'   \item{ParentChild_gz}{Record of taxonomic tree}
#'   \item{ParHat}{Parameter estimates and predictions}
#'   \item{g_i}{Associates every observation with a level of the taxonomic tree}
#'   \item{Y_ij}{Life-history parameters from FishBase}
#'   \item{Z_ik}{Taxonomy for each datum}
#'   \item{SR_obs}{Stock-recruit records from RAM Legacy stock-recruit database}
#'   \item{StockData}{Auxiliary information for every stock with stock-recruit information}
#'   \item{RAM}{Matrix specifying SEM parameters}
#'   \item{SEM_model}{Description of SEM path diagram}
#'   \item{tree}{Evolutionary tree in phylo format, for use in plotting, with order matching ParentChild_gz}
#'   ...
#' }
"FishBase_and_Morphometrics"

