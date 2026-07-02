#' Database of size, growth, maturity, and mortality parameters
#'
#' Output from `Fit_model` applied to database scraped from www.FishBase.org
#' using `rfishbase`
#'
#' @details
#' Length at maturity (Lmat) and total length (Loo) measurements are only included
#' when measured as Total Length (TL), and excluding measurements using
#' standard length (SL) or fork length (FL).
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
#' @details
#'
#' **Model output**
#'
#' \code{FishBase_and_Morphometrics$beta_gv} contains model output.
#' Each row corresponds to a taxon or ancestral node, and each column corresponds
#' to a predicted trait.  Note that log always refers to natural log.
#' Several columns are labeled base, and this corresponds to the base level for
#' each categorical trait.  Values are calculated as the mean when sampling from the
#' predictive distribution for a given trait, and this allows us to approximate the
#' back-transformed mean value.
#' Columns (in respective order) are defined as:
#'
#' | column | trait | source |
#' | --- | --- | --- |
#' | log(age_max)| log of maximum age (years) | FishBase |
#' | trophic level | trophic level (dimensionless), where 1 is primary producers, etc. | FishBase |
#' | log(aspect_ratio) | log of caudal fin height and length divided by area (dimensionless) | FishBase |
#' | log(fecundity) | Annual eggs produced (number/year) | FishBase |
#' | log(growth_coefficient) | von Bertalannffy growth coefficient (year^-1) | FishBase |
#' | temperature | average temperature from portion of population sampled (celcius) | FishBase |
#' | log(length_max) | log of maximum length (cm) | FishBase |
#' | log(length_infinity) | log of von Bertalanffy asymptotic maximum length (cm) | FishBase |
#' | log(length_maturity) | log length at 50% maturity (cm) | FishBase |
#' | log(age_maturity) | log age at 50% sexual maturity (years) | FishBase |
#' | log(natural mortality) | log of natural mortality rate M (year^-1) | FishBase |
#' | log(weight_infinity) | log of asymptotic maximum weight (grams) | FishBase |
#' | log(max_body_depth) | log maximum body depth (cm) | Morphometrics |
#' | log(max_body_width) | log maximum body width (cm) | Morphometrics |
#' | log(lower_jaw_length) | log length of lower jaw (cm) | Morphometrics |
#' | log(min_caudal_pedoncule_depth) | log depth of caudal pedoncule (connecting caudal fin to body) | Morphometrics |
#' | log(offspring_size) | log size of offspring (kg) | Morphometrics |
#' | base | Probability of nonguarder level for spawning type category | FishBase |
#' | spawning_typeguarders | Probability of guarder level for spawning type category | FishBase |
#' | spawning_typebearers | Probability of live bearer level for spawning type category | FishBase |
#' | base | Probability of demersal type for habitat category | FishBase |
#' | habitatbathymetric | Probability of bathmetric type for habitat category | FishBase |
#' | habitatbenthopelagic | Probability of benthopelagic type for habitat category | FishBase |
#' | habitatreefassociated | Probability of reef-associated type for habitat category | FishBase |
#' | habitatpelagic | Probability of pelagic type for habitat category | FishBase |
#' | base | Probability of generalist type for feeding mode category | FishBase |
#' | feeding_modemacrofauna | Probability of macrofauna type for feeding mode category | FishBase |
#' | feeding_modeplanktivorous_or_other | Probability of planktivorous or other type for feeding mode category | FishBase |
#' | base | Probability of fusiform or normal body type for body shape category | FishBase |
#' | body_shapeelongated | Probability of elongated type for body shape category | FishBase |
#' | body_shapeshort_and_or_deep | Probability of short and/or deep type for body shape category | FishBase |
#' | body_shapeeellike | Probability of eel-like type for body shape category | FishBase  |
#' | body_shapeother | Probability of other type for body shape category | FishBase |
#'
#'
#' @return
#' A tagged list containing data and predictions
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
#'
#'
"FishBase_and_Morphometrics"

