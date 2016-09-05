
#.onLoad <- function(libname, pkgname) {
#}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("###########################################################################################")
  packageStartupMessage("Loading package FishTraits, developed by James Thorson for the Northwest Fisheries Science Center")
  packageStartupMessage("For details and citation guidance, please see http://github.com/james-thorson/FishTraits/")
  packageStartupMessage("###########################################################################################")

  if( !"ThorsonUtilities" %in% utils::installed.packages()[,1] ){
    packageStartupMessage("Installing package: ThorsonUtilities...")
    devtools::install_github("james-thorson/utilities")
  }

  if( !"TMBhelper" %in% utils::installed.packages()[,1] ){
    packageStartupMessage("Installing package: TMBhelper...")
    devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")
  }
}

#' Load previous results
#'
#' \code{Load_previous_results} loads previous results distributed with package
#'
#' @param results_dir Directory containing object \code{Estimate_database.RData}, archiving previous results

#' @export
Load_previous_results = function( results_dir=system.file("extdata",package="FishTraits") ){
  # Load existing dataabase
  message("###########################################################################################")
  message("Loading results from model run: list 'Estimate_database'")
  load( file.path(results_dir,"Estimate_database.RData") )
  message("###########################################################################################")
  return( Estimate_database )
}
