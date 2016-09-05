
#.onLoad <- function(libname, pkgname) {
#}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("###########################################################################################")
  packageStartupMessage("Loading package FishTraits, developed by James Thorson for the Northwest Fisheries Science Center")
  packageStartupMessage("For details and citation guidance, please see http://github.com/james-thorson/FishTraits/")

  if( !"ThorsonUtilities" %in% utils::installed.packages()[,1] ){
    packageStartupMessage("Installing package: ThorsonUtilities...")
    devtools::install_github("james-thorson/utilities")
  }

  if( !"TMBhelper" %in% utils::installed.packages()[,1] ){
    packageStartupMessage("Installing package: TMBhelper...")
    devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")
  }

  # Load existing dataabase
  packageStartupMessage("Also loading results from model run: list 'Estimate_database'")
  results_path <- system.file("extdata", package="FishTraits")
  load( file.path(results_path,"Estimate_database.RData") )
  packageStartupMessage("###########################################################################################")
}
