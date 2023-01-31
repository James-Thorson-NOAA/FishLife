
# .onLoad <- function(libname, pkgname) {
# }

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("###########################################################################################")
  packageStartupMessage("Loading package FishLife, developed by James Thorson for the National Marine Fisheries Service")
  packageStartupMessage("For details and citation guidance, please see http://github.com/James-Thorson-NOAA/FishLife/")
  packageStartupMessage("###########################################################################################")

  # Check rfishbase version
  if( utils::packageVersion("rfishbase") >= numeric_version("4.0.0") ){
    packageStartupMessage("Please install an earlier version of rfishbase, e.g. using by running:")
    packageStartupMessage("  remotes::install_github( 'ropensci/rfishbase@fb-21.06', force=TRUE )")
  }
  
  #if( !"ThorsonUtilities" %in% utils::installed.packages()[,1] ){
  #  packageStartupMessage("Installing package: ThorsonUtilities...")
  #  devtools::install_github("james-thorson/utilities")
  #}

  #if( !"TMBhelper" %in% utils::installed.packages()[,1] ){
  #  packageStartupMessage("Installing package: TMBhelper...")
  #  devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")
  #}
}

