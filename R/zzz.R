
# .onLoad <- function(libname, pkgname) {
# }

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("###########################################################################################")
  packageStartupMessage("Loading package FishLife, developed by James Thorson for the National Marine Fisheries Service")
  packageStartupMessage("For details and citation guidance, please see http://github.com/James-Thorson-NOAA/FishLife/")
  packageStartupMessage("###########################################################################################")

  #if( !"ThorsonUtilities" %in% utils::installed.packages()[,1] ){
  #  packageStartupMessage("Installing package: ThorsonUtilities...")
  #  devtools::install_github("james-thorson/utilities")
  #}

  #if( !"TMBhelper" %in% utils::installed.packages()[,1] ){
  #  packageStartupMessage("Installing package: TMBhelper...")
  #  devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")
  #}
}

