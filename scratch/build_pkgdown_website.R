# Author: Kevin See
# Purpose: Build a website for this package using pkgdown
# Created: 5/17/2021
# Last Modified: 5/17/2021
# Notes: Based on instructions found here: https://pkgdown.r-lib.org/index.html

setwd(R'(C:\Users\James.Thorson\Desktop\Git\FishLife)')

#-----------------------------------------------------------------
# load needed libraries
library(pkgdown)


# Only needed once
if( FALSE ){
  # set up to automatically publish pkgdown site to GitHub
  # usethis::create_github_token()
  # gitcreds::gitcreds_set(url = "https://github.com")
  usethis::use_pkgdown_github_pages()

  # Run once to configure your package to use pkgdown
  usethis::use_pkgdown()

  # check that _pkgdown.yml looks good
  pkgdown::check_pkgdown()
}

# Run to build the website
pkgdown::build_site( examples=TRUE )
# build_articles( lazy = FALSE )

# to look at the site
pkgdown::preview_site()

#-----------------------------------------------------------------
# deploy site to gh-pages branch on GitHub
pkgdown::deploy_to_branch()
