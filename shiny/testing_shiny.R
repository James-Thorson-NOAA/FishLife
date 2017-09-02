

#Details of how to deploy the package at:
# http://shiny.rstudio.com/articles/shinyapps.html

#Details of other ways to host a Shiny App online:
# http://shiny.rstudio.com/tutorial/lesson7/

#Sign-in for hosting on shinyapps.io:
# http://www.shinyapps.io/

#Shiny cheat-sheet
# https://www.rstudio.com/wp-content/uploads/2016/01/shiny-cheatsheet.pdf

# Shiny examples
# https://github.com/rstudio/shiny-examples/tree/master/063-superzip-example

#renderUI examples
# https://www.youtube.com/watch?v=JUop-YfRAuw
# https://github.com/aagarw30/R-Shinyapp-Tutorial/blob/master/renderUI%20Example1/server.R

#Tutorial
# http://rstudio.github.io/shiny/tutorial/

#sliderInput animation example
# http://shiny.rstudio.com/gallery/image-output.html

#Button to turn off all check-boxes
# http://stackoverflow.com/questions/35043641/selecting-all-the-check-boxes-at-once-in-shiny

#Action buttons side-by-side
# http://stackoverflow.com/questions/20637248/shiny-4-small-textinput-boxes-side-by-side

# REAL VERSION -- DO NOT USE FOR TESTING
setwd( "C:/Users/James.Thorson/Desktop/Project_git/FishLife/shiny" )

# DEVELOPMENT VERSION -- USE FOR TESTING
setwd( "C:/Users/James.Thorson/Desktop/Project_git/FishLife/shiny/testing" )

##############
# Test locally
##############
library(shiny)
runApp()
#library(FishLife)
#Search_species( Class="actinopterygii", Order="acipenseriformes", Family="acipenseridae", Genus="acipenser", Species="baerii", add_ancestors=TRUE, ParentChild_gz=Estimate_database$ParentChild_gz )
#Search_species( Genus="acipenser", Species="baerii", add_ancestors=TRUE, ParentChild_gz=Estimate_database$ParentChild_gz )

##############
# Push to server
# Apparently, must ...
#  ... have all packages installed locally on the R version being used
#  .. currently working with Revolution Open R v3.2.2
##############
library(rsconnect)
deployApp( appName="FishLife" )
