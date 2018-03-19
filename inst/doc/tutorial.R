## ---- echo=TRUE, message=FALSE-------------------------------------------
# Install and load package
devtools::install_github("james-thorson/FishLife")
library( FishLife )

## ---- echo=TRUE, message=FALSE, fig.width=8, fig.height=8----------------
# Get basic plot for Lutjanus campechanus
Plot_taxa( Search_species(Genus="Lutjanus",Species="campechanus")$match_taxonomy, mfrow=c(2,2) )

## ---- echo=TRUE, message=FALSE, fig.width=8, fig.height=8----------------
# Get basic plot for Sebastes cortezi
Predict = Plot_taxa( Search_species(Genus="Sebastes",Species="cortezi")$match_taxonomy, mfrow=c(2,2) )

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(Predict[[1]]$Mean_pred, digits=3)

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(Predict[[1]]$Cov_pred, digits=3)

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(c(exp(Predict[[1]]$Mean_pred[-8]),Predict[[1]]$Mean_pred['Temperature']), digits=3)

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(c(exp(Predict[[1]]$Mean_pred[-8]+0.5*diag(Predict[[1]]$Cov_pred)[-8]),Predict[[1]]$Mean_pred['Temperature']), digits=3)

## ---- echo=TRUE, message=FALSE, fig.width=8, fig.height=8----------------
# Get basic plot for Family Scombridae 
Plot_taxa( Search_species(Family="Scombridae")$match_taxonomy, mfrow=c(2,2) )

## ---- echo=TRUE, message=FALSE, fig.width=8, fig.height=8----------------
# Compare two species
Taxa = c( Search_species(Genus="Oncorhynchus",Species="mykiss",add_ancestors=FALSE)$match_taxonomy,
  Search_species(Genus="Salmo",Species="Trutta",add_ancestors=FALSE)$match_taxonomy )
Plot_taxa( Taxa, mfrow=c(2,2) )

## ---- echo=TRUE, eval=FALSE, message=FALSE-------------------------------
#  # Load TMB
#  library( TMB )
#  
#  # Re-run results with a different model configuration
#  Estimate_database = Fit_model( N_factors=-3, N_obsfactors=-3, Use_REML=TRUE)

## ---- echo=TRUE, message=FALSE-------------------------------------------
Ynew_ij = matrix( c("Loo"=log(40),"K"=NA,"Winfinity"=NA,"tmax"=NA,"tm"=NA,"M"=NA,"Lm"=NA,"Temperature"=NA), nrow=1)

## ---- echo=TRUE, message=FALSE, warning=FALSE----------------------------
library(TMB)
Update = Update_prediction( Taxon=Search_species(Genus="Sebastes",Species="cortezi",add_ancestors=FALSE)$match_taxonomy, Ynew_ij=Ynew_ij)

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(Update$updateMean_j, digits=3)

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(Update$updateCov_j, digits=3)

