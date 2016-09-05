# FishTraits
Estimate fish traits for all marine fish species globally

# Example usage

### Load the package
```R
# Install and load package
devtools::install_github("james-thorson/FishTraits")
library( FishTraits )
Estimate_database = Load_previous_results()
```

### Get predictions for a given taxon
```R
# Get basic plot for Lutjanus campechanus (in database, so prediction is informed by species-specific data)
Plot_taxa( Search_species(Genus="Lutjanus",Species="campechanus")$match_taxonomy )

# Get basic plot for Sebastes cortezi (not in database, so uses predictive distribution for genus Sebastes)
Plot_taxa( Search_species(Genus="Sebastes",Species="cortezi")$match_taxonomy )

# Get basic plot for Family Scombridae 
Plot_taxa( Search_species(Family="Scombridae")$match_taxonomy )

# Compare two species
Taxa = c( Search_species(Genus="Oncorhynchus",Species="mykiss",add_ancestors=FALSE)$match_taxonomy,
  Search_species(Genus="Salmo",Species="Trutta",add_ancestors=FALSE)$match_taxonomy )
Plot_taxa( Taxa )

```

### Re-run the model
```R
# Load TMB
library( TMB )

# Re-run results with a different model configuration
Fit = Fit_model( N_factors=-3, N_obsfactors=-3, Use_REML=TRUE, RunDir=temp )

# Extract estimates and load into local memory in `Estimate_database` as default for plotting
Estimate_database = Fit[c("N_factors","N_obsfactors","Use_REML","Cov_gjj","ParentChild_gz","ParHat","g_i","Y_ij","Z_ik")]
```
