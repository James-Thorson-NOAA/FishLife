# FishLife
Estimate fish traits for all marine fish species globally

[![DOI](https://zenodo.org/badge/67250650.svg)](https://zenodo.org/badge/latestdoi/67250650)

# Visualize predictions

A graphical user interface (GUI) is available [online](https://james-thorson.shinyapps.io/FishLife/)

# Example usage

### Load the package
```R
# Install and load package
devtools::install_github("james-thorson/FishLife")
library( FishLife )
```

### Vignette available
Please see the `FishLife` vignette for details on how to extract predictions frmo the package, update predictions using new data, or replicate the analysis using a new data set.  
```R
vignette("tutorial","FishLife")
```

### Get predictions for a given taxon
I also show a few simple examples of life-history predictions using `FishLife`, as archived in the package.  
```R
# Get basic plot for Lutjanus campechanus (in database, so prediction is informed by species-specific data)
Plot_taxa( Search_species(Genus="Lutjanus",Species="campechanus")$match_taxonomy )

# Get basic plot for Sebastes cortezi (not in database, so uses predictive distribution for genus Sebastes)
Plot_taxa( Search_species(Genus="Sebastes",Species="cortezi")$match_taxonomy )

# Get basic plot for Family Scombridae 
Plot_taxa( Search_species(Family="Scombridae")$match_taxonomy )
```

Description of package
=============
### Please cite if using the software
* Thorson, J. T., S. B. Munch, J. M. Cope, and J. Gao. In press. Predicting life history parameters for all fishes worldwide. Ecological Applications.

Further reading
=============
### Evaluating accuracy of data and life-history predictions in FishBase
* Thorson, J. T., J. M. Cope, and W. S. Patrick. 2014. Assessing the quality of life history information in publicly available databases. Ecological Applications 24:217–226.

