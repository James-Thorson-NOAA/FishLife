# FishLife
Estimate growth, size, maturity, mortality, stock-recruit, and population-dynamics parameters for all fish species globally

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

### Use old database
By default `FishLife` uses the most-recent version published.  This currently includes both growth, size, maturity, and mortality parameters from FishBase, as well as stock-recruit parameters estimated using the RAM Legacy stock-recruit database.  To use earlier versions, use the `Database` argument in each function:

```R
# Get basic plot for Lutjanus campechanus (in database, so prediction is informed by species-specific data)
Plot_taxa( Search_species(Genus="Lutjanus",Species="campechanus")$match_taxonomy, Database="FishBase" )
```

or expliclty use the updated database using:

```R
Plot_taxa( Search_species(Genus="Lutjanus",Species="campechanus")$match_taxonomy, Database="FishBase_and_RAM" )`
```

Description of package
=============
### Please cite if using the software
* Thorson, J. T. In press.  Predicting recruitment density dependence and intrinsic growth rate for all fishes worldwide using a data-integrated life-history model.  Fish and Fisheries. 
* Thorson, J. T., S. B. Munch, J. M. Cope, and J. Gao. 2017. Predicting life history parameters for all fishes worldwide. Ecological Applications. 27(8): 2262–2276. http://onlinelibrary.wiley.com/doi/10.1002/eap.1606/full

Further reading
=============
### Evaluating accuracy of data and life-history predictions in FishBase
* Thorson, J. T., J. M. Cope, and W. S. Patrick. 2014. Assessing the quality of life history information in publicly available databases. Ecological Applications 24:217–226. http://onlinelibrary.wiley.com/doi/10.1890/12-1855.1/abstract

Applications for stock assessment
=============

* Uku, Pacific Islands, 2020 (link [here](https://www.researchgate.net/profile/Marc_Nadon/publication/341385433_Stock_assessment_of_uku_Aprion_virescens_in_Hawaii_2020/links/5ebd99bf92851c11a867bf18/Stock-assessment-of-uku-Aprion-virescens-in-Hawaii-2020.pdf))

