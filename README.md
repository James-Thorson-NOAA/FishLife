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

# Get basic plot and extract values for Family Scombridae 
( Predictions = Plot_taxa(Search_species(Family="Scombridae")$match_taxonomy) )
```

### Extract other values
You can also see the full set of parameters calculated for each taxon, either for internal use of anticipated to be useful for users:
```R
head(FishLife::FishBase_and_RAM$beta_gv)
```
These can similarly be extracted and plotted:
```
params = matrix( c("K","M", "G","ln_MASPS"), ncol=2, byrow=TRUE)
Plot_taxa( Search_species(Genus="Lutjanus",Species="campechanus")$match_taxonomy, params=params )
```
while other values (e.g., slope at the origin for the Beverton-Holt stock recruit curve) can then be calculated from the set of available parameters.  

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

Description of research
=============
Presentation of research program available [online](https://www.youtube.com/watch?v=efVXe0J80oU&feature=youtu.be)

Applications for stock assessment
=============

* Sea mullet, east coast Australia, Fisheries Queensland 2022 (link [here](https://era.daf.qld.gov.au/id/eprint/8600/1/sea_mullet_rtex_2022.pdf))
* Spanish mackerel, eastern Australia, Fisheries Queensland 2021 (link [here](http://era.daf.qld.gov.au/id/eprint/8226/25/Spanish%20mackerel%20EC%20stock%20assessment%20report%202021.pdf))
* Uku, Hawaii, PIFSC, 2020 (link [here](https://www.researchgate.net/profile/Marc_Nadon/publication/341385433_Stock_assessment_of_uku_Aprion_virescens_in_Hawaii_2020/links/5ebd99bf92851c11a867bf18/Stock-assessment-of-uku-Aprion-virescens-in-Hawaii-2020.pdf))
* Pollock, Eastern Bering Sea, AFSC, 2018 (link [here](https://archive.fisheries.noaa.gov/afsc/REFM/docs/2018/BSAI/2018EBSpollock.pdf))
* Black Marlin, Indian Ocean, IOTC, 2018 (link [here](https://www.iotc.org/sites/default/files/documents/2018/09/IOTC-2018-WPB16-15_-_BLM_JABBA_Final.pdf))
* Striped Marlin, Indian Ocean, IOTC, 2018 (link [here](https://www.iotc.org/sites/default/files/documents/2018/09/IOTC-2018-WPB16-16_-_MLS_JABBA_Final.pdf))
* Smoothound shark, South Africa, DAFF, 2018 (link [here](https://www.researchgate.net/publication/338491221_Assessment_of_smoothhound_shark_Mustelus_mustelus_in_South_Africa))
* Soupfin shark, South Africa, DAFF, 2018 (link [here](https://www.researchgate.net/publication/338491033_First_comprehensive_assessment_of_soupfin_shark_Galeorhinus_galeus_in_South_Africa))

Journal Arcticles using FishLife
=============

1.  Mora, P., Figueroa-Muñoz, G., Cubillos, L., Strange-Olate, P., 2022. A data-limited approach to determine the status of the artisanal fishery of sea silverside in southern Chile. Mar. Fish. Sci. MAFIS 35, 275–298.

3.	Omori, K.L., Tribuzio, C.G., Babcock, E.A., Hoenig, J.M., 2021. Methods for Identifying Species Complexes Using a Novel Suite of Multivariate Approaches and Multiple Data Sources: A Case Study With Gulf of Alaska Rockfish. Front. Mar. Sci. 1084.

2.	Pawluk, M., Fujiwara, M., Martinez-Andrade, F., 2022. Climate change linked to functional homogenization of a subtropical estuarine system. Ecol. Evol. 12, e8783. https://doi.org/10.1002/ece3.8783

3.	Pons, M., Cope, J.M., Kell, L.T., 2020. Comparing performance of catch-based and length-based stock assessment methods in data-limited fisheries. Can. J. Fish. Aquat. Sci. 77, 1026–1037. https://doi.org/10.1139/cjfas-2019-0276

4.	Rudd, M.B., Thorson, J.T., Sagarese, S.R., 2019. Ensemble models for data-poor assessment: accounting for uncertainty in life-history information. ICES J. Mar. Sci. 76, 870–883. https://doi.org/10.1093/icesjms/fsz012

5.	Thorson, J.T., 2020. Predicting recruitment density dependence and intrinsic growth rate for all fishes worldwide using a data-integrated life-history model. Fish Fish. 21, 237–251. https://doi.org/10.1111/faf.12427

6.	Thorson, J.T., Munch, S.B., Cope, J.M., Gao, J., 2017. Predicting life history parameters for all fishes worldwide. Ecol. Appl. 27, 2262–2276. https://doi.org/10.1002/eap.1606




