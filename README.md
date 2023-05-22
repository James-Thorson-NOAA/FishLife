# FishLife
A package for phylogenetic comparative methods (PCM) and phylogenetic trait imputation (PTI) using phylogenetic factor analysis and/or phylogenetic structural equation models.

The package also includes results for these analyses applied to all fishes fishes globally, estimating:
* life history parameters (growth, maturity, mortality);
* juvenile productivity (stock-recruit parameters);
* life-cycle characteristics (generation time and intrinsic growth rate);
* spawning, behavioral, reproductive, and foraging traits;
* morphometric characteristics;
as described below

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7590994.svg)](https://doi.org/10.5281/zenodo.7590994)

# Database of results
The FishLife package also includes results for three prior analyses, which can be used to extract life-history predictions for all fishes:
* [Thorson et al. 2023](https://doi.org/10.1111/2041-210X.14076):  This includes life-history parameters based on data from FishBase, morphometric information from [FishShapes]([url](http://onlinelibrary.wiley.com/doi/abs/10.1002/ecy.3829)), and spawning, behavioral, reproductive, and trophic traits.  Results are accessed using `data(FishBase, package="FishLife_and_Mophometrics")`, or `Plot_taxa(..., Database = FishLife::FishBase_and_Mophometrics)`
* [Thorson 2020]([url](https://doi.org/10.1111/faf.12427)):  This includes life-history parameters based on data from FishBase as well as stock-recruit information from the RAM Legacy database, and combines these to get life-cycle predictions.  Results are accessed using `data(FishBase, package="FishLife_and_RAM")`, or `Plot_taxa(..., Database = FishLife::FishBase_and_RAM)`
* [Thorson et al. 2017]([url](http://onlinelibrary.wiley.com/doi/10.1002/eap.1606/full)):  This includes life-history parameters based on data from FishBase, and is accessed using `data(FishBase, package="FishLife")`, or `Plot_taxa(..., Database = FishLife::FishBase)`

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
* Thorson, J.T., Maureaud, A.A., Frelat, R., Mérigot, B., Bigman, J.S., Friedman, S.T., Palomares, M.L.D., Pinsky, M.L., Price, S.A., Wainwright, P., 2023. Identifying direct and indirect associations among traits by merging phylogenetic comparative methods and structural equation models. Methods Ecol. Evol. n/a. https://doi.org/10.1111/2041-210X.14076

### Previous software versions and analytical descriptions
* Thorson, J.T., 2020. Predicting recruitment density dependence and intrinsic growth rate for all fishes worldwide using a data-integrated life-history model. Fish Fish. 21, 237–251. https://doi.org/10.1111/faf.12427 
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
 * Yang, W.-H., Martin, T.S., Moffitt, D., 2022. Stock assessment of Queensland east coast dusky flathead (Platycephalus fuscus), Australia, with data to December 2020.
 * Parker, D., Kikuchi, E., Mourato, B.L., 2022. Assessment of the South Atlantic swordfish (Xiphias gladius) stock using JABBA. Collect Vol Sci Pap ICCAT 79, 608–639. (link [here](https://www.researchgate.net/profile/Bruno-Mourato/publication/365748009_ASSESSMENT_OF_THE_SOUTH_ATLANTIC_SWORDFISH_XIPHIAS_GLADIUS_STOCK_USING_JABBA/links/638142f1c2cb154d29293a78/ASSESSMENT-OF-THE-SOUTH-ATLANTIC-SWORDFISH-XIPHIAS-GLADIUS-STOCK-USING-JABBA.pdf))

Journal Arcticles using FishLife
=============
1.	Auber, A., Waldock, C., Maire, A., Goberville, E., Albouy, C., Algar, A.C., McLean, M., Brind’Amour, A., Green, A.L., Tupper, M., Vigliola, L., Kaschner, K., Kesner-Reyes, K., Beger, M., Tjiputra, J., Toussaint, A., Violle, C., Mouquet, N., Thuiller, W., Mouillot, D., 2022. A functional vulnerability framework for biodiversity conservation. Nat. Commun. 13, 4774. https://doi.org/10.1038/s41467-022-32331-y

2.	Fitz, K.S., Montes Jr., H.R., Thompson, D.M., Pinsky, M.L., n.d. Isolation-by-distance and isolation-by-oceanography in Maroon Anemonefish (Amphiprion biaculeatus). Evol. Appl. n/a. https://doi.org/10.1111/eva.13448

3.	Fujiwara, M., Simpson, A., Torres-Ceron, M., Martinez-Andrade, F., 2022. Life-history traits and temporal patterns in the incidence of coastal fishes experiencing tropicalization. Ecosphere 13, e4188. https://doi.org/10.1002/ecs2.4188

4.	Hay, A., Riggins, C.L., Heard, T., Garoutte, C., Rodriguez, Y., Fillipone, F., Smith, K.K., Menchaca, N., Williamson, J., Perkin, J.S., 2022. Movement and mortality of invasive suckermouth armored catfish during a spearfishing control experiment. Biol. Invasions. https://doi.org/10.1007/s10530-022-02834-2

5.	Hirota, D.S., Haimovici, M., Sant’Ana, R., Mourato, B.L., Santos, E.K., Cardoso, L.G., 2022. Life history, population dynamics and stock assessment of the bycatch species Brazilian flathead (Percophis brasiliensis) in southern Brazil. Reg. Stud. Mar. Sci. 102597. https://doi.org/10.1016/j.rsma.2022.102597

6.	Mora, P., Figueroa-Muñoz, G., Cubillos, L., Strange-Olate, P., 2022. A data-limited approach to determine the status of the artisanal fishery of sea silverside in southern Chile. Mar. Fish. Sci. MAFIS 35, 275–298.

7.	Omori, K.L., Tribuzio, C.G., Babcock, E.A., Hoenig, J.M., 2021. Methods for Identifying Species Complexes Using a Novel Suite of Multivariate Approaches and Multiple Data Sources: A Case Study With Gulf of Alaska Rockfish. Front. Mar. Sci. 1084.

8.	Pawluk, M., Fujiwara, M., Martinez-Andrade, F., 2022. Climate change linked to functional homogenization of a subtropical estuarine system. Ecol. Evol. 12, e8783. https://doi.org/10.1002/ece3.8783

9.	Pons, M., Cope, J.M., Kell, L.T., 2020. Comparing performance of catch-based and length-based stock assessment methods in data-limited fisheries. Can. J. Fish. Aquat. Sci. 77, 1026–1037. https://doi.org/10.1139/cjfas-2019-0276

10.	Rudd, M.B., Thorson, J.T., Sagarese, S.R., 2019. Ensemble models for data-poor assessment: accounting for uncertainty in life-history information. ICES J. Mar. Sci. 76, 870–883. https://doi.org/10.1093/icesjms/fsz012

11.	Safaraliev, I.A., Popov, N.N., 2022. Qualitative Assessment of the Stock Status of Freshwater Bream Abramis brama (Cyprinidae) from the Ural Stock Based on the LB-SPR Method. J. Ichthyol. 62, 476–486. https://doi.org/10.1134/S0032945222030134

12.	Thorson, J.T., 2020. Predicting recruitment density dependence and intrinsic growth rate for all fishes worldwide using a data-integrated life-history model. Fish Fish. 21, 237–251. https://doi.org/10.1111/faf.12427

13.	Thorson, J.T., Munch, S.B., Cope, J.M., Gao, J., 2017. Predicting life history parameters for all fishes worldwide. Ecol. Appl. 27, 2262–2276. https://doi.org/10.1002/eap.1606



