# Sarcoscypha in Ukraine
This repository contains the code, data, and extended methodology for the article "Genus *Sarcoscypha* (Pezizomycetes, Ascomycota) in Ukraine: species composition, distribution, threats" by ..., published in *journalname* (year). The article is available at *DOI*.

The repository includes the following:
- R scripts for data analysis and visualisation
- Data files used in the analysis
- Metadata describing the data sources and processing steps
- Extended methodology for species distribution modelling

# Species Disdistribution Modelling
## General information
To explore environmental drivers for distribution of *Sarcoscypha* in Ukraine, we used species distribution modelling. Since morphological identification to a species level in the field is dubious, and both species are known to follow similar distributional pattern, we built a model at the genus level.

We used the MaxEnt method (Phillips et al., 2006) to model the distribution of *Sarcoscypha* in Ukraine, based on presence-background data. The model was built using 13 environmental variables, including bioclimatic variables, vegetation indices, and tree cover data. The model was evaluated using the Area Under the Receiver Operating Characteristic Curve (AUC) metric. We derived prediction uncertainty as the standard deviation of predictions among 100 models. 

We used the R programming language (R Core Team, 2024) and the following packages: *dismo* v. 1.3 (Hijmans et al., 2024), *spThin* v. 0.2 (Aiello-Lammens et al., 2015), *usdm* v. 2.1 (Naimi et al., 2014), *terra* v. 1.8 (Hijmans et al., 2025), *raster* v. 3.6 (Hijmans, 2025), *ggplot2* v. 3.5 (Wickham, 2016).

## Modelling procedure

We obtained *Sarcoscypha* occurrence data (regardless of the species) from the “Fungi of Ukraine” Facebook group (Rudenko, 2025), 29 obs.), iNaturalist (https://www.inaturalist.org/, 618 obs.), GBIF (GBIF.org, 2025), 738 obs.), UkrBIN (UkrBIN, 2017), 27 obs.), PlutoF (Nilsson et al., 2010), 64 obs.), as well as personal communications (24 obs.). 

iNaturalist search query was:

quality_grade=any&identifications=any&place_id=8860&taxon_id=49136&verifiable=true&spam=false

GBIF search query was:

{
  "and" : [
    "Country is Ukraine",
    "OccurrenceStatus is Present",
    "TaxonKey is Sarcoscypha (Fr.) Boud."
  ]
}

All raw occurrence data are stored in the [*data*](https://github.com/olehprylutskyi/Sarcoscypha/tree/main/SDM/data/raw_occ_data) folder of this repository. The data were filtered to remove duplicates, records with unknown coordinates, and records with known coordinate uncertainty greater than 200 m. The final dataset included 298 unique occurrence records.

To reduce both spatial bias and clustering in occurrence data, we performed a two-step spatial thinning procedure by (i) selecting only one occurrence point within each pixel of the covariate rasters; (ii) leaving only points separated from each other by no less than 5 km, using *thin()* function with 100 replications from *spThin* R package v. 0.2 (Aiello-Lammens et al., 2015). After that procedure, 216 presence points were left for further analysis.

Since we do not have true absence points, we employed the presence-background modelling method. For this purpose, we randomly generated 100 times more background points than we had presence points within Ukraine administrative boundary spatial polygon, covering the variation of environmental conditions in the area of modelling, using the *spsample()* function from *sp* R package v. 1.6 (Pebesma and Bivand, 2005). After removing duplicates, 20,916 background points were left.

As covariates for the *Sarcoscypha* distribution model, we used CHELSA climatological variables (Karger et al., 2017) at a spatial resolution of 30 arcsec (~1 km). We selected 19 bioclimatic variables for current conditions (bio1–bio19, means of 1981–2010), as well as additional climatological covariates: annual range (cmi_range) and mean monthly (cmi_mean) climate moisture index, annual range (hurs_range) and mean monthly (hurs_mean) near-surface relative humidity, frost change frequency (fcf), annual number of days with snow cover (scd), and growing season length (gsl).

We also incorporated Earth observation data as three season-constrained vegetation indices – Normalized Difference Vegetation Index (NDVI), Normalized Difference Water Index (NDWI) and Enhanced Vegetation Index (EVI) for summer (May–September) and winter (November–March) seasons – six indices in total, which were supposed to reflect vegetation dynamics. The data source was two MODIS-derived datasets: MOD13A1.061 Terra Vegetation Indices 16-Day Global 500m	(https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MOD13A1) and MODIS Combined 16-Day NDWI datasets (https://developers.google.com/earth-engine/datasets/catalog/MODIS_MCD43A4_006_NDWI), accessed through the Google Earth Engine (Gorelick et al., 2017) Data Catalogue. All indices were calculated using Google Earth Engine (Gorelick et al., 2017) by applying a “median” function for the temporal subset of all available images between 2010 and 2022, for respective months.

We also incorporated tree cover data from the Global Plant Functional Types (PFT) dataset, v2.0.8, by European Space Agency (Harper et al., 2023). Two datasets describing the percentage cover (0-100%) of a plant functional type at a spatial resolution of 300 m for two categories – broadleaved deciduous trees and ​​needle-leaved evergreen trees – were picked.

Since both MODIS-derived and tree cover data have finer spatial resolution than CHELSA data, we upscaled rasters of vegetation indices to the same spatial extent and resolution as CHELSA, using the *terra* R package v. 1.8 (Hijmans et al., 2025). All covariates were then cropped by the Ukraine administrative boundary and by a custom polygon mask of large inland water bodies using the *raster* R package v. 3.6 (Hijmans, 2025). The total number of environmental variables was 37.

To control collinearity among environmental variables, we first applied Variance Inflation Factor analysis (Akinwande et al., 2015) using the *vifstep()* function from the *usdm* R package v. 2.1 (Naimi et al., 2014), with a threshold of 10. Covariates with Variance Inflation Factor values above the threshold level were excluded. After that, all remaining variables were tested for Pearson correlation using the *cor()* function from the *stats* R package (R Core Team, 2024). For all pairs of variables exhibiting correlation above |0.7|, only one was left based on presumable relevance for spring cup fungi.

As a result, we selected 13 covariates: mean diurnal air temperature range (bio02), precipitation seasonality (bio15), mean monthly precipitation amount of the coldest quarter (bio19), frost change frequency (fcf), growing season length (gsl), mean (hurs_mean) and annual range (hurs_range) of monthly near-surface relative humidity, snow cover days (scd), Enhanced Vegetation Index for both summer and winter seasons (EVIs, EVIw), Normalized Difference Water Index for winter season (NDWIw), broadleaved deciduous tree cover (TREES_DB), and needleleaved evergreen tree cover (TREES_NE).

As a modelling method, we choose MaxEnt (Phillips et al., 2006), implemented in the *dismo* R package v. 1.3 (Hijmans et al., 2024). For the same set of covariates, we fitted 100 models, each with a random seed and independent splitting data into training and testing sets. Each splitting was performed using the kfold() function from the *dismo* R package, 80% randomly assigned for training and 20% for a testing subset. The Area under the Receiver Operating Curve (AUC) was used to assess the predictive performance of the models (Swets, 1988), obtained for each model with the evaluate function of the *dismo* R package.

After that, we calculated model prediction, variables' contribution, and predictive performance as mean values between 100 models' outputs. We derived prediction uncertainty as the standard deviation of predictions among 100 models. All raster calculations were performed using the raster R package.

Final visualisations were made using the *ggplot2* R package v. 3.5 (Wickham, 2016). All calculations were performed using R 4.2.2 (R Core Team, 2024).

![Predicted distribution of Sarcoscypha in Ukraine, prediction uncertainty, and key environmental drivers for its distribution](https://github.com/olehprylutskyi/Sarcoscypha/blob/main/SDM/outputs/publication_ready/Maxent_chelsa_plate.png)

## References
Aiello-Lammens, M. E., Boria, R. A., Radosavljevic, A., Vilela, B., & Anderson, R. P. (2015). spThin: An R package for spatial thinning of species occurrence records for use in ecological niche models. Ecography, 38(5), 541–545. https://doi.org/10.1111/ecog.01132

Akinwande, M. O., Dikko, H. G., & Samson, A. (2015). Variance Inflation Factor: As a Condition for the Inclusion of Suppressor Variable(s) in Regression Analysis. Open Journal of Statistics, 5(7), Article 7. https://doi.org/10.4236/ojs.2015.57075

GBIF.org. (2025). Genus Sarcoscypha. GBIF Occurrence Download [Dataset]. GBIF. https://doi.org/10.15468/dl.j3qx6j

Gorelick, N., Hancher, M., Dixon, M., Ilyushchenko, S., Thau, D., & Moore, R. (2017). Google Earth Engine: Planetary-scale geospatial analysis for everyone. Remote Sensing of Environment. https://doi.org/10.1016/j.rse.2017.06.031

Harper, K. L., Lamarche, C., Hartley, A., Peylin, P., Ottlé, C., Bastrikov, V., San Martín, R., Bohnenstengel, S. I., Kirches, G., Boettcher, M., Shevchuk, R., Brockmann, C., & Defourny, P. (2023). ESA Land Cover Climate Change Initiative (Land_Cover_cci): Global Plant Functional Types (PFT) Dataset, v2.0.8 [Application/xml]. NERC EDS Centre for Environmental Data Analysis. https://doi.org/10.5285/26A0F46C95EE4C29B5C650B129AAB788

Hijmans, R. J. (2025). Raster: Geographic Data Analysis and Modeling. R package (Version 3.6-32) [Computer software]. https://CRAN.R-project.org/package=raster

Hijmans, R. J., Bivand, R., Pebesma, E., & Sumner, M. D. (2025). terra: Spatial Data Analysis (Version 1.8-29) [Computer software]. https://CRAN.R-project.org/package=terra

Hijmans, R. J., Phillips, S., Leathwick, J., & Elith, J. (2024). dismo: Species Distribution Modeling (Version 1.3-16) [Computer software]. https://CRAN.R-project.org/package=dismo

Karger, D. N., Conrad, O., Böhner, J., Kawohl, T., Kreft, H., Soria-Auza, R. W., Zimmermann, N. E., Linder, H. P., & Kessler, M. (2017). Climatologies at high resolution for the earth’s land surface areas. Scientific Data, 4(1), Article 1. https://doi.org/10.1038/sdata.2017.122

Naimi, B., Hamm, N. A. S., Groen, T. A., Skidmore, A. K., & Toxopeus, A. G. (2014). Where is positional uncertainty a problem for species distribution modelling? Ecography, 37(2), 191–203. https://doi.org/10.1111/j.1600-0587.2013.00205.x

Nilsson, H., Kessy Abarenkov, Leho Tedersoo, Kai Vellak, Irja Saar, Vilmar Veldre, Erast Parmasto, Marko Prous, Anne Aan, Margus Ots, Olavi Kurina, Ivika Ostonen, Janno Jõgeva, Siim Halapuu, Kadri Põldmaa, Märt Toots, Jaak Truu, Karl-Henrik Larsson, & Urmas Kõljalg. (2010). PlutoF—a Web Based Workbench for Ecological and Taxonomic Research, with an Online Implementation for Fungal ITS Sequences. Evolutionary Bioinformatics, 189. https://doi.org/10.4137/EBO.S6271

Pebesma, E. J., & Bivand, R. S. (2005). Classes and methods for spatial data in R. R News, 5(2), 9–13.

Phillips, S. J., Anderson, R. P., & Schapire, R. E. (2006). Maximum entropy modeling of species geographic distributions. Ecological Modelling, 190(3–4), 231–259. https://doi.org/10.1016/j.ecolmodel.2005.03.026

R Core Team. (2024). R: A language and environment for statistical computing [Computer software]. Foundation for Statistical Computing. https://www.R-project.org/

Rudenko, Y. (2025). Fungi of Ukraine. https://www.facebook.com/groups/Hryby.Ukrayiny

UkrBIN. (2017). UkrBIN: Ukrainian Biodiversity Information Network [public project & web application]. UkrBIN, Database on Biodiversity Information. http://ukrbin.com/

Wickham, H. (2016). ggplot2: Elegant Graphics for Data Analysis (2nd ed. 2016). Springer International Publishing : Imprint: Springer. https://doi.org/10.1007/978-3-319-24277-4