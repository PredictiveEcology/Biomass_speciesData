---
title: "Biomass_speciesData"
author: "Eliot McIntire, Alex M Chubaty, Ceres Barros"
date: "11 May 2021"
output:
  html_document:
    keep_md: yes
editor_options:
  chunk_output_type: console
---



[![Gitter](https://badges.gitter.im/PredictiveEcology/LandR_Biomass.svg)](https://gitter.im/PredictiveEcology/LandR_Biomass?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)

# Overview

Download and pre-process species % cover data to prepare species cover rasters.
This module can currently access different data sources (kNN, CASFRI, Paul Pickell's dataset and data from forest inventories in Alberta), however only kNN datasets are freely available.
Remaining datasets can only be accessed by authorized Google users.

The module defaults to processing cover data for 5 species/genera *Abies sp.*, *Picea glauca*, *Picea mariana*, *Pinus sp.* and  *Populus tremuloides* in a random polygon in Alberta.

# Usage

After downloading the data that can be accessed for the species chosen, the module merges the datasets (if several are used) by filling in the gaps of higher quality data with lower quality data.
The parameter `types` lists the data sources to be used and its order determines the hierarchy of data filling from the lowest quality to highest quality. So, e.g., if using `types = c("KNN", "CASFRI")`, CASFRI is considered the highest quality source, and filled in with information from KNN.
If the user wishes to use another data source, they will need to provide an accompanying function prefixed `prepSpeciesLayers_*`, where * corresponds to the new data source's character string in `types`. For instance, the user would create a `prepSpeciesLayers_MyDataType` function and pass `types = c("MyDataType")` to the `types` parameter.
We suggest following the `LandR::prepSpeciesLayers_KNN` and `LandR::prepSpeciesLayers_CASFRI` functions as templates to build a custom function.



```r
if (!require("SpaDES.install")) {
  devtools::install_github("PredictiveEcology/SpaDES.install")
}
SpaDES.install::makeSureAllPackagesInstalled(modulePath = "../Biomass_speciesData")

# User may want to set some options -- see ?reproducibleOptions 
#    -- often will be set outside of project by user --
# options(reproducible.inputPaths = "E:/Data/LandR_related/") # to re-use datasets across projects
```


```r
library(raster)
library(SpaDES)

setPaths(inputPath = file.path(tempDir, "inputs"), 
         cachePath = file.path(tempDir, "cache"), 
         modulePath = "../", 
         outputPath = file.path(tempDir, "outputs"))

## do you want to hand-draw a map or use defaults?
# - note that large areas will take longer to compute
handDrawMap <- TRUE

if (handDrawMap) {
  dev()
  clearPlot()
  canadaMap <- Cache(getData, 'GADM', country = 'CAN', level = 1, path = "data/",
                     cacheRepo = getPaths()$cachePath, quick = FALSE)
  possibleStudyArea <- "../LandscapesInMotion/data/maps/Foothills_study_area.shp"
  Plot(canadaMap, speedup = 5, visualSqueeze = 0.9) # 5 seemed optimal
  if (file.exists(possibleStudyArea)) {
    LIM_SA <- shapefile(possibleStudyArea)
    Plot(LIM_SA, addTo = "canadaMap", col = "green")
  }
  
  ## hand-drawn study area
  if (!exists("studyAreaLarge")) {
    message("Since there is no object called 'studyAreaLarge', please draw a study area with 10 points")
    
    # use Cache here, so that if you restart, can use same random polygon, allowing Cache to be effective during
    #   the simInit and spades calls below
    severalrandompoints <- Cache(clickCoordinates, 10) 
    if (startsWith(attr(severalrandompoints, "tags"), "cache")) message("Taking studyAreaLarge from Cache")
    studyAreaLarge <- SpatialPolygons(list(
      Polygons(list(Polygon(severalrandompoints$coords)),
               ID = "handDrawnPoly")),
      proj4string = crs(canadaMap))
  }
} else {
  studyAreaLarge <- Cache(randomStudyArea, size = 1e8) # cache this so it creates a random one only once on a machine
}

modules <- list("Biomass_speciesData")
objects <- list("studyAreaLarge" = studyAreaLarge)

opts <- options(reproducible.useCache = TRUE)
mySim <- simInit(modules = modules, objects = objects, paths = getPaths())

mySimOut <- spades(mySim)
options(opts)
```

# Parameters

Summary of user-visible parameters:


```
## defineParameter: 'coverThresh' is not of specified type 'integer'.
```

```
## defineParameter: '.useCache' is not of specified type 'logical'.
```



|paramName            |paramClass |default |min |max |paramDesc                                                                                                                                                                                                                                                                                                                             |
|:--------------------|:----------|:-------|:---|:---|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|coverThresh          |integer    |10      |NA  |NA  |The minimum % cover a species needs to have (per pixel) in the study area to be considered present                                                                                                                                                                                                                                    |
|demoMode             |logical    |FALSE   |NA  |NA  |if TRUE, the module can be run with no input objects. Else, at a minimum, studyAreaLarge must be provided                                                                                                                                                                                                                             |
|sppEquivCol          |character  |Boreal  |NA  |NA  |The column in sim$specieEquivalency data.table to use as a naming convention                                                                                                                                                                                                                                                          |
|types                |character  |KNN     |NA  |NA  |The possible data sources. These must correspond to a function named paste0('prepSpeciesLayers_', types). Defaults to 'KNN', to get the Canadian Forestry Service, National Forest Inventory, kNN-derived species cover maps from 2001 - see https://open.canada.ca/data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990 for metadata |
|vegLeadingProportion |numeric    |0.8     |0   |1   |a number that define whether a species is leading for a given pixel                                                                                                                                                                                                                                                                   |
|.plotInitialTime     |numeric    |NA      |NA  |NA  |This describes the simulation time at which the first plot event should occur                                                                                                                                                                                                                                                         |
|.plotInterval        |numeric    |NA      |NA  |NA  |This describes the simulation time interval between plot events                                                                                                                                                                                                                                                                       |
|.saveInitialTime     |numeric    |NA      |NA  |NA  |This describes the simulation time at which the first save event should occur                                                                                                                                                                                                                                                         |
|.saveInterval        |numeric    |NA      |NA  |NA  |This describes the simulation time interval between save events                                                                                                                                                                                                                                                                       |
|.studyAreaName       |character  |NA      |NA  |NA  |Human-readable name for the study area used. If NA, a hash of studyArea will be used.                                                                                                                                                                                                                                                 |
|.useCache            |logical    |init    |NA  |NA  |Controls cache; caches the init event by default                                                                                                                                                                                                                                                                                      |
|.useParallel         |numeric    |16      |NA  |NA  |Used in reading csv file with fread. Will be passed to data.table::setDTthreads.                                                                                                                                                                                                                                                      |

# Data dependencies

## Input data

Summary of input objects:


```
## defineParameter: 'coverThresh' is not of specified type 'integer'.
```

```
## defineParameter: '.useCache' is not of specified type 'logical'.
```



|objectName         |objectClass              |desc                                                                                                                                                                                                                                                                                                                                                                 |sourceURL |
|:------------------|:------------------------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:---------|
|rasterToMatchLarge |RasterLayer              |a raster of the studyAreaLarge in the same resolution and projection as biomassMap                                                                                                                                                                                                                                                                                   |          |
|sppColorVect       |character                |A named vector of colors to use for plotting. The names must be in sim$speciesEquivalency[[sim$sppEquivCol]], and should also contain a color for 'Mixed'                                                                                                                                                                                                            |NA        |
|sppEquiv           |data.table               |table of species equivalencies. See LandR::sppEquivalencies_CA.                                                                                                                                                                                                                                                                                                      |          |
|studyAreaLarge     |SpatialPolygonsDataFrame |Polygon to use as the parametrisation study area. (studyAreaLarge is only used for parameter estimation, and can be larger than the actual study area used for LandR Biomass simulations). If not provided by the user, it will default to an area in Southwestern Alberta, Canada (which is the same as the default study area used for LandR Biomass simulations). |NA        |
|studyAreaReporting |SpatialPolygonsDataFrame |multipolygon (typically smaller/unbuffered than studyArea) to use for plotting/reporting. Defaults to an area in Southwestern Alberta, Canada.                                                                                                                                                                                                                       |NA        |

## Output data

Summary of the module outputs:


```
## defineParameter: 'coverThresh' is not of specified type 'integer'.
```

```
## defineParameter: '.useCache' is not of specified type 'logical'.
```



|objectName    |objectClass |desc                                                                               |
|:-------------|:-----------|:----------------------------------------------------------------------------------|
|speciesLayers |RasterStack |biomass percentage raster layers by species in Canada species map                  |
|treed         |data.table  |one logical column for each species, indicating whether there were non-zero values |
|numTreed      |numeric     |a named vector with number of pixels with non-zero cover values                    |
|nonZeroCover  |numeric     |A single value indicating how many pixels have non-zero cover                      |

# Links to other modules

Intended for use with the LandR-Biomass suite of models - if LandWeb project permissions are valid, then the full datasets used for the LandWeb can be used here.

## Getting help

- <https://gitter.im/PredictiveEcology/LandR_Biomass>
