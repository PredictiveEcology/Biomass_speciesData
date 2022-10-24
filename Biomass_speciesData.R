# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects, use sim$xxx, and are thus globally available
# to all modules. Functions can be used without sim$ as they are namespaced, like functions
# in R packages. If exact location is required, functions will be: sim$<moduleName>$FunctionName
defineModule(sim, list(
  name = "Biomass_speciesData",
  description = paste("Download and pre-process species % cover raster data, overlaying",
                      "lower quality data with higher quality data."),
  keywords = c("LandWeb", "LandR", "LandR Biomass", "species percent cover"),
  authors = c(
    person("Ceres", "Barros", email = "ceres.barros@ubc.ca", role = c("aut", "cre")),
    person(c("Eliot", "J", "B"), "McIntire", email = "eliot.mcintire@nrcan-rncan.gc.ca", role = c("aut")),
    person(c("Alex", "M."), "Chubaty", email = "achubaty@for-cast.ca", role = c("aut"))
  ),
  childModules = character(0),
  version = list(Biomass_speciesData = "1.0.0.9000"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "Biomass_speciesData.Rmd"),
  reqdPkgs = list("data.table", "magrittr", "pryr",
                  "raster", "reproducible (>= 1.2.6.9005)", "SpaDES.core", "SpaDES.tools",
                  # "curl", "httr", ## called directly by this module, but pulled in by LandR (Sep 6th 2022).
                  ## Excluded because loading is not necessary (just installation)
                  "PredictiveEcology/LandR@development (>= 1.0.9.9000)",
                  "PredictiveEcology/pemisc@development"),
  parameters = bindrows(
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description"),
    defineParameter("coverThresh", "integer", 10L, NA, NA,
                    paste("The minimum % cover a species needs to have (per pixel) in the study",
                          "area to be considered present")),
    defineParameter("dataYear", "numeric", 2001, NA, NA,
                    paste("Passed to `paste0('prepSpeciesLayers_', types)` function to fetch data",
                          "from that year (if applicable). Defaults to 2001 as the default kNN year.")),
    defineParameter("sppEquivCol", "character", "Boreal", NA, NA,
                    paste("The column in `sim$sppEquiv` data.table to group species by and use as a",
                          "naming convention. If different species in, e.g., the kNN data have the same",
                          "name in the chosen column, their data are merged into one species by summing",
                          "their % cover in each raster cell.")),
    defineParameter("types", "character", "KNN", NA, NA,
                    paste("The possible data sources. These must correspond to a function named",
                          "`paste0('prepSpeciesLayers_', types)`. Defaults to 'KNN'",
                          "to get the Canadian Forestry Service, National Forest Inventory,",
                          "kNN-derived species cover maps from year 'dataYear', using the",
                          "`LandR::prepSpeciesLayers_KNN` function (see https://open.canada.ca/",
                          "data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990 for details on these data).",
                          "Other currently available options are 'ONFRI', 'CASFRI', 'Pickell' and",
                          "'ForestInventory', which attempt to get proprietary data - the user must be granted",
                          "access first. A custom function can be used to retrieve any data, just as long as",
                          "it is accessible by the module (e.g., in the global environment) and is named as",
                          "`paste0('prepSpeciesLayers_', types)`.")),
    defineParameter("vegLeadingProportion", "numeric", 0.8, 0, 1,
                    "a number that defines whether a species is leading for a given pixel. Only used for plotting."),
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA,
                    "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", NA, NA, NA,
                    "This describes the simulation time interval between plot events"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA,
                    "This describes the simulation time at which the first save event should occur"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA,
                    "This describes the simulation time interval between save events"),
    defineParameter(".sslVerify", "numeric", curl::curl_options("^ssl_verifypeer$"), NA , NA,
                    paste("Passed to `httr::config(ssl_verifypeer = P(sim)$.sslVerify)` when downloading KNN",
                          "(NFI) datasets. Set to 0L if necessary to bypass checking the SSL certificate (this",
                          "may be necessary when NFI's FTP website SSL certificate is down/out-of-date).")),
    defineParameter(".studyAreaName", "character", NA, NA, NA,
                    "Human-readable name for the study area used. If NA, a hash of `studyAreaLarge` will be used."),
    defineParameter(".useCache", "character", "init", NA, NA,
                    desc = "Controls cache; caches the init event by default"),
    defineParameter(".useParallel", "numeric", parallel::detectCores(), NA, NA,
                    "Used in reading csv file with fread. Will be passed to `data.table::setDTthreads`.")
  ),
  inputObjects = bindrows(
    expectsInput("rasterToMatchLarge", "RasterLayer",
                 desc = paste("a raster of `studyAreaLarge` in the same resolution and projection the simulation's.",
                              "Defaults to the using the Canadian Forestry Service, National Forest Inventory,",
                              "kNN-derived stand biomass map."),
                 sourceURL = ""),
    expectsInput("sppColorVect", "character",
                 desc = paste("A named vector of colors to use for plotting.",
                              "The names must be in `sim$sppEquiv[[sim$sppEquivCol]]`,",
                              "and should also contain a color for 'Mixed'"),
                 sourceURL = NA),
    expectsInput("sppEquiv", "data.table",
                 desc = "table of species equivalencies. See `LandR::sppEquivalencies_CA`.",
                 sourceURL = ""),
    expectsInput("studyAreaLarge", "SpatialPolygonsDataFrame",
                 desc =  paste("Polygon to use as the parametrisation study area. Must be provided by the user.",
                               "Note that `studyAreaLarge` is only used for parameter estimation, and",
                               "can be larger than the actual study area used for LandR simulations (e.g,",
                               "larger than `studyArea` in LandR Biomass_core)."),
                 sourceURL = NA),
    expectsInput("studyAreaReporting", "SpatialPolygonsDataFrame",
                 desc = paste("multipolygon (typically smaller/unbuffered than `studyAreaLarge` and `studyArea`",
                              "in LandR Biomass_core) to use for plotting/reporting.",
                              "If not provided, will default to `studyAreaLarge`."),
                 sourceURL = NA)
  ),
  outputObjects = bindrows(
    createsOutput("speciesLayers", "RasterStack",
                  desc = "biomass percentage raster layers by species in Canada species map"),
    createsOutput("treed", "data.table",
                  desc = paste("Table with one logical column for each species, indicating whether",
                               "there were non-zero cover values in each pixel.")),
    createsOutput("numTreed", "numeric",
                  desc = paste("a named vector with number of pixels with non-zero cover values for",
                               "each species")),
    createsOutput("nonZeroCover", "numeric",
                  desc = "A single value indicating how many pixels have non-zero cover")
  )
))

## event types
#   - type `init` is required for initialiazation

doEvent.Biomass_speciesData <- function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "Biomass_speciesData", "initPlot",
                           eventPriority = 1)

      sim <- biomassDataInit(sim)
    },
    initPlot = {
      newDev <- if (!is.null(dev.list())) {
        devCur <- dev.cur()
        max(dev.list()) + 1
      } else {
        1
      }
      dev.set(newDev)
      plotVTM(speciesStack = raster::mask(sim$speciesLayers, sim$studyAreaReporting) %>%
                raster::stack(),
              vegLeadingProportion = P(sim)$vegLeadingProportion,
              sppEquiv = sim$sppEquiv,
              sppEquivCol = P(sim)$sppEquivCol,
              colors = sim$sppColorVect,
              title = "Initial Types")
      if (exists("devCur")) dev.set(devCur)
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

### template initialization
biomassDataInit <- function(sim) {
  cacheTags <- c(currentModule(sim), "function:biomassDataInit")
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": biomassInit() using dataPath '", dPath, "'.")

  if (!exists("speciesLayers", envir = envir(sim), inherits = FALSE))
    sim$speciesLayers <- list()

  for (type in P(sim)$types) {
    fnName <- paste0("prepSpeciesLayers_", type)
    whereIsFnName <- pryr::where(fnName)

    envirName <- attr(whereIsFnName, "name")
    if (is.null(envirName))
      envirName <- environmentName(whereIsFnName)

    message("#############################################")
    message(type, " -- Loading using ", fnName, " located in ", envirName)
    message("#############################################")
    if (!exists(fnName)) {
      stop(fnName, " does not exist. Please make it accessible in a package, as an object, ",
           " or in the .GlobalEnv")
    }

    fn <- get(fnName)
    httr::with_config(config = httr::config(ssl_verifypeer = P(sim)$.sslVerify), {
    speciesLayersNew <- Cache(fn,
                              destinationPath = dPath, # this is generic files (preProcess)
                              outputPath = outputPath(sim), # this will be the studyArea-specific files (postProcess)
                              studyArea = sim$studyAreaLarge,
                              studyAreaName = P(sim)$.studyAreaName,
                              rasterToMatch = sim$rasterToMatchLarge,
                              sppEquiv = sim$sppEquiv,
                              sppEquivCol = P(sim)$sppEquivCol,
                              thresh = P(sim)$coverThresh,
                              year = P(sim)$dataYear,
                              userTags = c(cacheTags, fnName, "prepSpeciesLayers"),
                              omitArgs = c("userTags"))
    })

    sim$speciesLayers <- if (length(sim$speciesLayers) > 0) {
      Cache(overlayStacks,
            highQualityStack = speciesLayersNew,
            lowQualityStack = sim$speciesLayers,
            destinationPath = outputPath(sim),
            userTags = c(cacheTags, "overlayStacks"),
            omitArgs = c("userTags"))
    } else {
      speciesLayersNew
    }
    rm(speciesLayersNew)
  }

  assertSpeciesLayers(sim$speciesLayers, P(sim)$coverThresh)

  species <- names(sim$speciesLayers)

  origFilenames <- vapply(layerNames(sim$speciesLayers),
                          function(r) filename(sim$speciesLayers[[r]]),
                          character(1))

  ## re-enforce study area mask (merged/summed layers are losing the mask)
  sim$speciesLayers <- raster::mask(sim$speciesLayers, sim$rasterToMatchLarge)

  ## make sure empty pixels inside study area have 0 cover, instead of NAs.
  ## this can happen when data has NAs instead of 0s and is not merged/overlayed (e.g. CASFRI)
  tempRas <- sim$rasterToMatchLarge
  tempRas[!is.na(tempRas[])] <- 0
  sim$speciesLayers <- raster::cover(sim$speciesLayers, tempRas)
  rm(tempRas)

  sim$speciesLayers <- if (inMemory(sim$speciesLayers)) {
    sim$speciesLayers
  } else {
    lapply(seq_along(layerNames(sim$speciesLayers)), function(r) {
      writeRaster(sim$speciesLayers[[r]], filename = origFilenames[r], overwrite = TRUE)
    })
  }
  sim$speciesLayers <- raster::stack(sim$speciesLayers) %>% setNames(species)

  singular <- length(P(sim)$types) == 1
  message("sim$speciesLayers is from ", paste(P(sim)$types, collapse = ", "),
          " overlaid in that sequence, higher quality last"[!singular])

  message("------------------")
  message("There are ", sum(!is.na(sim$speciesLayers[[1]][])), " pixels with trees in them")

  # Calculate number of pixels with species cover
  speciesLayersDT <- as.data.table(sim$speciesLayers[] > 0)
  speciesLayersDT[, pixelId := seq(NROW(speciesLayersDT))]
  sim$treed <- na.omit(speciesLayersDT)
  colNames <- names(sim$treed)[!names(sim$treed) %in% "pixelId"]
  sim$numTreed <- sim$treed[, append(
    lapply(.SD, sum),
    list(total = NROW(sim$treed))), .SDcols = colNames]

  # How many have zero cover
  bb <- speciesLayersDT[, apply(.SD, 1, any), .SDcols = 1:nlayers(sim$speciesLayers)]
  sim$nonZeroCover <- sum(na.omit(bb))
  message("There are ", sim$nonZeroCover, " pixels with non-zero tree cover in them.")

  return(invisible(sim))
}

.inputObjects <- function(sim) {
  cacheTags <- c(currentModule(sim), "function:.inputObjects")
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")

  if (!suppliedElsewhere("studyAreaLarge", sim)) {
    stop("Please provide a 'studyAreaLarge' polygon.
         If parameterisation is to be done on the same area as 'studyArea'
         provide the same polygon to 'studyAreaLarge'")
    # message("'studyAreaLarge' was not provided by user. Using the same as 'studyArea'")
    # sim <- objectSynonyms(sim, list(c("studyAreaLarge", "studyArea"))) # Jan 2021 we agreed to force user to provide a SA/SAL
  }

  if (is.na(P(sim)$.studyAreaName)) {
    params(sim)[[currentModule(sim)]][[".studyAreaName"]] <- reproducible::studyAreaName(sim$studyAreaLarge)
    message("The .studyAreaName is not supplied; derived name from sim$studyAreaLarge: ",
            params(sim)[[currentModule(sim)]][[".studyAreaName"]])
  }

  if (!suppliedElsewhere("studyAreaReporting", sim)) {
    message("'studyAreaReporting' was not provided by user. Using the same as 'studyAreaLarge'.")
    sim$studyAreaReporting <- sim$studyAreaLarge
  }

  needRTM <- FALSE
  if (is.null(sim$rasterToMatchLarge)) {
    if (!suppliedElsewhere("rasterToMatchLarge", sim)) {      ## if one is not provided, re do both (safer?)
      needRTM <- TRUE
      message("There is no rasterToMatchLarge supplied; will attempt to use rawBiomassMap")
    } else {
      stop("rasterToMatchLarge is going to be supplied, but ", currentModule(sim), " requires it ",
           "as part of its .inputObjects. Please make it accessible to ", currentModule(sim),
           " in the .inputObjects by passing it in as an object in simInit(objects = list(rasterToMatchLarge = aRaster)",
           " or in a module that gets loaded prior to ", currentModule(sim))
    }
  }

  if (needRTM) {
    ## if rawBiomassMap exists, it needs to match SALarge, if it doesn't make it
    if (!suppliedElsewhere("rawBiomassMap", sim) ||
        !compareRaster(sim$rawBiomassMap, sim$studyAreaLarge, stopiffalse = FALSE)) {
      rawBiomassMapURL <- paste0("http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                                 "canada-forests-attributes_attributs-forests-canada/",
                                 "2001-attributes_attributs-2001/",
                                 "NFI_MODIS250m_2001_kNN_Structure_Biomass_TotalLiveAboveGround_v1.tif")

      httr::with_config(config = httr::config(ssl_verifypeer = P(sim)$.sslVerify), {
      #necessary for KNN
      rawBiomassMapFilename <- basename(rawBiomassMapURL)
      rawBiomassMap <- Cache(prepInputs,
                             targetFile = rawBiomassMapFilename,
                             url = rawBiomassMapURL,
                             destinationPath = dPath,
                             studyArea = sim$studyAreaLarge,
                             rasterToMatch = NULL,
                             maskWithRTM = FALSE,
                             useSAcrs = FALSE,     ## never use SA CRS
                             method = "bilinear",
                             datatype = "INT2U",
                             filename2 = NULL,
                             userTags = c(cacheTags, "rawBiomassMap"),
                             omitArgs = c("destinationPath", "targetFile", "userTags", "stable"))
      })
    } else {
      rawBiomassMap <- Cache(postProcess,
                             x = sim$rawBiomassMap,
                             studyArea = sim$studyAreaLarge,
                             useSAcrs = FALSE,
                             maskWithRTM = FALSE,   ## mask with SA
                             method = "bilinear",
                             datatype = "INT2U",
                             filename2 = NULL,
                             overwrite = TRUE,
                             userTags = cacheTags,
                             omitArgs = c("destinationPath", "targetFile", "userTags", "stable"))
    }

    ## if we need rasterToMatchLarge, that means a) we don't have it, but b) we will have rawBiomassMap
    if (is.null(sim$rasterToMatchLarge))
      warning(paste0("rasterToMatchLarge is missing and will be created \n",
                     "from rawBiomassMap and studyAreaLarge.\n",
                     "If this is wrong, provide raster"))

    sim$rasterToMatchLarge <- rawBiomassMap
    RTMvals <- getValues(sim$rasterToMatchLarge)
    sim$rasterToMatchLarge[!is.na(RTMvals)] <- 1
    sim$rasterToMatchLarge <- Cache(writeOutputs, sim$rasterToMatchLarge,
                                    filename2 = .suffix(file.path(dPath, "rasterToMatchLarge.tif"),
                                                        paste0("_", P(sim)$.studyAreaName)),
                                    datatype = "INT2U", overwrite = TRUE,
                                    userTags = c(cacheTags, "rasterToMatchLarge"),
                                    omitArgs = c("userTags"))
  }

  if (!compareCRS(sim$studyAreaLarge, sim$rasterToMatchLarge)) {
    warning(paste0("studyAreaLarge and rasterToMatchLarge projections differ.\n",
                   "studyAreaLarge will be projected to match rasterToMatchLarge"))
    sim$studyAreaLarge <- spTransform(sim$studyAreaLarge, raster::crs(sim$rasterToMatchLarge))
    sim$studyAreaLarge <- fixErrors(sim$studyAreaLarge)
  }

  if (!suppliedElsewhere("sppEquiv", sim)) {
    if (!is.null(sim$sppColorVect))
      stop("If you provide sppColorVect, you MUST also provide sppEquiv")

    data("sppEquivalencies_CA", package = "LandR", envir = environment())
    sim$sppEquiv <- as.data.table(sppEquivalencies_CA)

    ## check spp column to use
    if (P(sim)$sppEquivCol == "Boreal") {
      message(paste("There is no 'sppEquiv' table supplied;",
                    "will attempt to use species listed under 'Boreal'",
                    "in the 'LandR::sppEquivalencies_CA' table"))
    } else {
      if (any(grepl(P(sim)$sppEquivCol, names(sim$sppEquiv)))) {
        message(paste("There is no 'sppEquiv' table supplied,",
                      "will attempt to use species listed under", P(sim)$sppEquivCol,
                      "in the 'LandR::sppEquivalencies_CA' table"))
      } else {
        stop("You changed 'sppEquivCol' without providing 'sppEquiv',",
             "and the column name can't be found in the default table ('LandR::sppEquivalencies_CA').",
             "Please provide conforming 'sppEquivCol', 'sppEquiv' and 'sppColorVect'")
      }
    }

    ## remove empty lines/NAs
    sim$sppEquiv <- sim$sppEquiv[!"", on = P(sim)$sppEquivCol]
    sim$sppEquiv <- na.omit(sim$sppEquiv, P(sim)$sppEquivCol)

    ## add default colors for species used in model
    sim$sppColorVect <- sppColors(sim$sppEquiv, P(sim)$sppEquivCol,
                                  newVals = "Mixed", palette = "Accent")
  } else {
    if (is.null(sim$sppColorVect))
      stop("If you provide 'sppEquiv' you MUST also provide 'sppColorVect'")
  }

  return(invisible(sim))
}
