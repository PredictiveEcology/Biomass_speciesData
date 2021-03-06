# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects, use sim$xxx, and are thus globally available
# to all modules. Functions can be used without sim$ as they are namespaced, like functions
# in R packages. If exact location is required, functions will be: sim$<moduleName>$FunctionName
defineModule(sim, list(
  name = "Biomass_speciesData",
  description = "Download and pre-process proprietary LandWeb data.",
  keywords = c("LandWeb", "LandR"),
  authors = c(
    person(c("Eliot", "J", "B"), "McIntire", email = "eliot.mcintire@canada.ca", role = c("aut", "cre")),
    person(c("Alex", "M."), "Chubaty", email = "achubaty@for-cast.ca", role = c("aut")),
    person("Ceres", "Barros", email = "cbarros@mail.ubc.ca", role = c("aut"))
  ),
  childModules = character(0),
  version = list(Biomass_speciesData = "1.0.0"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "Biomass_speciesData.Rmd"),
  reqdPkgs = list("data.table", "magrittr", "pryr",
                  "raster", "reproducible (>=1.0.0.9011)", "SpaDES.core", "SpaDES.tools",
                  "PredictiveEcology/LandR@development (>=0.0.11.9008)",
                  "PredictiveEcology/pemisc@development"),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description"),
    defineParameter("coverThresh", "integer", 10, NA, NA,
                    paste("The minimum % cover a species needs to have (per pixel) in the study",
                          "area to be considered present")),
    defineParameter("sppEquivCol", "character", "Boreal", NA, NA,
                    "The column in sim$specieEquivalency data.table to use as a naming convention"),
    defineParameter("types", "character", "KNN", NA, NA,
                    paste("The possible data sources. These must correspond to a function named",
                          "paste0('prepSpeciesLayers_', types). Defaults to 'KNN', to get the",
                          "Canadian Forestry Service, National Forest Inventory, kNN-derived species",
                          "cover maps from 2001 - see https://open.canada.ca/data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990",
                          "for metadata")),
    defineParameter("vegLeadingProportion", "numeric", 0.8, 0, 1,
                    "a number that define whether a species is leading for a given pixel"),
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA,
                    "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", NA, NA, NA,
                    "This describes the simulation time interval between plot events"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA,
                    "This describes the simulation time at which the first save event should occur"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA,
                    "This describes the simulation time interval between save events"),
    defineParameter(".studyAreaName", "character", NA, NA, NA,
                    "Human-readable name for the study area used. If NA, a hash of studyArea will be used."),
    defineParameter(".useCache", "logical", "init", NA, NA,
                    desc = "Controls cache; caches the init event by default"),
    defineParameter(".useParallel", "numeric", parallel::detectCores(), NA, NA,
                    "Used in reading csv file with fread. Will be passed to data.table::setDTthreads.")
  ),
  inputObjects = bindrows(
    expectsInput("rasterToMatchLarge", "RasterLayer",
                 desc = paste("a raster of the studyAreaLarge in the same resolution and projection as biomassMap"),
                 sourceURL = ""),
    expectsInput("sppColorVect", "character",
                 desc = paste("A named vector of colors to use for plotting.",
                              "The names must be in sim$speciesEquivalency[[sim$sppEquivCol]],",
                              "and should also contain a color for 'Mixed'"),
                 sourceURL = NA),
    expectsInput("sppEquiv", "data.table",
                 desc = "table of species equivalencies. See LandR::sppEquivalencies_CA.",
                 sourceURL = ""),
    expectsInput("studyAreaLarge", "SpatialPolygonsDataFrame",
                 desc =  paste("Polygon to use as the parametrisation study area.",
                               "(studyAreaLarge is only used for parameter estimation, and",
                               "can be larger than the actual study area used for LandR Biomass simulations).",
                               "If not provided by the user, it will default to an area in Southwestern Alberta,",
                               "Canada (which is the same as the default study area used for LandR Biomass simulations)."),
                 sourceURL = NA),
    expectsInput("studyAreaReporting", "SpatialPolygonsDataFrame",
                 desc = paste("multipolygon (typically smaller/unbuffered than studyArea) to use for plotting/reporting.",
                              "Defaults to an area in Southwestern Alberta, Canada."),
                 sourceURL = NA)
  ),
  outputObjects = bindrows(
    createsOutput("speciesLayers", "RasterStack",
                  desc = "biomass percentage raster layers by species in Canada species map"),
    createsOutput("treed", "data.table",
                  desc = "one logical column for each species, indicating whether there were non-zero values"),
    createsOutput("numTreed", "numeric",
                  desc = "a named vector with number of pixels with non-zero cover values"),
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
      devCur <- dev.cur()
      newDev <- if (!is.null(dev.list())) max(dev.list()) + 1 else 1
      quickPlot::dev(newDev)
      plotVTM(speciesStack = raster::mask(sim$speciesLayers, sim$studyAreaReporting) %>%
                raster::stack(),
              vegLeadingProportion = P(sim)$vegLeadingProportion,
              sppEquiv = sim$sppEquiv,
              sppEquivCol = P(sim)$sppEquivCol,
              colors = sim$sppColorVect,
              title = "Initial Types")
      quickPlot::dev(devCur)
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
      envirName <- whereIsFnName

    message("#############################################")
    message(type, " -- Loading using ", fnName, " located in ", envirName)
    message("#############################################")
    if (!exists(fnName)) {
      stop(fnName, " does not exist. Please make it accessible in a package, as an object, ",
           " or in the .GlobalEnv")
    }

    fn <- get(fnName)
    speciesLayersNew <- Cache(fn,
                              destinationPath = dPath, # this is generic files (preProcess)
                              outputPath = outputPath(sim), # this will be the studyArea-specific files (postProcess)
                              studyArea = sim$studyAreaLarge,
                              studyAreaName = P(sim)$.studyAreaName,
                              rasterToMatch = sim$rasterToMatchLarge,
                              sppEquiv = sim$sppEquiv,
                              sppEquivCol = P(sim)$sppEquivCol,
                              thresh = P(sim)$coverThresh,
                              userTags = c(cacheTags, fnName, "prepSpeciesLayers"),
                              omitArgs = c("userTags"))

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
  sim$speciesLayers <- cover(sim$speciesLayers, tempRas)
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
    message("'studyAreaLarge' was not provided by user. Using a polygon (6250000 m^2) in southwestern Alberta, Canada")
    sim$studyAreaLarge <- randomStudyArea(seed = 1234, size = (250^2)*100)
  }

  if (is.na(P(sim)$.studyAreaName)) {
    params(sim)[[currentModule(sim)]][[".studyAreaName"]] <- studyAreaName(sim$studyAreaLarge)
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
      httr::with_config(config = httr::config(ssl_verifypeer = 0L), { ## TODO: re-enable verify
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
                                    filename2 = file.path(cachePath(sim), "rasters", "rasterToMatchLarge.tif"),
                                    datatype = "INT2U", overwrite = TRUE,
                                    userTags = c(cacheTags, "rasterToMatchLarge"),
                                    omitArgs = c("userTags"))
  }

  if (!compareCRS(sim$studyAreaLarge, sim$rasterToMatchLarge)) {
    warning(paste0("studyAreaLarge and rasterToMatchLarge projections differ.\n",
                   "studyAreaLarge will be projected to match rasterToMatchLarge"))
    sim$studyAreaLarge <- spTransform(sim$studyAreaLarge, crs(sim$rasterToMatchLarge))
    sim$studyAreaLarge <- fixErrors(sim$studyAreaLarge)
  }

  if (!suppliedElsewhere("sppEquiv", sim)) {
    if (!is.null(sim$sppColorVect))
      stop("If you provide sppColorVect, you MUST also provide sppEquiv")

    data("sppEquivalencies_CA", package = "LandR", envir = environment())
    sim$sppEquiv <- as.data.table(sppEquivalencies_CA)
    ## By default, Abies_las is renamed to Abies_sp
    sim$sppEquiv[KNN == "Abie_Las", LandR := "Abie_sp"]

    ## check spp column to use
    if (P(sim)$sppEquivCol == "Boreal") {
      message(paste("There is no 'sppEquiv' table supplied;",
                    "will attempt to use species listed under 'Boreal'",
                    "in the 'LandR::sppEquivalencies_CA' table"))
    } else {
      if (grepl(P(sim)$sppEquivCol, names(sim$sppEquiv))) {
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
