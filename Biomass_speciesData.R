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
    person(c("Alex", "M."), "Chubaty", email = "achubaty@friresearch.ca", role = c("aut")),
    person("Ceres", "Barros", email = "cbarros@mail.ubc.ca", role = c("aut"))
  ),
  childModules = character(0),
  version = list(SpaDES.core = "0.2.3.9009", Biomass_speciesData = "1.0.0", LandR = "0.0.3.9000"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "Biomass_speciesData.Rmd"),
  reqdPkgs = list("RCurl", "XML", "data.table", "magrittr", "raster",
                  "reproducible", "SpaDES.core", "SpaDES.tools",
                  "PredictiveEcology/LandR@development",
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
    defineParameter(".useCache", "logical", "init", NA, NA,
                    desc = "Controls cache; caches the init event by default"),
    defineParameter(".useParallel", "numeric", parallel::detectCores(), NA, NA,
                    "Used in reading csv file with fread. Will be passed to data.table::setDTthreads.")
  ),
  inputObjects = bind_rows(
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
                               "can be larger than the actual study area of interest).",
                               "If not provided by the user, it will first default to 'studyArea',",
                               "if this object exists. If not, it will default to an area in",
                               "Southwestern Alberta, Canada (the same as the default used for 'studyArea')."),
                 sourceURL = NA),
    expectsInput("studyAreaReporting", "SpatialPolygonsDataFrame",
                 desc = paste("multipolygon (typically smaller/unbuffered than studyArea) to use for plotting/reporting.",
                              "Defaults to an area in Southwestern Alberta, Canada."),
                 sourceURL = NA)
  ),
  outputObjects = bind_rows(
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
      quickPlot::dev(2)
      plotVTM(speciesStack = raster::mask(sim$speciesLayers, sim$studyAreaReporting) %>% stack(),
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
                              rasterToMatch = sim$rasterToMatchLarge,
                              sppEquiv = sim$sppEquiv,
                              sppEquivCol = P(sim)$sppEquivCol,
                              thresh = P(sim)$coverThresh,
                              userTags = c(cacheTags, "prepSpeciesLayers"),
                              omitArgs = c("userTags"))
    sim$speciesLayers <- if (length(sim$speciesLayers) > 0) {
      Cache(overlayStacks,
            highQualityStack = speciesLayersNew,
            lowQualityStack = sim$speciesLayers,
            destinationPath = dPath,
            userTags = c(cacheTags, "overlayStacks"),
            omitArgs = c("userTags"))
    } else {
      speciesLayersNew
    }
    rm(speciesLayersNew)
  }

  assertSpeciesLayers(sim$speciesLayers, P(sim)$coverThresh)

  species <- names(sim$speciesLayers)

  origFilenames <- vapply(layerNames(sim$speciesLayers), function(r) filename(sim$speciesLayers[[r]]),
                          character(1))

  ## re-enforce study area mask (merged/summed layers are losing the mask)
  sim$speciesLayers <- raster::mask(sim$speciesLayers, sim$studyAreaLarge)

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
  message("There are ", sum(!is.na(sim$speciesLayers[[1]][])),
          " pixels with trees in them")

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
    if (suppliedElsewhere("studyArea", sim) && !is.null(sim$studyArea)) {
      message("'studyAreaLarge' was not provided by user. Using the same as 'studyArea'")
      sim$studyAreaLarge <- sim$studyArea
    } else {
      message("'studyAreaLarge' was not provided by user. Using a polygon (6250000 m^2) in southwestern Alberta, Canada")
      sim$studyAreaLarge <- randomStudyArea(seed = 1234, size = (250^2)*100)
    }
  }

  if (!suppliedElsewhere("studyAreaReporting", sim)) {
    if (suppliedElsewhere("studyArea", sim) && !is.null(sim$studyArea)) {
      message("'studyAreaReporting' was not provided by user. Using the same as 'studyArea'.")
      sim$studyAreaReporting <- sim$studyArea
    } else {
      message("'studyAreaReporting' was not provided by user. Using the same as 'studyAreaLarge'.")
      sim$studyAreaReporting <- sim$studyAreaLarge
    }
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
    if (!suppliedElsewhere("rawBiomassMap", sim)) {
      url <- paste0("http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                    "canada-forests-attributes_attributs-forests-canada/2001-attributes_attributs-2001/")
      fileURLs <- getURL(url, dirlistonly = TRUE)
      fileNames <- getHTMLLinks(fileURLs)
      rawBiomassMapFilename <- grep("Biomass_TotalLiveAboveGround.*.tif$", fileNames, value = TRUE)
      rawBiomassMapURL <- paste0(url, rawBiomassMapFilename)

      rawBiomassMap <- Cache(prepInputs,
                                 targetFile = rawBiomassMapFilename,
                                 url = rawBiomassMapURL,
                                 destinationPath = dPath,
                                 studyArea = sim$studyAreaLarge,   ## Ceres: makePixel table needs same no. pixels for this, RTM rawBiomassMap, LCC.. etc
                                 # studyArea = sim$studyArea,
                                 rasterToMatch = if (!needRTM) sim$rasterToMatchLarge else NULL,
                                 # maskWithRTM = TRUE,    ## if RTM not supplied no masking happens (is this intended?)
                                 maskWithRTM = if (!needRTM) TRUE else FALSE,
                                 ## TODO: if RTM is not needed use SA CRS? -> this is not correct
                                 # useSAcrs = if (!needRTM) TRUE else FALSE,
                                 useSAcrs = FALSE,     ## never use SA CRS
                                 method = "bilinear",
                                 datatype = "INT2U",
                                 filename2 = NULL,
                                 userTags = c(cacheTags, "rawBiomassMap"),
                                 omitArgs = c("destinationPath", "targetFile", "userTags", "stable"))
    }

    ## if we need rasterToMatchLarge, that means a) we don't have it, but b) we will have rawBiomassMap
    if (is.null(sim$rasterToMatchLarge))
      warning(paste0("rasterToMatchLarge is missing and will be created \n",
                     "from rawBiomassMap and studyAreaLarge.\n
              If this is wrong, provide raster"))

    sim$rasterToMatchLarge <- rawBiomassMap
    RTMvals <- getValues(sim$rasterToMatchLarge)
    sim$rasterToMatchLarge[!is.na(RTMvals)] <- 1
    sim$rasterToMatchLarge <- Cache(writeOutputs, sim$rasterToMatchLarge,
                                    filename2 = file.path(cachePath(sim), "rasters", "rasterToMatchLarge.tif"),
                                    datatype = "INT2U", overwrite = TRUE,
                                    userTags = c(cacheTags, "rasterToMatchLarge"),
                                    omitArgs = c("userTags"))

    ## this is old, and potentially not needed anymore
    if (FALSE) {
      studyArea <- sim$studyArea # temporary copy because it will be overwritten if it is suppliedElsewhere
      message("  Rasterizing the studyArea polygon map")
      if (!is(studyArea, "SpatialPolygonsDataFrame")) {
        dfData <- if (is.null(rownames(studyArea))) {
          polyID <- sapply(slot(studyArea, "polygons"), function(x) slot(x, "ID"))
          data.frame("field" = as.character(seq_along(length(studyArea))), row.names = polyID)
        } else {
          polyID <- sapply(slot(studyArea, "polygons"), function(x) slot(x, "ID"))
          data.frame("field" = rownames(studyArea), row.names = polyID)
        }
        studyArea <- SpatialPolygonsDataFrame(studyArea, data = dfData)
      }
      if (!identical(crs(studyArea), crs(sim$rasterToMatch))) {
        studyArea <- spTransform(studyArea, crs(sim$rasterToMatch))
        studyArea <- fixErrors(studyArea)

        ## TODO: OVERWRITE sim$studyArea here? what about SAlarge?
      }


      #TODO: review whether this is necessary (or will break LandWeb if removed) see Git Issue #22
      # layers provided by David Andison sometimes have LTHRC, sometimes LTHFC ... chose whichever
      LTHxC <- grep("(LTH.+C)", names(studyArea), value = TRUE)
      fieldName <- if (length(LTHxC)) {
        LTHxC
      } else {
        if (length(names(studyArea)) > 1) {
          ## study region may be a simple polygon
          names(studyArea)[1]
        } else NULL
      }

      sim$rasterToMatch <- crop(fasterizeFromSp(studyArea, sim$rasterToMatch, fieldName),
                                studyArea)
      sim$rasterToMatch <- Cache(writeRaster, sim$rasterToMatch,
                                 filename = file.path(dataPath(sim), "rasterToMatch.tif"),
                                 datatype = "INT2U", overwrite = TRUE,
                                 userTags = c(cacheTags, "rasterToMatch"),
                                 omitArgs = c("userTags"))
    }
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
