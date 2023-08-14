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
  version = list(Biomass_speciesData = "1.0.2"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "Biomass_speciesData.Rmd"),
  reqdPkgs = list("data.table", "gdalUtilities", ## LandR needs gdalUtilities to overlay rasters
                  "sf", "magrittr", "pryr", "RCurl", "reproducible (>= 2.0.2)", "terra",
                  "SpaDES.core (>= 2.0.2)", "XML",
                  # "curl", "httr", ## called directly by this module, but pulled in by LandR (Sep 6th 2022).
                  ## Excluded because loading is not necessary (just installation)
                  "CeresBarros/LandR@deprecatedArgs (HEAD)",
                  "PredictiveEcology/pemisc@development",
                  "PredictiveEcology/SpaDES.tools@development (>= 1.0.2)"),
  parameters = bindrows(
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description"),
    defineParameter("coverThresh", "integer", 10L, NA, NA,
                    desc = paste("The minimum % cover a species needs to have (per pixel) in the study",
                          "area to be considered present")),
    defineParameter("dataYear", "numeric", 2001, NA, NA,
                    paste("Passed to `paste0('prepSpeciesLayers_', types)` function to fetch data",
                          "from that year (if applicable). Defaults to 2001 as the default kNN year.")),
    defineParameter("sppEquivCol", "character", "Boreal", NA, NA,
                    desc = paste("The column in `sim$sppEquiv` data.table to group species by and use as a",
                          "naming convention. If different species in, e.g., the kNN data have the same",
                          "name in the chosen column, their data are merged into one species by summing",
                          "their % cover in each raster cell.")),
    defineParameter("types", "character", "KNN", NA, NA,
                    desc = paste("The possible data sources. These must correspond to a function named",
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
                    desc = "a number that defines whether a species is leading for a given pixel. Only used for plotting."),
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA,
                    desc = "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", NA, NA, NA,
                    desc = "This describes the simulation time interval between plot events"),
    defineParameter(".plots", "character", c("screen"), NA, NA,
                    desc = paste("Passed to `types` in `Plots` (see `?Plots`).",
                          "There are a few plots that are made within this module, if set.",
                          "Note that plots (or their data) saving will ONLY occur at `end(sim)`.",
                          "If `NA`, plotting is turned off completely (this includes plot saving).")),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA,
                    desc = "This describes the simulation time at which the first save event should occur"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA,
                    "This describes the simulation time interval between save events"),
    defineParameter(".sslVerify", "integer", as.integer(unname(curl::curl_options("^ssl_verifypeer$"))), NA_integer_, NA_integer_,
                    desc = paste("Passed to `httr::config(ssl_verifypeer = P(sim)$.sslVerify)` when downloading KNN",
                          "(NFI) datasets. Set to 0L if necessary to bypass checking the SSL certificate (this",
                          "may be necessary when NFI's website SSL certificate is not correctly configured).")),
    defineParameter(".studyAreaName", "character", NA, NA, NA,
                    "Human-readable name for the study area used. If NA, a hash of `studyAreaLarge` will be used."),
    defineParameter(".useCache", "character", "init", NA, NA,
                    desc = "Controls cache; caches the init event by default"),
    defineParameter(".useParallel", "numeric", parallel::detectCores(), NA, NA,
                    desc = "Used in reading csv file with fread. Will be passed to `data.table::setDTthreads`.")
  ),
  inputObjects = bindrows(
    expectsInput("rasterToMatchLarge", "SpatRaster",
                 desc = paste("a raster of `studyAreaLarge` in the same resolution and projection the simulation's.",
                              "Defaults to the using the Canadian Forestry Service, National Forest Inventory,",
                              "kNN-derived stand biomass map."),
                 sourceURL = ""),
    # expectsInput("rawBiomassMap", "SpatRaster",
    #              desc = paste("total biomass raster layer in study area. Only used to create `rasterToMatchLarge`",
    #                           "if necessary. Defaults to the Canadian Forestry Service, National Forest Inventory,",
    #                           "kNN-derived total aboveground biomass map from 2001 (in tonnes/ha), unless",
    #                           "'dataYear' != 2001. See ",
    #                           "https://open.canada.ca/data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990",
    #                           "for metadata."),
    #              sourceURL = ""), ## sourceURL varies by `dataYear`
    expectsInput("sppColorVect", "character",
                 desc = paste("A named vector of colors to use for plotting.",
                              "The names must be in `sim$sppEquiv[[sim$sppEquivCol]]`,",
                              "and should also contain a color for 'Mixed'"),
                 sourceURL = NA),
    expectsInput("sppEquiv", "data.table",
                 desc = "table of species equivalencies. See `LandR::sppEquivalencies_CA`.",
                 sourceURL = ""),
    expectsInput("sppNameVector", "character",
                 desc = paste("an optional vector of species names to be pulled from `sppEquiv`. Species names must match",
                              "`P(sim)$sppEquivCol` column in `sppEquiv`. If not provided, then species will be taken from",
                              "the entire `P(sim)$sppEquivCol` column in `sppEquiv`.",
                              "See `LandR::sppEquivalencies_CA`.")),
    expectsInput("studyAreaLarge", "sfc",
                 desc =  paste("Polygon to use as the parametrisation study area. Must be provided by the user.",
                               "Note that `studyAreaLarge` is only used for parameter estimation, and",
                               "can be larger than the actual study area used for LandR simulations (e.g,",
                               "larger than `studyArea` in LandR Biomass_core)."),
                 sourceURL = NA),
    expectsInput("studyAreaReporting", "sfc",
                 desc = paste("multipolygon (typically smaller/unbuffered than `studyAreaLarge` and `studyArea`",
                              "in LandR Biomass_core) to use for plotting/reporting.",
                              "If not provided, will default to `studyAreaLarge`."),
                 sourceURL = NA)
  ),
  outputObjects = bindrows(
    createsOutput("speciesLayers", "SpatRaster",
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
      ## TODO: use Plots() here to allow saving of the maps to png etc.
      if (anyPlotting(P(sim)$.plots) && any(P(sim)$.plots == "screen")) {
        newDev <- if (!is.null(dev.list())) {
          devCur <- dev.cur()
          max(dev.list()) + 1
        } else {
          1
        }
        dev.set(newDev)

        plotVTM(speciesStack = mask(sim$speciesLayers, sim$studyAreaReporting),
                vegLeadingProportion = P(sim)$vegLeadingProportion,
                sppEquiv = sim$sppEquiv,
                sppEquivCol = P(sim)$sppEquivCol,
                colors = sim$sppColorVect,
                title = "Initial Types")
        if (exists("devCur")) dev.set(devCur)
      }
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
                                userTags = c(cacheTags, fnName, "prepSpeciesLayers", P(sim)$.studyAreaName),
                                omitArgs = c("userTags"))
    })

    sim$speciesLayers <- if (length(sim$speciesLayers) > 0) {
      Cache(overlayStacks,
            highQualityStack = speciesLayersNew,
            lowQualityStack = sim$speciesLayers,
            destinationPath = outputPath(sim),
            userTags = c(cacheTags, "overlayStacks", P(sim)$.studyAreaName),
            omitArgs = c("userTags"))
    } else {
      speciesLayersNew
    }
    rm(speciesLayersNew)
  }

  assertSpeciesLayers(sim$speciesLayers, P(sim)$coverThresh)

  species <- names(sim$speciesLayers)

  origFilenames <- vapply(names(sim$speciesLayers),
                          function(r) Filenames(sim$speciesLayers[[r]], allowMultiple = FALSE),
                          character(1))

  ## re-enforce study area mask (merged/summed layers are losing the mask)
  sim$speciesLayers <- maskTo(sim$speciesLayers, sim$rasterToMatchLarge)

  ## make sure empty pixels inside study area have 0 cover, instead of NAs.
  ## this can happen when data has NAs instead of 0s and is not merged/overlayed (e.g. CASFRI)
  tempRas <- sim$rasterToMatchLarge
  tempRas[!is.na(tempRas[])] <- 0
  sim$speciesLayers <- cover(sim$speciesLayers, tempRas)
  rm(tempRas)

  ## filter out species with no data, or too little cover (some prepSpeciesLayers_*/overlay are not doing this)
  layersWdata <- vapply(names(sim$speciesLayers), function(nn) {
    xx <- sim$speciesLayers[[nn]]
    if (maxFn(xx) < P(sim)$coverThresh) FALSE else TRUE
  }, logical(1))
  sppKeep <- names(sim$speciesLayers)[layersWdata]
  if (sum(!layersWdata) > 0) {
    if (length(sppKeep)) {
      message("removing ", sum(!layersWdata), " species because they had <", P(sim)$coverThreshresh,
              " % cover in the study area\n",
              "  These species are retained (and could be further culled manually, if desired):\n",
              paste(sppKeep, collapse = " "))
    } else {
      message("no pixels for ", paste(names(layersWdata), collapse = " "),
              " were found with >=", thresh, " % cover in the study area.",
              "\n  No species layers were retained. Try lowering the threshold",
              " to retain species with low % cover")
    }
  }
  sim$speciesLayers <- sim$speciesLayers[[sppKeep]]
  species <- sppKeep

  ## speciesLayers brick/stack may have filename but layers do not...
  if (nzchar(Filenames(sim$speciesLayers, allowMultiple = FALSE)) && !all(nzchar(origFilenames))) {
    sim$speciesLayers[] <- sim$speciesLayers[] ## bring to memory
  }

  sim$speciesLayers <- if (inMemory(sim$speciesLayers)) {
    sim$speciesLayers
  } else {
    lapply(seq_along(names(sim$speciesLayers)), function(r) {
      writeRaster(sim$speciesLayers[[r]], filename = origFilenames[r], overwrite = TRUE)
    })
  }

  if (is(sim$speciesLayers, "list")) {
    sim$speciesLayers <- .stack(sim$speciesLayers)
  }

  setNames(sim$speciesLayers, species)

  singular <- length(P(sim)$types) == 1
  message("sim$speciesLayers is from ", paste(P(sim)$types, collapse = ", "),
          " overlaid in that sequence, higher quality last"[!singular])

  message("------------------")
  sumCoverPerPixel <- as.vector(values(sum(sim$speciesLayers, na.rm = TRUE)))
  message("There are ", sum(!is.na(sumCoverPerPixel)), " pixels with non-NA tree cover in them")
  message("There are ", sum(sumCoverPerPixel > 0, na.rm = TRUE), " pixels with non-zero tree cover in them")

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

  needRTML <- FALSE
  if (is.null(sim$rasterToMatchLarge)) {
    if (!suppliedElsewhere("rasterToMatchLarge", sim)) {      ## if one is not provided, re do both (safer?)
      needRTML <- TRUE
      message("There is no rasterToMatchLarge supplied; will attempt to use rawBiomassMap, if it is there; otherwise will try KNN")
    } else {
      stop("rasterToMatchLarge is going to be supplied, but ", currentModule(sim), " requires it ",
           "as part of its .inputObjects. Please make it accessible to ", currentModule(sim),
           " in the .inputObjects by passing it in as an object in simInit(objects = list(rasterToMatchLarge = aRaster)",
           " or in a module that gets loaded prior to ", currentModule(sim))
    }
  }

  if (needRTML) {
    ## if rawBiomassMap exists, it needs to match SALarge, if it doesn't make it
    if (is.null(sim$rawBiomassMap)) {
      if (P(sim)$dataYear == 2001) {
        biomassURL <- paste0("http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                             "canada-forests-attributes_attributs-forests-canada/",
                             "2001-attributes_attributs-2001/",
                             "NFI_MODIS250m_2001_kNN_Structure_Biomass_TotalLiveAboveGround_v1.tif")
      } else {
        if (P(sim)$dataYear == 2011) {
          biomassURL <- paste0("http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                               "canada-forests-attributes_attributs-forests-canada/2011-attributes_attributs-2011/",
                               "NFI_MODIS250m_2011_kNN_Structure_Biomass_TotalLiveAboveGround_v1.tif")
        } else {
          stop("'P(sim)$dataYear' must be 2001 OR 2011")
        }
      }

      httr::with_config(config = httr::config(ssl_verifypeer = P(sim)$.sslVerify), {
        rawBiomassMap <- prepRawBiomassMap(url = biomassURL,
                                         studyAreaName = P(sim)$.studyAreaName,
                                         cacheTags = cacheTags,
                                         cropTo = sim$studyAreaLarge,
                                         maskTo = sim$studyAreaLarge,
                                         projectTo = NA,  ## don't project to SA
                                         destinationPath = dPath)
    })
  } else {
    rawBiomassMap <- sim$rawBiomassMap
    if (!.compareCRS(sim$rawBiomassMap, sim$studyAreaLarge)) {
      ## note that extents may never align if the resolution and projection do not allow for it
      rawBiomassMap <- Cache(postProcess,
                             rawBiomassMap,
                             method = "bilinear",
                             cropTo = sim$studyAreaLarge,
                             maskTo = sim$studyAreaLarge,
                             projectTo = NA,  ## don't project to SA
                               overwrite = TRUE)
      }
    }
  }

  RTMs <- prepRasterToMatch(studyArea = sim$studyAreaLarge,
                            studyAreaLarge = sim$studyAreaLarge,
                            rasterToMatch = NULL,
                            rasterToMatchLarge = if (needRTML) NULL else sim$rasterToMatchLarge,
                            destinationPath = dPath,
                            templateRas = rawBiomassMap,
                            studyAreaName = P(sim)$.studyAreaName,
                            cacheTags = cacheTags)
  sim$rasterToMatchLarge <- RTMs$rasterToMatchLarge
  rm(RTMs)


  if (st_crs(sim$studyAreaLarge) != st_crs(sim$rasterToMatchLarge)) {
    warning(paste0("studyAreaLarge and rasterToMatchLarge projections differ.\n",
                   "studyAreaLarge will be projected to match rasterToMatchLarge"))
    sim$studyAreaLarge <- projectTo(sim$studyAreaLarge, sim$rasterToMatchLarge)
  }

  ## Species equivalencies table and associated columns ----------------------------
  ## make sppEquiv table and associated columns, vectors
  ## do not use suppliedElsewhere here as we need the tables to exist (or not)
  ## already (rather than potentially being supplied by a downstream module)
  ## the function checks whether the tables exist internally.
  ## check parameter consistency across modules
  paramCheckOtherMods(sim, "sppEquivCol", ifSetButDifferent = "error")
  paramCheckOtherMods(sim, "vegLeadingProportion", ifSetButDifferent = "error")

  sppOuts <- sppHarmonize(sim$sppEquiv, sim$sppNameVector, P(sim)$sppEquivCol,
                          sim$sppColorVect, P(sim)$vegLeadingProportion, sim$studyAreaLarge)
  ## the following may, or may not change inputs
  sim$sppEquiv <- sppOuts$sppEquiv
  sim$sppNameVector <- sppOuts$sppNameVector
  P(sim, module = currentModule(sim))$sppEquivCol <- sppOuts$sppEquivCol
  sim$sppColorVect <- sppOuts$sppColorVect

  return(invisible(sim))
}
