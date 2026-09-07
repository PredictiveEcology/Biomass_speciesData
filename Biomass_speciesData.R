# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects, use sim$xxx, and are thus globally available
# to all modules. Functions can be used without sim$ as they are namespaced, like functions
# in R packages. If exact location is required, functions will be: sim$<moduleName>$FunctionName
defineModule(sim, list(
  name = "Biomass_speciesData",
  description = paste("Download and pre-process species percent cover raster data, overlaying",
                      "lower quality data with higher quality data."),
  keywords = c("LandWeb", "LandR", "LandR Biomass", "species percent cover"),
  authors = c(
    person("Ceres", "Barros", email = "ceres.barros@ubc.ca", role = c("aut", "cre")),
    person(c("Eliot", "J", "B"), "McIntire", email = "eliot.mcintire@nrcan-rncan.gc.ca", role = c("aut")),
    person(c("Alex", "M."), "Chubaty", email = "achubaty@for-cast.ca", role = c("aut"))
  ),
  childModules = character(0),
  version = list(Biomass_speciesData = "1.0.6"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "Biomass_speciesData.Rmd"),
  loadOrder = list(before = c("Biomass_borealDataPrep", "Biomass_core")),
  reqdPkgs = list("data.table", "RCurl",
                  "sf", "terra", "XML",
                  "reproducible (>= 2.1.0)",
                  "SpaDES.core (>= 2.1.4)", "SpaDES.tools (>= 1.0.2)",
                  "PredictiveEcology/LandR@development (>= 1.1.5.9063)",
                  "PredictiveEcology/pemisc@development"),
  parameters = bindrows(
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description"),
    defineParameter("coverThresh", "integer", 10L, NA, NA,
                    desc = paste("The minimum % cover a species needs to have (per pixel) in the study",
                                 "area to be considered present")),
    defineParameter("dataYear", "numeric", 2020, NA, NA,
                    paste("Passed to `paste0('prepSpeciesLayers_', types)` function to fetch data",
                          "from that year (if applicable). Defaults to 2020 as the default SCANFI year.")),
    defineParameter("sppEquivCol", "character", "LandR", NA, NA,
                    desc = paste("The column in `sim$sppEquiv` data.table to group species by and use as a",
                                 "naming convention. If different species in, e.g., the SCANFI data have the same",
                                 "name in the chosen column, their data are merged into one species by summing",
                                 "their percent cover in each raster cell.")),
    defineParameter("types", "character", "SCANFI", NA, NA,
                    desc = paste(
                      "The possible data sources. These must correspond to a function named",
                      "`paste0('prepSpeciesLayers_', types)`. Defaults to 'SCANFI'.",
                      "Other currently available options are:",
                      "'CASFRI', 'ForestInventory', 'KNN', 'MBFRI', 'NTEMS', 'ONFRI', 'Pickell'.",
                      "All non-KNN datasets attempt to get proprietary data; the user must be granted access first.",
                      "A custom function can be used to retrieve any data, just as long as it is",
                      "accessible by the module (e.g., in the global environment) and is named as",
                      "`paste0('prepSpeciesLayers_', types)`."
                    )),
    defineParameter("vegLeadingProportion", "numeric", 0.8, 0, 1,
                    desc = "a number that defines whether a species is leading for a given pixel. Only used for plotting."),
    defineParameter(".plotInitialTime", "numeric", start(sim), NA, NA,
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
                    desc = paste(
                      "Passed to `httr::config(ssl_verifypeer = P(sim)$.sslVerify)` when downloading NFI datasets.",
                      "Set to 0L if necessary to bypass checking the SSL certificate",
                      "(this may be necessary when NFI's website SSL certificate is not correctly configured)."
                    )),
    defineParameter(".studyAreaName", "character", NA, NA, NA,
                    paste("Human-readable name for the study area used.",
                          "If `NA`, a hash of `studyArea_biomassParam` will be used.")),
    defineParameter(".useCache", "character", "init", NA, NA,
                    desc = "Controls cache; caches the init event by default"),
    defineParameter(".useParallel", "numeric", 2L, NA, NA,
                    desc = "Used in reading csv file with fread. Will be passed to `data.table::setDTthreads`.")
  ),
  inputObjects = bindrows(
    expectsInput("rasterToMatch", "SpatRaster",
                 desc = paste("conditionally used as template raster if studyArea_rasterToMatch_biomassParam",
                              "and rasterTomatch_biomassParam are not supplied")),
    expectsInput("rasterToMatch_biomassParam", "SpatRaster",
                 desc = paste(
                   "a raster of `studyArea_biomassParam` in the same resolution and projection as",
                   "the simulation's."
                 )),
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
    expectsInput("studyArea", "sf",
                 desc = paste("`sf` polygon or terra `SpatVector` to use as the study area - `nrow` must be one")),
    expectsInput("studyArea_biomassParam", "sf",
                 desc =  paste(
                   "Polygon to use as the parametrisation study area. Must be provided by the user.",
                   "Note that `studyArea_biomassParam` is only used for parameter estimation, and",
                   "can be larger than the actual study area used for LandR simulations",
                   "(e.g., larger than `studyArea` in LandR `Biomass_core`)."),
                 sourceURL = NA),
    expectsInput("studyAreaReporting", "sf",
                 desc = paste("multipolygon (typically smaller/unbuffered than `studyArea_biomassParam` and `studyArea`",
                              "in LandR `Biomass_core`) to use for plotting/reporting.",
                              "If not provided, will default to `studyArea_biomassParam`."),
                 sourceURL = NA)
  ),
  outputObjects = bindrows(
    createsOutput("speciesLayers", "SpatRaster",
                  desc = "biomass percentage raster layers by species in Canada species map")
  )
))

## event types
#   - type `init` is required for initialiazation

doEvent.Biomass_speciesData <- function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      sim <- scheduleEvent(sim, start(sim), "Biomass_speciesData", "initPlot", eventPriority = .first())

      sim <- biomassDataInit(sim)
    },
    initPlot = {
      ## TODO: use Plots() here to allow saving of the maps to png etc.
      if (anyPlotting(P(sim)$.plots)) {
        # browser()
        plt <- plotVTM(
          speciesStack = mask(sim$speciesLayers, sim$studyAreaReporting),
          vegLeadingProportion = P(sim)$vegLeadingProportion,
          sppEquiv = sim$sppEquiv,
          sppEquivCol = P(sim)$sppEquivCol,
          colors = sim$sppColorVect,
          title = "Initial Types"
        )
        Plots(plt, filename = "Vegetation Type Map and Histogram")
      }
    },
    warning(paste(
      "Undefined event type: '",
      current(sim)[1, "eventType", with = FALSE],
      "' in module '",
      current(sim)[1, "moduleName", with = FALSE],
      "'",
      sep = ""
    ))
  )
  return(invisible(sim))
}

### template initialization
biomassDataInit <- function(sim) {
  cacheTags <- c(
    currentModule(sim),
    "otherFunctions:biomassDataInit",
    P(sim)$.studyAreaName,
    P(sim)$dataYear
  )
  dPath <- asPath(inputPath(sim), 1)
  message(currentModule(sim), ": biomassInit() using dataPath '", dPath, "'.")

  if (!exists("speciesLayers", envir = envir(sim), inherits = FALSE)) {
    sim$speciesLayers <- list()
  }

  for (type in P(sim)$types) {
    fnName <- paste0("prepSpeciesLayers_", type)
    whereIsFnName <- where(fnName)

    envirName <- attr(whereIsFnName, "name")
    if (is.null(envirName)) {
      envirName <- environmentName(whereIsFnName)
    }
    if (!is.character(envirName)) {
      ## this is from "box" package; slightly different
      envirName <- environmentName(envirName)
    }

    message("#############################################")
    message(type, " -- Loading using ", fnName, " located in ", envirName)
    message("#############################################")
    if (!exists(fnName)) {
      stop(fnName, " does not exist. Please make it accessible in a package, as an object, ",
           " or in the .GlobalEnv")
    }
    fn <- get(fnName)
    httr::with_config(config = httr::config(ssl_verifypeer = P(sim)$.sslVerify), {
      speciesLayersNew <- fn(
        destinationPath = dPath, # this is generic files (preProcess)
        outputPath = outputPath(sim), # this will be the studyArea-specific files (postProcess)
        to = sim$studyArea_biomassParam,
        studyAreaName = P(sim)$.studyAreaName,
        projectTo = sim$rasterToMatch_biomassParam,
        sppEquiv = sim$sppEquiv,
        sppEquivCol = P(sim)$sppEquivCol,
        thresh = P(sim)$coverThresh,
        year = P(sim)$dataYear,     # read from `...` by prepSpeciesLayers_KNN()
        dataYear = P(sim)$dataYear, # formal of prepSpeciesLayers_SCANFI() / _NTEMS(); was never reaching them
        .functionName = fnName,
        userTags = c(cacheTags, fnName, "prepSpeciesLayers"),
        omitArgs = c("userTags")
      ) ## |> Cache() ## Cache retrieving RasterList as SpatRaster, so subsequent use breaks
    })

    sim$speciesLayers <- if (length(sim$speciesLayers) > 0) {
      overlayStacks(
        highQualityStack = speciesLayersNew,
        lowQualityStack = sim$speciesLayers,
        destinationPath = outputPath(sim)
      ) |>
        Cache(userTags = c(cacheTags, "overlayStacks"), omitArgs = c("userTags"))
    } else {
      speciesLayersNew
    }
    rm(speciesLayersNew)
  }

  assertSpeciesLayers(sim$speciesLayers, P(sim)$coverThresh)

  species <- names(sim$speciesLayers)

  origFilenames <- vapply(
    names(sim$speciesLayers),
    function(r) Filenames(sim$speciesLayers[[r]], allowMultiple = FALSE),
    character(1)
  )

  ## re-enforce study area mask (merged/summed layers are losing the mask)
  sim$speciesLayers <- maskTo(sim$speciesLayers, sim$rasterToMatch_biomassParam)

  ## make sure empty pixels inside study area have 0 cover, instead of NAs.
  ## this can happen when data has NAs instead of 0s and is not merged/overlayed (e.g. CASFRI)
  tempRas <- sim$rasterToMatch_biomassParam
  tempRas[!is.na(tempRas[])] <- 0
  sim$speciesLayers <- raster::cover(sim$speciesLayers, tempRas)
  names(sim$speciesLayers) <- species
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
    message("Writing speciesLayers to disk...")
    out <- lapply(seq_along(names(sim$speciesLayers)), function(r) {
      writeRaster(sim$speciesLayers[[r]], filename = origFilenames[r], overwrite = TRUE)
    })
    message("      ... Done")
    out
  }

  if (is(sim$speciesLayers, "list")) {
    sim$speciesLayers <- .stack(sim$speciesLayers)
  }

  sim$speciesLayers <- setNames(sim$speciesLayers, species)

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
  cacheTags <- c(currentModule(sim), "otherFunctions:.inputObjects")
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")

  rtm_res <- 240 ## SCANFI is 30m resolution and would be aggregated to this

  ## Study area(s) ------------------------------------------------
  if (!suppliedElsewhere("studyArea", sim)) {
    sim$studyArea <- LandR::randomStudyArea(seed = 1234, size = (rtm_res^2) * 100)
  }

  if (!suppliedElsewhere("studyArea_biomassParam", sim)) {
    if (is.null(sim$studyAreaLarge)) {
      sim$studyArea_biomassParam <- sim$studyArea
    } else {
      warning("please replace studyAreaLarge with studyArea_biomassParam")
      sim$studyArea_biomassParam <- sim$studyAreaLarge
    }
  }

  if (is.na(P(sim)$.studyAreaName)) {
    params(sim)[[currentModule(sim)]][[".studyAreaName"]] <- reproducible::studyAreaName(sim$studyArea_biomassParam)
    message("The .studyAreaName is not supplied; derived name from sim$studyArea_biomassParam: ",
            params(sim)[[currentModule(sim)]][[".studyAreaName"]])
  }

  if (!suppliedElsewhere("studyAreaReporting", sim)) {
    message("'studyAreaReporting' was not provided by user. Using the same as 'studyArea_biomassParam'.")
    sim$studyAreaReporting <- sim$studyArea_biomassParam
  }

  if (!suppliedElsewhere("rasterToMatch", sim)) {
    studyArea <- sim$studyArea
    if (!inherits(studyArea, "SpatVector")) {
      studyArea <- terra::vect(studyArea)
    }
    if (terra::is.lonlat(studyArea)) {
      ## use SCANFI projection - LandR requires projected rasters for dispersal
      studyArea <- terra::project(
        studyArea,
        paste("+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77",
              "+x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
      )
    }
    sim$rasterToMatch <- terra::rast(studyArea, res = c(rtm_res, rtm_res), vals = 1) |>
      terra::mask(mask = studyArea)
  }

  if (!suppliedElsewhere("rasterToMatch_biomassParam", sim)) {
    if (!is.null(sim$rasterToMatchLarge)) {
      warning("please use rasterToMatch_biomassParam in place of rasterToMatchLarge")
      sim$rasterToMatch_biomassParam <- sim$rasterToMatchLarge
    # } else if (!terra::compareGeom(sim$studyArea_biomassParam, sim$studyArea)) {
    } else if (!LandR::.compareCRS(sim$studyArea_biomassParam, sim$studyArea)) {
      ## SA_BP was supplied but not RTM_BP
      sim$rasterToMatch_biomassParam <- terra::rast(
        sim$studyArea_biomassParam,
        resolution = terra::res(sim$rasterToMatch),
        vals = 1
      ) |>
        reproducible::postProcess(maskTo = sim$studyArea_biomassParam)
    } else {
      sim$rasterToMatch_biomassParam <- sim$rasterToMatch
    }
  }

  if (st_crs(sim$studyArea_biomassParam) != st_crs(sim$rasterToMatch_biomassParam)) {
    warning(paste0(
      "studyArea_biomassParam and rasterToMatch_biomassParam projections differ.\n",
      "studyArea_biomassParam will be projected to match rasterToMatch_biomassParam"
    ))
    sim$studyArea_biomassParam <- reproducible::projectTo(sim$studyArea_biomassParam, sim$rasterToMatch_biomassParam)
  }

  ## Species equivalencies table and associated columns ----------------------------
  ## make sppEquiv table and associated columns, vectors
  ## do not use suppliedElsewhere here as we need the tables to exist (or not)
  ## already (rather than potentially being supplied by a downstream module)
  ## the function checks whether the tables exist internally.
  ## check parameter consistency across modules
  paramCheckOtherMods(sim, "sppEquivCol", ifSetButDifferent = "error")
  paramCheckOtherMods(sim, "vegLeadingProportion", ifSetButDifferent = "error")

  sppOuts <- LandR::sppHarmonize(
    sppEquiv = sim$sppEquiv,
    sppNameVector = sim$sppNameVector,
    sppEquivCol = P(sim)$sppEquivCol,
    sppColorVect = sim$sppColorVect,
    vegLeadingProportion = P(sim)$vegLeadingProportion,
    studyArea = sim$studyArea_biomassParam,
    dPath = dPath
  )
  ## the following may, or may not change inputs
  sim$sppEquiv <- sppOuts$sppEquiv
  sim$sppNameVector <- sppOuts$sppNameVector
  P(sim, module = currentModule(sim))$sppEquivCol <- sppOuts$sppEquivCol
  sim$sppColorVect <- sppOuts$sppColorVect

  return(invisible(sim))
}


where <- function(name, env = parent.frame()) {
  while (!identical(env, emptyenv())) {
    if (exists(name, envir = env, inherits = FALSE)) {
      return(env)
    }
    env <- parent.env(env)
  }
  stop(sprintf("Can't find '%s'", name), call. = FALSE)
}