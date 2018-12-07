# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects, use sim$xxx, and are thus globally available
# to all modules. Functions can be used without sim$ as they are namespaced, like functions
# in R packages. If exact location is required, functions will be: sim$<moduleName>$FunctionName
defineModule(sim, list(
  name = "BiomassSpeciesData",
  description = "Download and pre-process proprietary LandWeb data.",
  keywords = c("LandWeb", "LandR"),
  authors = c(person(c("Eliot", "J", "B"), "McIntire", email = "eliot.mcintire@canada.ca", role = c("aut", "cre")),
              person(c("Alex", "M."), "Chubaty", email = "achubaty@friresearch.ca", role = c("aut")),
              person("Ceres", "Barros", email = "cbarros@mail.ubc.ca", role = c("aut"))),
  childModules = character(0),
  version = list(SpaDES.core = "0.1.2", BiomassSpeciesData = "0.0.1"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "BiomassSpeciesData.Rmd"),
  reqdPkgs = list("data.table", "googledrive", "gdalUtils", "magrittr", "pryr", "raster", ## TODO: is gdalUtils actually used?
                  "reproducible", "SpaDES.core", "SpaDES.tools",
                  "PredictiveEcology/pemisc@development"),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description"),
    defineParameter("types", "character", "KNN", NA, NA,
                    "The possible data sources. These must correspond to a function named paste0('prepSpeciesLayers_', type)"),
    defineParameter("speciesEquivalencyColumn", "character", "LandR", NA, NA,
                    "The column in sim$specieEquivalency data.table to use as a naming convention"),
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA,
                    "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", NA, NA, NA,
                    "This describes the simulation time interval between plot events"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA,
                    "This describes the simulation time at which the first save event should occur"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA,
                    "This describes the simulation time interval between save events"),
    defineParameter(".useCache", "logical", TRUE, NA, NA,
                    paste("Should this entire module be run with caching activated?",
                          "This is generally intended for data-type modules, where stochasticity and time are not relevant")),
    defineParameter(".useParallel", "numeric", parallel::detectCores(), NA, NA,
                    "Used in reading csv file with fread. Will be passed to data.table::setDTthreads")
  ),
  inputObjects = bind_rows(
    # expectsInput("CASFRIRas", "RasterStack",
    #              desc = "biomass percentage raster layers by species in Canada species map, created by Pickell et al., UBC, resolution 100m x 100m from LandSat and kNN based on CASFRI.",
    #              sourceURL = "https://drive.google.com/file/d/1y0ofr2H0c_IEMIpx19xf3_VTBheY0C9h/view?usp=sharing"),
    # expectsInput("Pickell", "RasterStack",
    #              desc = paste("biomass percentage raster layers by species in Canada species map,",
    #                           "created by Pickell et al., UBC, resolution 100m x 100m from LandSat and kNN based on CASFRI."),
    #              sourceURL = "https://drive.google.com/open?id=1M_L-7ovDpJLyY8dDOxG3xQTyzPx2HSg4"),
    expectsInput("speciesLayers", "RasterStack",
                 desc = "biomass percentage raster layers by species in Canada species map",
                 sourceURL = "http://tree.pfc.forestry.ca/kNN-Species.tar"),
    expectsInput("rasterToMatch", "RasterLayer",
                 desc = "Raster layer of study area used for cropping, masking and projecting.
                 Defaults to the kNN biomass map masked with `studyArea`",
                 sourceURL = "http://tree.pfc.forestry.ca/kNN-StructureBiomass.tar"),
    expectsInput("speciesEquivalency", "data.table",
                 desc = "table of species equivalencies. See pemisc::sppEquivalencies_CA.",
                 sourceURL = ""),
    expectsInput("sppMerge", c("list"),
                 desc = paste("list of kNN species layers that should be merged.",
                              "If none, create an empty list. Defaults to merging",
                              "of Pinus contorta and P. banksiana into Pinus sp."),
                 sourceURL = ""),
    expectsInput("sppNameVector", "character",
                 desc = "vector of species to select", sourceURL = ""),
    expectsInput("studyArea", "SpatialPolygonsDataFrame",
                 desc =  paste("Multipolygon to use as the study area,",
                               "Defaults to a square shapefile in Southwestern Alberta, Canada."),
                 sourceURL = NA),
    expectsInput("studyAreaLarge", "SpatialPolygonsDataFrame",
                 desc = paste("multipolygon (larger area than studyArea) to use for parameter estimation.",
                              "Defaults to a square shapefile in Southwestern Alberta, Canada."),
                 sourceURL = NA)
  ),
  outputObjects = bind_rows(
    createsOutput("speciesLayers", "RasterStack",
                  desc = "biomass percentage raster layers by species in Canada species map"),
    createsOutput("sppNameVector",  "character",
                  desc = "vector of species names to select.")
  )
  ))

## event types
#   - type `init` is required for initialiazation

doEvent.BiomassSpeciesData <- function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      sim <- biomassDataInit(sim)
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

### template initialization
biomassDataInit <- function(sim) {
  cacheTags <- c(currentModule(sim), "function:biomassDataInit")
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)))

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
                              studyArea = sim$studyArea,
                              rasterToMatch = sim$rasterToMatch,
                              sppNameVector = sim$sppNameVector,
                              speciesEquivalency = sim$speciesEquivalency,
                              speciesEquivalencyColumn = P(sim)$speciesEquivalencyColumn,
                              sppMerge = sim$sppMerge,
                              userTags = cacheTags)
    sim$speciesLayers <- if (length(sim$speciesLayers) > 0) {
      overlayStacks(highQualityStack = speciesLayersNew,
                    lowQualityStack = sim$speciesLayers,
                    destinationPath = dPath)
    } else {
      speciesLayersNew
    }
    rm(speciesLayersNew)

  }

  singular <- length(P(sim)$types) == 1
  message("sim$speciesLayers is from ", paste(P(sim)$types, collapse = ", "),
          "overlaid in that sequence, higher quality last"[!singular])

  return(invisible(sim))
}

.inputObjects <- function(sim) {
  dPath <- asPath(dataPath(sim), 1)
  cacheTags <- c(currentModule(sim), "function:.inputObjects")

  if (!suppliedElsewhere("studyArea", sim)) {
    message("'studyArea' was not provided by user. Using a polygon in southwestern Alberta, Canada,")

    sim$studyArea <- randomStudyArea(seed = 1234)
  }

  if (!suppliedElsewhere("studyAreaLarge", sim)) {
    message("'studyAreaLarge' was not provided by user. Using the same as 'studyArea'.")
    sim$studyAreaLarge <- sim$studyArea
  }

  if (is.null(sim$rasterToMatch)) {
    if (!suppliedElsewhere("rasterToMatch", sim)) {
      message("There is no rasterToMatch supplied; will attempt to use biomassMap")

      biomassMap <- Cache(prepInputs,
                          targetFile = asPath(basename(biomassMapFilename)),
                          archive = asPath(c("kNN-StructureBiomass.tar",
                                             "NFI_MODIS250m_kNN_Structure_Biomass_TotalLiveAboveGround_v0.zip")),
                          url = extractURL("rasterToMatch"),
                          destinationPath = dPath,
                          studyArea = sim$studyAreaLarge,   ## TODO: should this be studyAreaLarge? in RTM below it is...
                          useSAcrs = TRUE,
                          method = "bilinear",
                          datatype = "INT2U",
                          filename2 = TRUE, overwrite = TRUE,
                          userTags = cacheTags)

      sim$rasterToMatch <- biomassMap
      message("  Rasterizing the studyAreaLarge polygon map")
      #TODO: check whether this LandWeb centric stuf is necessary. Does rasterToMatch need FRI? see Issue #10
      # Layers provided by David Andison sometimes have LTHRC, sometimes LTHFC ... chose whichever
      LTHxC <- grep("(LTH.+C)",names(sim$studyAreaLarge), value = TRUE)
      fieldName <- if (length(LTHxC)) {
        LTHxC
      } else {
        if (length(names(sim$studyAreaLarge)) > 1) {   ## study region may be a simple polygon
          names(sim$studyAreaLarge)[1]
        } else NULL
      }

      sim$rasterToMatch <- crop(fasterizeFromSp(sim$studyAreaLarge, sim$rasterToMatch, fieldName),
                                sim$studyAreaLarge)
      sim$rasterToMatch <- Cache(writeRaster, sim$rasterToMatch,
                                 filename = file.path(dataPath(sim), "rasterToMatch.tif"),
                                 datatype = "INT2U", overwrite = TRUE)
    } else {
      stop("rasterToMatch is going to be supplied, but ", currentModule(sim), " requires it ",
           "as part of its .inputObjects. Please make it accessible to ", currentModule(sim),
           " in the .inputObjects by passing it in as an object in simInit(objects = list(rasterToMatch = aRaster)",
           " or in a module that gets loaded prior to ", currentModule(sim))
    }
  }

  if (!suppliedElsewhere("sppNameVector", sim)) {
    ## default to 6 species (see below)
    sim$sppNameVector <- c("Abie_sp", "Pice_gla", "Pice_mar", "Pinu_ban", "Pinu_con", "Popu_tre")
  }

  if (!suppliedElsewhere("sppMerge", sim)) {
    ## two merged into one (Pinu_ban, Pinu_con to Pinu_sp)
    sim$sppMerge <- list(Pinu_sp = c("Pinu_Ban", "Pinu_Con"))
  }

  if (!suppliedElsewhere("speciesEquivalency", sim)) {
    data("sppEquivalencies_CA", package = "pemisc", envir = environment())
    sim$speciesEquivalency <- as.data.table(sppEquivalencies_CA)

    ## By default, Abies_las is renamed to Abies_sp
    sim$speciesEquivalency[KNN == "Abie_Las", LandR := "Abie_sp"]

    ## add default colors for species used in model
    defaultCols <- RColorBrewer::brewer.pal(6, "Accent")
    LandRNames <- c("Pice_mar", "Pice_gla", "Popu_tre", "Pinu_sp", "Abie_sp")
    sim$speciesEquivalency[LandR %in% LandRNames, cols := defaultCols[-4]]
    sim$speciesEquivalency[EN_generic_full == "Mixed", cols := defaultCols[4]]
  }

  return(invisible(sim))
}

prepSpeciesLayers_CASFRI <- function(destinationPath, outputPath,
                                     url = "https://drive.google.com/file/d/1y0ofr2H0c_IEMIpx19xf3_VTBheY0C9h/view?usp=sharing",
                                     studyArea, rasterToMatch,
                                     sppNameVector,
                                     speciesEquivalency,
                                     speciesEquivalencyColumn,
                                     sppMerge) {
  CASFRItiffFile <- asPath(file.path(destinationPath, "Landweb_CASFRI_GIDs.tif"))
  CASFRIattrFile <- asPath(file.path(destinationPath, "Landweb_CASFRI_GIDs_attributes3.csv"))
  CASFRIheaderFile <- asPath(file.path(destinationPath,"Landweb_CASFRI_GIDs_README.txt"))

  message("  Loading CASFRI layers...")
  CASFRIRas <- Cache(prepInputs,
                     #targetFile = asPath("Landweb_CASFRI_GIDs.tif"),
                     targetFile = basename(CASFRItiffFile),
                     archive = asPath("CASFRI for Landweb.zip"),
                     url = url,
                     alsoExtract = c(CASFRItiffFile, CASFRIattrFile, CASFRIheaderFile),
                     destinationPath = destinationPath,
                     fun = "raster::raster",
                     studyArea = studyArea,
                     rasterToMatch = rasterToMatch,
                     method = "bilinear", ## ignore warning re: ngb (#5)
                     datatype = "INT4U",
                     filename2 = NULL, #TRUE,
                     overwrite = TRUE,
                     userTags =  c("CASFRIRas", "stable"))

  message("Load CASFRI data and headers, and convert to long format, and define species groups")
  #if (P(sim)$.useParallel > 1) data.table::setDTthreads(P(sim)$.useParallel)

  loadedCASFRI <- Cache(loadCASFRI,
                        CASFRIRas = CASFRIRas,
                        attrFile = CASFRIattrFile,
                        headerFile = CASFRIheaderFile, ## TODO: this isn't used internally
                        sppNameVector = pemisc::equivalentName(sppNameVector, speciesEquivalency, "CASFRI"),
                        speciesEquivalency = speciesEquivalency,
                        sppEndNamesCol = speciesEquivalencyColumn,
                        sppMerge = sppMerge,
                        userTags = c("function:loadCASFRI", "BigDataTable",
                                     "speciesLayers", "KNN"))

  #uniqueKeepSp <- unique(loadedCASFRI$keepSpecies$spGroup) ## TODO: the names don't match at all (#6)
  #if (!all(names(sim$speciesLayers) %in% uniqueKeepSp))
  #  warning("some kNN species not in CASFRI layers.")

  message('Make stack from CASFRI data and headers')
  CASFRISpStack <- #Cache(CASFRItoSpRasts,
    CASFRItoSpRasts(CASFRIRas = CASFRIRas,
                         #loadedCASFRI = loadedCASFRI,
                         CASFRIattrLong = loadedCASFRI$CASFRIattrLong,
                         keepSpecies = loadedCASFRI$keepSpecies,
                         CASFRIdt = loadedCASFRI$CASFRIdt,
                         #speciesLandR = pemisc::equivalentName(sppNameVector, speciesEquivalency, "LandR"), # don't want this because the spMerges are not here - use keepSpecies object
                         destinationPath = outputPath#,
                         #userTags = c("function:CASFRItoSpRasts", "CASFRIstack")
                         )

  return(CASFRISpStack)
}

prepSpeciesLayers_KNN <- function(destinationPath, outputPath,
                                  url = "http://tree.pfc.forestry.ca/kNN-Species.tar",
                                  studyArea, rasterToMatch,
                                  sppNameVector,
                                  speciesEquivalency,
                                  speciesEquivalencyColumn,
                                  sppMerge) {
  loadkNNSpeciesLayers(
    dPath = destinationPath,
    rasterToMatch = rasterToMatch,
    studyArea = studyArea,
    sppNameVector = sppNameVector,
    speciesEquivalency = speciesEquivalency,
    sppMerge = sppMerge,
    knnNamesCol = "KNN",
    sppEndNamesCol = speciesEquivalencyColumn,
    thresh = 10,
    url = url,
    userTags = c("speciesLayers", "KNN"))
}


prepSpeciesLayers_Pickell <- function(destinationPath, outputPath,
                                      url = "https://drive.google.com/open?id=1M_L-7ovDpJLyY8dDOxG3xQTyzPx2HSg4",
                                      studyArea, rasterToMatch,
                                      sppNameVector,
                                      speciesEquivalency,
                                      speciesEquivalencyColumn,
                                      sppMerge) {

  speciesLayers <- Cache(prepInputs,
                         targetFile = asPath("SPP_1990_100m_NAD83_LCC_BYTE_VEG_NO_TIES_FILLED_FINAL.dat"),
                         url = url,
                         archive = asPath("SPP_1990_100m_NAD83_LCC_BYTE_VEG_NO_TIES_FILLED_FINAL.zip"),
                         alsoExtract = asPath("SPP_1990_100m_NAD83_LCC_BYTE_VEG_NO_TIES_FILLED_FINAL.hdr"),
                         destinationPath = destinationPath,
                         fun = "raster::raster",
                         studyArea = studyArea,
                         rasterToMatch = rasterToMatch,
                         method = "bilinear", ## ignore warning re: ngb (#5)
                         datatype = "INT2U",
                         filename2 = NULL,
                         overwrite = TRUE,
                         userTags = c("speciesLayers", "KNN", "Pickell", "stable"))

  makePickellStack(PickellRaster = speciesLayers,
                   sppNameVector = sppNameVector,
                   speciesEquivalency = speciesEquivalency,
                   speciesEquivalencyColumn = speciesEquivalencyColumn,
                   sppMerge = sppMerge,
                   destinationPath = destinationPath)

}


prepSpeciesLayers_ForestInventory <- function(destinationPath, outputPath,
                                      url = "https://drive.google.com/file/d/1JnKeXrw0U9LmrZpixCDooIm62qiv4_G1/view?usp=sharing",
                                      studyArea, rasterToMatch,
                                      sppNameVector,
                                      speciesEquivalency,
                                      speciesEquivalencyColumn,
                                      sppMerge) {

  # This includes LandType because it will use that at the bottom of this function to
  #  remove NAs
  CClayerNames <- c("Pine", "Black Spruce", "Deciduous", "Fir", "White Spruce", "LandType")
  CClayerNamesWDots <- gsub(" ", ".", CClayerNames)
  CClayerNamesLandR <- equivalentName(CClayerNamesWDots, speciesEquivalency, speciesEquivalencyColumn)
  CClayerNamesFiles <- paste0(gsub(" ", "", CClayerNames), "1.tif")

  options(map.useParallel = FALSE)
  ml <- mapAdd(rasterToMatch, isRasterToMatch = TRUE, layerName = "rasterToMatch",
               #useSAcrs = TRUE, #poly = TRUE,
               #      columnNameForLabels = "NSN",
               filename2 = NULL)
  ml <- mapAdd(studyArea, map = ml, isStudyArea = TRUE, layerName = "studyArea",
                     useSAcrs = TRUE, #poly = TRUE,
               #      columnNameForLabels = "NSN",
               filename2 = NULL)

  ml <- mapAdd(map = ml, url = url, layerName = CClayerNames, CC = TRUE,
               destinationPath = destinationPath,
               #studyArea = studyArea,
               #rasterToMatch = rasterToMatch,
               targetFile = CClayerNamesFiles, filename2 = NULL, ## TODO: check this for file creation sadness
               alsoExtract = "similar", leaflet = FALSE, method = "ngb")

  ccs <- ml@metadata[CC == TRUE & !(layerName == "LandType"), ]
  CCs <- maps(ml, layerName = ccs$layerName)
  CCstack <- raster::stack(CCs)
  #CCstack[NA_ids] <- NA
  CCstack[CCstack[] < 0] <- 0
  CCstack[CCstack[] > 10] <- 10
  CCstack <- CCstack * 10 # convert back to percent
  NA_ids <- which(is.na(ml$LandType[]) | ml$LandType[] == 5)
  CCstack[NA_ids] <- NA

  names(CCstack) <- equivalentName(names(CCstack), speciesEquivalency, "LandWeb")

  CCstack
}
