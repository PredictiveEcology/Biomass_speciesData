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
  reqdPkgs = list("data.table", "googledrive", "gdalUtils", "magrittr", "raster", ## TODO: is gdalUtils actually used?
                  "reproducible", "SpaDES.core", "SpaDES.tools",
                  "PredictiveEcology/pemisc@development"),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description"),
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
    expectsInput("CASFRIRas", "RasterStack",
                 desc = "biomass percentage raster layers by species in Canada species map, created by Pickell et al., UBC, resolution 100m x 100m from LandSat and kNN based on CASFRI.",
                 sourceURL = "https://drive.google.com/file/d/1y0ofr2H0c_IEMIpx19xf3_VTBheY0C9h/view?usp=sharing"),
    expectsInput("Pickell", "RasterStack",
                 desc = paste("biomass percentage raster layers by species in Canada species map,",
                              "created by Pickell et al., UBC, resolution 100m x 100m from LandSat and kNN based on CASFRI."),
                 sourceURL = "https://drive.google.com/open?id=1M_L-7ovDpJLyY8dDOxG3xQTyzPx2HSg4"),
    expectsInput("rasterToMatch", "RasterLayer",
                 desc = "Raster layer of study area used for cropping, masking and projecting.
                 Defaults to the kNN biomass map masked with `studyArea`",
                 sourceURL = "http://tree.pfc.forestry.ca/kNN-StructureBiomass.tar"),
    expectsInput("speciesEquivalency", "data.table",
                 desc = "table of species equivalencies. See pemisc::sppEquivalencies_CA.",
                 sourceURL = ""),
    expectsInput("speciesLayers", "RasterStack",
                 desc = "biomass percentage raster layers by species in Canada species map",
                 sourceURL = "http://tree.pfc.forestry.ca/kNN-Species.tar"),
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
  dPath <- asPath(dataPath(sim))

  message("Load CASFRI data and headers, and convert to long format, and define species groups")
  if (P(sim)$.useParallel > 1) data.table::setDTthreads(P(sim)$.useParallel)

  loadedCASFRI <- Cache(loadCASFRI,
                        CASFRIRas = sim$CASFRIRas,
                        attrFile = mod$CASFRIattrFile,
                        headerFile = mod$CASFRIheaderFile, ## TODO: this isn't used internally
                        sppNameVector = sim$sppNameVector,
                        speciesEquivalency = sim$speciesEquivalency,
                        sppEndNamesCol = "LandR",
                        sppMerge = sim$sppMerge,
                        userTags = c(cacheTags, "function:loadCASFRI", "BigDataTable"))

  message("Make stack of species layers from Pickell's layer")

  ## check if all species found in kNN database are in CASFRI
  uniqueKeepSp <- unique(loadedCASFRI$keepSpecies$spGroup) ## TODO: the names don't match at all (#6)
  if (!all(names(sim$speciesLayers) %in% uniqueKeepSp))
    warning("some kNN species not in CASFRI layers.")

  ## TODO: weird bug when useCache = "overwrite"
  if (getOption("reproducible.useCache") == "overwrite")
    opt <- options(reproducible.useCache = TRUE) ## TODO: temporary workaround

  PickellSpStack <- Cache(makePickellStack,
                          PickellRaster = sim$Pickell,
                          uniqueKeepSp = uniqueKeepSp,
                          speciesKnn = names(sim$speciesLayers),
                          destinationPath = dPath,
                          userTags = c(cacheTags, "function:makePickellStack", "PickellStack"))

  if (exists("opt", inherits = FALSE)) options(opt) ## TODO: temporary workaround with above

  crs(PickellSpStack) <- crs(sim$rasterToMatch) # bug in writeRaster

  message('Make stack from CASFRI data and headers')
  CASFRISpStack <- Cache(CASFRItoSpRasts,
                         CASFRIRas = sim$CASFRIRas,
                         loadedCASFRI = loadedCASFRI,
                         speciesKnn = names(sim$speciesLayers),
                         destinationPath = dPath,
                         userTags = c(cacheTags, "function:CASFRItoSpRasts", "CASFRIstack"))

  message("Overlay Pickell and CASFRI stacks")
  outStack <- Cache(overlayStacks,
                    highQualityStack = CASFRISpStack,
                    lowQualityStack = PickellSpStack,
                    outputFilenameSuffix = "CASFRI_Pickell",
                    destinationPath = dPath,
                    userTags = c(cacheTags, "function:overlayStacks", "Pickell_CASFRI"))

  crs(outStack) <- crs(sim$rasterToMatch) # bug in writeRaster

  message("Overlay Pickell_CASFRI with open data set stacks")
  speciesLayers2 <- Cache(overlayStacks,
                          highQualityStack = outStack,
                          lowQualityStack = sim$speciesLayers,
                          outputFilenameSuffix = "CASFRI_Pickell_kNN",
                          destinationPath = dPath,
                          userTags = c(cacheTags, "function:overlayStacks", "CASFRI_Pickell_kNN"))
  crs(speciesLayers2) <- crs(sim$rasterToMatch)

  ## replace species layers
  sim$speciesLayers <- speciesLayers2
  message("Using overlaid datasets from CASFRI, Pickell and CFS kNN")

  return(invisible(sim))
}

.inputObjects <- function(sim) {
  dPath <- asPath(dataPath(sim), 1)
  cacheTags <- c(currentModule(sim), "function:.inputObjects")

  if (!suppliedElsewhere("studyAreaLarge", sim)) {
    message("'studyAreaLarge' was not provided by user. Using a polygon in southwestern Alberta, Canada,")

    sim$studyAreaLarge <- randomStudyArea(seed = 1234)
  }

  if (!suppliedElsewhere("studyArea", sim)) {
    message("'studyArea' was not provided by user. Using the same as 'studyAreaLarge'.")
    sim$studyArea <- sim$studyAreaLarge
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
                          studyArea = sim$studyArea,   ## TODO: should this be studyAreaLarge? in RTM below it is...
                          useSAcrs = TRUE,
                          method = "bilinear",
                          datatype = "INT2U",
                          filename2 = TRUE, overwrite = TRUE,
                          userTags = cacheTags)

      sim$rasterToMatch <- biomassMap
      message("  Rasterizing the studyAreaLarge polygon map")

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
    data("sppEquivalencies_CA", package = "pemisc")
    sim$speciesEquivalency <- as.data.table(sppEquivalencies_CA)

    ## By default, Abies_las is renamed to Abies_sp
    sim$speciesEquivalency[KNN == "Abie_Las", LandR := "Abie_sp"]

    ## add default colors for species used in model
    defaultCols <- RColorBrewer::brewer.pal(6, "Accent")
    LandRNames <- c("Pice_mar", "Pice_gla", "Popu_tre", "Pinu_sp", "Abie_sp")
    sim$speciesEquivalency[LandR == LandRNames, cols := defaultCols[-4]]
    sim$speciesEquivalency[EN_generic_full == "Mixed", cols := defaultCols[4]]
  }

  if (!suppliedElsewhere("speciesLayers")) {
    speciesLayersList <- Cache(loadkNNSpeciesLayers,
                               dPath = dPath,
                               rasterToMatch = sim$rasterToMatch,
                               studyArea = sim$studyAreaLarge,
                               sppNameVector = sim$sppNameVector,
                               speciesEquivalency = sim$speciesEquivalency,
                               sppMerge = sim$sppMerge,
                               knnNamesCol = "KNN",
                               sppEndNamesCol = "LandR",
                               thresh = 10,
                               url = extractURL("speciesLayers"),
                               userTags = c(cacheTags, "speciesLayers"))

    sim$speciesLayers <- speciesLayersList$speciesLayers

    ## update the list of original species to use;
    ## as kNN is the lowest resolution, if species weren't found they will be excluded
    sim$sppNameVector <- speciesLayersList$sppNameVector
  }

  if (!suppliedElsewhere("Pickell") | !suppliedElsewhere("CASFRIRas")) {
    if (!suppliedElsewhere("Pickell")) {
      message("  Loading Pickell et al. layers...")
      sim$Pickell <- Cache(prepInputs,
                           targetFile = asPath("SPP_1990_100m_NAD83_LCC_BYTE_VEG_NO_TIES_FILLED_FINAL.dat"),
                           url = extractURL("Pickell"),
                           archive = asPath("SPP_1990_100m_NAD83_LCC_BYTE_VEG_NO_TIES_FILLED_FINAL.zip"),
                           alsoExtract = asPath("SPP_1990_100m_NAD83_LCC_BYTE_VEG_NO_TIES_FILLED_FINAL.hdr"),
                           destinationPath = dPath,
                           fun = "raster::raster",
                           studyArea = sim$studyArea,
                           rasterToMatch = sim$rasterToMatch,
                           method = "bilinear", ## ignore warning re: ngb (#5)
                           datatype = "INT2U",
                           filename2 = TRUE,
                           overwrite = TRUE,
                           userTags = c(cacheTags, "Pickell", "stable"))
    }

    if (!suppliedElsewhere("CASFRIRas")) {
      mod$CASFRItiffFile <- asPath(file.path(dPath, "Landweb_CASFRI_GIDs.tif"))
      mod$CASFRIattrFile <- asPath(file.path(dPath, "Landweb_CASFRI_GIDs_attributes3.csv"))
      mod$CASFRIheaderFile <- asPath(file.path(dPath,"Landweb_CASFRI_GIDs_README.txt"))

      message("  Loading CASFRI layers...")
      sim$CASFRIRas <- Cache(prepInputs,
                             targetFile = asPath("Landweb_CASFRI_GIDs.tif"),
                             archive = asPath("CASFRI for Landweb.zip"),
                             url = extractURL("CASFRIRas"),
                             alsoExtract = c(mod$CASFRItiffFile, mod$CASFRIattrFile, mod$CASFRIheaderFile),
                             destinationPath = dPath,
                             fun = "raster::raster",
                             studyArea = sim$studyArea,
                             rasterToMatch = sim$rasterToMatch,
                             method = "bilinear", ## ignore warning re: ngb (#5)
                             datatype = "INT4U",
                             filename2 = TRUE,
                             overwrite = TRUE,
                             userTags =  c(cacheTags, "CASFRIRas", "stable"))
    }
  }
  return(invisible(sim))
}

