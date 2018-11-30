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
    defineParameter(".crsUsed", "character",
                    paste("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95",
                          "+x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"),
                    NA, NA, "CRS to be used. Defaults to the biomassMap projection"), ## TODO: remove .crsUsed
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
    defineParameter(".useParallel", "numeric", parallel::detectCores(),
                    "Used in reading csv file with fread. Will be passed to data.table::setDTthreads")
  ),
  inputObjects = bind_rows(
    expectsInput("biomassMap", "RasterLayer",
                 desc = "total biomass raster layer in study area, default is Canada national biomass map",
                 sourceURL = "http://tree.pfc.forestry.ca/kNN-StructureBiomass.tar"), ## TODO: use rasterToMatch instead of biomassMap
    expectsInput("CASFRIRas", "RasterStack",
                 desc = "biomass percentage raster layers by species in Canada species map, created by Pickell et al., UBC, resolution 100m x 100m from LandSat and kNN based on CASFRI.",
                 sourceURL = "https://drive.google.com/file/d/1y0ofr2H0c_IEMIpx19xf3_VTBheY0C9h/view?usp=sharing"),
    expectsInput("Pickell", "RasterStack",
                 desc = paste("biomass percentage raster layers by species in Canada species map,",
                              "created by Pickell et al., UBC, resolution 100m x 100m from LandSat and kNN based on CASFRI."),
                 sourceURL = "https://drive.google.com/open?id=1M_L-7ovDpJLyY8dDOxG3xQTyzPx2HSg4"),
    expectsInput("speciesEquivalency", c("data.table"),
                 desc = "table of species equivalencies. See pemisc::sppEquivalencies_CA for further information",
                 sourceURL = ""),
    expectsInput("speciesLayers", "RasterStack",
                 desc = "biomass percentage raster layers by species in Canada species map",
                 sourceURL = "http://tree.pfc.forestry.ca/kNN-Species.tar"),
    expectsInput("sppMerge", c("list"),
                 desc = paste("list of kNN species layers that should be merged.",
                              "If none, create an empty list. Defaults to merging",
                              "of Pinus contorta and P. banksiana into Pinus sp."),
                 sourceURL = ""),
    expectsInput("sppNameVector", c("character"),
                 desc = "vector of species to select", sourceURL = ""),
    expectsInput("studyArea", "SpatialPolygonsDataFrame",
                 desc = paste("multipolygon to use as the study area,",
                              "with attribute LTHFC describing the fire return interval.",
                              "Defaults to a square shapefile in Southwestern Alberta, Canada."),
                 sourceURL = ""), ## TODO: I don't believe fire return interval is needed for this module!
    expectsInput("studyAreaLarge", "SpatialPolygonsDataFrame",
                 desc = paste("multipolygon (larger area than studyArea) to use for parameter estimation,",
                              "with attribute LTHFC describing the fire return interval.",
                              "Defaults to a square shapefile in Southwestern Alberta, Canada."),
                 sourceURL = "") ## TODO: I don't believe fire return interval is needed for this module!
  ),
  outputObjects = bind_rows(
    createsOutput("speciesLayers", "RasterStack",
                  desc = "biomass percentage raster layers by species in Canada species map"),
    createsOutput("sppNameVector",  c("character", "matrix"), ## TODO: no longer a matrix?
                  desc = "vector or matrix of species to select. If matrix, should have two columns of raw and 'end' species names")
  )
))

## event types
#   - type `init` is required for initialiazation

doEvent.BiomassSpeciesData = function(sim, eventTime, eventType) {
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
  ## load Pickell Pickell et al. and CASFRI

  ## TODO: remove this whole session cache file block?
  if (!exists("sessionCacheFile")) {
    sessionCacheFile <<- tempfile()
  }
  isKnownUser <- (!grepl("shiny", Sys.info()["user"]))
  .cacheVal <<- if (isKnownUser) {
    oauthFilePath <- file.path(modulePath(sim), "..", ".httr-oauth")
    options(httr_oauth_cache = oauthFilePath)
    oauthFilePath
  } else if (grepl("VIC-A", Sys.info()["nodename"])) {
    sessionCacheFile
  } else {
    FALSE
  }

  ## download a small file to confirm access to private Google Drive files
  aaa <- testthat::capture_error({
    googledrive::drive_auth(use_oob = TRUE, verbose = TRUE, cache = .cacheVal)
    file_url <- "https://drive.google.com/file/d/1sJoZajgHtsrOTNOE3LL8MtnTASzY0mo7/view?usp=sharing"
    googledrive::drive_download(googledrive::as_id(file_url), path = tempfile(),
                                overwrite = TRUE, verbose = TRUE)
  })

  if (is.null(aaa)) { # means got the file
    dPath <- asPath(dataPath(sim))
    cacheTags <- c(currentModule(sim), "event:init", "stable")

    message("  Loading CASFRI and Pickell et al. layers")
    ## Ceres: i keep having problems with cached prepInputs here. outside Cache removed ## TODO: restore it
    Pickell <- prepInputs(targetFile = asPath("SPP_1990_100m_NAD83_LCC_BYTE_VEG_NO_TIES_FILLED_FINAL.dat"),
                          url = extractURL("Pickell"),
                          archive = asPath("SPP_1990_100m_NAD83_LCC_BYTE_VEG_NO_TIES_FILLED_FINAL.zip"),
                          alsoExtract = asPath("SPP_1990_100m_NAD83_LCC_BYTE_VEG_NO_TIES_FILLED_FINAL.hdr"),
                          destinationPath = dPath,
                          fun = "raster::raster",
                          studyArea = sim$studyArea,
                          rasterToMatch = sim$biomassMap, ## TODO: use sim$rasterToMatch instead of biomassMap
                          method = "bilinear",
                          datatype = "INT2U",
                          filename2 = TRUE,
                          overwrite = TRUE,
                          userTags = c(cacheTags, "function:prepInputs", "Pickell"))

    CASFRItiffFile <- asPath(file.path(dPath, "Landweb_CASFRI_GIDs.tif"))
    CASFRIattrFile <- asPath(file.path(dPath, "Landweb_CASFRI_GIDs_attributes3.csv"))
    CASFRIheaderFile <- asPath(file.path(dPath,"Landweb_CASFRI_GIDs_README.txt"))

    CASFRIRas <- prepInputs(targetFile = asPath("Landweb_CASFRI_GIDs.tif"), ## TODO: restore use of Cache
                            archive = asPath("CASFRI for Landweb.zip"),
                            url = extractURL("CASFRIRas"),
                            alsoExtract = c(CASFRItiffFile, CASFRIattrFile, CASFRIheaderFile),
                            destinationPath = dPath,
                            fun = "raster::raster",
                            studyArea = sim$studyArea,
                            rasterToMatch = sim$biomassMap, ## TODO: use sim$rasterToMatch instead of biomassMap
                            method = "bilinear",
                            datatype = "INT4U",
                            filename2 = TRUE,
                            overwrite = TRUE,
                            userTags =  c(cacheTags, "function:prepInputs", "CASFRIRas"))

    message("Load CASFRI data and headers, and convert to long format, and define species groups")
    if (P(sim)$.useParallel > 1) data.table::setDTthreads(P(sim)$.useParallel)

    loadedCASFRI <- Cache(loadCASFRI, ## TODO: restore use of Cache
                          CASFRIRas = CASFRIRas,
                          attrFile = CASFRIattrFile,
                          headerFile = CASFRIheaderFile,
                          sppNameVector = sim$sppNameVector,
                          speciesEquivalency = sim$speciesEquivalency,
                          sppEndNamesCol = "LandR",
                          sppMerge = sim$sppMerge,
                          userTags = c("function:loadCASFRI", "BigDataTable"))

    message("Make stack of species layers from Pickell's layer")

    ## check if all species found in kNN database are in CASFRI
    uniqueKeepSp <- unique(loadedCASFRI$keepSpecies$spGroup)
    if (!all(names(sim$speciesLayers) %in% uniqueKeepSp))
      warning("some kNN species not in CASFRI layers.")

    ## TODO: weird bug when useCache = "overwrite"
    if (getOption("reproducible.useCache") == "overwrite")
      opt <- options(reproducible.useCache = TRUE) ## TODO: temporary workaround

    PickellSpStack <- Cache(makePickellStack,
                            PickellRaster = Pickell,
                            uniqueKeepSp = uniqueKeepSp,
                            speciesKnn = names(sim$speciesLayers),
                            destinationPath = dPath,
                            userTags = c(cacheTags, "function:makePickellStack", "PickellStack"))

    if (exists("opt", inherits = FALSE)) options(opt) ## TODO: temporary workaround with above

    crs(PickellSpStack) <- crs(sim$biomassMap) # bug in writeRaster

    message('Make stack from CASFRI data and headers')
    CASFRISpStack <- Cache(CASFRItoSpRasts,
                           CASFRIRas = CASFRIRas,
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

    crs(outStack) <- crs(sim$biomassMap) # bug in writeRaster ## TODO: use sim$rasterToMatch instead of biomassMap

    message("Overlay Pickell_CASFRI with open data set stacks")
    specieslayers2 <- Cache(overlayStacks,
                            highQualityStack = outStack,
                            lowQualityStack = sim$speciesLayers,
                            outputFilenameSuffix = "CASFRI_Pickell_kNN",
                            destinationPath = dPath,
                            userTags = c(cacheTags, "function:overlayStacks", "CASFRI_Pickell_kNN"))
    crs(specieslayers2) <- crs(sim$biomassMap) ## TODO: use sim$rasterToMatch instead of biomassMap

    ## replace species layers
    sim$speciesLayers <- specieslayers2
    message("Using overlaid datasets from CASFRI, Pickell and CFS kNN")
  }

  return(invisible(sim))
}

.inputObjects <- function(sim) {
  dPath <- asPath(dataPath(sim), 1)
  cacheTags = c(currentModule(sim), "function:.inputObjects")

  if (!suppliedElsewhere("studyAreaLarge", sim)) {
    message("'studyAreaLarge' was not provided by user. Using a polygon in Southwestern Alberta, Canada.")

    canadaMap <- Cache(getData, 'GADM', country = 'CAN', level = 1, path = asPath(dPath),
                       cacheRepo = getPaths()$cachePath, quick = FALSE)
    smallPolygonCoords = list(
      coords = data.frame(x = c(-115.9022, -114.9815, -114.3677, -113.4470, -113.5084,
                                -114.4291, -115.3498, -116.4547, -117.1298, -117.3140),
                          y = c(50.45516, 50.45516, 50.51654, 50.51654, 51.62139,
                                52.72624, 52.54210, 52.48072, 52.11243, 51.25310))
    )

    sim$studyAreaLarge <- SpatialPolygons(list(Polygons(list(
      Polygon(smallPolygonCoords$coords)), ID = "swAB_polygon")),
      proj4string = crs(canadaMap))
  }

  if (!suppliedElsewhere("studyArea", sim)) {
    message("'studyArea' was not provided by user. Using the same as 'studyAreaLarge'.")
    sim$studyArea <- sim$studyAreaLarge
  }

  ## check projection
  ## TODO: check more than just projection; remove .crsUsed, use rasterToMatch
  if (!identical(as.character(P(sim)$.crsUsed),
                 as.character(crs(sim$studyAreaLarge)))) {
    sim$studyAreaLarge <- spTransform(sim$studyAreaLarge, P(sim)$.crsUsed) #faster without Cache
  }

  if (!identical(as.character(P(sim)$.crsUsed),
                 as.character(crs(sim$studyArea)))) {
    sim$studyArea <- spTransform(sim$studyArea, P(sim)$.crsUsed) #faster without Cache
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
  }

  if (!suppliedElsewhere("biomassMap", sim)) { ## TODO: use rasterToMatch
    biomassMapFilename <- file.path(dPath, "NFI_MODIS250m_kNN_Structure_Biomass_TotalLiveAboveGround_v0.tif")
    sim$biomassMap <- Cache(prepInputs,
                            targetFile = asPath(basename(biomassMapFilename)),
                            archive = asPath(c("kNN-StructureBiomass.tar",
                                               "NFI_MODIS250m_kNN_Structure_Biomass_TotalLiveAboveGround_v0.zip")),
                            url = extractURL("biomassMap", sim),
                            destinationPath = dPath,
                            studyArea = sim$studyArea,
                            useSAcrs = TRUE,
                            method = "bilinear",
                            datatype = "INT2U",
                            filename2 = TRUE,
                            overwrite = TRUE,
                            userTags = c(cacheTags, "biomassMap"))
  }

  if (!suppliedElsewhere("speciesLayers")) { ## TODO: rule should check for sppNameVector
    specieslayersList <- Cache(loadkNNSpeciesLayers,
                               dPath = dPath,
                               rasterToMatch = sim$biomassMap, ## TODO: use rasterToMatch
                               studyArea = sim$studyAreaLarge,
                               sppNameVector = sim$sppNameVector,
                               speciesEquivalency = sim$speciesEquivalency,
                               sppMerge = sim$sppMerge,
                               knnNamesCol = "KNN",
                               sppEndNamesCol = "LandR",
                               thresh = 10, ## TODO: do we use this? not used in ProprietaryData
                               url = extractURL("speciesLayers"),
                               userTags = c(cacheTags, "speciesLayers"))

    sim$speciesLayers <- specieslayersList$speciesLayers

    ## update the list of original species to use;
    ## as kNN is the lowest resolution, if species weren't found they will be excluded
    sim$sppNameVector <- specieslayersList$sppNameVector
  }

  return(invisible(sim))
}
