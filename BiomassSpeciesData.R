# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects, use sim$xxx, and are thus globally available
# to all modules. Functions can be used without sim$ as they are namespaced, like functions
# in R packages. If exact location is required, functions will be: sim$<moduleName>$FunctionName
defineModule(sim, list(
  name = "BiomassSpeciesData",
  description = "Download and pre-process proprietary LandWeb data.",
  keywords = c("LandWeb"),
  authors = c(person(c("Eliot", "J", "B"), "McIntire", email = "eliot.mcintire@canada.ca", role = c("aut", "cre")),
              person(c("Alex", "M."), "Chubaty", email = "alexander.chubaty@canada.ca", role = c("aut")),
              person("Ceres", "Barros", email = "cbarros@mail.ubc.ca", role = c("aut"))),
  childModules = character(0),
  version = list(SpaDES.core = "0.1.2", BiomassSpeciesData = "0.0.1"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "BiomassSpeciesData.Rmd"),
  reqdPkgs = list("googledrive", "data.table", "raster", "magrittr",
                  "PredictiveEcology/SpaDES.core@development",
                  "PredictiveEcology/SpaDES.tools@development",
                  "PredictiveEcology/reproducible@development",
                  "ygc2l/webDatabases"),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description"),
    defineParameter("crsUsed", "character", "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0",
                    NA, NA, "CRS to be used. Defaults to the biomassMap projection"),
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", NA, NA, NA, "This describes the simulation time interval between plot events"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first save event should occur"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA, "This describes the simulation time interval between save events"),
    defineParameter(".useCache", "logical", TRUE, NA, NA, "Should this entire module be run with caching activated? This is generally intended for data-type modules, where stochasticity and time are not relevant"),
    defineParameter(name = "useParallel", class = "numeric", default = parallel::detectCores(),
                    desc = "Used in reading csv file with fread. Will be passed to data.table::setDTthreads")
  ),
  inputObjects = bind_rows(
    expectsInput(objectName = "speciesList", objectClass = c("character", "matrix"),
                 desc = "vector or matrix of species to select. If matrix, should have two columns of raw and 'end' species names", sourceURL = ""),
    expectsInput(objectName = "biomassMap", objectClass = "RasterLayer",
                 desc = "total biomass raster layer in study area, default is Canada national biomass map",
                 sourceURL = "http://tree.pfc.forestry.ca/kNN-StructureBiomass.tar"),
    expectsInput(objectName = "specieslayers", objectClass = "RasterStack",
                 desc = "biomass percentage raster layers by species in Canada species map",
                 sourceURL = "http://tree.pfc.forestry.ca/kNN-Species.tar"),
    expectsInput(objectName = "shpStudySubRegion", objectClass = "SpatialPolygonsDataFrame",
                 desc = "this shape file contains two informaton: Sub study area with fire return interval attribute. 
                 Defaults to a square shapefile in Southwestern Alberta, Canada", sourceURL = ""),
    expectsInput(objectName = "shpStudyRegionFull", objectClass = "SpatialPolygonsDataFrame",
                 desc = "this shape file contains two informaton: Full study area with fire return interval attribute.
                 Defaults to a square shapefile in Southwestern Alberta, Canada", sourceURL = ""), # i guess this is study area and fire return interval
    expectsInput(objectName = "Pickell", objectClass = "RasterStack",
                 desc = "biomass percentage raster layers by species in Canada species map, created by Pickell et al., UBC, resolution 100m x 100m from LandSat and kNN based on CASFRI.",
                 sourceURL = "https://drive.google.com/open?id=1M_L-7ovDpJLyY8dDOxG3xQTyzPx2HSg4"),
    expectsInput(objectName = "CASFRIRas", objectClass = "RasterStack",
                 desc = "biomass percentage raster layers by species in Canada species map, created by Pickell et al., UBC, resolution 100m x 100m from LandSat and kNN based on CASFRI.",
                 sourceURL = "https://drive.google.com/file/d/1y0ofr2H0c_IEMIpx19xf3_VTBheY0C9h/view?usp=sharing")
    ),
  outputObjects = bind_rows(
    createsOutput(objectName = "specieslayers", objectClass = "RasterStack",
                  desc = "biomass percentage raster layers by species in Canada species map"),
    createsOutput(objectName = "speciesList", objectClass =  c("character", "matrix"),
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
  
  aaa <- testthat::capture_error({
    googledrive::drive_auth(use_oob = TRUE, verbose = TRUE, cache = .cacheVal)
    file_url <- "https://drive.google.com/file/d/1sJoZajgHtsrOTNOE3LL8MtnTASzY0mo7/view?usp=sharing"
    googledrive::drive_download(googledrive::as_id(file_url), path = tempfile(),
                                overwrite = TRUE, verbose = TRUE)
  })
  
  if (is.null(aaa)) { # means got the file
    dPath <- dataPath(sim)
    cacheTags <- c(currentModule(sim), "event:init", "function:prepInputs")
    
    message("  Loading CASFRI and Pickell et al. layers")
    Pickell <- Cache(prepInputs,
                     targetFile = asPath("SPP_1990_100m_NAD83_LCC_BYTE_VEG_NO_TIES_FILLED_FINAL.dat"),
                     archive = asPath("SPP_1990_100m_NAD83_LCC_BYTE_VEG_NO_TIES_FILLED_FINAL.zip"),
                     url = extractURL(objectName = "Pickell"),
                     alsoExtract = asPath("SPP_1990_100m_NAD83_LCC_BYTE_VEG_NO_TIES_FILLED_FINAL.hdr"),
                     destinationPath = asPath(dPath),
                     fun = "raster::raster",
                     studyArea = sim$shpStudySubRegion,
                     rasterToMatch = sim$biomassMap,
                     method = "bilinear",
                     datatype = "INT2U",
                     filename2 = TRUE,
                     userTags = c(cacheTags, "Pickell"))#, notOlderThan = Sys.time())
    
    CASFRITifFile <- asPath(file.path(dPath, "Landweb_CASFRI_GIDs.tif"))
    CASFRIattrFile <- asPath(file.path(dPath, "Landweb_CASFRI_GIDs_attributes3.csv"))
    CASFRIheaderFile <- asPath(file.path(dPath,"Landweb_CASFRI_GIDs_README.txt"))
    
    CASFRIRas <- Cache(prepInputs,
                       targetFile = asPath("Landweb_CASFRI_GIDs.tif"),
                       archive = asPath("CASFRI for Landweb.zip"),
                       url = extractURL(objectName = "CASFRIRas"),
                       alsoExtract = c(CASFRITifFile, CASFRIattrFile, CASFRIheaderFile),
                       destinationPath = asPath(dPath),
                       fun = "raster::raster",
                       studyArea = sim$shpStudySubRegion,
                       rasterToMatch = sim$biomassMap,
                       method = "bilinear",
                       datatype = "INT4U",
                       filename2 = TRUE,
                       userTags =  c(cacheTags, "CASFRIRas"))
    
    message("Load CASFRI data and headers, and convert to long format, and define species groups")
    if (P(sim)$useParallel > 1) data.table::setDTthreads(P(sim)$useParallel)
    loadedCASFRI <- Cache(loadCASFRI, CASFRIRas, CASFRIattrFile, CASFRIheaderFile,
                          speciesList,
                          # destinationPath = asPath(dPath),
                          # debugCache = "complete",
                          userTags = c("stable", "BigDataTable"))
    
    message("Make stack of species layers from Pickell's layer")
    
    ## check if all kNN species are in CASFRI and if there are case issues
    uniqueKeepSp <- unique(loadedCASFRI$keepSpecies$spGroup)
    if(all(tolower(sim$speciesList[,2]) %in% tolower(uniqueKeepSp))) {
      if(!all(sim$speciesList[,2] %in% uniqueKeepSp))
        stop("end species names in kNN do not match CASFRI species codes.
             \nUpper/lower case issues")
    } else warning("some kNN species not in CASFRI")
    
    PickellSpStack <- Cache(makePickellStack, #paths = lapply(paths(sim), basename),   # paths was throwing an error in cache 
                            PickellRaster = Pickell, uniqueKeepSp, speciesList = sim$speciesList,
                            destinationPath = dPath, userTags = c("stable", "PickellStack"))
    
    crs(PickellSpStack) <- crs(sim$biomassMap) # bug in writeRaster
    
    message('Make stack from CASFRI data and headers')
    CASFRISpStack <- Cache(CASFRItoSpRasts, CASFRIRas, loadedCASFRI, speciesList = sim$speciesList,
                           destinationPath = dPath, userTags = c("stable", "CASFRIstk"))
    
    
    message("Overlay Pickell and CASFRI stacks")
    outStack <- Cache(overlayStacks, 
                      highQualityStack = CASFRISpStack, lowQualityStack = PickellSpStack, 
                      outputFilenameSuffix = "CASFRI_Pickell", destinationPath = dPath,
                      userTags = c("stable", "Pickell_CASFRI")) 
    
    crs(outStack) <- crs(sim$biomassMap) # bug in writeRaster
    
    message("Overlay Pickell_CASFRI with open data set stacks")
    specieslayers2 <- Cache(overlayStacks, 
                            highQualityStack = outStack, lowQualityStack = sim$specieslayers,
                            outputFilenameSuffix = "CASFRI_Pickell_KNN",
                            destinationPath = dPath, 
                            userTags = c("stable", "CASFRI_Pickell_KNN"))
    crs(specieslayers2) <- crs(sim$biomassMap)
    
    ## replace species layers
    sim$specieslayers <- specieslayers2
    message("Using overlaid datasets from CASFRI, Pickell and CFS kNN")
    }
  
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  dPath <- dataPath(sim)
  cacheTags = c(currentModule(sim), "function:.inputObjects")
  
  if (!suppliedElsewhere("shpStudyRegionFull", sim)) {
    message("'shpStudyRegionFull' was not provided by user. Using a polygon in Southwestern Alberta, Canada")
    
    canadaMap <- Cache(getData, 'GADM', country = 'CAN', level = 1, path = asPath(dPath),
                       cacheRepo = getPaths()$cachePath, quick = FALSE) 
    smallPolygonCoords = list(coords = data.frame(x = c(-115.9022,-114.9815,-114.3677,-113.4470,-113.5084,-114.4291,-115.3498,-116.4547,-117.1298,-117.3140), 
                                                  y = c(50.45516,50.45516,50.51654,50.51654,51.62139,52.72624,52.54210,52.48072,52.11243,51.25310)))
    
    sim$shpStudyRegionFull <- SpatialPolygons(list(Polygons(list(Polygon(smallPolygonCoords$coords)), ID = "swAB_polygon")),
                                              proj4string = crs(canadaMap))
  }
  
  if (!suppliedElsewhere("shpStudySubRegion", sim)) {
    message("'shpStudySubRegion' was not provided by user. Using the same as 'shpStudyRegionFull'")
    sim$shpStudySubRegion <- sim$shpStudyRegionFull
  }
  
  ## check projection
  if (!identical(P(sim)$crsUsed, crs(sim$shpStudyRegionFull))) {
    sim$shpStudyRegionFull <- spTransform(sim$shpStudyRegionFull, P(sim)$crsUsed) #faster without Cache
  }
  
  if (!identical(P(sim)$crsUsed, crs(sim$shpStudySubRegion))) {
    sim$shpStudySubRegion <- spTransform(sim$shpStudySubRegion, P(sim)$crsUsed) #faster without Cache
  }
  
  if (!suppliedElsewhere("speciesList", sim)) {
    ## default to 6 species, one changing name, and two merged into one
    sim$speciesList <- as.matrix(data.frame(speciesNamesRaw = c("Abie_Las", "Pice_Gla", "Pice_Mar", "Pinu_Ban", "Pinu_Con", "Popu_Tre"),
                                            speciesNamesEnd =  c("Abie_sp", "Pice_gla", "Pice_mar", "Pinu_sp", "Pinu_sp", "Popu_tre")))
  }
  
  if (!suppliedElsewhere("biomassMap", sim)) {
    biomassMapFilename <- file.path(dPath, "NFI_MODIS250m_kNN_Structure_Biomass_TotalLiveAboveGround_v0.tif")
    sim$biomassMap <- Cache(prepInputs,
                            targetFile = biomassMapFilename,
                            url = extractURL(objectName = "biomassMap"),
                            archive = asPath(c("kNN-StructureBiomass.tar",
                                               "NFI_MODIS250m_kNN_Structure_Biomass_TotalLiveAboveGround_v0.zip")),
                            destinationPath = asPath(dPath),
                            studyArea = sim$shpStudySubRegion,
                            useSAcrs = TRUE,
                            method = "bilinear",
                            datatype = "INT2U",
                            filename2 = TRUE,
                            userTags = c(cacheTags, "biomassMap"))
  }
  
  if (!suppliedElsewhere("specieslayers")) {
    specieslayersList <- Cache(loadkNNSpeciesLayers,
                               dataPath = asPath(dPath), 
                               rasterToMatch = sim$biomassMap, 
                               studyArea = sim$shpStudyRegionFull,
                               speciesList = sim$speciesList,
                               thresh = 10,
                               url = extractURL("specieslayers"),
                               cachePath = cachePath(sim),
                               userTags = c(cacheTags, "specieslayers"))
    
    sim$specieslayers <- specieslayersList$specieslayers
    sim$speciesList <- specieslayersList$speciesList
    
  }
  
  return(invisible(sim))
} 
