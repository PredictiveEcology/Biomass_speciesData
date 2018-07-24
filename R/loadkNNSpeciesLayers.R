## ------------------------------------------------------------------
## Function to load kNN species layers from online data repository
## ------------------------------------------------------------------

## dataPath: directory to data folder
## rasterToMatch: passed to prepInputs
## studyArea: passed to prepInputs
## species is either a character vector of species names to download,
##    or a two-column matrix with the species names to download and final names, with column names = c("speciesnamesRaw", "speciesNamesEnd") 
##    should two raw species names share the same final name, their biomass data will be considered as the "same species"
## thresh: is the minimum number of pixels where the species must have biomass > 0 to be considered present in the study area. 
##    Defaults to 1
## url: is the source url for the data, passed to prepInputs.

loadkNNSpeciesLayers <- function(dataPath, rasterToMatch, studyArea, 
                                 speciesList = "all", thresh = 1, url, cachePath, ...) {
  require(magrittr)
  
  ## get all kNN species
  allSpp <- Cache(untar, tarfile = file.path(dataPath, "kNN-Species.tar"), list = TRUE) 
  allSpp <- allSpp %>%
    grep(".zip", ., value = TRUE) %>%
    sub("_v0.zip", "", .) %>%
    sub(".*Species_", "", .) 
  
  ## check if species is a vector/matrix
  if (class(speciesList) == "character") {
    if (speciesList == "all") {
      ## get all species layers from .tar
      speciesList <- allSpp
    }
    
    ## check for missing species
    if(any(!speciesList %in% allSpp)) {
      warning("Some species not present in kNN database. 
              /n  Check if this is correct")
      speciesList <- speciesList[speciesList %in% allSpp]
    }
    
    ## make a matrix of raw and final species names
    speciesList <-  matrix(data = rep(speciesList, 2),
                           nrow = length(speciesList), ncol = 2, byrow = FALSE)
    colnames(speciesList) = c("speciesnamesRaw", "speciesNamesEnd")
    
  } else if(class(speciesList) == "matrix") {
    ## check column names
    if(!setequal(colnames(speciesList), c("speciesnamesRaw", "speciesNamesEnd")))
      stop("names(species) must be c('speciesnamesRaw', 'speciesNamesEnd'), for raw species names and final species names respectively")
    
    ## check for missing species
    if(any(!speciesList[,1] %in% allSpp)) {
      warning("Some species not present in kNN database. 
              /n  Check if this is correct")
      speciesList <- speciesList[speciesList[, 1] %in% allSpp,]
    }
  } else stop("species must be a character vector or a two-column matrix")
  
  suffix <- if (basename(cachePath) == "cache") paste0(as.character(ncell(rasterToMatch)),"px") else
    basename(cachePath)
  suffix <- paste0("_", suffix)
  
  loadFun <- function(sp) {
    targetFile <- paste0("NFI_MODIS250m_kNN_Species_", sp, "_v0.tif")
    postProcessedFilename <- .suffix(targetFile, suffix = suffix)
    
    species1 <- prepInputs(
      targetFile = targetFile,
      url = url,
      archive = asPath(c("kNN-Species.tar", paste0("NFI_MODIS250m_kNN_Species_", sp, "_v0.zip"))),
      destinationPath = asPath(dataPath),
      fun = "raster::raster",
      studyArea = studyArea,
      rasterToMatch = rasterToMatch,
      method = "bilinear",
      datatype = "INT2U",
      filename2 = postProcessedFilename)
    
    names(species1) <- sp
    return(species1)
  }
  
  species1 <- Cache(lapply,
                    speciesList[, "speciesnamesRaw"],
                    loadFun,
                    userTags = "kNN_SppLoad")
  
  names(species1) <- speciesList[, "speciesnamesRaw"]
  
  ## Sum species that share same final name
  if(any(duplicated(speciesList[, 2]))) {
    dubs <- speciesList[duplicated(speciesList[, 2]), 2]   ## get the duplicated final names
    
    ## make a list of species that will be summed (those with duplicated final names)
    spp2sum <- lapply(dubs, FUN = function(x) {
      speciesList[speciesList[, 2] %in% x, 1]
    })
    
    names(spp2sum) = dubs     
    
    for(i in 1:length(spp2sum)) {
      sumSpecies <- spp2sum[[i]]
      newLayerName <- names(spp2sum)[i]
      
      fname <- .suffix(file.path(dataPath, paste0("KNN", newLayerName, ".tif")), suffix)
      a <- Cache(sumRastersBySpecies,
                 speciesLayers = species1[sumSpecies], 
                 newLayerName = newLayerName,
                 filenameToSave = asPath(fname),
                 ...)
      a <- raster(fname) ## ensure a gets a filename
      
      ## replace spp rasters by the summed one
      species1[sumSpecies] <- NULL
      species1[[newLayerName]] <- a
    }
  }
  
  ## Rename species layers - note: merged species were renamed already
  nameReplace <- as.matrix(speciesList[,2])
  rownames(nameReplace) = speciesList[, 1]
  
  toReplace <- names(species1)[names(species1) %in% rownames(nameReplace)]
  names(species1)[names(species1) %in% toReplace] <- nameReplace[toReplace, 1]
  
  ## remove layers that have less data than thresh (i.e. spp absent in study area)
  ## count no. of pixels that have biomass
  layerData <- Cache(sapply, X = species1, function(x) sum(x[] > 0, na.rm = TRUE))
  
  ## remove layers that had < thresh pixels with biomass
  species1[layerData < thresh] <- NULL
  
  ## return stack and final species matrix
  list(specieslayers = stack(species1), speciesList = speciesList)
  }


## ------------------------------------------------------------------
## Function to sum rasters of species layers
## ------------------------------------------------------------------

## speciesLayers: stack of species layers rasters
## layersToSum: names/indices of layers to be summed - optional
## filenameToSave: file path to save output raster
## newLayerName: name of the output raster layer

sumRastersBySpecies <- function(speciesLayers, layersToSum,
                                filenameToSave, newLayerName) {
  ras_out <- raster::calc(raster::stack(speciesLayers[layersToSum]), sum) 
  names(ras_out) <- newLayerName 
  writeRaster(ras_out, filename = filenameToSave, datatype = "INT2U", overwrite = TRUE) 
  ras_out # Work around for Cache 
}
