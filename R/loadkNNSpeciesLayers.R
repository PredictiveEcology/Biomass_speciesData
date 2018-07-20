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

loadkNNSpeciesLayers <- function(dataPath, rasterToMatch, studyArea, 
                                 species = "all", thresh = 1, cachePath, ...) {
  require(magrittr)
  
  ## check if species is a vector/matrix
  if (class(species) == "character") {
    if (species == "all") {
      ## get all species layers from .tar
      species <- untar(tarfile = file.path(dataPath, "kNN-Species.tar"), list = TRUE) %>%
        grep(".zip", ., value = TRUE) %>%
        sub("_v0.zip", "", .) %>%
        sub(".*Species_", "", .)
    }
    
    ## make a matrix of raw and final species names
    species <-  matrix(data = rep(species, 2),
                       nrow = length(species), ncol = 2, byrow = FALSE)
    colnames(species) = c("speciesnamesRaw", "speciesNamesEnd")
    
  } else if(class(species) == "matrix") {
    ## check column names
    if(!setequal(colnames(species), c("speciesnamesRaw", "speciesNamesEnd")))
      stop("names(species) must be c('speciesnamesRaw', 'speciesNamesEnd'), for raw species names and final species names respectively")
  } else stop("species must be a character vector or a two-column matrix")
  
  species1 <- list()
  a11 <- 1
  suffix <- if (basename(cachePath) == "cache") paste0(as.character(ncell(rasterToMatch)),"px") else
    basename(cachePath)
  suffix <- paste0("_", suffix)
  
  for (sp in species[, "speciesnamesRaw"]) {
    targetFile <- paste0("NFI_MODIS250m_kNN_Species_", sp, "_v0.tif")
    postProcessedFilename <- .suffix(targetFile, suffix = suffix)
    # postProcessedFilename <- basename(.suffix(targetFile, suffix = suffix))   ## remove root appended by suffix
    
    species1[[sp]] <- prepInputs(
      targetFile = targetFile,
      url = extractURL(objectName = "specieslayers", module = "BiomassSpeciesData"),
      archive = asPath(c("kNN-Species.tar", paste0("NFI_MODIS250m_kNN_Species_", sp, "_v0.zip"))),
      #alsoExtract = if (sp == speciesnamesRaw[1]) paste0("NFI_MODIS250m_kNN_Species_", speciesnamesRaw[-1], "_v0.tif"),
      destinationPath = asPath(dataPath),
      fun = "raster::raster",
      studyArea = studyArea,
      rasterToMatch = rasterToMatch,
      method = "bilinear",
      datatype = "INT2U",
      filename2 = postProcessedFilename
    )
  }
  
  ## Sum species that share same final name
  if(any(duplicated(species[, 2]))) {
    dubs <- species[duplicated(species[, 2]), 2]   ## get the duplicated final names
    
    ## make a list of species that will be summed (those with duplicated final names)
    spp2sum <- lapply(dubs, FUN = function(x) {
      species[species[, 2] %in% dubs, 1]
    })
    
    names(spp2sum) = dubs     
    
    for(i in 1:length(spp2sum)) {
      sumSpecies <- spp2sum[[i]]
      newLayerName <- names(spp2sum)[i]
      
      fname <- .suffix(file.path(dataPath, "KNNPinu_sp.tif"), suffix)
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
    
    ## Rename species whose raw and final names differ
    nameReplace <- species[!species[, 2] %in% dubs,, drop = FALSE] %>%
      .[which(.[, 1] != .[,2]),, drop = FALSE]
    rownames(nameReplace) = nameReplace[, 1]
    toReplace <- names(species1)[names(species1) %in% nameReplace[,1]]
    names(species1)[names(species1) %in% toReplace] <- nameReplace[toReplace, 2]
  }

  ## remove layers that have less data than thresh (i.e. spp absent in study area)
  ## count no. of pixels that have biomass
  layerData <- Cache(sapply, X = species1, function(x) sum(x[] > 0, na.rm = TRUE))
  
  ## remove layers that had < thresh pixels with biomass
  species1[layerData < thresh] <- NULL

  ## return stack
  stack(species1)
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
