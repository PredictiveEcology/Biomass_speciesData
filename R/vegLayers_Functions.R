loadCASFRI <- function(CASFRIRas, attrFile, headerFile, speciesList) {
  CASFRIattr <- Cache(fread, attrFile)
  
  ## WORKAROUND: maunally extracted the column names because the text file is
  ## very inconsistent (blank lines and uses mix of tabs and spaces as delimiters)
  CASFRIheader <- c("GID", "CAS_ID", "SPECIES_1", "SPECIES_PER_1", "SPECIES_2",
                    "SPECIES_PER_2", "SPECIES_3", "SPECIES_PER_3", "SPECIES_4",
                    "SPECIES_PER_4", "SPECIES_5", "SPECIES_PER_5",
                    "LYR_CROWN_CLOSURE_LOWER", "LYR_CROWN_CLOSURE_UPPER",
                    "LYR_HEIGHT_LOWER", "LYR_HEIGHT_UPPER", "WETLAND_TYPE",
                    "SOIL_MOIST_REG", "AGE", "NAT_NON_VEG", "NON_FOR_ANTH",
                    "NON_FOR_VEG", "HAS_DST", "HAS_NFL", "HAS_LYR")
  #CASFRIheader <- fread(headerFile, skip = 14, nrows = 50, header = FALSE,
  #                      sep = "\t") ## doesn't work
  
  setnames(CASFRIattr, CASFRIheader)
  set(CASFRIattr, , grep(CASFRIheader, pattern = "^SPECIES|^GID|^AGE", invert = TRUE), NULL)
  #setnames(CASFRIattr, CASFRIheader$V1)
  #set(CASFRIattr, , grep(CASFRIheader$V1, pattern = "^SPECIES|^GID|^AGE", invert = TRUE), NULL)
  setkey(CASFRIattr, "GID")
  NAVals <- c("XXXX MISS", "UNDEF", "XXXX ERRC")
  for (i in 1:5) {
    set(CASFRIattr, which(CASFRIattr[[paste0("SPECIES_", i)]] %in% NAVals),
        paste0("SPECIES_", i), NA_character_)
    set(CASFRIattr, which(CASFRIattr[[paste0("SPECIES_PER_", i)]] %in% NAVals),
        paste0("SPECIES_", i), NA_character_)
  }
  for (i in 1:1) {
    CASFRIattr <- CASFRIattr[which(CASFRIattr[[paste0("SPECIES_PER_", i)]] > 15), ]
  }
  for (i in 2:5) {
    set(CASFRIattr, which(CASFRIattr[[paste0("SPECIES_PER_", i)]] <= 15),
        paste0("SPECIES_", i), NA_character_)
  }
  
  keepSpecies <- whSpecies(CASFRIattr, speciesList) # select species not based on abundance but on user inputs
  CASFRIattrLong <- melt(CASFRIattr, id.vars = c("GID"),
                         measure.vars = paste0("SPECIES_", 1:5))
  CA2 <- melt(CASFRIattr, id.vars = c("GID"),
              measure.vars = c(paste0("SPECIES_PER_", 1:5)))
  CASFRIattrLong[, pct := CA2$value]
  rm(CA2)
  CASFRIattrLong <- na.omit(CASFRIattrLong)
  
  CASFRIdt <- CASFRIRas[] %>% data.table(GID = ., rastInd = 1:ncell(CASFRIRas))
  CASFRIdt <- CASFRIdt[, isNA := is.na(GID)]
  CASFRIdt <- CASFRIdt[isNA == FALSE]
  setkey(CASFRIdt, GID)
  set(CASFRIdt, , "isNA", NULL)
  
  return(list(keepSpecies = keepSpecies, CASFRIattrLong = CASFRIattrLong,
              CASFRIdt = CASFRIdt))
}

whSpecies <- function(CASFRIattr, speciesList) {
  keepSpecies <- data.table(keepSpecies = unique(CASFRIattr$SPECIES_1))
  
  ## make compatible names
  keepSpecies$keepSpecies2 <- substring(sub(" ", "_",  keepSpecies$keepSpecies), 1, 8) %>%
    sub("_spp", "_sp", .)
  
  ## species groups according to user-supplied list
  ## make compatible species names first
  kNNnames <- lapply(strsplit(speciesList[,1], "_"), function(x) {
    x[1] <- substring(x[1], 1, 4)
    x[2] <- paste0(toupper(substring(x[2], 1, 1)), substring(x[2], 2, 3))
    x
  })
  kNNnames <- sapply(kNNnames, function(x) paste(x, collapse = "_"))
  
  rownames(speciesList) = tolower(kNNnames)
  matchNames <- tolower(keepSpecies[tolower(keepSpecies2) %in% rownames(speciesList), keepSpecies2])
  keepSpecies[tolower(keepSpecies2) %in% rownames(speciesList), spGroup := speciesList[matchNames,2]]
  
  ## Because some sister species are usually poorly distinguished in the CASFRI data,
  ## some need to be pooled.
  ## add Picea engelmannii x glauca hybrid if one of the others is in the list
  setkey(keepSpecies, keepSpecies)
  if(any(speciesList[,2] %in% c("Pice_eng", "Pice_gla")))
    keepSpecies["Pice hybr", spGroup := "Pice_sp"]
  
  ## add other Populus to Populus sp is there is Populus in the list
  if(any(grep("Popu", speciesList[,2])))
    keepSpecies[c("Popu spp.", "Popu balb", "Popu balt", "Popu hybr", "Popu delt"), spGroup := "Popu_tre"]
  
  ## add other Betula to Betula sp is there is Populus in the list
  if(any(grep("Betu", speciesList[,2])))
    keepSpecies[c("Betu papy", "Betu neoa" , "Betu spp."), spGroup := "Popu_tre"]
  
  ## Finally, filter species
  ## (note that there might be more species than in the original data)
  keepSpecies <- keepSpecies[!is.na(spGroup)]
  
  keepSpecies
}

makePickellStack <- function(PickellRaster, uniqueKeepSp, speciesList, destinationPath) {
  PickellRaster[] <- PickellRaster[]
  PickellRaster[PickellRaster[] %in% c(230, 220, 255)] <- NA_integer_ # water, non veg
  PickellStack <- list()
  
  rasterOptions(maxmemory = 1e9)
  
  ## species in Pickel's data
  PickellSpp <- c("Pice_mar", "Pice_gla", "Pinu_sp", "Popu_tre")
  
  ## selected spp absent from Pickell's data
  NA_Sp <- unique(speciesList[,2][!speciesList[,2] %in% PickellSpp])
  
  ## selected spp present in Pickell's data
  OK_Sp <- unique(speciesList[,2][speciesList[,2] %in% PickellSpp])
  
  ## All NA_Sp species codes should be in CASFRI spp list
  if (length(NA_Sp) > 1)
    warning(cat("Not all species selected are in Pickell's data. Check if this is correct:\n",
                paste0(NA_Sp, collapse = ", ")))
  
  ## empty rasters for NA_sp
  for(N in NA_Sp){  
    message("  running ", N, ", assigning NA because not enough data")
    PickellStack[[N]] <- raster(PickellRaster) %>% setValues(x =  ., values = NA_integer_)
    PickellStack[[N]] <- Cache(writeRaster, PickellStack[[N]],
                               filename = asPath(file.path(destinationPath, paste0("Pickell", N, ".tif"))),
                               overwrite = TRUE, datatype = "INT2U")
  }
  
  ## converting existing species codes into percentages
  for(N in lapply(OK_Sp, grep, uniqueKeepSp, value = TRUE)) {
    message("  converting Pickell's codes to pct cover raster, for ", N)
    
    if(N == "Pice_gla") {
      PickellStack[[N]] <- raster(PickellRaster) %>% setValues(NA_integer_)
      PickellStack[[N]][PickellRaster[] %in% c(41, 42, 43)] <- 60
      PickellStack[[N]][PickellRaster[] %in% c(44)] <- 80
      PickellStack[[N]][PickellRaster[] %in% c(14, 34)] <- 40
      PickellStack[[N]] <- Cache(writeRaster, PickellStack[[N]] ,
                                 filename = asPath(file.path(destinationPath, paste0("Pickell", N, ".tif"))),
                                 overwrite = TRUE, datatype = "INT1U")
    }
    if(N == "Pice_mar") {
      PickellStack[[N]] <- raster(PickellRaster) %>% setValues(NA_integer_)
      PickellStack[[N]][PickellRaster[] %in% c(23, 26)] <- 60
      PickellStack[[N]][PickellRaster[] %in% c(22)] <- 80
      PickellStack[[N]][PickellRaster[] %in% c(32, 42)] <- 40
      PickellStack[[N]] <- Cache(writeRaster, PickellStack[[N]],
                                 filename = asPath(file.path(destinationPath, paste0("Pickell", N, ".tif"))),
                                 overwrite = TRUE, datatype = "INT1U")
    }
    if(N == "Pinu_sp") {
      PickellStack[[N]] <- raster(PickellRaster) %>% setValues(NA_integer_)
      PickellStack[[N]][PickellRaster[] %in% c(31, 32, 34)] <- 60
      PickellStack[[N]][PickellRaster[] %in% c(33)] <- 80
      PickellStack[[N]][PickellRaster[] %in% c(23, 43)] <- 40
      PickellStack[[N]] <- Cache(writeRaster, PickellStack[[N]],
                                 filename = asPath(file.path(destinationPath, paste0("Pickell", N, ".tif"))),
                                 overwrite = TRUE, datatype = "INT1U")
    }
    if(N == "Popu_tre") {
      PickellStack[[N]] <- raster(PickellRaster) %>% setValues(NA_integer_)
      PickellStack[[N]][PickellRaster[] %in% c(14)] <- 60
      PickellStack[[N]][PickellRaster[] %in% c(11)] <- 80
      PickellStack[[N]][PickellRaster[] %in% c(31, 41)] <- 40
      PickellStack[[N]] <- Cache(writeRaster, PickellStack[[N]],
                                 filename = asPath(file.path(destinationPath, paste0("Pickell", N, ".tif"))),
                                 overwrite = TRUE, datatype = "INT2U")
    }
  }
  
  stack(PickellStack)
}

## ---------------------------------------------------------------------------------

CASFRItoSpRasts <- function(CASFRIRas, loadedCASFRI, speciesList, destinationPath) {
  spRasts <- list()
  spRas <- raster(CASFRIRas) %>% setValues(., NA_integer_)
  
  ## selected spp absent from CASFRI data
  NA_Sp <- unique(speciesList[,2][!speciesList[,2] %in% unique(loadedCASFRI$keepSpecies$spGroup)])
  
  ## All NA_Sp species codes should be in CASFRI spp list
  if (length(NA_Sp) > 1)
    warning(cat("Not all species selected are in loadedCASFRI. Check if this is correct:\n",
                paste0(NA_Sp, collapse = ", ")))
  
  ## empty rasters for NA_sp
  for(sp in NA_Sp){  
    message("  running ", sp, ", assigning NA because not enough data")
    spRasts[[sp]] <- spRas
    spRasts[[sp]] <- Cache(writeRaster, spRasts[[sp]],
                           filename = asPath(file.path(destinationPath, paste0("CASFRI", sp,".tif"))),
                           overwrite = TRUE, datatype = "INT2U")
  }
  
  sppTODO <- intersect(unique(loadedCASFRI$keepSpecies$spGroup), speciesList[,2])
  
  for (sp in sppTODO) {
    spRasts[[sp]] <- spRas
    message("starting ", sp)
    aa2 <- loadedCASFRI$CASFRIattrLong[
      value %in% loadedCASFRI$keepSpecies[spGroup == sp, keepSpecies]][
        , min(100L, sum(pct)), by = GID]
    setkey(aa2, GID)
    cc <- aa2[loadedCASFRI$CASFRIdt] %>% na.omit()
    rm(aa2)
    spRasts[[sp]][cc$rastInd] <- cc$V1
    message("  ", sp, " writing to disk")
    
    startCRS <- crs(spRasts[[sp]])
    spRasts[[sp]] <- writeRaster(spRasts[[sp]],
                                 filename = asPath(file.path(destinationPath, paste0("CASFRI", sp,".tif"))),
                                 datatype = "INT1U", overwrite = TRUE)
    
    if (is(spRasts[[sp]], "Raster")) {
      # Rasters need to have their disk-backed value assigned, but not shapefiles
      # This is a bug in writeRaster was spotted with crs of rastTmp became
      # +proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs
      # should have stayed at
      # +proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0
      if (!identical(startCRS, crs(spRasts[[sp]])))
        crs(spRasts[[sp]]) <- startCRS
    }
    message("  ", sp, " done")
  }
  
  stack(spRasts)
}

## ---------------------------------------------------------------------------------


overlayStacks <- function(highQualityStack, lowQualityStack, outputFilenameSuffix = "overlay",
                          destinationPath) {
  ## check if HQ resolution > LQ resolutions
  hqLarger <- ncell(lowQualityStack) * prod(res(lowQualityStack)) <
    ncell(highQualityStack) * prod(res(highQualityStack))
  
  ## make table of species layers in HQ and LQ
  dt1 <- data.table(SPP = layerNames(highQualityStack), HQ = NA)
  dt2 <- data.table(SPP = layerNames(lowQualityStack), LQ = NA)
  setkey(dt1, SPP); setkey(dt2, SPP)
  dtj <- merge(dt1, dt2, all = TRUE)
  
  ## check which layers have species info in HQ and LQ
  dtj[, HQ := any(!is.na(highQualityStack[[SPP]][])), by = 1:nrow(dtj)]
  dtj[, LQ := any(!is.na(lowQualityStack[[SPP]][])), by = 1:nrow(dtj)]
  
  stackRas <- list()
  for(x in 1:nrow(dtj)) { 
    stackRas[[x]] <- dtj[x, overlay.fun(SPP, HQ, LQ, 
                                        HQStack = highQualityStack,
                                        LQStack = lowQualityStack,
                                        fileSuff = outputFilenameSuffix,
                                        destPath = destinationPath)]   ## this is not working unless the function is created in here maybe try piping
  }
  names(stackRas) = dtj$SPP
  
  stack(stackRas)
}

## ---------------------------------------------------------------------------------
## Overlaying function to use in overlayStacks -------------------------------------
## Function to be applied to each row of a data.table containing information 
##    of whether the species layer exists in the HQ and LQ data. 
##    Only overlays if data exists in both layers, otherwise returns the layer with data
## SPP: data.table column of species layer name 
## HQ: data.table column of whether SPP is present in HQ layers
## LQ: data.table column of whether SPP is present in LQ layers
## HQStack: high quality list/stack of rasters (will be used preferencially)
## LQStack: high quality list/stack of rasters (will be used to fill NAs in HQStack)
## fileSuff: file suffix to save raster if there was overlaying
## destPath: directory for saved rasters

overlay.fun <- function(SPP, HQ, LQ, HQStack, LQStack,
                        fileSuff, destPath) {
  
  ## if HQ & LQ have data, pool
  if (HQ & LQ) {
    ## check equality of raster attributes and correct if necessary
    if (!all(
      isTRUE(all.equal(extent(LQStack), extent(HQStack))),
      isTRUE(all.equal(crs(LQStack), crs(HQStack))),
      isTRUE(all.equal(res(LQStack), res(HQStack))))) {
      message("  ", SPP, " extents, or resolution, or projection did not match; ",
              "using gdalwarp to make them overlap")
      LQRastName <- basename(tempfile(fileext = ".tif"))
      if (!nzchar(filename(LQStack[[SPP]]))) {
        LQCurName <- basename(tempfile(fileext = ".tif"))
        LQStack[[SPP]][] <- as.integer(LQStack[[SPP]][])
        LQStack[[SPP]] <- writeRaster(LQStack[[SPP]], filename = LQCurName,
                                      datatype = "INT2U")
      }
      
      LQRastInHQcrs <- projectExtent(LQStack, crs = crs(HQStack))
      # project LQ raster into HQ dimensions
      gdalwarp(overwrite = TRUE,
               dstalpha = TRUE,
               s_srs = as.character(crs(LQStack[[SPP]])),
               t_srs = as.character(crs(HQStack[[SPP]])),
               multi = TRUE, of = "GTiff",
               tr = res(HQStack),
               te = c(xmin(LQRastInHQcrs), ymin(LQRastInHQcrs),
                      xmax(LQRastInHQcrs), ymax(LQRastInHQcrs)),
               filename(LQStack[[SPP]]), ot = "Byte",
               LQRastName)
      
      LQRast <- raster(LQRastName)
      LQRast[] <- LQRast[]
      unlink(LQRastName)
      
      try(unlink(LQCurName), silent = TRUE)
      
      if (hqLarger) {
        tmpHQName <- basename(tempfile(fileext = ".tif"))
        
        gdalwarp(overwrite = TRUE,
                 dstalpha = TRUE,
                 s_srs = as.character(crs(HQStack[[SPP]])),
                 t_srs = as.character(crs(HQStack[[SPP]])),
                 multi = TRUE, of = "GTiff",
                 tr = res(HQStack),
                 te = c(xmin(LQRastInHQcrs), ymin(LQRastInHQcrs),
                        xmax(LQRastInHQcrs), ymax(LQRastInHQcrs)),
                 filename(HQStack[[SPP]]), ot = "Byte", tmpHQName)
        HQRast <- raster(tmpHQName)
        HQRast[] <- HQRast[]
        HQRast[HQRast[] == 255] <- NA_integer_
        unlink(tmpHQName)
      } else {
        HQRast <- HQStack[[SPP]]
      }
    } else {
      LQRast <- LQStack[[SPP]]
      HQRast <- HQStack[[SPP]]
    }
    
    message("  Writing new, overlaid ", SPP, " raster to disk.")
    if (!compareRaster(LQRast, HQRast))
      stop("Stacks not identical, something is wrong with overlayStacks function.")
    
    NAs <- is.na(HQRast[])
    
    ## complete missing HQ data with LQ data
    HQRast[NAs] <- LQRast[][NAs]
    HQRast <- writeRaster(HQRast, datatype = "INT1U",
                          filename = file.path(destPath,
                                               paste0(SPP, "_", fileSuff, ".tif")),
                          overwrite = TRUE)
    names(HQRast) <- SPP
    return(HQRast)
  } else {
    
    ## if only HQ/LQ exist return one of them
    ## if none have data return one of the empty to keep all layers
    if (HQ) {
      HQRast <- HQStack[[SPP]]
      names(HQRast) <- SPP
      return(HQRast)
    } else if(LQ) {
      LQRast <- LQStack[[SPP]]
      names(LQRast) <- SPP
      return(LQRast)
    } else {
      HQRast <- HQStack[[SPP]]
      names(HQRast) <- SPP
      return(HQRast)
    }
  }
}

## ---------------------------------------------------------------------------------

gdalwarp2 <- function(rasterWithDiskBacked, dstfilename, ...) {
  dstfilenameTmp <- .suffix(dstfilename, "_tmp")
  gdalwarp(srcfile = basename(filename(rasterWithDiskBacked)),
           dstfile = basename(dstfilenameTmp), ...)
  
  rr <- raster(dstfilenameTmp)
  rr[] <- rr[]
  if (is.integer(rr[])) {
    dt <- if (maxValue(rr) > 65534) {
      "INT4U"
    } else {
      "INT2U"
    }
  }
  rr[rr[] == 255] <- NA_integer_
  rr <- writeRaster(rr, filename = dstfilename,
                    datatype = dt,  overwrite = TRUE)
  unlink(dstfilenameTmp)
  rr
}


