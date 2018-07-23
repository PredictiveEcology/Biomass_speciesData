loadCASFRI <- function(CASFRIRas, attrFile, headerFile) {
  CASFRIattr <- fread(attrFile)

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

  keepSpecies <- whSpecies(CASFRIattr, topN = 16) # 16 most abundant species
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

whSpecies <- function(CASFRIattr, topN = 16) {
  spAbund <- CASFRIattr[, .N, by = "SPECIES_1"] %>% setkeyv("N") #%>% print()
  spAbund2 <- CASFRIattr[, .N, by = "SPECIES_2"] %>% setkeyv("N") #%>% print()
  setorder(spAbund, -N)
  setorder(spAbund2, N)
  keepSpecies <- data.table(keepSpecies = spAbund$SPECIES_1[1:topN])
  set(keepSpecies, , "spGroup", keepSpecies$keepSpecies)
  setkey(keepSpecies, keepSpecies)
  keepSpecies <- keepSpecies[!"Pseu menz"]
  keepSpecies[c("Pice glau", "Pice enge", "Pice hybr", "Pice spp."), spGroup := "Pice_gla"]
  keepSpecies["Pice mari", spGroup := "Pice_mar"]
  keepSpecies["Betu papy", spGroup := "Betu_pap"]
  keepSpecies[c("Abie bals", "Abie lasi"), spGroup := "Abie_sp"]
  keepSpecies[c("Lari lari"), spGroup := "Lari_lar"]
  keepSpecies[c("Pinu cont", "Pinu conl"), spGroup := "Pinu_sp"]
  keepSpecies[c("Pinu bank", "Pinu spp."), spGroup := "Pinu_sp"]
  keepSpecies[c("Popu trem", "Popu balb"), spGroup := "Popu_tre"]
  keepSpecies
}

makePickellStack <- function(PickellRaster, uniqueKeepSp, species, destinationPath) {
  PickellRaster[] <- PickellRaster[]
  PickellRaster[PickellRaster[] %in% c(230, 220, 255)] <- NA_integer_ # water, non veg
  PickellStack <- list()
  
  rasterOptions(maxmemory = 1e9)
  
  ## species in Pickel's data
  PickellSpp <- c("Pice_mar", "Pice_gla", "Pinu_sp", "Popu_tre")
  
  ## selected spp absent from Pickell's data
  NA_Sp <- species[,2][!species[,2] %in% PickellSpp]
  
  ## selected spp present in Pickell's data
  OK_Sp <- species[,2][species[,2] %in% PickellSpp]
    
  ## All NA_Sp species codes should be in CASFRI spp list
  if (!(all(NA_Sp %in% uniqueKeepSp)))
    stop("Codes in loadedCASFRI have changed: expecting ", NA_Sp[!(NA_Sp %in% uniqueKeepSp)])
  
  ## empty rasters for NA_sp
  for (N in lapply(NA_Sp, grep, uniqueKeepSp, value = TRUE)) {
    message("  running ", N, ", assigning NA because not enough data")
    PickellStack[[N]] <- raster(PickellRaster) %>% setValues(x =  ., values = NA_integer_)
    PickellStack[[N]] <- Cache(writeRaster, PickellStack[[N]],
                               filename = asPath(file.path(destinationPath, paste0("Pickell", N, ".tif"))),
                               overwrite = TRUE, datatype = "INT2U")
  }

  ## converting species codes into percentages
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

## TODO filter CASFRI by species

CASFRItoSpRasts <- function(CASFRIRas, loadedCASFRI, species, destinationPath) {
  spRasts <- list()
  spRas <- raster(CASFRIRas) %>% setValues(., NA_integer_)
  
  sppTODO <- setequal(unique(loadedCASFRI$keepSpecies$spGroup) %in% species)
  
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
    rastTmp <- writeRaster(
      spRasts[[sp]],
      filename = asPath(file.path(destinationPath, paste0("CASFRI",sp,".tif"))),
      datatype = "INT1U", overwrite = TRUE
    )
    if (is(rastTmp, "Raster")) { # Rasters need to have their disk-backed value assigned, but not shapefiles
      # This is a bug in writeRaster was spotted with crs of rastTmp became
      # +proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs
      # should have stayed at
      # +proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0
      if (!identical(crs(rastTmp), crs(spRasts[[sp]])))
        crs(rastTmp) <- crs(spRasts[[sp]])

      spRasts[[sp]] <- rastTmp
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
  dt1 <- data.table(SPP = layerNames(highQualityStack), HQ = TRUE)
  dt2 <- data.table(SPP = layerNames(lowQualityStack), LQ = TRUE)
  dtj <- dt1[dt2, on = .(SPP)]
  
  browser()   ## TODO overlaying still needs to be checked
  stack(dtj[, overlay.fun(SPP, HQ, LQ),
            by = 1:nrow(dtj)])
  
}

## ---------------------------------------------------------------------------------
## Overlaying function to use in overlayStacks -------------------------------------
## Function to be applied to each row of a data.table containing information 
##    of whether the species layer exists in the HQ and LQ data. 
## SPP: data.table column of species layer name 
## HQ: data.table column of whether SPP is present in HQ layers
## LQ: : data.table column of whether SPP is present in LQ layers

overlay.fun <- function(SPP, HQ, LQ) {
  ## if HQ & LQ have data, pool
  if (HQ & LQ) {
    ## check equality of raster attributes and correct if necessary
    if (!all(
      isTRUE(all.equal(extent(lowQualityStack), extent(highQualityStack))),
      isTRUE(all.equal(crs(lowQualityStack), crs(highQualityStack))),
      isTRUE(all.equal(res(lowQualityStack), res(highQualityStack))))) {
      message("  ", SPP, " extents, or resolution, or projection did not match; ",
              "using gdalwarp to make them overlap")
      LQRastName <- basename(tempfile(fileext = ".tif"))
      if (!nzchar(filename(lowQualityStack[[SPP]]))) {
        LQCurName <- basename(tempfile(fileext = ".tif"))
        lowQualityStack[[SPP]][] <- as.integer(lowQualityStack[[SPP]][])
        lowQualityStack[[SPP]] <- writeRaster(lowQualityStack[[SPP]], filename = LQCurName,
                                                 datatype = "INT2U")
      }
      
      LQRastInHQcrs <- projectExtent(lowQualityStack, crs = crs(highQualityStack))
      # project LQ raster into HQ dimensions
      gdalwarp(overwrite = TRUE,
               dstalpha = TRUE,
               s_srs = as.character(crs(lowQualityStack[[SPP]])),
               t_srs = as.character(crs(highQualityStack[[SPP]])),
               multi = TRUE, of = "GTiff",
               tr = res(highQualityStack),
               te = c(xmin(LQRastInHQcrs), ymin(LQRastInHQcrs),
                      xmax(LQRastInHQcrs), ymax(LQRastInHQcrs)),
               filename(lowQualityStack[[SPP]]), ot = "Byte",
               LQRastName)
      
      LQRast <- raster(LQRastName)
      LQRast[] <- LQRast[]
      unlink(LQRastName)
      
      try(unlink(LQCurName), silent = TRUE)
      
      if (hqLarger) {
        tmpHQName <- basename(tempfile(fileext = ".tif"))
        
        gdalwarp(overwrite = TRUE,
                 dstalpha = TRUE,
                 s_srs = as.character(crs(highQualityStack[[SPP]])),
                 t_srs = as.character(crs(highQualityStack[[SPP]])),
                 multi = TRUE, of = "GTiff",
                 tr = res(highQualityStack),
                 te = c(xmin(LQRastInHQcrs), ymin(LQRastInHQcrs),
                        xmax(LQRastInHQcrs), ymax(LQRastInHQcrs)),
                 filename(highQualityStack[[SPP]]), ot = "Byte", tmpHQName)
        HQRast <- raster(tmpHQName)
        HQRast[] <- HQRast[]
        HQRast[HQRast[] == 255] <- NA_integer_
        unlink(tmpHQName)
      } else {
        HQRast <- highQualityStack[[SPP]]
      }
    } else {
      LQRast <- lowQualityStack[[SPP]]
      HQRast <- highQualityStack[[SPP]]
    }
    
    message("  Writing new, overlaid ", SPP, " raster to disk.")
    if (!compareRaster(LQRast, HQRast))
      stop("Stacks not identical, something is wrong with overlayStacks function.")
    
    NAs <- is.na(HQRast[])
    
    ## complete missing HQ data with LQ data
    HQRast[NAs] <- LQRast[][NAs]
    HQRast <- writeRaster(HQRast, datatype = "INT1U",
                          filename = file.path(destinationPath,
                                               paste0(SPP, "_", outputFilenameSuffix, ".tif")),
                          overwrite = TRUE)
    names(HQRast) <- SPP
    return(HQRast)
  } else {
    
    ## if only HQ/LQ exist return one of them
    if (HQ) {
      HQRast <- highQualityStack[[SPP]]
      names(HQRast) <- SPP
      return(HQRast)
    }
    
    if(LQ) {
      LQRast <- lowQualityStack[[SPP]]
      names(LQRast) <- SPP
      return(LQRast)
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


