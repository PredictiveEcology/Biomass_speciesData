loadCASFRI <- function(CASFRIRas, attrFile, headerFile, sppNameVector, speciesEquivalency,
                       sppEndNamesCol, sppMerge) {
  CASFRIattr <- Cache(fread, asPath(attrFile))
  if (length(names(CASFRIattr)) < dim(CASFRIattr)[2]) {
    ## Cache brings a table with partial headers when it's called for the second
    ## time after names are changed (data.table issue); do it again, overwriting Cache.
    CASFRIattr <- Cache(fread, asPath(attrFile), useCache = "overwrite")
  }
  ## WORKAROUND: manually extracted the column names because the text file is
  ## very inconsistent (blank lines and uses mix of tabs and spaces as delimiters)
  CASFRIheader <- c("GID", "CAS_ID", "SPECIES_1", "SPECIES_PER_1", "SPECIES_2",
                    "SPECIES_PER_2", "SPECIES_3", "SPECIES_PER_3", "SPECIES_4",
                    "SPECIES_PER_4", "SPECIES_5", "SPECIES_PER_5",
                    "LYR_CROWN_CLOSURE_LOWER", "LYR_CROWN_CLOSURE_UPPER",
                    "LYR_HEIGHT_LOWER", "LYR_HEIGHT_UPPER", "WETLAND_TYPE",
                    "SOIL_MOIST_REG", "AGE", "NAT_NON_VEG", "NON_FOR_ANTH",
                    "NON_FOR_VEG", "HAS_DST", "HAS_NFL", "HAS_LYR")
  #CASFRIheader <- fread(headerFile, skip = 14, nrows = 50, header = FALSE, sep = "\t") ## doesn't work

  setnames(CASFRIattr, CASFRIheader)
  set(CASFRIattr, NULL, grep(CASFRIheader, pattern = "^SPECIES|^GID|^AGE", invert = TRUE), NULL)
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

  ## select species not based on abundance but on user inputs:
  keepSpecies <- whSpecies(CASFRIattr, sppNameVector, speciesEquivalency, sppEndNamesCol, sppMerge)

  CASFRIattrLong <- melt(CASFRIattr, id.vars = c("GID"),
                         measure.vars = paste0("SPECIES_", 1:5))
  CA2 <- melt(CASFRIattr, id.vars = c("GID"),
              measure.vars = c(paste0("SPECIES_PER_", 1:5)))
  CASFRIattrLong[, pct := CA2$value]
  rm(CA2)
  CASFRIattrLong <- na.omit(CASFRIattrLong)

  CASFRIdt <- data.table(GID = CASFRIRas[], rastInd = 1:ncell(CASFRIRas))
  CASFRIdt <- CASFRIdt[, isNA := is.na(GID)]
  CASFRIdt <- CASFRIdt[isNA == FALSE]
  setkey(CASFRIdt, GID)
  set(CASFRIdt, NULL, "isNA", NULL)

  return(list(keepSpecies = keepSpecies,
              CASFRIattrLong = CASFRIattrLong,
              CASFRIdt = CASFRIdt))
}

whSpecies <- function(CASFRIattr, sppNameVector, speciesEquivalency, sppEndNamesCol, sppMerge) {
  keepSpecies <- na.omit(data.table(keepSpecies = unique(CASFRIattr$SPECIES_1)))

  ## convert to LandR format - remove the "." in "spp."
  keepSpecies$keepSpecies2 <- sub(".", "", keepSpecies$keepSpecies, fixed = TRUE) %>%
    equivalentName(., speciesEquivalency, sppEndNamesCol)

  ## species groups according to user-supplied list
  sppMerge2 <- data.table(toMerge = unlist(sppMerge, use.names = FALSE),
                          endName = rep(names(sppMerge), times = sapply(sppMerge, length)))
  sppMerge2$toMerge <- equivalentName(sppMerge2$toMerge, speciesEquivalency, sppEndNamesCol)
  keepSpecies[sppMerge2, on = "keepSpecies2==toMerge"]
  keepSpecies <- sppMerge2[keepSpecies, on = "toMerge==keepSpecies2"]
  setnames(keepSpecies, c("toMerge", "endName"), c("keepSpecies2", "spGroup"))

  ## Because some sister species are usually poorly distinguished in the CASFRI data,
  ## some need to be pooled.
  ## add Picea engelmannii x glauca hybrid if one of the others is in the list
  setkey(keepSpecies, keepSpecies)
  if (any(sppNameVector %in% c("Pice_eng", "Pice_gla")))
    keepSpecies["Pice hybr", spGroup := "Pice_sp"]

  ## add other Populus to Populus sp if there is Populus in the list
  if (any(grep("Popu", sppNameVector)))
    keepSpecies[grep("Popu spp.", "Popu balb", "Popu balt", "Popu hybr", "Popu delt"),
                spGroup := "Popu_tre"]

  ## add other Betula to Betula sp if there is Betula in the list
  if (any(grep("Betu", sppNameVector)))
    keepSpecies[c("Betu papy", "Betu neoa" , "Betu spp."),
                spGroup := "Popu_tre"]

  ## fill empty groups
  keepSpecies[is.na(spGroup), spGroup := keepSpecies2]

  ## Finally, filter species
  ## (note that there might be more species than in the original data)
  keepSpecies <- keepSpecies[keepSpecies2 %in% sppNameVector]

  keepSpecies
}

makePickellStack <- function(PickellRaster, uniqueKeepSp, speciesKnn, destinationPath) {
  ## bring to memory and replace water, non veg by NAs
  PickellRaster[] <- PickellRaster[]
  PickellRaster[PickellRaster[] %in% c(230, 220, 255)] <- NA_integer_

  ## create list and template raster
  spRasts <- list()
  spRas <- raster(PickellRaster) %>% setValues(., NA_integer_)

  rasterOptions(maxmemory = 1e9)

  ## species in Pickel's data
  PickellSpp <- c("Pice_mar", "Pice_gla", "Pinu_sp", "Popu_tre")

  ## selected Knn spp absent from Pickell's data
  NA_Sp <- setdiff(speciesKnn, PickellSpp)

  ## selected Knn spp present in Pickell's data
  sppTODO <- intersect(speciesKnn, PickellSpp)

  ## All NA_Sp species codes should be in CASFRI spp list
  if (length(NA_Sp))
    warning(cat("Not all selected species are in Pickell's data. Check if this is correct:\n",
                paste(paste0(NA_Sp, collapse = ", "), "absent\n")))

  ## empty rasters for NA_sp
  for (sp in NA_Sp) {
    message("  running ", sp, ". Assigning NA, because absent from Pickell's data")
    spRasts[[sp]] <- spRas
    spRasts[[sp]] <- Cache(writeRaster,   spRasts[[sp]],
                          filename = asPath(file.path(destinationPath, paste0("Pickell", sp, ".tif"))),
                          overwrite = TRUE, datatype = "INT2U")
  }

  ## converting existing species codes into percentages
  for (sp in lapply(sppTODO, grep, uniqueKeepSp, value = TRUE)) {
    message("  converting Pickell's codes to pct cover raster, for ", sp)

    if (sp == "Pice_gla") {
      spRasts[[sp]] <- spRas
      spRasts[[sp]][PickellRaster[] %in% c(41, 42, 43)] <- 60
      spRasts[[sp]][PickellRaster[] %in% c(44)] <- 80
      spRasts[[sp]][PickellRaster[] %in% c(14, 34)] <- 40
      spRasts[[sp]] <- Cache(writeRaster,   spRasts[[sp]] ,
                            filename = asPath(file.path(destinationPath, paste0("Pickell", sp, ".tif"))),
                            overwrite = TRUE, datatype = "INT1U")
    }
    if (sp == "Pice_mar") {
      spRasts[[sp]] <- spRas
      spRasts[[sp]][PickellRaster[] %in% c(23, 26)] <- 60
      spRasts[[sp]][PickellRaster[] %in% c(22)] <- 80
      spRasts[[sp]][PickellRaster[] %in% c(32, 42)] <- 40
      spRasts[[sp]] <- Cache(writeRaster,   spRasts[[sp]],
                            filename = asPath(file.path(destinationPath, paste0("Pickell", sp, ".tif"))),
                            overwrite = TRUE, datatype = "INT1U")
    }
    if (sp == "Pinu_sp") {
      spRasts[[sp]] <- spRas
      spRasts[[sp]][PickellRaster[] %in% c(31, 32, 34)] <- 60
      spRasts[[sp]][PickellRaster[] %in% c(33)] <- 80
      spRasts[[sp]][PickellRaster[] %in% c(23, 43)] <- 40
      spRasts[[sp]] <- Cache(writeRaster,   spRasts[[sp]],
                            filename = asPath(file.path(destinationPath, paste0("Pickell", sp, ".tif"))),
                            overwrite = TRUE, datatype = "INT1U")
    }
    if (sp == "Popu_tre") {
      spRasts[[sp]] <- spRas
      spRasts[[sp]][PickellRaster[] %in% c(14)] <- 60
      spRasts[[sp]][PickellRaster[] %in% c(11)] <- 80
      spRasts[[sp]][PickellRaster[] %in% c(31, 41)] <- 40
      spRasts[[sp]] <- Cache(writeRaster,   spRasts[[sp]],
                            filename = asPath(file.path(destinationPath, paste0("Pickell", sp, ".tif"))),
                            overwrite = TRUE, datatype = "INT2U")
    }
  }
  raster::stack(spRasts)
}

## ---------------------------------------------------------------------------------

CASFRItoSpRasts <- function(CASFRIRas, loadedCASFRI, speciesKnn, destinationPath) {
  ## create list and template raster
  spRasts <- list()
  spRas <- raster(CASFRIRas) %>% setValues(., NA_integer_)

  ## selected spp absent from CASFRI data
  NA_Sp <- setdiff(speciesKnn, unique(loadedCASFRI$keepSpecies$spGroup))

  ## All NA_Sp species codes should be in CASFRI spp list
  if (length(NA_Sp))
    warning(cat("Not all selected species are in loadedCASFRI. Check if this is correct:\n",
                paste(paste0(NA_Sp, collapse = ", "), "absent\n")))

  ## empty rasters for NA_sp
  for (sp in NA_Sp) {
    message("  running ", sp, ". Assigning NA, because absent from CASFRI")
    spRasts[[sp]] <- spRas
    spRasts[[sp]] <- Cache(writeRaster, spRasts[[sp]],
                           filename = asPath(file.path(destinationPath, paste0("CASFRI", sp,".tif"))),
                           overwrite = TRUE, datatype = "INT2U")
  }

  sppTODO <- intersect(unique(loadedCASFRI$keepSpecies$spGroup), speciesKnn)

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

  raster::stack(spRasts)
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


