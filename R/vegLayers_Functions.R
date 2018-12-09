loadCASFRI <- function(CASFRIRas, attrFile, headerFile, sppNameVector, sppEquiv,
                       sppEquivCol, sppMerge) {

  # The ones we want
  sppEquiv <- sppEquiv[!is.na(sppEquiv[[sppEquivCol]]),]
  if (missing(sppNameVector))
    sppNameVector <- unique(sppEquiv[[sppEquivCol]])
  sppNameVectorCASFRI <- equivalentName(sppNameVector, sppEquiv,  column = "CASFRI", multi = TRUE)


  # CASFRI stuff
  CASFRIheader <- fread(headerFile, skip = 14, nrows = 49, header = FALSE, sep = "", fill = TRUE)
  header <- apply(CASFRIheader, 1, function(x) sub(pattern = "(\t+| ).*$", "", x))
  CASFRIheader <- header[nchar(header) != 0]

  wantedColumns <- grep(CASFRIheader, pattern = "^SPECIES|^GID|^AGE")

  CASFRIattr <- fread(asPath(attrFile), select = wantedColumns)

  setnames(CASFRIattr, CASFRIheader[wantedColumns])
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

  #keepSpecies <- whSpecies(CASFRIattr, sppNameVector, sppEquiv, sppEquivCol, sppMerge)

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

  CASFRIattrLong <- CASFRIattrLong[value %in% sppNameVectorCASFRI]
  return(list(#keepSpecies = keepSpecies,
              CASFRIattrLong = CASFRIattrLong,
              CASFRIdt = CASFRIdt))
}

# whSpecies <- function(CASFRIattr, sppNameVector, sppEquiv, sppEquivCol, sppMerges) {
#   browser()
#   keepSpecies <- na.omit(data.table(CASFRI = unique(CASFRIattr$SPECIES_1)))
#
#   ## convert to LandR format - remove the "." in "spp."
#   keepSpecies[, (sppEquivCol) := list(equivalentName(
#     sub(".", "", CASFRI, fixed = TRUE), sppEquiv, sppEquivCol))]
#
#   ## species groups according to user-supplied list
#   sppMerge2 <- data.table(toMerge = unlist(sppMerge, use.names = FALSE),
#                           endName = rep(names(sppMerge), times = sapply(sppMerge, length)))
#   sppMerge2$toMerge <- equivalentName(sppMerge2$toMerge, sppEquiv, sppEquivCol)
#   #keepSpecies[sppMerge2, on = paste0(sppEquivCol,"==toMerge")]
#   keepSpecies <- sppMerge2[keepSpecies, on = paste0("toMerge==", sppEquivCol)]
#   setnames(keepSpecies, c("toMerge", "endName"), c(sppEquivCol, "spGroup"))
#
#   ## Because some sister species are usually poorly distinguished in the CASFRI data,
#   ## some need to be pooled.
#   ## add Picea engelmannii x glauca hybrid if one of the others is in the list
#   setkey(keepSpecies, CASFRI)
#   if (any(sppNameVector %in% c("Pice_eng", "Pice_gla")))
#     keepSpecies["Pice hybr", spGroup := equivalentName("Popu_sp", sppEquiv, sppEquivCol)]
#
#   ## add other Populus to Populus sp if there is Populus in the list
#   if (any(grep("Popu", sppNameVector)))
#     keepSpecies[CASFRI %in% c("Popu spp.", "Popu balb", "Popu balt", "Popu hybr", "Popu delt", "Popu trem"),
#                 spGroup := equivalentName("Popu_tre", sppEquiv, sppEquivCol)]
#
#   ## add other Betula to Betula sp if there is Betula in the list
#   if (any(grep("Betu", sppNameVector)))
#     keepSpecies[c("Betu papy", "Betu neoa" , "Betu spp."),
#                 spGroup := equivalentName("Popu_tre", sppEquiv, sppEquivCol)]
#
#   ## fill empty groups
#   keepSpecies[is.na(spGroup), spGroup := get(sppEquivCol)]
#
#   ## Finally, filter species
#   ## (note that there might be more species than in the original data)
#   keepSpecies <- keepSpecies[get(sppEquivCol) %in%
#                                equivalentName(sppNameVector, sppEquiv, sppEquivCol)]
#
#   keepSpecies
# }
#
makePickellStack <- function(PickellRaster, sppNameVector,
                             sppEquiv, sppEquivCol,
                             sppMerge,
                             destinationPath) {
  ## bring to memory and replace water, non veg by NAs
  sppEquiv <- speciesEquivalency[!is.na(speciesEquivalency[[sppEndNamesCol]]),]

  # Take this from the speciesEquivalency table; user cannot supply manually
  sppNameVector <- unique(sppEquiv[[sppEndNamesCol]])
  names(sppNameVector) <- sppNameVector

  PickellSpp <- c("Pice_mar", "Pice_gla", "Pinu_sp", "Popu_tre")
  names(PickellSpp) <- PickellSpp

  # Pick the full LandR dataset, which should be broad. We will change to sppEndNamesCol
  #   below
  sppOfInterest <- equivalentName(sppNameVector, sppEquiv, "LandR", multi = TRUE)
  sppInPickell <- lapply(PickellSpp, function(sp)
    equivalentName(sp, sppEquiv, "LandR", multi = TRUE)
  )
  # Check that each of the layers that Pickell did are actually desired in speciesEquivalency
  needPickel <- sapply(sppInPickell, function(sp) {
    any(sp %in% sppOfInterest)
  })
  # These are the ones in Pickell data set that we want according to speciesEquivalency
  PickellSpp <- equivalentName(PickellSpp[needPickel], sppEquiv, sppEndNamesCol)

  PickellRaster[] <- PickellRaster[]
  PickellRaster[PickellRaster[] %in% c(230, 220, 255)] <- NA_integer_

  ## create list and template raster
  spRasts <- list()
  spRas <- raster(PickellRaster) %>% setValues(., NA_integer_)

  rasterOptions(maxmemory = 1e9)

  # keepSpecies <- na.omit(data.table(unique(PickellSpp)))
  # names(keepSpecies) <- sppEndNamesCol
  # ## species groups according to user-supplied list
  # sppMerge2 <- data.table(toMerge = unlist(sppMerge, use.names = FALSE),
  #                         endName = rep(names(sppMerge), times = sapply(sppMerge, length)))
  # sppMerge2$toMerge <- equivalentName(sppMerge2$toMerge, speciesEquivalency, sppEndNamesCol)
  # #keepSpecies[sppMerge2, on = paste0(sppEndNamesCol,"==toMerge")]
  # keepSpecies <- sppMerge2[keepSpecies, on = paste0("endName==", sppEndNamesCol)]
  #
  # sppNameVectorMerged <- unique(keepSpecies$endName)
  #
  # ## selected Knn spp absent from Pickell's data
  # NA_Sp <- setdiff(sppNameVectorMerged, PickellSpp)
  #
  # ## selected Knn spp present in Pickell's data
  # sppTODO <- intersect(sppNameVectorMerged, PickellSpp)

  ## All NA_Sp species codes should be in CASFRI spp list
  # if (length(NA_Sp))
  #   warning(cat("Not all selected species are in Pickell's data. Check if this is correct:\n",
  #               paste(paste0(NA_Sp, collapse = ", "), "absent\n")))
  #
  # ## empty rasters for NA_sp
  # for (sp in NA_Sp) {
  #   message("  running ", sp, ". Assigning NA, because absent from Pickell's data")
  #   spRasts[[sp]] <- spRas
  #   spRasts[[sp]] <- Cache(writeRaster,   spRasts[[sp]],
  #                         filename = asPath(file.path(destinationPath, paste0("Pickell", sp, ".tif"))),
  #                         overwrite = TRUE, datatype = "INT2U")
  # }

  ## converting existing species codes into percentages
  for (sp in PickellSpp) {
    message("  converting Pickell's codes to pct cover raster, for ", sp)

    if (sp == equivalentName("Pice_gla", sppEquiv, sppEquivCol)) {
      spRasts[[sp]] <- spRas
      spRasts[[sp]][PickellRaster[] %in% c(41, 42, 43)] <- 60
      spRasts[[sp]][PickellRaster[] %in% c(44)] <- 80
      spRasts[[sp]][PickellRaster[] %in% c(14, 34)] <- 40
      spRasts[[sp]] <- Cache(writeRaster,   spRasts[[sp]] ,
                            filename = asPath(file.path(destinationPath, paste0("Pickell", sp, ".tif"))),
                            overwrite = TRUE, datatype = "INT1U")
    }
    if (sp == equivalentName("Pice_mar", sppEquiv, sppEquivCol)) {
      spRasts[[sp]] <- spRas
      spRasts[[sp]][PickellRaster[] %in% c(23, 26)] <- 60
      spRasts[[sp]][PickellRaster[] %in% c(22)] <- 80
      spRasts[[sp]][PickellRaster[] %in% c(32, 42)] <- 40
      spRasts[[sp]] <- Cache(writeRaster,   spRasts[[sp]],
                            filename = asPath(file.path(destinationPath, paste0("Pickell", sp, ".tif"))),
                            overwrite = TRUE, datatype = "INT1U")
    }
    if (sp == equivalentName("Pinu_sp", sppEquiv, sppEquivCol)) {
      spRasts[[sp]] <- spRas
      spRasts[[sp]][PickellRaster[] %in% c(31, 32, 34)] <- 60
      spRasts[[sp]][PickellRaster[] %in% c(33)] <- 80
      spRasts[[sp]][PickellRaster[] %in% c(23, 43)] <- 40
      spRasts[[sp]] <- Cache(writeRaster,   spRasts[[sp]],
                            filename = asPath(file.path(destinationPath, paste0("Pickell", sp, ".tif"))),
                            overwrite = TRUE, datatype = "INT1U")
    }
    if (sp == equivalentName("Popu_tre", sppEquiv, sppEquivCol)) {
      spRasts[[sp]] <- spRas
      spRasts[[sp]][PickellRaster[] %in% c(14)] <- 60
      spRasts[[sp]][PickellRaster[] %in% c(11)] <- 80
      spRasts[[sp]][PickellRaster[] %in% c(31, 41)] <- 40
      spRasts[[sp]] <- Cache(writeRaster,   spRasts[[sp]],
                            filename = asPath(file.path(destinationPath, paste0("Pickell", sp, ".tif"))),
                            overwrite = TRUE, datatype = "INT2U")
    }
  }
  ## species in Pickel's data
  raster::stack(spRasts)
}

## ---------------------------------------------------------------------------------

CASFRItoSpRasts <- function(CASFRIRas, CASFRIattrLong, CASFRIdt,
                            sppEquiv, sppEquivCol, destinationPath) {
  # The ones we want
  sppEquiv <- sppEquiv[!is.na(sppEquiv[[sppEquivCol]]),]

  # Take this from the sppEquiv table; user cannot supply manually
  sppNameVector <- unique(sppEquiv[[sppEquivCol]])
  names(sppNameVector) <- sppNameVector

  # This
  sppListMergesCASFRI <-lapply(sppNameVector, function(x)
    equivalentName(x, sppEquiv,  column = "CASFRI", multi = TRUE)
  )

  ## create list and template raster
  spRasts <- list()
  spRas <- raster(CASFRIRas) %>% setValues(., NA_integer_)

  ## NOT SURE IF THESE LINES ABOUT NA are relevant -- ELiot Dec 7
  ## selected spp absent from CASFRI data
  NA_Sp <- which(is.na(sppListMergesCASFRI))#setdiff(speciesLandR, unique(keepSpecies$spGroup))

  ## All NA_Sp species codes should be in CASFRI spp list
  if (length(NA_Sp))
    warning(cat("Not all selected species are in loadedCASFRI. Check if this is correct:\n",
                paste(paste0(keepSpecies$CASFRI[NA_Sp], collapse = ", "), "absent\n")))

  ## empty rasters for NA_sp
  for (sp in NA_Sp) {
    message("  running ", sp, ". Assigning NA, because absent from CASFRI")
    spRasts[[sp]] <- spRas
    spRasts[[sp]] <- Cache(writeRaster, spRasts[[sp]],
                           filename = asPath(file.path(destinationPath, paste0("CASFRI", sp,".tif"))),
                           overwrite = TRUE, datatype = "INT2U")
  }

  sppTODO <- unique(names(sppListMergesCASFRI))

  for (sp in sppTODO) {
    spCASFRI <- sppListMergesCASFRI[[sp]]
    spRasts[[sp]] <- spRas
    message("starting ", sp)
    if (length(spCASFRI) > 1)
      message("  Merging ", paste(spCASFRI, collapse = ", "), "; becoming: ", sp)
    aa2 <- CASFRIattrLong[
      value %in% spCASFRI][
        , min(100L, sum(pct)), by = GID]
    setkey(aa2, GID)
    cc <- aa2[CASFRIdt] %>% na.omit()
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


