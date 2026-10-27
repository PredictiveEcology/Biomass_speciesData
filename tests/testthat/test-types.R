## `types` names the data source: the module calls `paste0("prepSpeciesLayers_", type)`, which can
## be a function the user supplies (e.g. in the global environment), not only LandR's SCANFI one.

test_that("a user-supplied prepSpeciesLayers_<type> function is called and its layers used", {
  skip_if_not_installed("terra")
  skip_if_not_installed("LandR")
  rtm <- terra::rast(terra::ext(0, 2400, 0, 2400), resolution = 240, crs = "EPSG:3978", vals = 1)
  sa <- sf::st_as_sf(terra::as.polygons(terra::ext(rtm), crs = terra::crs(rtm)))
  cover <- c(Pice_mar = 50, Pinu_ban = 30)

  received <- NULL
  assign("prepSpeciesLayers_MyData", function(...) {
    received <<- list(...)
    lyrs <- terra::rast(lapply(cover, function(x) rtm * x))
    names(lyrs) <- names(cover)
    lyrs
  }, envir = globalenv())
  on.exit(rm("prepSpeciesLayers_MyData", envir = globalenv()))

  mySim <- SpaDES.core::simInit(
    times = list(start = 0, end = 1),
    params = list(Biomass_speciesData = list(types = "MyData", coverThresh = 10L, dataYear = 2015,
                                             .plots = NA, .studyAreaName = "test",
                                             .useCache = FALSE)),
    modules = moduleName,
    objects = list(studyArea = sa, studyArea_biomassParam = sa,
                   rasterToMatch = rtm, rasterToMatch_biomassParam = rtm,
                   sppEquiv = LandR::sppEquivalencies_CA[LandR %in% names(cover)],
                   sppNameVector = names(cover)),
    paths = testPaths
  )
  out <- SpaDES.core::spades(mySim, debug = FALSE)

  expect_false(is.null(received))
  expect_true(all(c("destinationPath", "outputPath", "to", "cropTo", "projectTo", "studyAreaName",
                    "sppEquiv", "sppEquivCol", "thresh", "year", "dataYear") %in% names(received)))
  expect_identical(received$dataYear, 2015)
  expect_identical(received$thresh, 10L)
  expect_setequal(names(out$speciesLayers), names(cover))
})
