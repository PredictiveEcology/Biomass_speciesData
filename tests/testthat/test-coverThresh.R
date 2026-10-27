## Regression tests for the species-removal messages in biomassDataInit(). When every species
## layer was below `coverThresh` the module stopped on an undefined `thresh` while building the
## message, before reaching the "no tree species" state (speciesLayers = NULL). When only some
## were below, the message printed `P(sim)$coverThreshresh`, which is NULL.

## run the module's init with a stand-in data source giving these maximum covers per species
runWithCover <- function(cover) {
  rtm <- terra::rast(terra::ext(0, 2400, 0, 2400), resolution = 240, crs = "EPSG:3978", vals = 1)
  sa <- sf::st_as_sf(terra::as.polygons(terra::ext(rtm), crs = terra::crs(rtm)))

  assign("prepSpeciesLayers_TestCover", function(...) {
    lyrs <- terra::rast(lapply(cover, function(x) rtm * x))
    names(lyrs) <- names(cover)
    lyrs
  }, envir = globalenv())
  on.exit(rm("prepSpeciesLayers_TestCover", envir = globalenv()))

  mySim <- SpaDES.core::simInit(
    times = list(start = 0, end = 1),
    params = list(Biomass_speciesData = list(types = "TestCover", coverThresh = 10L,
                                             .plots = NA, .studyAreaName = "test",
                                             .useCache = FALSE)),
    modules = moduleName,
    objects = list(studyArea = sa, studyArea_biomassParam = sa,
                   rasterToMatch = rtm, rasterToMatch_biomassParam = rtm,
                   sppEquiv = LandR::sppEquivalencies_CA[LandR %in% names(cover)],
                   sppNameVector = names(cover)),
    paths = testPaths
  )
  SpaDES.core::spades(mySim, debug = FALSE)
}

test_that("no species above coverThresh gives NULL speciesLayers, not an error", {
  skip_if_not_installed("terra")
  skip_if_not_installed("LandR")
  expect_message(out <- runWithCover(c(Pice_mar = 2, Pinu_ban = 5)),
                 "were found with >=10 % cover")
  expect_null(out$speciesLayers)
})

test_that("species below coverThresh are removed and the message gives the threshold", {
  skip_if_not_installed("terra")
  skip_if_not_installed("LandR")
  expect_message(out <- runWithCover(c(Pice_mar = 50, Pinu_ban = 5)),
                 "removing 1 species because they had <10 % cover")
  expect_identical(names(out$speciesLayers), "Pice_mar")
})
