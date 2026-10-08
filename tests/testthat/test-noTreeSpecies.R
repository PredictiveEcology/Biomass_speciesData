## Regression test for the "no tree species" path (non-forested ELFs): biomassDataInit()
## must carry speciesLayers = NULL through instead of stopping. (A zero-layer SpatRaster was
## tried first: terra cannot wrap(), unwrap() or write one, so it did not survive Cache.)

test_that("NULL speciesLayers survives the assertion; the sppKeep guard is still needed", {
  skip_if_not_installed("terra")

  ## (a) the assertion no longer stops on "no layers"
  expect_silent(LandR::assertSpeciesLayers(NULL, 10L))

  ## the reason the module guards instead of calling `[[sppKeep]]`: terra errors on character(0)
  template <- terra::rast(terra::ext(0, 100, 0, 100), resolution = 10, crs = "EPSG:3978")
  template[] <- 1
  expect_error(template[[character(0)]])

  ## (c) the normal path is unchanged: subsetting by kept names still works
  twoLyr <- c(template, template)
  names(twoLyr) <- c("Pice_mar", "Pinu_ban")
  expect_equal(names(twoLyr[["Pice_mar"]]), "Pice_mar")
})
