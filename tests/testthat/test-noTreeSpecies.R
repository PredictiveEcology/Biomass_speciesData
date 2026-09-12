## Regression test for the "no tree species" path (non-forested ELFs): biomassDataInit()
## must carry a zero-layer speciesLayers through instead of stopping.

test_that("empty speciesLayers survives the assertion and the sppKeep subsetting", {
  skip_if_not_installed("terra")

  template <- terra::rast(terra::ext(0, 100, 0, 100), resolution = 10, crs = "EPSG:3978")
  template[] <- 1

  empty <- LandR:::.emptySpatRaster(template)

  ## contract shape
  expect_equal(terra::nlyr(empty), 0L)
  expect_equal(names(empty), character(0))
  expect_true(terra::compareGeom(empty, template))

  ## (a) the zero-layer assertion no longer stops
  expect_silent(LandR::assertSpeciesLayers(empty, 10L))

  ## the reason the module guards instead of calling `[[sppKeep]]`: terra errors on character(0)
  expect_error(template[[character(0)]])

  ## (c) the normal path is unchanged: subsetting by kept names still works
  twoLyr <- c(template, template)
  names(twoLyr) <- c("Pice_mar", "Pinu_ban")
  expect_equal(names(twoLyr[["Pice_mar"]]), "Pice_mar")
})
