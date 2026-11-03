## Regression test: speciesLayers must have exactly the geometry of rasterToMatch_biomassParam,
## because biomassDataInit() calls cover(speciesLayers, <rasterToMatch_biomassParam>). With only
## `to = studyArea_biomassParam` (a polygon) and `projectTo = rasterToMatch_biomassParam`, the
## output is cropped to the polygon's bounding box, which loses the rows and columns where the
## grid reaches past the polygon: "[cover] raster dimensions do not match" (BC part of ELF 12.1).

test_that("cropTo = rasterToMatch_biomassParam keeps the grid when it reaches past the polygon", {
  skip_if_not_installed("terra")
  skip_if_not_installed("reproducible")
  ## a 240 m grid with whole rows and columns beyond the polygon's bounding box
  rtm <- terra::rast(terra::ext(0, 2400, 0, 2400), resolution = 240, crs = "EPSG:3978", vals = 1)
  poly <- sf::st_as_sf(terra::vect("POLYGON ((500 500, 1900 550, 1850 1900, 550 1850, 500 500))",
                                   crs = "EPSG:3978"))
  ## the source data on another grid and CRS
  src <- terra::rast(terra::ext(-1000, 3500, -1000, 3500), resolution = 30, crs = "EPSG:3978", vals = 1)
  src <- terra::project(src, "EPSG:3005")

  without <- reproducible::postProcessTo(src, to = poly, projectTo = rtm, verbose = -2)
  expect_false(isTRUE(terra::compareGeom(without, rtm, stopOnError = FALSE)))

  with <- reproducible::postProcessTo(src, to = poly, cropTo = rtm, projectTo = rtm, verbose = -2)
  expect_true(terra::compareGeom(with, rtm))
  ## still masked to the polygon
  expect_true(all(is.na(terra::values(terra::mask(with, terra::vect(poly), inverse = TRUE)))))
})
