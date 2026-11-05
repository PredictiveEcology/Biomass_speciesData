## When `rasterToMatch_biomassParam` was not supplied, .inputObjects built it over
## `studyArea_biomassParam` only when the two study areas had different CRSs; otherwise it used
## `rasterToMatch`, and the species layers are cropped to it, so a larger `studyArea_biomassParam`
## had no effect. The default is now the `rasterToMatch` grid extended over
## `studyArea_biomassParam` (same rule as Biomass_borealDataPrep, which crops it back to `rasterToMatch`).

rtmBiomassParamBranch <- function() {
  exprs <- parse(testthat::test_path("..", "..", "Biomass_speciesData.R"), keep.source = FALSE)
  def <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                  identical(x[[2]], as.name(".inputObjects")), exprs)
  stopifnot(length(def) == 1L)
  bod <- as.list(body(eval(def[[1]][[3]])))
  branch <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("if")) &&
                     identical(deparse(x[[2]]), "!suppliedElsewhere(\"rasterToMatch_biomassParam\", sim)"),
                   bod)
  stopifnot(length(branch) == 1L)
  branch[[1]]
}

squarePolyRTM <- function(xmin, xmax, ymin, ymax, crs = "EPSG:3978") {
  sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(xmin, ymin), c(xmax, ymin), c(xmax, ymax), c(xmin, ymax), c(xmin, ymin)
  ))), crs = crs))
}

runRtmBiomassParamDefault <- function(studyArea, studyArea_biomassParam, rasterToMatch) {
  e <- new.env()
  e$sim <- new.env()
  e$sim$studyArea <- studyArea
  e$sim$studyArea_biomassParam <- studyArea_biomassParam
  e$sim$rasterToMatch <- rasterToMatch
  e$suppliedElsewhere <- function(...) FALSE
  eval(rtmBiomassParamBranch(), envir = e)
  e$sim$rasterToMatch_biomassParam
}

## 240 m grid with an origin that is not a multiple of 240 from 0, as from rast(studyArea, res = 240)
sa <- squarePolyRTM(1010, 3410, 1010, 3410)
rtm <- terra::rast(terra::vect(sa), res = c(240, 240), vals = 1) |> terra::mask(mask = terra::vect(sa))

test_that("without a larger studyArea_biomassParam, the default is rasterToMatch", {
  out <- runRtmBiomassParamDefault(sa, sa, rtm)
  expect_identical(out, rtm)
})

test_that("a larger studyArea_biomassParam gives a larger default on the rasterToMatch grid", {
  sa_bp <- squarePolyRTM(-5000, 8000, -4000, 9000)
  out <- runRtmBiomassParamDefault(sa, sa_bp, rtm)

  expect_true(terra::ncell(out) > terra::ncell(rtm))
  expect_equal(terra::res(out), terra::res(rtm))
  expect_true(terra::same.crs(out, rtm))
  expect_equal(terra::origin(out), terra::origin(rtm))
  e_out <- as.vector(terra::ext(out))
  e_bp <- as.vector(terra::ext(terra::vect(sa_bp)))
  expect_true(all(e_out[c("xmin", "ymin")] <= e_bp[c("xmin", "ymin")]) &&
                all(e_out[c("xmax", "ymax")] >= e_bp[c("xmax", "ymax")]))
  ## masked like the default rasterToMatch: cells touching studyArea_biomassParam are kept
  ## (a triangle around studyArea, so some cells of its extent are outside it)
  sa_bp_tri <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(-5000, -5000), c(15000, -5000), c(-5000, 15000), c(-5000, -5000)
  ))), crs = "EPSG:3978"))
  out_tri <- runRtmBiomassParamDefault(sa, sa_bp_tri, rtm)
  expect_true(anyNA(terra::values(out_tri)))
  expect_equal(is.na(as.vector(terra::values(out_tri))),
               is.na(as.vector(terra::values(terra::rasterize(terra::vect(sa_bp_tri), out_tri, touches = TRUE)))))

  ## the downsize step crops it back to exactly the rasterToMatch grid
  idx <- terra::setValues(out, seq(terra::ncell(out)))
  cropped <- terra::crop(idx, rtm) |> terra::mask(rtm)
  expect_true(terra::compareGeom(cropped, rtm, stopOnError = FALSE))
  expect_equal(is.na(as.vector(terra::values(cropped))), is.na(as.vector(terra::values(rtm))))
})

test_that("a studyArea_biomassParam in another CRS is projected to the rasterToMatch CRS", {
  sa_bp <- squarePolyRTM(-5000, 8000, -4000, 9000) |> sf::st_transform(4326)
  out <- runRtmBiomassParamDefault(sa, sa_bp, rtm)
  expect_true(terra::same.crs(out, rtm))
  expect_equal(terra::res(out), terra::res(rtm))
  expect_equal(terra::origin(out), terra::origin(rtm))
  expect_true(terra::ncell(out) > terra::ncell(rtm))
})

## The reverse: rasterToMatch_biomassParam supplied without studyArea_biomassParam. The polygon
## defaulted to studyArea, so the masks used the smaller area.
saBiomassParamBranch <- function() {
  exprs <- parse(testthat::test_path("..", "..", "Biomass_speciesData.R"), keep.source = FALSE)
  def <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                  identical(x[[2]], as.name(".inputObjects")), exprs)
  bod <- as.list(body(eval(def[[1]][[3]])))
  branch <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("if")) &&
                     identical(deparse(x[[2]]), "!suppliedElsewhere(\"studyArea_biomassParam\", sim)"),
                   bod)
  stopifnot(length(branch) == 1L)
  branch[[1]]
}

runSaBiomassParamDefault <- function(studyArea, rasterToMatch_biomassParam = NULL) {
  e <- new.env()
  e$sim <- new.env()
  e$sim$studyArea <- studyArea
  e$sim$rasterToMatch_biomassParam <- rasterToMatch_biomassParam
  supplied <- if (is.null(rasterToMatch_biomassParam)) character(0) else "rasterToMatch_biomassParam"
  e$suppliedElsewhere <- function(object, ...) object %in% supplied
  eval(saBiomassParamBranch(), envir = e)
  e$sim$studyArea_biomassParam
}

test_that("without either biomassParam object, studyArea_biomassParam is studyArea", {
  expect_identical(runSaBiomassParamDefault(sa), sa)
})

test_that("a supplied rasterToMatch_biomassParam sets studyArea_biomassParam, not studyArea", {
  sa_bp <- squarePolyRTM(-5000, 8000, -4000, 9000)
  rtm_bp <- terra::rast(terra::vect(sa_bp), res = c(240, 240), vals = 1) |>
    terra::mask(mask = terra::vect(sa_bp))
  out <- runSaBiomassParamDefault(sa, rtm_bp)
  expect_s4_class(out, "SpatVector")
  expect_equal(nrow(out), 1L)
  expect_true(terra::same.crs(out, rtm_bp))
  ## covers the raster's non-NA cells, so it is larger than studyArea
  expect_equal(as.vector(terra::ext(out)), as.vector(terra::ext(rtm_bp)))
  expect_true(terra::expanse(out) > as.numeric(sf::st_area(sa)))
})
