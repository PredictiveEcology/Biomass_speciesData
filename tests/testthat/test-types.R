## Only SCANFI is available. The older data sources (kNN, CASFRI, ForestInventory, MBFRI, ONFRI,
## Pickell, NTEMS) were removed: they could not be run from the module (they need `studyArea`
## and `rasterToMatch`, which it does not pass, and NTEMS was never implemented).

test_that("a types other than SCANFI stops with a clear error", {
  for (tp in list("KNN", c("KNN", "SCANFI"))) {
    expect_error(
      SpaDES.core::simInit(
        times = list(start = 0, end = 1),
        params = list(Biomass_speciesData = list(types = tp)),
        modules = moduleName,
        paths = testPaths
      ),
      "`types` must be \"SCANFI\".*were removed"
    )
  }
})
