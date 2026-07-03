Known issues: <https://github.com/PredictiveEcology/Biomass_speciesData/issues>

# Biomass_speciesData 1.0.5

* switched the default species-layer data source from kNN to SCANFI (default year 2001 to 2020), updating input/output metadata, download logic, and references accordingly
* reworked study-area inputs: renamed `studyAreaLarge` to `studyArea_biomassParam` and `rasterToMatchLarge` to `rasterToMatch_biomassParam`, changed `studyArea` to a `SpatVector`, now require `rasterToMatch` as an input and work in metres throughout, and simplified `.inputObjects`
* changed default `sppEquivCol` to "LandR"
* more robust CRS handling: use `.compareCRS()` instead of `terra::compareGeom()` since the study area may not be a terra object
* plotting: `plotVTM` now uses `Plots()`, with adjusted `.plotInitialTime` handling
* maintenance: raised the minimum `LandR` version, piped caching through `Cache()`, removed the `pryr` dependency, and refreshed CI (`render-module-rmd.yaml`) with rebuilt module `.Rmd`

# Biomass_speciesData 1.0.4 (2024-06-06)

* bumped module to 1.0.4 and raised minimum `SpaDES.core` (>= 2.1.4) and `LandR` versions
* removed stale comments and old code
* modernized GitHub Actions CI (new PredictiveEcology actions and fixed `Require`-based package installation) and refreshed the manual/`.Rmd`, including handling `.sslVerify` as an integer and moving `spades.inputPath` to `reproducible.destinationPath`
