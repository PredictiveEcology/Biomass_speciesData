Known issues: <https://github.com/PredictiveEcology/Biomass_speciesData/issues>

# Biomass_speciesData (development version)

* The message for an unset `.studyAreaName` comes from `reproducible::studyAreaName(notSupplied = ".studyAreaName")` (PredictiveEcology/reproducible#638), so it reads the same in every module that uses it: "`.studyAreaName` not supplied; using a hash of `<object>`: <hash>". With an older reproducible the name is the same and there is no message.
* `speciesLayers` are now cropped to `rasterToMatch_biomassParam`, not to the bounding box of `studyArea_biomassParam`. When the grid reached past that box (here, a 240 m grid aggregated from 30 m around an irregular polygon) the layers lost a row and `biomassDataInit` stopped with "[cover] raster dimensions do not match".
* `reqdPkgs` now lists `curl`, `httr` and `raster`, which the module's code uses.

* `vegLeadingProportion` now defaults to `LandR::leadingSpeciesProp()` (option
  `LandR.leadingSpeciesProp`, which takes `LandR.mixedwoodProp`, 0.75, unless set), so the
  leading-species threshold is set once for every module and LandR function instead of being
  hard-coded per module. **The default changes from 0.8 to 0.75**, which changes vegetation type
  maps. Requires LandR >= 1.2.0.9024 (PredictiveEcology/LandR#234).

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
