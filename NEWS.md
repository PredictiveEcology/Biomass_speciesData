# Biomass_speciesData (development version)

* `.inputObjects` stops, before any download, when another module has a different `dataYear` (`paramCheckOtherMods(ifSetButDifferent = "error")`). `tests/testthat/test-template.R`, the unmodified SpaDES template (it tested events this module does not have, from a hard-coded path), is replaced by `test-dataYear.R`.
