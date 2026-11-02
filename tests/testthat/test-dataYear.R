## dataYear is shared with Biomass_core, Biomass_borealDataPrep, burnSummaries and fireSense_dataPrepPredict;
## a module with another value, set or default, is an error, before any download.
## (Replaces test-template.R, the unmodified SpaDES template, which tested events this module does not have.)

modRoot <- normalizePath(testthat::test_path("..", ".."), winslash = "/")

stubVegModule <- function(dataYear, env = parent.frame()) {
  mp <- withr::local_tempdir(.local_envir = env)
  dir.create(file.path(mp, "vegStub"))
  writeLines(sprintf('
defineModule(sim, list(name = "vegStub", description = "", keywords = "",
  authors = person("A", "B", email = "a@b.com", role = c("aut", "cre")),
  childModules = character(0), version = list(vegStub = "0.0.1"),
  timeframe = as.POSIXlt(c(NA, NA)), timeunit = "year", citation = list(), documentation = list(),
  reqdPkgs = list(),
  parameters = rbind(defineParameter("dataYear", "numeric", %s, NA, NA, "")),
  inputObjects = bindrows(), outputObjects = bindrows()))
doEvent.vegStub <- function(sim, eventTime, eventType, debug = FALSE) invisible(sim)
', dataYear), file.path(mp, "vegStub", "vegStub.R"))
  mp
}

test_that("the default dataYear is 2020", {
  md <- SpaDES.core::moduleMetadata(module = basename(modRoot), path = dirname(modRoot))
  expect_identical(md$parameters$default[[match("dataYear", md$parameters$paramName)]], 2020)
})

test_that("another module's different dataYear stops .inputObjects", {
  withr::local_options(spades.useRequire = FALSE, spades.moduleCodeChecks = FALSE)
  mp <- stubVegModule(2011)
  expect_error(suppressMessages(SpaDES.core::simInit(
    times = list(start = 2020, end = 2020), modules = c(basename(modRoot), "vegStub"),
    paths = list(modulePath = c(dirname(modRoot), mp), inputPath = withr::local_tempdir(),
                 outputPath = withr::local_tempdir(), cachePath = withr::local_tempdir()))),
    "multiple values for dataYear")
})
