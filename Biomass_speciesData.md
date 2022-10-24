---
title: "LandR _Biomass_speciesData_ Manual"
date: "Last updated: 2022-10-24"
output:
  bookdown::html_document2:
    toc: true
    toc_float: true
    toc_depth: 4
    theme: sandstone
    number_sections: false
    df_print: paged
    keep_md: yes
editor_options:
  chunk_output_type: console
  markdown: 
    wrap: 80
bibliography: citations/references_Biomass_speciesData.bib
citation-style: citations/ecology-letters.csl
link-citations: true
always_allow_html: true
---

<!-- the following are text references used in captions for LaTeX compatibility -->

(ref:Biomass-speciesData) *Biomass_speciesData*

(ref:percent) %

(ref:Abie-bal) *Abies balsamea*

(ref:Abie-las) *A. lasiocarpa*

(ref:Pinu-con) *Pinus contorta*

(ref:Pice-sp) *Picea spp.*





[![module-version-Badge](D:/GitHub/LandR-Manual/modules/Biomass_speciesData/figures/moduleVersionBadge.png)](https://github.com/PredictiveEcology/Biomass_speciesData/commit/48c20226689f4b7c4ed10032153e96f52e42e860)

[![Issues-badge](D:/GitHub/LandR-Manual/modules/Biomass_speciesData/figures/issuesBadge.png)](https://github.com/PredictiveEcology/Biomass_speciesData/issues)

#### Authors:

Ceres Barros <ceres.barros@ubc.ca> [aut, cre], Eliot J B McIntire <eliot.mcintire@nrcan-rncan.gc.ca> [aut], Alex M. Chubaty <achubaty@for-cast.ca> [aut]
<!-- ideally separate authors with new lines, '\n' not working -->

**This documentation is work in progress. Potential discrepancies and omissions
may exist for the time being. If you find any, contact us using the "Get help"
link above.**

## Module Overview

### Quick links

-   [General functioning](#bsppdata-general-functioning)

-   [List of input objects](#bsppdata-inputs-list)

-   [List of parameters](#bsppdata-params-list)

-   [List of outputs](#bsppdata-outputs-list)

-   [Simulation flow and module events](#bsppdata-sim-flow)

### Module summary

LandR *Biomass_speciesData* (hereafter *Biomass_speciesData*) downloads and
pre-processes species percent ((ref:percent) cover) data layers used by other
LandR data modules (e.g., *Biomass_borealDataPrep*) and by the LandR forest
simulation module *Biomass_core*.

### Links to other modules {#bsppdata-links-modules}

*Biomass_speciesData* is intended to be used with any LandR modules that require
species (ref:percent) cover raster layers (see examples below). See
[here](https://rpubs.com/PredictiveEcology/LandR_Module_Ecosystem) for all
available modules in the LandR ecosystem and select *Biomass_speciesData* from
the drop-down menu to see potential linkages.

-   [*Biomass_borealDataPrep*](https://github.com/PredictiveEcology/Biomass_borealDataPrep):
prepares all parameters and inputs (including initial landscape conditions)
that *Biomass_core* needs to run a realistic simulation. Default
values/inputs produced are relevant for boreal forests of Western Canada.
Used downstream from *Biomass_speciesData*;

-   [*Biomass_core*](https://github.com/PredictiveEcology/Biomass_core): core
forest dynamics simulation module. Used downstream from
*Biomass_speciesData*.

## Module manual

### General functioning {#bsppdata-general-functioning}

*Biomass_speciesData* accesses and processes species (ref:percent) cover data
for the parametrisation and initialisation of LandR *Biomass_core*. This module
ensures 1) that all data use the same geospatial geometries and 2) that these
are correctly re-projected to the study area used for parametrisation
(`studyAreaLarge` polygon), and 3) attempts to sequentially fill-in and replace
the lowest quality data with higher quality data when several data sources are
used. It's primary output is a `RasterStack` of species (ref:percent) cover,
with each layer corresponding to a species.

Currently, the module can access the Canadian National Forest Inventory (NFI)
forest attributes kNN dataset [the default; @BeaudoinEtAl2017], the Common
Attribute Schema for Forest Resource Inventories dataset [CASFRI; @Cosco2011],
the Ontario Forest Resource Inventory (ONFRI), a dataset specific to Alberta
compiled by Paul Pickell, and other Alberta forest inventory datasets. However,
**only the NFI kNN data are freely available** and access to the other datasets
must be granted by module developers and data owners, and requires a Google
account. Nevertheless, the module is flexible enough that any user can use it to
process additional datasets, provided that an adequate R function is passed to
the module (see `types` parameter details in the [list of
parameters](#bsppdata-params-list))

When multiple data sources are used, the module will replace lower quality data
with higher quality data following the order specified in the `types` parameter.

When multiple species of a given data source are to be grouped, (ref:percent)
cover is summed across species of the same group within each pixel. Please see
the `sppEquiv` in the [list of input objects](#bsppdata-inputs-list) for
information on how to define species groups.

The module can also exclude species (ref:percent) cover layers if they don't
have a minimum (ref:percent) cover value in at least one pixel. The user should
still inspect where species is deemed present (e.g., in how many pixels in
total), as it is possible that some datasets only have a few pixels where the
species is present, but with reported high (ref:percent) cover. In this case,
the user may choose to exclude these species *a posteriori*. The summary plot
automatically shown by *Biomass_speciesData* can help diagnose whether certain
species are present in very few pixels (see Fig.
\@ref(fig:fig-Biomass-speciesDataOutPlots)).

### List of input objects {#bsppdata-inputs-list}

Below is the full list of input objects that *Biomass_speciesData* requires
(Table \@ref(tab:moduleInputs2-Biomass-speciesData)). Of these, the only input
that **must** be provided (i.e., *Biomass_speciesData* does not have a default
for) is `studyAreaLarge`.

Of the inputs in Table \@ref(tab:moduleInputs2-Biomass-speciesData), the
following are particularly important and deserve special attention:

-   `studyAreaLarge` -- the polygon defining the area for which species cover
data are desired. It can be larger (but never smaller) that the study area
used in the simulation of forest dynamics (i.e., `studyArea` object in
*Biomass_core*), in which case it should fully cover it.

-   `sppEquiv` -- a table of correspondences between different species naming
conventions. This table is used across several LandR modules, including
*Biomass_core*. It is particularly important here because it will determine
whether and which species (and their cover layers) are merged. For instance,
if the user wishes to simulate a generic *Picea spp.* that includes, *Picea
glauca*, *Picea mariana* and *Picea engelmannii*, they will need to provide
these three species names in the data column (e.g., `KNN` if obtaining
forest attribute kNN data layers from the National Forest Inventory), but
the same name (e.g., "Pice_Spp") in the column chosen for the naming
convention used throughout the simulation (defined by the `sppEquivCol`
parameter). See Table \@ref(tab:mergingSpp-Biomass-speciesData) for an
example.

<table>
<caption>(\#tab:mergingSpp-Biomass-speciesData)Example of species merging for simulation. Here the user wants to model (ref:Abie-bal), (ref:Abie-las) and (ref:Pinu-con) as separate species, but all (ref:Pice-sp) as a genus-level group. For this, all six species are separately identified in the 'KNN' column, so that their (ref:percent) cover layers can be obtained, but in the 'Boreal' column (which defines the naming convention used in the simulation in this example) all (ref:Pice-sp) have the same name. (ref:Biomass-speciesData) will merge their (ref:percent) cover data into a single layer by summing their cover per pixel.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> Species </th>
   <th style="text-align:left;"> KNN </th>
   <th style="text-align:left;"> Boreal </th>
   <th style="text-align:left;"> Modelled as </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;font-style: italic;"> Abies balsamea </td>
   <td style="text-align:left;"> Abie_Bal </td>
   <td style="text-align:left;"> Abie_Bal </td>
   <td style="text-align:left;font-style: italic;"> Abies balsamea </td>
  </tr>
  <tr>
   <td style="text-align:left;font-style: italic;"> Abies lasiocarpa </td>
   <td style="text-align:left;"> Abie_Las </td>
   <td style="text-align:left;"> Abie_Las </td>
   <td style="text-align:left;font-style: italic;"> Abies lasiocarpa </td>
  </tr>
  <tr>
   <td style="text-align:left;font-style: italic;"> Picea engelmannii </td>
   <td style="text-align:left;"> Pice_Eng </td>
   <td style="text-align:left;"> Pice_Spp </td>
   <td style="text-align:left;font-style: italic;"> Picea spp. </td>
  </tr>
  <tr>
   <td style="text-align:left;font-style: italic;"> Picea glauca </td>
   <td style="text-align:left;"> Pice_Gla </td>
   <td style="text-align:left;"> Pice_Spp </td>
   <td style="text-align:left;font-style: italic;"> Picea spp. </td>
  </tr>
  <tr>
   <td style="text-align:left;font-style: italic;"> Picea mariana </td>
   <td style="text-align:left;"> Pice_Mar </td>
   <td style="text-align:left;"> Pice_Spp </td>
   <td style="text-align:left;font-style: italic;"> Picea spp. </td>
  </tr>
  <tr>
   <td style="text-align:left;font-style: italic;"> Pinus contorta </td>
   <td style="text-align:left;"> Pinu_Con </td>
   <td style="text-align:left;"> Pinu_Con </td>
   <td style="text-align:left;font-style: italic;"> Pinus contorta </td>
  </tr>
</tbody>
</table>

\newpage
\blandscape

<table>
<caption>(\#tab:moduleInputs2-Biomass-speciesData)List of (ref:Biomass-speciesData) input objects and their description.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> objectName </th>
   <th style="text-align:left;"> objectClass </th>
   <th style="text-align:left;"> desc </th>
   <th style="text-align:left;"> sourceURL </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> rasterToMatchLarge </td>
   <td style="text-align:left;"> RasterLayer </td>
   <td style="text-align:left;"> a raster of `studyAreaLarge` in the same resolution and projection the simulation's. Defaults to the using the Canadian Forestry Service, National Forest Inventory, kNN-derived stand biomass map. </td>
   <td style="text-align:left;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sppColorVect </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> A named vector of colors to use for plotting. The names must be in `sim$sppEquiv[[sim$sppEquivCol]]`, and should also contain a color for 'Mixed' </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sppEquiv </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> table of species equivalencies. See `LandR::sppEquivalencies_CA`. </td>
   <td style="text-align:left;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> studyAreaLarge </td>
   <td style="text-align:left;"> SpatialPolygonsDataFrame </td>
   <td style="text-align:left;"> Polygon to use as the parametrisation study area. Must be provided by the user. Note that `studyAreaLarge` is only used for parameter estimation, and can be larger than the actual study area used for LandR simulations (e.g, larger than `studyArea` in LandR Biomass_core). </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> studyAreaReporting </td>
   <td style="text-align:left;"> SpatialPolygonsDataFrame </td>
   <td style="text-align:left;"> multipolygon (typically smaller/unbuffered than `studyAreaLarge` and `studyArea` in LandR Biomass_core) to use for plotting/reporting. If not provided, will default to `studyAreaLarge`. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
</tbody>
</table>

\elandscape

### List of parameters {#bsppdata-params-list}

Table \@ref(tab:moduleParams2-Biomass-speciesData) lists all parameters used in
*Biomass_speciesData* and their detailed information. All these parameters have
default values specified in the module's metadata.

Of these parameters, the following are particularly important:

-   `coverThresh` -- integer. Defines a minimum (ref:percent) cover value (from
0-100) that the species must have in at least one pixel to be considered
present in the study area, otherwise it is excluded from the final stack of
species layers (`speciesLayers`). Note that this will affect what species
have data for an eventual simulation and the user will need to adjust
simulation parameters accordingly (e.g., species in trait tables will need
to match the species in `speciesLayers`).

-   `types` -- character. Which (ref:percent) cover data sources are to be used
(see [General functioning](#bsppdata-general-functioning)). Several data
sources can be passed, in which case the module will overlay the lower
quality layers with higher quality ones following the order of data sources
in `types`. For instance, if
`types == c("KNN", "CASFRI", "ForestInventory")`, KNN is assumed to be the
lowest quality data set and ForestInventory the highest, hence values in KNN
layers are replaced with overlapping values from CASFRI layers and values
from KNN and CASFRI layers are replaced with overlapping values of
ForestInventory layers.

\newpage
\blandscape

<table>
<caption>(\#tab:moduleParams2-Biomass-speciesData)List of (ref:Biomass-speciesData) parameters and their description.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> paramName </th>
   <th style="text-align:left;"> paramClass </th>
   <th style="text-align:left;"> default </th>
   <th style="text-align:left;"> min </th>
   <th style="text-align:left;"> max </th>
   <th style="text-align:left;"> paramDesc </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> coverThresh </td>
   <td style="text-align:left;"> integer </td>
   <td style="text-align:left;"> 10 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> The minimum % cover a species needs to have (per pixel) in the study area to be considered present </td>
  </tr>
  <tr>
   <td style="text-align:left;"> dataYear </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 2001 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Passed to `paste0('prepSpeciesLayers_', types)` function to fetch data from that year (if applicable). Defaults to 2001 as the default kNN year. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sppEquivCol </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> Boreal </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> The column in `sim$sppEquiv` data.table to group species by and use as a naming convention. If different species in, e.g., the kNN data have the same name in the chosen column, their data are merged into one species by summing their % cover in each raster cell. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> types </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> KNN </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> The possible data sources. These must correspond to a function named `paste0('prepSpeciesLayers_', types)`. Defaults to 'KNN' to get the Canadian Forestry Service, National Forest Inventory, kNN-derived species cover maps from year 'dataYear', using the `LandR::prepSpeciesLayers_KNN` function (see https://open.canada.ca/ data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990 for details on these data). Other currently available options are 'ONFRI', 'CASFRI', 'Pickell' and 'ForestInventory', which attempt to get proprietary data - the user must be granted access first. A custom function can be used to retrieve any data, just as long as it is accessible by the module (e.g., in the global environment) and is named as `paste0('prepSpeciesLayers_', types)`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> vegLeadingProportion </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 0.8 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> a number that defines whether a species is leading for a given pixel. Only used for plotting. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .plotInitialTime </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> This describes the simulation time at which the first plot event should occur </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .plotInterval </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> This describes the simulation time interval between plot events </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .saveInitialTime </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> This describes the simulation time at which the first save event should occur </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .saveInterval </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> This describes the simulation time interval between save events </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .sslVerify </td>
   <td style="text-align:left;"> integer </td>
   <td style="text-align:left;"> 64 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Passed to `httr::config(ssl_verifypeer = P(sim)$.sslVerify)` when downloading KNN (NFI) datasets. Set to 0L if necessary to bypass checking the SSL certificate (this may be necessary when NFI's FTP website SSL certificate is down/out-of-date). </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .studyAreaName </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Human-readable name for the study area used. If NA, a hash of `studyAreaLarge` will be used. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .useCache </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> init </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Controls cache; caches the init event by default </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .useParallel </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 16 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Used in reading csv file with fread. Will be passed to `data.table::setDTthreads`. </td>
  </tr>
</tbody>
</table>

\elandscape

### List of outputs {#bsppdata-outputs-list}

The module produces the outputs in Table
\@ref(tab:moduleOutputs-Biomass-speciesData), and automatically saves the
processed species cover layers in the output path defined in
`getPaths(sim)$outputPath`.

<table>
<caption>(\#tab:moduleOutputs-Biomass-speciesData)List of (ref:Biomass-speciesData) output objects and their description.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> objectName </th>
   <th style="text-align:left;"> objectClass </th>
   <th style="text-align:left;"> desc </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> speciesLayers </td>
   <td style="text-align:left;"> RasterStack </td>
   <td style="text-align:left;"> biomass percentage raster layers by species in Canada species map </td>
  </tr>
  <tr>
   <td style="text-align:left;"> treed </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> Table with one logical column for each species, indicating whether there were non-zero cover values in each pixel. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> numTreed </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> a named vector with number of pixels with non-zero cover values for each species </td>
  </tr>
  <tr>
   <td style="text-align:left;"> nonZeroCover </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> A single value indicating how many pixels have non-zero cover </td>
  </tr>
</tbody>
</table>

### Simulation flow and module events {#bsppdata-sim-flow}

*Biomass_speciesData* initialises itself and prepares all inputs provided that
it has internet access to download the raw data layers, or that these layers
have been previously downloaded and stored in the folder specified by
`options("reproducible.destinationPath")`[^biomass_speciesdata-1].

The module defaults to processing cover data fo all species listed in the
`Boreal` column of the default `sppEquiv` input `data.table` object, for which
there are available (ref:percent) cover layers in the kNN dataset (Table
\@ref(tab:defaultSppBiomass-speciesData); see `?LandR::sppEquivalencies_CA` for
more information):

<table>
<caption>(\#tab:defaultSppBiomass-speciesData)List of species cover data downloaded by default by (ref:Biomass-speciesData).</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> Species </th>
   <th style="text-align:left;"> Generic name </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;font-style: italic;"> Abies balsamea </td>
   <td style="text-align:left;"> Balsam Fir </td>
  </tr>
  <tr>
   <td style="text-align:left;font-style: italic;"> Abies lasiocarpa </td>
   <td style="text-align:left;"> Fir </td>
  </tr>
  <tr>
   <td style="text-align:left;font-style: italic;"> Acer negundo </td>
   <td style="text-align:left;"> Boxelder maple </td>
  </tr>
  <tr>
   <td style="text-align:left;font-style: italic;"> Acer pensylvanicum </td>
   <td style="text-align:left;"> Striped maple </td>
  </tr>
  <tr>
   <td style="text-align:left;font-style: italic;"> Acer saccharum </td>
   <td style="text-align:left;"> Sugar maple </td>
  </tr>
  <tr>
   <td style="text-align:left;font-style: italic;"> Acer spicatum </td>
   <td style="text-align:left;"> Mountain maple </td>
  </tr>
  <tr>
   <td style="text-align:left;font-style: italic;"> Acer spp. </td>
   <td style="text-align:left;"> Maple </td>
  </tr>
  <tr>
   <td style="text-align:left;font-style: italic;"> Alnus spp </td>
   <td style="text-align:left;"> Alder </td>
  </tr>
  <tr>
   <td style="text-align:left;font-style: italic;"> Betula alleghaniensis </td>
   <td style="text-align:left;"> Swamp birch </td>
  </tr>
  <tr>
   <td style="text-align:left;font-style: italic;"> Betula papyrifera </td>
   <td style="text-align:left;"> Paper birch </td>
  </tr>
  <tr>
   <td style="text-align:left;font-style: italic;"> Betula populifolia </td>
   <td style="text-align:left;"> Gray birch </td>
  </tr>
  <tr>
   <td style="text-align:left;font-style: italic;"> Betula spp. </td>
   <td style="text-align:left;"> Birch </td>
  </tr>
  <tr>
   <td style="text-align:left;font-style: italic;"> Fagus grandifolia </td>
   <td style="text-align:left;"> American beech </td>
  </tr>
  <tr>
   <td style="text-align:left;font-style: italic;"> Fraxinus americana </td>
   <td style="text-align:left;"> American ash </td>
  </tr>
  <tr>
   <td style="text-align:left;font-style: italic;"> Fraxinus nigra </td>
   <td style="text-align:left;"> Black ash </td>
  </tr>
  <tr>
   <td style="text-align:left;font-style: italic;"> Larix laricina </td>
   <td style="text-align:left;"> Tamarack </td>
  </tr>
  <tr>
   <td style="text-align:left;font-style: italic;"> Larix lyallii </td>
   <td style="text-align:left;"> Alpine larch </td>
  </tr>
  <tr>
   <td style="text-align:left;font-style: italic;"> Larix occidentalis </td>
   <td style="text-align:left;"> Western larch </td>
  </tr>
  <tr>
   <td style="text-align:left;font-style: italic;"> Larix spp. </td>
   <td style="text-align:left;"> Larch </td>
  </tr>
  <tr>
   <td style="text-align:left;font-style: italic;"> Picea engelmannii </td>
   <td style="text-align:left;"> Engelmann's spruce </td>
  </tr>
  <tr>
   <td style="text-align:left;font-style: italic;"> Picea glauca </td>
   <td style="text-align:left;"> White.Spruce </td>
  </tr>
  <tr>
   <td style="text-align:left;font-style: italic;"> Picea mariana </td>
   <td style="text-align:left;"> Black.Spruce </td>
  </tr>
  <tr>
   <td style="text-align:left;font-style: italic;"> Picea spp. </td>
   <td style="text-align:left;"> Spruce </td>
  </tr>
  <tr>
   <td style="text-align:left;font-style: italic;"> Pinus albicaulis </td>
   <td style="text-align:left;"> Whitebark pine </td>
  </tr>
  <tr>
   <td style="text-align:left;font-style: italic;"> Pinus banksiana </td>
   <td style="text-align:left;"> Jack pine </td>
  </tr>
  <tr>
   <td style="text-align:left;font-style: italic;"> Pinus contorta </td>
   <td style="text-align:left;"> Lodgepole pine </td>
  </tr>
  <tr>
   <td style="text-align:left;font-style: italic;"> Pinus monticola </td>
   <td style="text-align:left;"> Western white pine </td>
  </tr>
  <tr>
   <td style="text-align:left;font-style: italic;"> Pinus resinosa </td>
   <td style="text-align:left;"> Red pine </td>
  </tr>
  <tr>
   <td style="text-align:left;font-style: italic;"> Pinus spp. </td>
   <td style="text-align:left;"> Pine </td>
  </tr>
  <tr>
   <td style="text-align:left;font-style: italic;"> Populus balsamifera v. balsamifera </td>
   <td style="text-align:left;"> Balsam poplar </td>
  </tr>
  <tr>
   <td style="text-align:left;font-style: italic;"> Populus trichocarpa </td>
   <td style="text-align:left;"> Black cottonwood </td>
  </tr>
  <tr>
   <td style="text-align:left;font-style: italic;"> Populus grandidentata </td>
   <td style="text-align:left;"> White poplar </td>
  </tr>
  <tr>
   <td style="text-align:left;font-style: italic;"> Populus spp. </td>
   <td style="text-align:left;"> Poplar </td>
  </tr>
  <tr>
   <td style="text-align:left;font-style: italic;"> Populus tremuloides </td>
   <td style="text-align:left;"> Trembling poplar </td>
  </tr>
  <tr>
   <td style="text-align:left;font-style: italic;"> Tsuga canadensis </td>
   <td style="text-align:left;"> Eastern hemlock </td>
  </tr>
  <tr>
   <td style="text-align:left;font-style: italic;"> Tsuga spp. </td>
   <td style="text-align:left;"> Hemlock </td>
  </tr>
</tbody>
</table>

*Biomass_speciesData* only runs two events, the `init` event where all species
cover layers are processed and a plotting event (`initPlot`) that plots the
final layers.

The general flow of *Biomass_speciesData* processes is:

1.  Download (if necessary) and spatial processing of species cover layers from
the first data source listed in the `types` parameter. Spatial processing
consists in sub-setting the data to the area defined by `studyAreaLarge` and
ensuring that the spatial projection and resolution match those of
`rasterToMatchLarge`. After spatial processing, species layers that have no
pixels with values $\ge$ `coverThresh` are excluded.

2.  If more than one data source is listed in `types`, the second set of species
cover layers is downloaded and processed as above.

3.  The second set of layers is assumed to be the highest quality dataset and
used to replaced overlapping pixel values on the first (including for
species whose layers may have been initially excluded after applying the
`coverThresh` filter).

4.  Steps 2 and 3 are repeated for remaining data sources listed in `types`.

5.  Final layers are saved to disk and plotted (`initPlot` event). A summary of
number of pixels with forest cover are calculated (`treed`and `numTreed`
output objects; see [list of outputs](#bsppdata-outputs-list)).

[^biomass_speciesdata-1]: Raw data layers downloaded by the module are saved in
\`dataPath(sim)\`, which can be controlled via
\`options(reproducible.destinationPath = ...)\`.

## Usage example {#bsppdata-example}

This module can be run stand-alone, but it only compiles species (ref:percent)
cover data into layers used by other modules.

### Load `SpaDES` and other packages.

### Set up R libraries {#bsppdata-example-libs}


```r
options(repos = c(CRAN = "https://cloud.r-project.org"))
# tempDir <- tempdir()
tempDir <- "C:/Users/cbarros/AppData/Local/Temp/Biomass_sppData-example"

pkgPath <- file.path(tempDir, "packages", version$platform,
                     paste0(version$major, ".", strsplit(version$minor, "[.]")[[1]][1]))
dir.create(pkgPath, recursive = TRUE)
.libPaths(pkgPath, include.site = FALSE)

if (!require(Require, lib.loc = pkgPath)) {
  install.packages("Require")
  library(Require, lib.loc = pkgPath)
}

setLinuxBinaryRepo()
```

### Get the module and module dependencies {#bsppdataexample-pkg-mods}


```r
Require("PredictiveEcology/SpaDES.project@6d7de6ee12fc967c7c60de44f1aa3b04e6eeb5db", 
        require = FALSE, upgrade = FALSE, standAlone = TRUE)

paths <- list(inputPath = normPath(file.path(tempDir, "inputs")), 
              cachePath = normPath(file.path(tempDir, "cache")), 
              modulePath = normPath(file.path(tempDir, "modules")), 
              outputPath = normPath(file.path(tempDir, "outputs")))

SpaDES.project::getModule(modulePath = paths$modulePath,
                          c("PredictiveEcology/Biomass_speciesData@master"),
                          overwrite = TRUE)

## make sure all necessary packages are installed:
outs <- SpaDES.project::packagesInModules(modulePath = paths$modulePath)
Require(c(unname(unlist(outs)), "SpaDES"),
        require = FALSE, standAlone = TRUE)

## load necessary packages
Require(c("SpaDES", "LandR", "reproducible"), upgrade = FALSE, install = FALSE)
```

### Setup simulation

For this demonstration we are using all default parameter values, except
`coverThresh` , which is lowered to 5(ref:percent). The species layers (the
major output of interest) are saved automatically, so there is no need to tell
`spades` what to save using the `outputs` argument (see
`?SpaDES.core::outputs`).

We pass the global parameter `.plotInitialTime = 1` in the `simInitAndSpades`
function to activate plotting.


```r
# User may want to set some options -- see ?reproducibleOptions 
#    -- e.g., often the path to the 'inputs' folder will be set outside of project by user:
# options(reproducible.inputPaths = "E:/Data/LandR_related/") # to re-use datasets across projects
studyAreaLarge <- Cache(randomStudyArea, size = 1e7,
                        cacheRepo = paths$cachePath) # cache this so it creates a random one only once on a machine

# Pick the species you want to work with -- here we use the naming convention in "Boreal" column of LandR::sppEquivalencies_CA (default)
speciesNameConvention <- "Boreal"
speciesToUse <- c("Pice_Gla", "Popu_Tre", "Pinu_Con")

sppEquiv <- LandR::sppEquivalencies_CA[get(speciesNameConvention) %in% speciesToUse]
# Assign a colour convention for graphics for each species
sppColorVect <- LandR::sppColors(sppEquiv, speciesNameConvention,
                                 newVals = "Mixed", palette = "Set1")

## Usage example
modules <- list("Biomass_speciesData")
objects <- list("studyAreaLarge" = studyAreaLarge,
                "sppEquiv" = sppEquiv,
                "sppColorVect" = sppColorVect)
params <- list("Biomass_speciesData" = list("coverThresh" = 5L))
```

### Run module

Note that because this is a data module (i.e., only attempts to prepare data for
the simulation) we are not iterating it and so both the start and end times are
set to `1` here.


```r
opts <- options(reproducible.useCache = TRUE,
                spades.inputPath = paths$inputPath)

mySimOut <- simInitAndSpades(times = list(start = 1, end = 1),
                             modules = modules, 
                             parameters = params,
                             objects = objects, 
                             paths = paths,
                             .plotInitialTime = 1)
options(opts)
```

Here are some of outputs of *Biomass_speciesData* (dominant species) in a
randomly generated study area within Canada.

<div class="figure" style="text-align: center">
<img src="D:/GitHub/LandR-Manual/modules/Biomass_speciesData/figures/testRunFigure.png" alt="(ref:Biomass-speciesData) automatically generates a plot of species dominance and number of presences in the study area when `.plotInitialTime=1` is passed as an argument." width="70%" />
<p class="caption">(\#fig:fig-Biomass-speciesDataOutPlots)(ref:Biomass-speciesData) automatically generates a plot of species dominance and number of presences in the study area when `.plotInitialTime=1` is passed as an argument.</p>
</div>

## References {#bsppdata-refs}
