---
title: "LandR _Biomass_speciesData_ Manual"
subtitle: "v.1.0.0"
date: "Last updated: 2022-03-01"
output:
  bookdown::html_document2:
    toc: true
    toc_float: true
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

# LandR *Biomass_speciesData* Module

<!-- the following are text references used in captions for LaTeX compatibility -->

(ref:Biomass-speciesData) *Biomass_speciesData*



[![made-with-Markdown](https://img.shields.io/badge/Made%20with-Markdown-1f425f.png)](http://commonmark.org)
[![Generic
badge](https://img.shields.io/badge/Get%20help-Report%20issues-%3CCOLOR%3E.png)](https://github.com/PredictiveEcology/Biomass_speciesData/issues)

<!-- if knitting to pdf remember to add the pandoc_args: ["--extract-media", "."] option to yml in order to get the badge images -->

**This documentation is work in progress. Potential discrepancies and omissions
may exist for the time being. If you find any, do contact us using the link
above\^\^**

#### Authors:

Eliot J B McIntire <eliot.mcintire@nrcan-rncan.gc.ca> [aut, cre], Alex M. Chubaty <achubaty@for-cast.ca> [aut], Ceres Barros <cbarros@mail.ubc.ca> [aut]
<!-- ideally separate authors with new lines, '\n' not working -->

## Module Overview

### Module summary

This module downloads and pre-process species % cover data layers to be passed
to other LandR data modules (e.g., *Biomass_borealDataPrep*) or to the LandR
forest simulation module *Biomass_core*.

### Module inputs and parameters at a glance

Below are the full list of input objects (Table
\@ref(tab:moduleInputs-Biomass-speciesData)) and parameters (Table
\@ref(tab:moduleParams-Biomass-speciesData)) that *Biomass_speciesData* expects.
Of these, the only input that **must** be provided (i.e., *Biomass_speciesData*
does not have a default for) is `studyAreaLarge`.

Raw data layers downloaded by the module are saved in `dataPath(sim)`, which can
be controlled via `options(reproducible.destinationPath = ...)`.

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>List of (ref:Biomass-speciesData) input objects and their description.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> objectName </th>
   <th style="text-align:left;"> desc </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> rasterToMatchLarge </td>
   <td style="text-align:left;"> a raster of `studyAreaLarge` in the same resolution and projection the simulation's. Defaults to the using the Canadian Forestry Service, National Forest Inventory, kNN-derived stand biomass map. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sppColorVect </td>
   <td style="text-align:left;"> A named vector of colors to use for plotting. The names must be in sim$sppEquiv[[sim$sppEquivCol]], and should also contain a color for 'Mixed' </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sppEquiv </td>
   <td style="text-align:left;"> table of species equivalencies. See `LandR::sppEquivalencies_CA`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> studyAreaLarge </td>
   <td style="text-align:left;"> Polygon to use as the parametrisation study area. Must be provided by the user. Note that `studyAreaLarge` is only used for parameter estimation, and can be larger than the actual study area used for LandR simulations (e.g, larger than `studyArea` in LandR Biomass_core). </td>
  </tr>
  <tr>
   <td style="text-align:left;"> studyAreaReporting </td>
   <td style="text-align:left;"> multipolygon (typically smaller/unbuffered than `studyAreaLarge` and `studyArea` in LandR Biomass_core) to use for plotting/reporting. If not provided, will default to `studyAreaLarge`. </td>
  </tr>
</tbody>
</table>

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>List of (ref:Biomass-speciesData) parameters and their description.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> paramName </th>
   <th style="text-align:left;"> paramDesc </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> coverThresh </td>
   <td style="text-align:left;"> The minimum % cover a species needs to have (per pixel) in the study area to be considered present </td>
  </tr>
  <tr>
   <td style="text-align:left;"> dataYear </td>
   <td style="text-align:left;"> Passed to `paste0('prepSpeciesLayers_', types)` function to fetch data from that year (if applicable). Defaults to 2001 as the default kNN year. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sppEquivCol </td>
   <td style="text-align:left;"> The column in `sim$sppEquiv` data.table to group species by and use as a naming convention. If different species in, e.g., the kNN data have the same name in the chosen column, their data are merged into one species by summing their % cover in each raster cell. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> types </td>
   <td style="text-align:left;"> The possible data sources. These must correspond to a function named paste0('prepSpeciesLayers_', types). Defaults to 'KNN' to get the Canadian Forestry Service, National Forest Inventory, kNN-derived species cover maps from year 'dataYear', using the `LandR::prepSpeciesLayers_KNN` function (see https://open.canada.ca/ data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990 for details on these data). Other currently available options are 'ONFRI', 'CASFRI', 'Pickell' and 'ForestInventory', which attempt to get proprietary data - the user must be granted access first. A custom function can be used to retrieve any data, just as long as it is accessible by the module (e.g., in the global environment) and is named as paste0('prepSpeciesLayers_', types). </td>
  </tr>
  <tr>
   <td style="text-align:left;"> vegLeadingProportion </td>
   <td style="text-align:left;"> a number that defines whether a species is leading for a given pixel. Only used for plotting. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .plotInitialTime </td>
   <td style="text-align:left;"> This describes the simulation time at which the first plot event should occur </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .plotInterval </td>
   <td style="text-align:left;"> This describes the simulation time interval between plot events </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .saveInitialTime </td>
   <td style="text-align:left;"> This describes the simulation time at which the first save event should occur </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .saveInterval </td>
   <td style="text-align:left;"> This describes the simulation time interval between save events </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .studyAreaName </td>
   <td style="text-align:left;"> Human-readable name for the study area used. If NA, a hash of `studyAreaLarge` will be used. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .useCache </td>
   <td style="text-align:left;"> Controls cache; caches the init event by default </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .useParallel </td>
   <td style="text-align:left;"> Used in reading csv file with fread. Will be passed to data.table::setDTthreads. </td>
  </tr>
</tbody>
</table>

### Events

*Biomass_speciesData* only runs two events:

-   Module "initiation" (`init` event), during which all species % cover layers
    are downloaded and processed.
-   Plotting of the processed species cover layers (`initPlot` event).

### Module outputs

The module produces the following outputs (Table
\@ref(tab:moduleOutputs-Biomass-speciesData)):

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>List of (ref:Biomass-speciesData) output objects and their description.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> objectName </th>
   <th style="text-align:left;"> desc </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> speciesLayers </td>
   <td style="text-align:left;"> biomass percentage raster layers by species in Canada species map </td>
  </tr>
  <tr>
   <td style="text-align:left;"> treed </td>
   <td style="text-align:left;"> Table with one logical column for each species, indicating whether there were non-zero cover values in each pixel. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> numTreed </td>
   <td style="text-align:left;"> a named vector with number of pixels with non-zero cover values for each species </td>
  </tr>
  <tr>
   <td style="text-align:left;"> nonZeroCover </td>
   <td style="text-align:left;"> A single value indicating how many pixels have non-zero cover </td>
  </tr>
</tbody>
</table>

and automatically saves the processed species cover layers in the output path
defined in `getPaths(sim)$outputPath`.

### Links to other modules

Intended to be used with other LandR data modules (e.g.,
*Biomass_borealDataPrep*) that require species cover data and the LandR forest
simulation *Biomass_core* module. You can see all *potential* module linkages
within the LandR ecosystem
[here](https://rpubs.com/PredictiveEcology/LandR_Module_Ecosystem). Select
*Biomass_speciesData* from the drop-down menu to see linkages.

### Getting help

-   <https://github.com/PredictiveEcology/Biomass_speciesData/issues>

## Module manual

### Detailed description

This module accesses and processes species percent cover (% cover) data for the
parametrisation and initialization of LandR *Biomass_core*. This module ensures
1) all data use the same geospatial geometries and 2) that these are correctly
re-projected to studyAreaLarge, and 3) attempts to sequentially fill-in and
replace the lowest quality data with higher quality data when several data
sources are used. It's primary output is a `RasterStack` of species % cover,
with each layer corresponding to a species.

Currently, the module can access the Canadian Forest Inventory forest attributes
kNN dataset [the default; @BeaudoinEtAl2017], the Common Attribute Schema for
Forest Resource Inventories [CASFRI; @Cosco2011] dataset, the Ontario Forest
Resource Inventory (ONFRI), a dataset specific to Alberta compiled by Paul
Pickell, and other Alberta forest inventory datasets. However, **only the NFI
kNN data are freely available** -- access to the other datasets must be granted
by module developers and data owners, and a Google account is required.
Nevertheless, the module is flexible enough that any user can use it to process
additional datasets, provided that an adequate R function is passed to the
module (see `types` parameter details in [Parameters])

When multiple data sources are used, the module will use replace lower quality
data with higher quality data following the order specified by the parameter
`types` (see [Parameters]).

The module can also exclude species % cover layers if they don't have a minimum
% cover value in at least one pixel. This means that the user should still
inspect in how many pixels the species is deemed present, as it is possible that
some data have only a few pixels with high % cover for a given species. In this
case, the user may choose to exclude these species *a posteriori*. The summary
plot automatically shown by *Biomass_speciesData* can help diagnose whether
certain species are present in very few pixels (see Fig.
\@ref(fig:fig-Biomass-speciesDataOutPlots)).

### Initialization, inputs and parameters

*Biomass_speciesData* initializes itself and prepares all inputs provided that
it has internet access to download the raw data layers (or that these layers
have been previously downloaded and stored in the folder specified by
`options("reproducible.destinationPath")`).

The module defaults to processing cover data fo all species listed in the
`Boreal` column of the default `sppEquiv` input `data.table` object, for which
there are available % cover layers in the kNN dataset (Table
\@ref(tab:defaultSppBiomass-speciesData); see `?LandR::sppEquivalencies_CA` for
more information):


```
## Error in set(x, j = name, value = value): RHS of assignment to existing column 'Species' is zero length but not NULL. If you intend to delete the column use NULL. Otherwise, the RHS must have length > 0; e.g., NA_integer_. If you are trying to change the column type to be an empty list column then, as with all column type changes, provide a full length RHS vector such as vector('list',nrow(DT)); i.e., 'plonk' in the new column.
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>List of species cover data downloaded by default by  (ref:Biomass-speciesData).</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> Species </th>
   <th style="text-align:left;"> Generic name </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Abies balsamea </td>
   <td style="text-align:left;"> Balsam Fir </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Abies lasiocarpa </td>
   <td style="text-align:left;"> Fir </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Acer negundo </td>
   <td style="text-align:left;"> Boxelder maple </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Acer pensylvanicum </td>
   <td style="text-align:left;"> Striped maple </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Acer saccharum </td>
   <td style="text-align:left;"> Sugar maple </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Acer spicatum </td>
   <td style="text-align:left;"> Mountain maple </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Acer spp. </td>
   <td style="text-align:left;"> Maple </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Alnus spp </td>
   <td style="text-align:left;"> Alder </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Betula alleghaniensis </td>
   <td style="text-align:left;"> Swamp birch </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Betula papyrifera </td>
   <td style="text-align:left;"> Paper birch </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Betula populifolia </td>
   <td style="text-align:left;"> Gray birch </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Betula spp. </td>
   <td style="text-align:left;"> Birch </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Fagus grandifolia </td>
   <td style="text-align:left;"> American beech </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Fraxinus americana </td>
   <td style="text-align:left;"> American ash </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Fraxinus nigra </td>
   <td style="text-align:left;"> Black ash </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Larix laricina </td>
   <td style="text-align:left;"> Tamarack </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Larix lyallii </td>
   <td style="text-align:left;"> Alpine larch </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Larix occidentalis </td>
   <td style="text-align:left;"> Western larch </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Larix spp. </td>
   <td style="text-align:left;"> Larch </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Picea engelmannii </td>
   <td style="text-align:left;"> Engelmann's spruce </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Picea glauca </td>
   <td style="text-align:left;"> White.Spruce </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Picea mariana </td>
   <td style="text-align:left;"> Black.Spruce </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Picea spp. </td>
   <td style="text-align:left;"> Spruce </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Pinus albicaulis </td>
   <td style="text-align:left;"> Whitebark pine </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Pinus banksiana </td>
   <td style="text-align:left;"> Jack pine </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Pinus contorta </td>
   <td style="text-align:left;"> Lodgepole pine </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Pinus monticola </td>
   <td style="text-align:left;"> Western white pine </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Pinus resinosa </td>
   <td style="text-align:left;"> Red pine </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Pinus spp. </td>
   <td style="text-align:left;"> Pine </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Populus balsamifera v. balsamifera </td>
   <td style="text-align:left;"> Balsam poplar </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Populus trichocarpa </td>
   <td style="text-align:left;"> Black cottonwood </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Populus grandidentata </td>
   <td style="text-align:left;"> White poplar </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Populus spp. </td>
   <td style="text-align:left;"> Poplar </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Populus tremuloides </td>
   <td style="text-align:left;"> Trembling poplar </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Tsuga canadensis </td>
   <td style="text-align:left;"> Eastern hemlock </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Tsuga spp. </td>
   <td style="text-align:left;"> Hemlock </td>
  </tr>
</tbody>
</table>

#### Input objects

*Biomass_speciesData* requires the following input data layers

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>List of (ref:Biomass-speciesData) input objects and their description.</caption>
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
   <td style="text-align:left;"> A named vector of colors to use for plotting. The names must be in sim$sppEquiv[[sim$sppEquivCol]], and should also contain a color for 'Mixed' </td>
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

Of the inputs in Table \@ref(tab:moduleInputs2-Biomass-speciesData), the
following are particularly important and deserve special attention:

-   `studyAreaLarge` -- the polygon defining the area for which species cover
    data area desired. It can be larger (but never smaller) that the study area
    used in the simulation of forest dynamics (i.e., `studyArea` object in
    *Biomass_core*).

-   `sppEquiv` -- a table of correspondences between different species naming
    conventions. This table is used across several LandR modules, including
    *Biomass_core*. It is particularly important here because it will determine
    whether and how species (and their cover layers) are merged, if this is
    desired by the user. For instance, if the user wishes to simulate a generic
    *Picea spp.* that includes, *Picea glauca*, *Picea mariana* and *Picea
    engelmannii*, they will need to provide these three species names in the
    data column (e.g., `KNN` if obtaining forest attribute kNN data layers from
    the Canadian Forest Inventory), but the same name (e.g., "Pice_Spp") in the
    coumn chosen for the naming convention used throughout the simulation (the
    `sppEquivCol` parameter); see Table
    \@ref(tab:mergingSpp-Biomass-speciesData) for an example)

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>Example of species merging for simulation. Here the user wants to model _Abies balsamea_, _A. lasiocarpa_ and _Pinus contorta_ as separate species, but all _Picea_ species as a generic _Picea spp._. For this, all six species are identified in the 'KNN' column, so that their % cover layers can be obtained, but in the 'Boreal' column (which defines the naming convention used in the simulation in this example) all _Picea_ species have the same name. (ref:Biomass-speciesData) will merge their % cover data into a single layer by summing their cover per pixel.</caption>
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
   <td style="text-align:left;"> *Abies balsamea* </td>
   <td style="text-align:left;"> Abie_Bal </td>
   <td style="text-align:left;"> Abie_Bal </td>
   <td style="text-align:left;"> *Abies balsamea* </td>
  </tr>
  <tr>
   <td style="text-align:left;"> *Abies lasiocarpa* </td>
   <td style="text-align:left;"> Abie_Las </td>
   <td style="text-align:left;"> Abie_Las </td>
   <td style="text-align:left;"> *Abies lasiocarpa* </td>
  </tr>
  <tr>
   <td style="text-align:left;"> *Picea engelmannii* </td>
   <td style="text-align:left;"> Pice_Eng </td>
   <td style="text-align:left;"> Pice_Spp </td>
   <td style="text-align:left;"> *Picea spp.* </td>
  </tr>
  <tr>
   <td style="text-align:left;"> *Picea glauca* </td>
   <td style="text-align:left;"> Pice_Gla </td>
   <td style="text-align:left;"> Pice_Spp </td>
   <td style="text-align:left;"> *Picea spp.* </td>
  </tr>
  <tr>
   <td style="text-align:left;"> *Picea mariana* </td>
   <td style="text-align:left;"> Pice_Mar </td>
   <td style="text-align:left;"> Pice_Spp </td>
   <td style="text-align:left;"> *Picea spp.* </td>
  </tr>
  <tr>
   <td style="text-align:left;"> *Pinus contorta* </td>
   <td style="text-align:left;"> Pinu_Con </td>
   <td style="text-align:left;"> Pinu_Con </td>
   <td style="text-align:left;"> *Pinus contorta* </td>
  </tr>
</tbody>
</table>

#### Parameters

Table \@ref(tab:moduleParams2-Biomass-speciesData) lists all parameters used in
*Biomass_speciesData* and their detailed information.

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>List of (ref:Biomass-speciesData) parameters and their description.</caption>
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
   <td style="text-align:left;"> The possible data sources. These must correspond to a function named paste0('prepSpeciesLayers_', types). Defaults to 'KNN' to get the Canadian Forestry Service, National Forest Inventory, kNN-derived species cover maps from year 'dataYear', using the `LandR::prepSpeciesLayers_KNN` function (see https://open.canada.ca/ data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990 for details on these data). Other currently available options are 'ONFRI', 'CASFRI', 'Pickell' and 'ForestInventory', which attempt to get proprietary data - the user must be granted access first. A custom function can be used to retrieve any data, just as long as it is accessible by the module (e.g., in the global environment) and is named as paste0('prepSpeciesLayers_', types). </td>
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
   <td style="text-align:left;"> .studyAreaName </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Human-readable name for the study area used. If NA, a hash of `studyAreaLarge` will be used. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .useCache </td>
   <td style="text-align:left;"> logical </td>
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
   <td style="text-align:left;"> Used in reading csv file with fread. Will be passed to data.table::setDTthreads. </td>
  </tr>
</tbody>
</table>

Of the parameters listed in Table \@ref(tab:moduleParams2-Biomass-speciesData),
the following are particularly important:

-   `coverThresh` -- integer. Defines a minimum % cover value (from 0-100) that
    the species must have in at least one pixel to be considered present in the
    study area, otherwise it is excluded from the final stack of species layers.
    Note that this will affect what species have data for an eventual simulation
    and the user will need to adjust simulation parameters (e.g., species in
    trait tables will need to match the species in the cover layers)
    accordingly.

-   `types` -- character. Which % cover data sources are to be used (see
    [Detailed description]). Several data sources can be passed, in which case
    the module will overlay the lower quality layers with higher quality ones
    **following the order of data sources specified by `types`** -- i.e., if
    `types == c("KNN", "CASFRI", "ForestInventory")`, *KNN* is assumed to be the
    lowest quality data set and *ForestInventory* the highest: values in *KNN*
    layers are replaced with overlapping values from *CASFRI* layers and values
    from *KNN* and *CASFRI* layers are replaced with overlapping values of
    *ForestInventory* layers.

### Simulation flow

The general flow of *Biomass_speciesData* processes is:

1.  Download (if necessary) of and spatial processing of species cover layers
    from the first data source listed in the `types` parameter. Spatial
    processing consists in sub-setting the data to the area defined by
    `studyAreaLarge` and ensuring that the spatial projection and resolution
    match those of `rasterToMatchLarge`. After spatial processing, species
    layers that have no pixels with values >= to the `coverThresh` parameter are
    excluded.

2.  If more than one data source is listed in `types`, the second set of species
    cover layers is downloaded and processed as above.

3.  The second set of layers is assumed to be the highest quality dataset and
    used to replaced overlapping pixel values on the first (including for
    species whose layers may have been initially excluded after applying the
    `coverThresh` filter).

4.  Steps 2 and 3 are repeated for remaining data sources listed in `types`.

5.  Final layers are saved to disk and plotted. A summary of number of pixels
    with forest cover are calculated (`treed`and `numTreed` output objects; see
    [Module outputs]).

## Usage example

### Load `SpaDES` and other packages.


```r
if (!require(Require)) {
  install.packages("Require")
  library(Require)
}

Require(c("PredictiveEcology/SpaDES.install",
          "SpaDES", "PredictiveEcology/SpaDES.core@development",
          "PredictiveEcology/LandR"), 
        install_githubArgs = list(dependencies = TRUE))
```

### Get module, necessary packages and set up folder directories


```r
tempDir <- tempdir()
spadesModulesDirectory <- file.path(tempDir, "modules")

paths <- list(inputPath = normPath(file.path(tempDir, "inputs")), 
              cachePath = normPath(file.path(tempDir, "cache")), 
              modulePath = normPath(spadesModulesDirectory), 
              outputPath = normPath(file.path(tempDir, "outputs")))

getModule("PredictiveEcology/Biomass_speciesData", modulePath = spadesModulesDirectory, overwrite = TRUE)

## make sure all necessary packages are installed:
makeSureAllPackagesInstalled(spadesModulesDirectory)
```

### Setup simulation

For this demonstration we are using all default parameter values, except
`coverThresh` , which is lowered to 5%. The species layers (the major output of
interest) are saved automatically, so there is no need to tell `spades` what to
save using the `outputs` argument (see `?SpaDES.core::outputs`).

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
                reproducible.inputPaths = paths$inputPath)

mySimOut <- simInitAndSpades(times = list(start = 1, end = 1),
                             modules = modules, 
                             parameters = params,
                             objects = objects, 
                             paths = paths,
                             .plotInitialTime = 1)
options(opts)
```

Here are some of the data retrieved by *Biomass_speciesData* for a randomly
generated study area within Canada.

<img src="figures/testRunFigure.png" title="(ref:Biomass-speciesData) automatically generates a plot of species dominance and number of presences in the study area, when `.plotInitialTime = 1` is passed as an argument." alt="(ref:Biomass-speciesData) automatically generates a plot of species dominance and number of presences in the study area, when `.plotInitialTime = 1` is passed as an argument." width="50%" />

## References
