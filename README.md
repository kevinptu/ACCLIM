# ACCLIM (v2)

Assimilation Capacity Coordinates with Light Interception Model

Satellite remote sensing of photosynthesis and evapotranspiration, developed by Kevin Tu.

## Code

- `gee/Acclim_v2_shapefile.js`: current Google Earth Engine application. Includes selectable NDVI²/calibrated EVI2 FAPAR, photosynthesis and evapotranspiration, carbon production, fitted extractable soil-water storage, temperature diagnostics, and optional Terra/Aqua LST comparisons.
- `r/gpp_acclimate_v2_GRD_K2.Rmd`: R Markdown implementation/reference, preserved as supplied. Its parameters and features are not fully synchronized with the GEE application; do not assume numerical equivalence.

## Running the Earth Engine app

1. Copy the JavaScript into the Google Earth Engine Code Editor using an Earth Engine-enabled account/project.
2. Review the constants, dataset access and asset folder settings. The preloaded geometry list currently references `projects/ee-landflux/assets/ecoexplorer`; replace this for another account or use drawn geometry/latitude and longitude.
3. Run the script, select a location, dates and FAPAR method, then click **Run ACCLIM**.
4. Optional chart series are selected from cached data. Satellite temperature observations load separately from model results.

Publishing or updating an Earth Engine App is a separate step; pushing to GitHub does not update a deployed app. Private assets and external uploader services are not bundled here.

## Running the R reference

Open the R Markdown file in RStudio, install its required packages, and update its input-data path. The script currently references a local FLUXNET CSV under `GPP_Project/`. Input datasets are not included; obtain them under their applicable access and use terms.

## Interpretation and development status

Soil moisture denotes effective extractable water. Reported effective depths derive from fitted storage and estimated available water capacity; they are model diagnostics. Satellite LST represents clear-sky overpass observations, whereas modeled temperatures use daytime energy-balance quantities.

The initial repository snapshot preserves the current local source files. JavaScript syntax was checked; a live Earth Engine run and end-to-end R execution were not performed for this snapshot. Existing source attribution and comments are preserved. No open-source license has been selected.
