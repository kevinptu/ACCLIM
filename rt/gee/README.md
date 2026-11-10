# ACCLIM-RT Earth Engine app

Paste the entire generated **ACCLIM_RT_GEE.js** into a new Earth Engine Code Editor script, save it, and run. Publish it using the Code Editor Apps dialog when ready. The script contains the optical constants, numerical fitter, data request, and UI. It requires no Python service or private model asset at runtime. A FeatureCollection supplied as a location must be accessible to the app's users.

The Drive reference is G:/My Drive/Model/ACCLIM/ACCLIM_RT_GEE_web_app_script.txt. The GitHub copy is maintained on the ACCLIM-RT branch. For edits, check the Drive copy for user changes before rebuilding or replacing it.

## Inputs and outputs

Specify coordinates, click the map, or select the center of a FeatureCollection asset. Enter an inclusive date range and a box side (default 10 km). The projected box supports latitudes -80 to 84 degrees.

The app uses Aqua MYD09A1, native 500 m MODIS pixels, and the same QA masks as rt/sample_modis.py. It excludes cloud, shadow, non-land, cirrus, internal cloud, snow, adjacent-cloud flags, poor red/NIR band QA, invalid reflectance, and invalid geometry. Actual DayOfYear values determine the fit dates. Each 8-day composite retains a selected acquisition for each pixel. Consequently, the app covers the acquisition dates retained in this product; a daily source would be needed to include every original clear-sky observation.

Each date with at least four valid pixels gets an independent broad-background fit. Dates with fewer pixels appear in the summary CSV with status insufficient_pixels. Request failures are reported in the summary and final status. The app continues to subsequent composites after a failed retrieval. Spatial coverage uses all native box pixels as its denominator, including pixels screened out or acquired on other days.

Outputs include:
- Median and mean nadir fcover, median LAI, joint red/NIR RMSE, common background red/NIR, pixel count, distinct geometry count, spatial coverage, and near-best background count.
- Cover ranges across backgrounds within the user-specified absolute RMSE tolerance (default .001 BRF). These are sensitivity ranges with other parameters fixed.
- Fcover and RMSE time series.
- A selectable per-date RMSE-versus-median-fcover scatterplot.
- Web-app CSV downloads for the date summary, selected-date pixels, and selected-date background/RMSE/fcover profile.

All background fits minimize mean squared red and NIR residuals with equal band weights. RMSE is sqrt(sum of both-band squared errors / (2 * number of pixels)). Red and NIR background values are constant across pixels within a date and may vary between dates.

## Exact optimization

The broad grid uses red .01-.40 and NIR from red to .50 at .01 increments (1,220 pairs). LAI is 0-8 at .1 increments. Fcover = 1-exp(-.5*LAI). The leaf properties, spherical Verhoef leaf angle distribution, hotspot .05, and spectral means match whole_pixel.py.

For each distinct measured solar/view/azimuth geometry, calculate the wavelength-resolved FourSAIL terms once per LAI, calculate band reflectance for 50 background levels, and reuse those values for all pixels sharing that geometry. Geometry values are used without rounding. Each pixel independently minimizes the same red/NIR residual as the Python implementation.

The geometry lookup table is discarded after its pixels have been scored. Composites are requested sequentially. Computation yields between geometries so the UI can update and cancel. Retrieved samples are reused during the session for an identical location/box/date selection. Changing that selection clears the sample cache. Source observations and completed date diagnostics remain in browser memory for the current selection, so memory use grows with box size and duration.

The soil-template subset is omitted. The broad grid and fitting objective match the previous analysis.

## Relation to ACCLIM FAPAR and FIPAR

The estimated quantity is projected nadir foliage cover under the specified LAI, leaf-angle, and clumping assumptions. The date and pixel exports make it available for research comparisons with ACCLIM FAPAR/FIPAR. Applying it to carbon or water calculations requires choosing the mapping and temporal/spatial alignment. This script provides the retrieval and diagnostics as a standalone research app.

## Build and verification

Use the same prosail 2.0.5 environment as the Python reference:

    python rt/gee/build.py
    node rt/gee/core.test.cjs

Edit core.js for numerical fitting and app.js for Earth Engine sampling/UI. build.py embeds PROSPECT-D leaf reflectance/transmittance at 1 nm and the spherical LIDF in the complete GEE script. Rebuild after changing either source. Leaf parameters are fixed at N=1.5, Cab=40, Car=8, Cbrown=0, Cw=.015, Cm=.009.

reference_brf.json contains independent Python FourSAIL values for six geometries, five LAIs, and four backgrounds. Tests include bare background, exact hotspot, near hotspot, nadir, and oblique geometries. They also check dates, duplicates, zero cover, and geometry reuse.

Local validation on 2026-11-10:

| Site/date | Pixels | Geometries | Median fcover | Joint RMSE | Local JavaScript seconds |
|---|---:|---:|---:|---:|---:|
| Tonzi 2012-05-01 | 451 | 105 | .451188364 | .011656387 | 1.98 |
| Tonzi 2012-08-30 | 453 | 112 | .139292024 | .030475659 | 2.04 |
| Jornada 2024-06-08 | 261 | 60 | .095162582 | .029434000 | 1.30 |
| Norunda 2024-06-02 | 322 | 78 | .572585068 | .011942522 | 1.82 |

These reproduce the saved Python broad-grid fits to floating-point precision. Maximum BRF difference in the independent numerical cases was 3.33e-16. Timing includes fitting and diagnostic construction in local Node, excludes Earth Engine retrieval and chart rendering, and will vary in browsers.

The live Earth Engine sampler was checked against all 451 May 1 Tonzi pixels and their geometry. CSV URL generation was also tested against Earth Engine. A controller smoke test used live Earth Engine with mock UI widgets to verify both composite requests, acquisition-date filtering, chart data construction, and repeat-run reuse of samples. Publishing the script as an app and checking the rendered UI remain deployment steps.

Equations and numerical reference:
- Verhoef et al. 2007, doi:10.1109/TGRS.2007.895844.
- https://github.com/jgomezdans/prosail (FourSAIL.py and PROSPECT-D).
- https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MYD09A1
