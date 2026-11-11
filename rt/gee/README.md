# ACCLIM-RT Earth Engine app

Paste the entire generated **ACCLIM_RT_GEE.js** into a new Earth Engine Code Editor script, save it, and run. Publish it using the Code Editor Apps dialog when ready. The script opens the full ACCLIM model. Choose **FAPAR = PROSAIL fcover**, **FIPAR = PROSAIL fcover**, or both in the existing menus. The retrieval runs for these selections and uses a configurable box centered on the model geometry (default 10 km). Selecting PROSAIL for FIPAR clears the shared-FAPAR override so the water selection takes effect. The PROSAIL retrieval diagnostics button opens the per-date plots and pixel exports. It requires no Python service or private model asset at runtime. A FeatureCollection supplied as a location must be accessible to the app's users.

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

For a selected PROSAIL fraction the model uses the date-specific spatial median of projected nadir foliage cover directly. A selected FAPAR equals cover even when fexp differs from 1. FIPAR feeds canopy energy and water partitioning, including the gridded net-radiation calculation. The model keeps its selected weather and satellite sources for its other inputs.

Daily cover uses linear interpolation between bracketing acquisitions up to 32 days apart. The retrieval includes 32 days on either side of the request for bracketing. Weekly periods take the mean with complete daily coverage. Monthly periods take the median of available daily estimates. rt_observed_days, rt_interpolated_days and rt_missing_pct describe coverage. Periods without sufficient coverage have unavailable physiology. Repeated experiments can reuse the same fitted dates through the existing input cache.


The FIPAR menu also offers **FIPAR = m_cover × NDVI + b_cover (PROSAIL fit)**. The model fits pixel fcover against NDVI from the same screened Aqua red/NIR observation within the selected box and dates. Ordinary least squares gives each acquisition date equal total weight. Observations outside the run dates serve only the direct-cover interpolation. A zero red-plus-NIR denominator excludes the pixel from calibration. The fit requires variation in NDVI and follows the same four-pixel minimum per acquisition as the RT retrieval.

The model applies the fitted slope and intercept to its selected satellite product's timestep NDVI and clips the result to [0, 1]. When a different sensor supplies the model NDVI, this transfers the Aqua calibration to that sensor. The model rows carry m_cover, b_cover, cover_fit_dates, and cover_fit_pixels.

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

## Green-leaf and chlorophyll absorption

Select **FAPAR = PROSAIL green-leaf absorption** or **FAPAR = PROSAIL chlorophyll absorption**. Each option uses the selected fraction directly in photosynthesis, seasonal FAPAR tracking and the carbon/water coupling. Both values are included in model rows and the model CSV. FIPAR follows its own selection unless shared FAPAR is enabled.

The absorption calculation uses the spatial mean over the fitted pixel LAI distribution. It weights each occupied LAI by its pixel count and uses the same fitted shared background and fixed leaf parameters as the red/NIR retrieval. Background reflectance is flat over PAR at the fitted red BRF. Chlorophyll content is fixed at 40 micrograms per square centimeter and carotenoids at 8. Seasonal changes arise through the fitted canopy state and illumination.

At build time absorption.py uses Python FourSAIL hemispherical fluxes at 1 nm spacing from 400 through 700 nm. At each wavelength, foliage absorption equals incoming irradiance minus upward hemispherical flux at canopy top minus background absorption. Downward flux at the background includes multiple reflections. Chlorophyll absorption is the foliage absorption multiplied by the chlorophyll share of PROSPECT-D optical depth at that wavelength. Carotenoids, water and dry matter contribute to total optical depth. Brown pigments and anthocyanins have zero concentration in the fixed leaf model. This constituent partition assumes homogeneous mixing within the leaf.

Direct and diffuse fractions use the respective PROSAIL library incident spectra converted from energy to photons by wavelength weighting, with trapezoidal endpoint weights. The embedded lookup spans the retrieval LAI grid, backgrounds 0–0.4 at 0.1 intervals, and solar zenith cosines from near the horizon to nadir. The app interpolates in background reflectance and zenith cosine. Solar zenith changes over the day. Sixteen midpoint intervals on either side of noon use incident PAR proportional to positive solar zenith cosine. Diffuse illumination is isotropic. The assumed fraction of incident PAR photons that is diffuse is adjustable (default 0.5); this is a sensitivity parameter.

Response curves are interpolated between valid observations up to 32 days apart before evaluating the day's solar illumination. Period fractions are ratios of total absorbed to incident photons. Daily shortwave irradiance supplies the relative photon weights with ACCLIM's constant shortwave-to-PAR conversion. When only period weather is available its mean shortwave is held constant across those days. Monthly rows use the available days; daily and weekly rows require complete coverage. A zero incident-photon denominator yields an unavailable fraction.

Outputs include FAPAR_leaf, FAPAR_chl, rt_leaf_direct, rt_chl_direct, rt_leaf_diffuse, rt_chl_diffuse, rt_diffuse_par_fraction, rt_abs_observed_days, rt_abs_interpolated_days, and rt_abs_missing_pct. The direct and diffuse columns are the two illumination limits before mixing.

The two chlorophyll-fit FAPAR choices regress modeled pixel chlorophyll absorption on matched Aqua NDVI or EVI2 over acquisition dates inside the run interval. Each date receives equal total weight. The targets use each pixel's fitted LAI, that date's shared background and solar geometry, and the selected diffuse PAR fraction. EVI2 uses 2.5 × (NIR − red)/(1 + NIR + 2.4 × red), clipped to [0, 1] as in the model forcing. The selected line predicts FAPAR from the model satellite's timestep index, clipped to [0, 1]. An index prediction may therefore be available between direct absorption estimates. The output includes the selected m_abs_ndvi/b_abs_ndvi or m_abs_evi2/b_abs_evi2, abs_fit_dates and abs_fit_pixels. FAPAR_chl remains the direct RT diagnostic.
