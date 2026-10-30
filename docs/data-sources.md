# Weather, reflectance, and temporal processing

The weather selector controls model cadence: TerraClimate runs monthly; ERA5-Land and GRIDMET run daily. Satellite cadence does not determine the model timestep. Charts retain the native model timestep, including GPP, NPP, water fluxes, soil moisture, physiology and biomass production. Tower observations are averaged from their daily record statistics over each model interval. The monthly and 16-day CSVs remain explicit reporting products; they do not set chart cadence.

## Weather variables

| Source | Precipitation | Shortwave radiation (24-hour mean) | Temperature | Actual vapor pressure |
|---|---|---|---|---|
| TerraClimate | `pr`, mm/month | `srad × 0.1`, W/m² | `tmmn`, `tmmx × 0.1`, °C; mean of extrema | `vap × 0.001`, kPa |
| ERA5-Land DAILY_AGGR | `total_precipitation_sum × 1000`, mm/day | `surface_solar_radiation_downwards_sum / 86400`, W/m² | `temperature_2m_min`, `temperature_2m_max`, `temperature_2m`, K minus 273.15 | Saturation vapor pressure at `dewpoint_temperature_2m` |
| GRIDMET | `pr`, mm/day | `srad`, W/m² | `tmmn`, `tmmx`, K minus 273.15; mean of extrema | `es(Tmin) × rmax / 100` |

ERA5-Land uses `surface_pressure / 1000` in kPa. TerraClimate and GRIDMET retain the existing pressure constant. Weather is sampled at its catalog resolution: approximately 4.6 km, 11.1 km, and 4.6 km respectively. Selecting a 30 m satellite source does not downscale weather.

All sources now supply `ea` to the physiological humidity calculations. RHmin = ea/es(Tmax), RHmax = ea/es(Tmin), VPDmin = es(Tmin)−ea, and VPDmax = es(Tmax)−ea. RH is bounded to [0,1]; negative VPD is bounded to zero. The previous gridded shortcut ea=es(Tmin) is no longer used by the model core.

Supplied product values, including supplied gap-filled values, are retained. Negative precipitation values are treated as unavailable rather than turned into valid zero rainfall. Original observation-gap percentages are unknown for gridded products that do not expose them; `weather_source_missing_pct` remains blank rather than implying zero original missing observations.

## Reflectance sources

* Landsat 8/32-day C02 T1 L2 composites: already scaled `red` and `nir`, with the catalog's compositing masks. These are latest-valid-pixel composites, not temporal means.
* Combined HLS L30/S30: `B4` red, `B5` (L30) or `B8A` (S30) NIR. Earth Engine supplies scaled reflectance. Exclude cloud, adjacent cloud/shadow, shadow, snow and high aerosol using Fmask.
* Daily MOD09GQ/MYD09GQ: red/NIR × 0.0001; highest-quality paired bands from `QC_250m`; cloud-state QA joined by date from MOD09GA/MYD09GA `state_1km`.
* Eight-day MOD09Q1/MYD09Q1: red/NIR × 0.0001; paired-band quality from `QA`, cloud state from `State`.
* Sixteen-day MOD13Q1/MYD13Q1: red/NIR × 0.0001; good SummaryQA and no adjacent/mixed cloud, snow, or shadow in DetailedQA. Both products are 250 m; MYD13A1 is not used.

MOD09 cloud-state screening accepts only state 0, explicitly excluding state 3 (not set/assumed clear). Shadow, cirrus, internal cloud, adjacent cloud and snow flags are also excluded. Daily cloud-state QA has 1 km resolution even though reflectance is 250 m; categorical QA uses nearest-neighbor resampling.

MOD13 pixels are assigned their `DayOfYear` observation date before spatial reduction. Other single-day products use the observation date. Products without selected-pixel dates use composite midpoints, including the shortened final MODIS composite of a year. Landsat retains its documented overlapping year-end window. Same-day samples are averaged. For polygons, this is a mean of valid area samples; no claim of complete spatial coverage is made.

Paired red and NIR are interpolated linearly between valid dates. NDVI and EVI2 are then derived, and the existing albedo proxy `(red + nir)/2` is retained. Maximum bracketing spans are configurable in `SATELLITE_SOURCES`: 32 days for daily/HLS/8-day, 48 days for 16-day, 64 days for 32-day. No one-sided extrapolation is used. Observations outside the requested window are retrieved for bracketing. These are research settings, not guarantees of actual observation frequency.

Monthly TerraClimate inputs use means of the aligned daily satellite series. A month requires supported vegetation values for all its days. Thermal products and thermal comparison methods are unchanged.

## Bucket-only gap filling

Bucket precipitation is normalized to mm/second before interpolation, then converted back to a period total. Potential evaporative demand is interpolated in mm/second; surface-water evaporation and supporting bucket fields are interpolated separately. Existing supplied values always take precedence. Only two-sided interpolation is used; unresolved leading/trailing gaps remain unavailable.

The bucket continues through interpolated periods without resetting storage. Interpolated periods update storage but are excluded from the capacity-fitting objective. Capacity fitting retains the duration-weighted objective (seconds) and full initial bucket, with at least 365 contiguous represented days. The separate ET-effective bucket remains a diagnostic fitted to the original humidity-constrained ET reference.

Interpolated bucket inputs do not create GPP, NPP, temperatures, or other physiological model outputs for missing forcing periods. Acclimation advances by elapsed time between usable physiological periods. Tower weather remains monthly and never substitutes gridded weather; valid supplied precipitation can still support a bucket even when another tower weather field is unavailable.

## Flags and exports

* `satellite_status`: observed, interpolated, or unavailable. For monthly rows, interpolated means at least one day used interpolation.
* `satellite_missing_pct`: calendar days without a direct valid satellite observation, including days supported by interpolation. Regular composite spacing contributes to this percentage; it is not exclusively a cloud percentage.
* `satellite_unavailable_pct`: days unsupported even after permitted interpolation.
* `satellite_left_date`, `satellite_right_date`, `satellite_gap_days`, `satellite_date_basis`: interpolation and temporal provenance.
* `weather_missing_pct`: percentage of the seven canonical weather fields unavailable at the model timestep, including precipitation.
* `precipitation_missing_pct`: whether supplied precipitation was unavailable before bucket interpolation.
* `bucket_forcing_missing_pct`: 100 for a period with any missing bucket input before interpolation, otherwise 0.
* `bucket_record_missing_pct`: period-duration-weighted percentage of the complete requested record with missing bucket forcing before interpolation.
* `bucket_forcing_status`: supplied, interpolated, or unavailable. Supplied includes the source product's own gap-filled values.

Native model CSV exports all native model periods and their diagnostics. Daily CSV combines native daily model diagnostics and carbon products for daily weather; for monthly weather it allocates native monthly carbon totals uniformly over their represented days, labeled `allocated_from_weather_period`. This preserves monthly mass; it does not reconstruct daily physiology or interpolate between months. Monthly model charts/CSV use daylight-duration-weighted daytime flux means for daily runs, and sum precipitation over represented days. Monthly and annual carbon CSVs include complete totals, separately labeled available-period sums, requested days, missing percentages, and calendar coverage percentages. Missing days never contribute fabricated zero carbon.

Daily GPP mass = mean daytime GPP × daylength × 3600 × 12×10⁻⁶, in gC/m²/day. Monthly and annual carbon totals sum daily masses; NPP retains the configured CUE and cover allocation. No daily physiology is reconstructed from monthly values when a daily weather source is selected.

The interactive app samples satellite data in bounded batches and materializes weather and physiological output once. It performs the same analytic bucket equation and capacity search numerically in the browser, then adds the remaining diagnostics in Earth Engine. This avoids the server memory limit reached by expanding multi-year capacity searches into an Earth Engine computation graph. The numerical implementation is checked against the server implementation, including a 365-day record with a two-day forcing gap. These intermediate values are held for the current run; this does not create an Earth Engine asset or permanent data store.

## Catalog references

* [TerraClimate](https://developers.google.com/earth-engine/datasets/catalog/IDAHO_EPSCOR_TERRACLIMATE)
* [ERA5-Land](https://developers.google.com/earth-engine/datasets/catalog/ECMWF_ERA5_LAND_DAILY_AGGR)
* [GRIDMET](https://developers.google.com/earth-engine/datasets/catalog/IDAHO_EPSCOR_GRIDMET)
* [Landsat 8-day](https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_COMPOSITES_C02_T1_L2_8DAY), [32-day](https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_COMPOSITES_C02_T1_L2_32DAY)
* [HLS L30](https://developers.google.com/earth-engine/datasets/catalog/NASA_HLS_HLSL30_v002), [HLS S30](https://developers.google.com/earth-engine/datasets/catalog/NASA_HLS_HLSS30_v002)
* [MOD09GQ](https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MOD09GQ), [MOD09GA QA](https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MOD09GA)
* [MOD09Q1](https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MOD09Q1)
* [MOD13Q1](https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MOD13Q1), [MYD13Q1](https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MYD13Q1)

## Shared equations and timestep duration

The weather adapter defines each model interval. All sources use the same physiological and bucket equations; TerraClimate supplies monthly intervals and ERA5-Land/GRIDMET supply daily intervals. Satellite reflectances are aligned using the interpolation rules above.

Each row carries `period_seconds`, `daylight_seconds`, interval start/end timestamps, and `elapsed_seconds`. Daylight seconds sum astronomical daylight across the interval; monthly daylength is the mean across its days. Daytime radiation and PAR are converted using `period_seconds / daylight_seconds`.

* Acclimation weight: `exp(-elapsed_seconds / (tc * 86400))`, with `tc_sm` used for the soil-humidity response. Both time constants remain in days.
* Carbon mass: daytime GPP in micromol CO2/m2/s times `daylight_seconds * 12e-6`, yielding `GPP_period_gC_m2`; NPP uses the existing CUE.
* Water mass: daytime LE in W/m2 times `daylight_seconds / 2450000`, yielding mm/period. Equivalent mm/day values divide by `period_seconds / 86400`.
* Bucket forcing rates are in mm/second and integrated over interval seconds. Fitting weights and missing-data percentages use interval duration.

Initialization uses a fixed five-day window rather than five records. Daily inputs retain the original five-day convention; monthly inputs use the first monthly value over that window. Low-light periods hold the physiological state without holding the elapsed-time clock.

Daily and monthly predictions need not be identical: evaluating nonlinear equations with monthly mean weather differs from aggregating daily evaluations. The shared duration contract removes inconsistent unit conversions, not this aggregation effect. Subdaily weather ingestion is not yet offered in the interface.

## Selectable tower-weather timesteps

The importer retains full-day daily meteorology from the supplied half-hourly/hourly CSV. After preparing a comparison, choose **Tower weather timestep**: daily (default for newly processed data), weekly, or monthly. The copy payload passes that cadence to ACCLIM. **Use Tower Weather Data** remains optional; when enabled, its chosen cadence overrides the gridded-weather selector.

Daily temperature extrema come from each day's supplied records. Weekly/monthly Tmin and Tmax are means of daily extrema; mean temperature, shortwave, actual vapor pressure, pressure and LE_CORR + H_CORR are full-day means. Precipitation is summed. Weekly blocks start January 1, span seven days, and end at December 31 with a shorter final block. Boundary weeks/months retain their full calendar period even when the selected comparison begins inside it; carbon reporting is restricted to the requested dates.

Every required weather variable must have all expected source intervals within a model period. Missing or absent days make the affected period unavailable. Supplied gap-filled fields remain preferred; no gridded substitute is introduced. Existing bucket-only bridging still applies, without inventing missing carbon predictions.

Saved JSON retains daily weather and monthly quality/context summaries, regardless of the selected transfer cadence. To use daily or weekly weather with an older monthly-only JSON, process the original CSV/ZIP again through the standard importer. Changing the selector thereafter does not require parsing the CSV again.

The daily copy payload is larger and the model performs more updates. In a synthetic one-year fixture, complete transfers were approximately 148 KB daily, 70 KB weekly and 59 KB monthly (including tower GPP statistics). A one-month Earth Engine smoke test took approximately 7 seconds daily and 4 seconds weekly/monthly; these are illustrative, not performance guarantees.

NPP charts show mean daily carbon flux at each native weather period. Biomass production charts similarly show mean lb dry matter/acre/day; cumulative aboveground production is the end-of-period accumulation. Their source daily allocations, annual totals, and explicit 16-day CSV remain available. Monthly model exports from weekly runs allocate interval totals by calendar-day overlap rather than rerunning monthly physiology.

## Standalone LST versus NDVI

Select **Plot MODIS LST vs NDVI (no model run)**, choose a calendar year and square box side (default 10 km), then press **Plot LST vs NDVI**. Existing coordinate, map-drawn geometry and asset controls select the location. Polygon/asset geometries use their centroid. The box is constructed in a local UTM metric grid (polar stereographic outside UTM coverage); pixel centers inside the box determine the number of sampled LST cells. This mode bypasses all ACCLIM forcing, physiology, bucket, and carbon calculations and ignores model weather/date selections.

The analysis pairs [MYD09A1 500 m reflectance](https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MYD09A1) and the user-selected [MYD11A2 1 km daytime LST](https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MYD11A2) on identical 8-day composite start timestamps. MYD11A1 is daily and is not used here. NDVI is calculated from red and NIR at 500 m, then valid NDVI values are area-averaged onto the native LST projection. It is not computed from averaged reflectances. The native nominal 1 km grid is approximately 926.6 m, so a 10 km square can contain more than 100 pixel centers.

Reflectance QA requires clear cloud state, no shadow/cirrus/internal cloud/snow/adjacent-cloud flags, highest-quality red/NIR bands, and valid nonnegative reflectance up to 1. The LST quality selector defaults to valid retrievals with reported LST error at most 2 K: mandatory QA 0 or 1, data-quality bits 0, and LST-error bits 0 or 1. The original strict option remains available (mandatory QA 0, data-quality bits 0, LST-error bits 0). Both reject missing/cloud retrieval flags and nonzero data-quality flags. The default accepts the product’s “other quality” retrieval class, so it trades stricter screening for coverage; it does not imply the same certainty as the strict subset. Raw `LST_QC` and the selected policy are retained with each sampled row. Temperatures use the product's 0.02 K scale, converted to Celsius. Valid NDVI subpixels are averaged even with partial coverage; `NDVI_valid_fraction` records their area fraction. No temporal gap filling is performed. MYD11A2 QA describes the composite; it does not allow re-screening individual daily observations within its supplied mean.

Each plotted point is one valid LST pixel in one composite interval. The annual scatter is colored by month; the interval selector isolates an 8-day spatial scatter. Composite start dates label intervals; the reflectance-selected observation and LST mean need not represent the same day. The last interval is clipped at the calendar-year boundary for labeling. Partial years (2002 or current-year availability) report only supplied intervals.

Sampling is requested in quarterly batches. Box sides are bounded to 0.5–50 km. For more than 20,000 valid pairs the interface offers interval views rather than drawing the entire annual scatter; no pixel pairs are subsampled. Reducing the box enables the annual view. Model CSV export state is cleared when starting this diagnostic; this first version does not add a separate diagnostic export.

Validation: a 2014 Tonzi 10 km test returned 117 native grid cells and 4,448 valid pairs under strict QA (4,780 under the default ≤2 K policy) over 46 paired composite intervals. At Jornada (32.5849, -106.6032), a 2014 10 km test returned zero pairs under strict QA and 5,294 under the default policy. Common QA codes 17 and 81 mark produced LST as “other quality,” with good data-quality bits and reported error ≤1 K and ≤2 K respectively. The strict mandatory-QA filter excluded these retrievals. Empty-result messages now distinguish absent product intervals, no grid-cell centers, and rejection by the paired quality masks. An independent synthetic Earth Engine check verified mean-of-NDVI (0.35, versus a different NDVI-of-mean-reflectance), partial coverage, and empty pre-Aqua intervals.

The standalone diagnostic includes NIR versus red surface reflectance in place of the former histogram. Red is the x-axis and NIR the y-axis, in reflectance fractions. Both bands are area-averaged over the same quality-screened 500 m subpixels within each native LST cell, retaining the same valid LST–NDVI pixel–interval pairs. The plot follows the annual/8-day scatter selector and uses matching month colors. The two full-year time series still show NDVI and daytime LST means and spatial 5th/95th percentiles per composite (linear-interpolated percentiles at `(n-1)*p`). Empty intervals remain gaps. Contributing valid pixels can change over time.

Monthly LST–NDVI lines use unweighted ordinary least squares over all valid pixel–interval pairs in each calendar month. The displayed domain extends from the largest observed NDVI toward zero, stopping at LST = 60 °C if that boundary is encountered first. Positive slopes are also clipped wherever they cross 60 °C. Months need at least two distinct NDVI values. Each line matches its month’s point color; individual interval views retain the corresponding pooled monthly fit. These fits are descriptive; no independence assumption or confidence interval is imposed on spatially/temporally repeated samples.

## Standalone soil-reference search

### Current default: pooled geometric edge and PVI-intercept comparison

The method selector now defaults to **Pooled edge + triangle geometry + PVI intercept**. The previous month-balanced quantile method remains available for comparison. Both use the same spatial search controls described below. Neither changes model parameters.

The geometric method pools all screened observations across dates, without an EVI2 cutoff, equal-month weighting, or a requirement for each date to reproduce the line. It retains the central 98% of observed red values and generates up to 27 candidate edges from 12/18/24 equal-width red bins, 1%/2%/5% NIR quantiles, and full/darker-two-thirds/brighter-two-thirds sets of boundary anchors. Bins need 20 observations and a fit needs four anchors. A median-pairwise-slope fit with median intercept is iterated twice using detrended within-bin quantiles. This gives comparable influence to brightness intervals rather than densely populated vegetation clouds. Candidate slope must be positive and ≤20.

Each candidate is evaluated against the whole pooled cloud, not only its originating anchors:

- Near-edge tolerance is `max(0.002, tolerance/2)` NIR reflectance. At least four populated bins need three near-edge observations each, at least 20 near-edge pairs must exist overall, and their 5th–95th percentile red span must cover at least half of the pooled red range.
- Where metadata exists, near-edge observations must span at least two composite dates and three native analysis-grid pixels. Different dates may support different portions of the edge. These are coverage checks, not independent sample-size claims.
- The mean fraction below the line by more than the edge tolerance, balanced across populated red bins, must not exceed 10%. This discourages lines cutting through the cloud.
- The fitted edge evaluated at supported dark/bright red endpoints supplies two candidate reference regions. The median red/NIR coordinates of observations in the top 5% of above-line NIR residuals supply a candidate vegetation region. A barycentric triangle test measures enclosure, with 5% coordinate tolerance. Vegetation height must exceed 0.02 NIR reflectance and bin-balanced enclosure must exceed 35%. These are exploratory defaults, not fixed universal tasseled-cap proportions.
- Nonphysical fitted reference reflectance is rejected. A geometrically plausible edge still does not establish exposed-soil identity. A vegetation–shadow cloud can pass some of these tests; missing soil remains a possible interpretation.

Candidate geometric score is `3 × below_fraction + max(0, 1 − supported_span) + (1 − enclosure) + max(0, 4 − supported_bins)`. Only candidates passing the above support checks enter the plausible-geometry comparison. If none passes, the best geometric exploratory fit is displayed rather than reporting that soil does not exist. Numerical/insufficient-data cases still return no fit.

### Regression-based EVI2min and NDVImin

For each plausible edge, draw **8,192 independent uniform red/NIR pairs**, using the same fixed seed (173) and initial draws for every geometry. Four visible controls set red and NIR minima/maxima (default 0–1 for each). Retain nonnegative PVI and positive red + NIR. At least 100 retained points and nonzero PVI variance are required. These are **simulated reflectances**, not observed satellite points or a claim that all combinations occur in nature.

Calculate `PVI = (NIR − a Red − b) / sqrt(1 + a²)`. Ordinary least-squares intercepts of EVI2 versus PVI and NDVI versus PVI are the new primary EVI2min and NDVImin estimates. Report slopes, R², retained draw count and bounds. Also refit using the lower 50% and 75% of simulated PVI values to show intercept sensitivity to PVI fitting range. Compare the selected EVI2 intercept with a simulation using the observed central red/NIR ranges. Sampling-bound and PVI-range dependence are warnings, not automatic rejection criteria.

Plausible candidates are ranked by geometric score plus `(absolute deviation from the plausible-candidate median EVI2 intercept + maximum EVI2 intercept shift across PVI ranges) / tolerance`. This prioritizes geometric support and intercept stability, never a predetermined intercept value. The research tolerance defaults to 0.01 and is editable. All candidate geometries remain inspectable. Broad intercept variation is reported explicitly rather than hidden by the selected candidate. The reported range across plausible geometries is a sensitivity range, not a confidence interval. Missing date/pixel metadata in imported diagnostic data is also reported.

The previous exact-line medians remain comparison diagnostics only in this mode. They are not the new primary estimates. Do not compare new intercepts with old exact medians as though they were the same estimator.

### Screening and outputs

Both search methods now retain only MYD09A1 `StateQA` land/water bits 3–5 equal to 1 (land). Water, ephemeral water and mixed shorelines are excluded **before** averaging native 500 m red/NIR reflectance. Existing cloud, cloud-shadow, cirrus, snow and band-quality masks remain active. This extra mask applies to the soil-reference search, not the ACCLIM model or existing LST–NDVI diagnostic. Classification can miss subpixel water; the QA mask is not perfect water identification. See the [MYD09A1 catalog bit definitions](https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MYD09A1).

The results panel provides box comparisons and a geometry selector showing candidate EVI2 intercepts. For each geometry it plots the pooled red/NIR cloud, near-edge support, candidate dark/bright/vegetation regions and fitted edge. It lists the contributing composite dates and pixel count. Separate simulated EVI2–PVI and NDVI–PVI plots show regression lines. Geometry, intercepts and warnings remain in `UI._soilSearch`; no model-baseline transfer is automatic.

Validation: unit tests cover sparse edges assembled from different dates, synthetic line recovery, one-date/one-pixel failures, missing metadata, deterministic simulations, OLS checks, degenerate ranges, controls and all three diagnostic charts for public/private scripts. Earth Engine synthetic tests verified that a half-water analysis pixel contains only land reflectance after masking, with valid fraction 0.5. Fresh 2014 10 km tests returned 4,752 Tonzi and 5,281 Jornada observations. Their selected edges had support on 14/30 dates and 69/106 pixels respectively. Under the default 0–1 simulation bounds, plausible-geometry EVI2 intercept ranges were approximately 0.021–0.179 at Tonzi and 0.090–0.094 at Jornada. Tonzi is reported as sensitive, not as a uniquely identified soil baseline. These checks verify behavior, not field validity.

### Spatial controls and previous-method reference

The paragraphs below describe the spatial controls shared by both methods and the retained **previous month-balanced** fitting/endpoint method. Its EVI2 ceiling, equal-month weighting and exact-line medians do not apply to the new geometric default.

Choose **Search for soil reference (no model run)** under Analysis. The original target uses the existing coordinate/map/asset controls. Choose a year, reference-box side (2–10 km; default 10), and extent: one box, 3 × 3 (default), or 5 × 5. Boxes are evaluated center-first, then outward in non-overlapping rings in a local metric projection. The default extent is 30 × 30 km; the maximum is 50 × 50 km. Each box is fitted separately rather than pooling different backgrounds. The full requested extent is evaluated, even if the center passes, to allow comparison among references.

Optional reference latitude/longitude shifts the search center without moving the original target or changing model coordinates. Leave both blank to expand around the target. Candidate distances are geodesic center-to-original-target distances. Original target and candidate boxes are added to the map; supported candidates are green and unsupported boxes amber. New runs/reset remove those diagnostic layers.

The search uses quality-screened Aqua MYD09A1 reflectance only, with the existing clear-cloud, shadow, cirrus, snow and band-quality checks. It does not require LST or use the LST-quality selector. Valid 500 m red/NIR reflectances are averaged on the product's native projection scaled by two (nominal 1 km); partial coverage is retained and recorded. This differs from the LST-paired diagnostic's sampling/selection. Sampling is sequential in quarterly batches per box. Large searches take longer; no model forcing, GPP, ET, or carbon calculations run.

### Candidate fitting and screening

Defaults are exploratory research settings, not validated universal soil-identification rules:

- Maximum observation EVI2 for fitting: **0.25**, adjustable from 0.05 to 0.8. EVI2 is computed from each averaged red/NIR pair.
- Fit the **5th conditional percentile** of NIR against red, balancing total weight among eligible months. The JavaScript implementation profiles the weighted-quantile intercept and minimizes convex pinball loss. Numerical slopes at the search bracket boundaries (±20) are rejected.
- Trim the outer 1% of selected red values on each side. Require at least **100 fitting observations**, **3 months with at least 30 observations each**, and a red span of **0.04** (adjustable). Report the actual supported span after excluding insufficiently sampled months.
- Default maximum NIR departure is **0.01 reflectance** (adjustable). Compare the 2%/10% quantile lines and each leave-one-month-out line to the selected fit at both ends of its red range.
- Six equal-width red bins test boundary shape. At least four must contain 20 observations; the largest absolute binned 5th-percentile residual must meet the same tolerance.
- Repeat fitting with the EVI2 ceiling decreased/increased by 0.05. Insufficient support or an endpoint shift beyond tolerance fails this check. This guards against identifying the EVI2 selection boundary itself as soil.
- Reject nonpositive slopes and implausible fitted reflectance/EVI2. These screens establish numerical support only; exposed-soil identity is **not** automatically established.

### EVI2min and transfer

For `NIR = a R + b`, signed `PVI = (NIR − a R − b) / sqrt(1 + a²)`. At PVI = 0:

`EVI2soil(R) = 2.5 ((a − 1) R + b) / ((a + 2.4) R + b + 1)`.

The analogous bare-soil NDVI is `NDVIsoil(R) = ((a − 1) R + b) / ((a + 1) R + b)`. Both indices are evaluated on the same 1,001-point uniform red grid. NDVImin is the median of those NDVI values; its minimum/maximum show brightness dependence. A nonpositive NIR + Red denominator makes NDVI unavailable and prevents support of that candidate. NDVImin and its range appear alongside EVI2min in candidate details, the box comparison and supported-candidate summary. The brightness plot displays both indices. No model NDVI scaling parameter is changed.

Evaluate this exact expression at 1,001 uniformly spaced red values over the supported fit range. Report the median as candidate EVI2min and the minimum/maximum as brightness dependence. This avoids deriving a linear intercept from arbitrary random red/NIR draws; it is not a confidence interval or an estimate weighted by the target's soil brightness distribution.

Results list every box with distance, candidate EVI2min where estimable, and support/failure reasons. A selector displays all QA-valid pairs (including observations not used in fitting), the fitted boundary, and EVI2 along PVI = 0. The summary compares the supported candidates' estimates and identifies the nearest candidate for initial viewing. Even an unsupported fitted line can be inspected, but it is explicitly marked unsupported. The summary never transfers a rejected fit.

**No model baseline is changed automatically.** The user must verify comparable soil/substrate, litter, terrain and exposure before transferring a reference to the target. Multiple supported lines may still represent different soils; their median is a descriptive comparison, not an automatic recommended calibration. If none passes, the tool reports no supported reference. The diagnostic state is retained in `UI._soilSearch` until the next run/reset; this version does not add an export or an apply-to-model button.

Validation: the Jornada chart CSV reproduces slope 1.162162/intercept 0.032946 and passes default checks, with exact-line median EVI2min 0.093225. The Sierra conifer chart CSV fails default support screening. A fresh Earth Engine 2014 Jornada 10 km reflectance-only sample returned 5,281 pairs and a supported candidate near EVI2min 0.0945. Earth Engine checks also verified a 25-box shifted layout, target-distance preservation, a populated adjacent box, and an empty pre-Aqua period. Unit checks cover the public/private scripts, synthetic linear recovery, positive contamination, curved/insufficient samples, controls, asynchronous errors/cancellation and model bypass. The tests establish implementation behavior, not field validation of soil transfer.
