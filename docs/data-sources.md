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

Reflectance QA requires clear cloud state, no shadow/cirrus/internal cloud/snow/adjacent-cloud flags, highest-quality red/NIR bands, and valid nonnegative reflectance up to 1. LST requires good mandatory/data QA and reported LST error at most 1 K. Temperatures use the product's 0.02 K scale, converted to Celsius. Valid NDVI subpixels are averaged even with partial coverage; `NDVI_valid_fraction` records their area fraction. No temporal gap filling is performed. MYD11A2 QA describes the composite; it does not allow re-screening individual daily observations within its supplied mean.

Each plotted point is one valid LST pixel in one composite interval. The annual scatter is colored by month; the interval selector isolates an 8-day spatial scatter. Composite start dates label intervals; the reflectance-selected observation and LST mean need not represent the same day. The last interval is clipped at the calendar-year boundary for labeling. Partial years (2002 or current-year availability) report only supplied intervals.

Sampling is requested in quarterly batches. Box sides are bounded to 0.5–50 km. For more than 20,000 valid pairs the interface offers interval views rather than drawing the entire annual scatter; no pixel pairs are subsampled. Reducing the box enables the annual view. Model CSV export state is cleared when starting this diagnostic; this first version does not add a separate diagnostic export.

Validation: a 2014 Tonzi 10 km test returned 117 native grid cells and 4,448 valid pairs over 46 paired composite intervals. An independent synthetic Earth Engine check verified mean-of-NDVI (0.35, versus a different NDVI-of-mean-reflectance), partial coverage, and empty pre-Aqua intervals.
