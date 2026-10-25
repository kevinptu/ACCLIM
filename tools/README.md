# FLUXNET Shuttle to ACCLIM

The [browser importer](https://kevinptu.github.io/ACCLIM/tower-csv-importer.html) accepts Quanterra CSVs and FLUXNET FLUXMET/FULLSET/SUBSET HH/HR CSVs or ZIP packages. Use the official Shuttle download website linked there, select the downloaded ZIP, review fields, then copy the prepared summary into ACCLIM. CSV/ZIP processing runs locally in a browser worker. There is no upload, server-side storage, or automatic download from the importer. “Save summary JSON” explicitly saves a local file.

## Calculation and missing data

- Use all finite GPP when the chosen incoming-shortwave variable is greater than zero. No quality or missing-percentage cutoff is applied. Keep negative and gap-filled GPP values.
- Monthly missing GPP percentage is `100 × (daylight records − valid GPP records) / daylight records`, separately for each GPP method. A month with no valid GPP has a blank mean, not zero; known daylight with all missing GPP reports 100%.
- Records with missing radiation cannot be classified as daytime. Report their percentage separately, with all supplied records as denominator. Entire absent timestamps are counted separately between file bounds. Neither is silently treated as night or included in the known-daylight denominator.
- Means pool subdaily observations; they are not unweighted averages of daily means. Daily sums/counts retain enough information to recalculate after an inclusive date override in ACCLIM.
- Support uniform 30- or 60-minute intervals. Reject duplicate/out-of-order timestamps, mixed interval lengths, and daily-only inputs. Compute local midpoint from interval start/end; preserve the source local clock without DST or UTC conversion.
- Recognize blanks, nonfinite values, -999, and FLUXNET -9999/-9999.0/-9999.9 variants. Report radiation and NEE QC distributions as supplied; NEE QC is not a direct GPP gap-filling flag. Consult the source package for flag meanings.
- Default to GPP_NT_VUT_REF with GPP_DT_VUT_REF as an additional estimate where available. Both are daytime averages here; NT/DT refer to partitioning methods. CUT fields can be explicitly selected. Do not silently substitute missing methods.
- Subdaily FLUXNET GPP fields use µmol CO₂ m⁻² s⁻¹. Daily/monthly aggregated products are not substitutes for daytime means.

The importer reads site coordinates from matching BIF/BADM rows or an optional Shuttle snapshot. Missing metadata can be entered manually. Citation/product ID may require the Shuttle snapshot or manual entry; the source license/readme remains in the downloaded archive. Retain the dataset citation, DOI, and required acknowledgments for downstream use. Archive members are read in memory, never extracted to filesystem paths. Large files require sufficient browser memory (600 MB compressed input; 1.2 GB selected CSV). The hosted page is the supported browser entry point; opening the HTML directly with `file://` can block its worker.

## Batch workflow (stage 2)

Install Node.js 18+ and Python 3.11+, then install the pinned Shuttle dependency in a Python virtual environment:

```powershell
python -m venv .venv
.venv\Scripts\python -m pip install -r tools/requirements-shuttle.txt
.venv\Scripts\python tools/shuttle_to_acclim.py --sites US-Var --hubs ameriflux --output-dir outputs/shuttle
```

For an already downloaded package, no Shuttle installation/network is needed:

```powershell
python tools/shuttle_to_acclim.py --input path/to/package.zip --snapshot path/to/shuttle_snapshot.csv --output-dir outputs/shuttle
```

Use `--other none` for just the primary GPP. Use `--gpp`, `--other`, and `--radiation` to select exact column names. Use `--archive-file` if a ZIP has multiple candidate subdaily CSVs. `--snapshot` reuses a catalog; `--hubs` restricts online discovery. Explicit site IDs are required; there is no download-all default. Duplicate catalog matches fail so the user can choose the intended product.

The batch command uses official Shuttle discovery/download functions and the **same JavaScript parser** as the browser. Unlike temporary browser processing, it intentionally saves downloaded archives, catalog snapshots and per-site JSON summaries to the output folder. Snapshot metadata can contain public site-team contact details; it is not copied into the ACCLIM summary. Read/copy the JSON into ACCLIM's Tower CSV data input, load settings and run. Source metadata, product filename and selected fields accompany the data; the summary does not include raw subdaily rows.

## Validation

`node tests/tower_import.cjs` tests daylight/missing handling, two GPP methods, ZIP selection, metadata, date filtering, backward compatibility and Quanterra input. Set `ACCLIM_PRIVATE_SCRIPT` to also validate a local private GEE script.

Validated with the downloaded AmeriFlux US-Var v1.3_r1 2000–2025 archive: 312 calendar months. An independent streaming Python calculation matched both GPP monthly means, valid counts and known-daylight counts for every month. October 2000: primary mean 4.575752309, 22.23903177% missing (514 valid / 661 daylight). January 2024: primary mean 4.478522672, 0% missing (684/684). The official Shuttle discovery → download → conversion path also completed for US-Var. Other hubs' packages remain to be tested.

Public/private Earth Engine scripts accept version 3 summaries and older version 2 imports. A coverage month selector beneath the GPP chart reports missing percentages, unknown radiation, and absent intervals. The GEE scripts still require separate republication; pushing this repository updates only the hosted importer.

ZIP decoding uses vendored fflate 0.8.2 under its accompanying MIT license.

Validation data citation: Siyan Ma, Liukang Xu, Joseph Verfaillie, Dennis Baldocchi (2026), AmeriFlux FLUXNET-1F US-Var Vaira Ranch- Ione, Ver. v1.3_r1 (2000–2025), AmeriFlux AMP. https://doi.org/10.17190/AMF/1993904. Funding for the AmeriFlux data service was provided by the U.S. Department of Energy Office of Science. FLUXNET data products were produced and harmonized by eddy covariance regional networks and data processing centers, including AmeriFlux, ChinaFlux, European Fluxes Database, ICOS, JapanFlux, KoFlux, OzFlux, SAEON, and TERN. These products also include a modified version of ERA5 hourly data provided by the Copernicus Climate Change Service.

### Reopen saved summaries and reduce paste size

Choose a saved version 3 ACCLIM summary `.json` using the importer Data file control. The importer restores the comparison and recomputes its monthly preview without the original CSV or ZIP. Save summary JSON keeps the full record and quality details. Copy ACCLIM data omits redundant monthly previews and QC histograms, while preserving daily sums and coverage counts. Use Copy from / Copy through to restrict the copied date range if a long record makes the Earth Engine text box slow. The preview and saved JSON still cover the full record. Existing GEE scripts accept this smaller version 3 payload; no app republication is needed for this change.

### Optional tower weather runs

Reprocess the source CSV/ZIP to add meteorology to the summary; old saved JSON still supports GPP comparisons but cannot supply weather. The importer offers weather-column selectors. Defaults prefer TA_F_MDS / VPD_F_MDS when present, then consolidated _F; pressure and precipitation prefer PA_F/P_F, then PA/P. The selected shortwave column also supplies weather. Quanterra names and their unit row are supported. No second gap-filling is performed. Consolidated FLUXNET _F products may contain reanalysis-derived filling; choose fields deliberately and inspect the retained monthly source-QC histograms. Those flags retain the source product's definitions, and are not a common measured/gap-filled classification.

In the GEE inputs, load the summary, use its coordinates, and enable **Use Tower Weather Data** before running. Unchecked is the original gridded-weather path. Run the two modes separately to compare them. Both public and private scripts require republication for this option.

The app retains its monthly cadence. Required weather inputs are air temperature, humidity (VPD or RH), incoming shortwave, pressure, precipitation and corrected latent/sensible heat fluxes. Full-day interval means produce monthly SWR, actual vapour pressure, pressure and available energy (LE_CORR + H_CORR); Tmin/Tmax are averages of complete daily minima/maxima; Tdaily is the full-day temperature mean. Rainfall is summed, not averaged or scaled to fill missing hours. FLUXNET VPD hPa converts to kPa; Quanterra VPD is already kPa. Actual vapour pressure is calculated at each interval using the model's Magnus constants and VPD or RH, then averaged. The model retains its Tmin/Tmax-derived daytime temperature and astronomical daylength, but uses tower actual vapour pressure to calculate VPD/RH, rather than the default Tmin-as-dewpoint approximation. The sum LE_CORR + H_CORR replaces the model radiation/available-energy input. Both fluxes must be present at the same interval; missing components are never treated as zero. This is turbulent available energy, rather than a direct net-radiation measurement. Weather summaries made with NETRAD must be regenerated from CSV/ZIP (weather schema version 2). CO2 remains the existing prescribed series; satellite vegetation inputs and parameter choices remain in place.

Strict missing-weather handling: every required field must have all expected 30/60-minute intervals in a calendar month. Partial boundary months or a missing required field make the model month unavailable, without any gridded substitution. This rule is separate from GPP, whose available daytime values are still averaged without a coverage cutoff. The preview reports missing percentages by weather variable. The selected copy date range includes full-calendar-month weather statistics, even for a partial-month model request.

Runs restart acclimation and soil-water states at each contiguous weather segment; no state or daily-carbon interpolation bridges missing months. The EVI2/FAPAR fit and MODIS lookup use the full requested period in both modes, so weather gaps do not change the vegetation calibration. Missing model months remain blank in charts and exports; incomplete annual totals remain unavailable. Annual temperature used in biomass allocation also comes from tower weather and requires a full year of complete temperature data. Outputs include weather source/status and selected source fields. GEE computation at a live site must be checked after republishing; automated tests cover aggregation, units, missing periods, summary transport and segmentation.

Weather units: https://fluxnet.org/data/aboutdata/data-variables/

### Tower LE_CORR overlay

The Evapotranspiration chart offers **Add Tower LE_CORR**, independently of the Use Tower Weather Data checkbox. Reprocess the source CSV/ZIP to include daily LE_CORR sums/counts; older JSON without those fields leaves this series unavailable. It uses all available 24-hour intervals, including night, negative fluxes and supplied gap-filled values. The selected-date mean in W/m² is multiplied by 86400 / 2450000 to give equivalent mm/day (the same latent heat constant as the model). No daylight multiplier is applied and no coverage cutoff is imposed on this observation series. Its missing percentage uses supplied intervals; absent timestamps and partial-month coverage remain separately reported. The complete-weather rule for model forcing is unchanged. The monthly export includes TowerLE_CORR_mm_day, valid count and missing percentage.
