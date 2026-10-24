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
