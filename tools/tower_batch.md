# Run a multi-site tower comparison

From the repository root:

    python tools/tower_batch.py path/to/ACCLIM_five_site_years.zip --output comparison.csv

The Python launcher requires authenticated earthengine-api credentials and Node with @google/earthengine available. Use --project to select your Earth Engine project. Pass --params with a JSON object to test different ACCLIM parameters across the entire batch.

Each subdaily FLUXMET member is parsed separately by the same parser as the web importer. Its site-specific BIF metadata supplies the coordinates. The runner executes the public ACCLIM model with monthly TerraClimate forcing and MODIS 16-day reflectances. Each site has an independent model state. The combined CSV identifies each site and date, including when two sites share a calendar year.

Modeled and tower GPP columns are mean daytime rates in µmol CO₂ m⁻² s⁻¹. DT is GPP_DT_VUT_REF and NT is GPP_NT_CUT_REF; both use the importer's incoming shortwave > 0 filter. The GPP missing percentage is for DT among intervals with known positive shortwave. Tower ET uses the full-day mean of finite LE_CORR, converted to mm/day. Empty cells represent unavailable values.

To inspect a ZIP without running Earth Engine, use:

    node tools/tower_batch.cjs path/to/ACCLIM_five_site_years.zip

The hosted importer opens this same ZIP with a site-year selector. Multi-site execution uses the batch command above.

## Reuse sampled inputs

The first run saves JSON bundles beside the ZIP in <archive_stem>_inputs. Subsequent runs, including a new process on another day, automatically reuse the matching samples. Each bundle contains native sampled gridded weather, screened satellite observations (including padding for interpolation), and GE/VI diagnostics for the selected periods. Use --cache-dir to choose another directory, or --refresh-inputs to fetch fresh samples when source archives or Earth Engine products have been revised.

Parameter overrides recompute alignment, model forcing, fluxes and bucket fits. The cache identity includes sampling schema version, geometry, dates, weather source, satellite source and timestep. A different selection fetches a separate bundle. Optional SoilGrids and 10-year historical baseline data are requested separately. Model execution still requires Earth Engine.

The Earth Engine web app uses the same input bundles in memory for up to five distinct selections per app session. Leave **Reuse sampled weather and satellite inputs** checked for parameter experiments. Clear it for a fresh pull. Reloading the web app clears its in-memory bundles.
