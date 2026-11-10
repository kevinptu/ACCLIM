# Run a multi-site tower comparison

From the repository root:

    python tools/tower_batch.py path/to/ACCLIM_five_site_years.zip --output comparison.csv

The Python launcher requires authenticated earthengine-api credentials and Node with @google/earthengine available. Use --project to select your Earth Engine project. Pass --params with a JSON object to test different ACCLIM parameters across the entire batch.

Each subdaily FLUXMET member is parsed separately by the same parser as the web importer. Its site-specific BIF metadata supplies the coordinates. The runner executes the public ACCLIM model with monthly TerraClimate forcing and MODIS 16-day reflectances. Each site has an independent model state. The combined CSV identifies each site and date, including when two sites share a calendar year.

Modeled and tower GPP columns are mean daytime rates in µmol CO₂ m⁻² s⁻¹. DT is GPP_DT_VUT_REF and NT is GPP_NT_CUT_REF; both use the importer's incoming shortwave > 0 filter. The GPP missing percentage is for DT among intervals with known positive shortwave. Tower ET uses the full-day mean of finite LE_CORR, converted to mm/day. Empty cells represent unavailable values.

To inspect a ZIP without running Earth Engine, use:

    node tools/tower_batch.cjs path/to/ACCLIM_five_site_years.zip

The hosted importer opens this same ZIP with a site-year selector. Multi-site execution uses the batch command above.
