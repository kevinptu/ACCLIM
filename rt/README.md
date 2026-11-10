# ACCLIM-RT

This spatial experiment asks whether differences in canopy patch cover can explain pixel-to-pixel variation in directional red and NIR reflectance within one area on one observation date. Each site has its own soil and vegetated patch spectra. The working assumption is that cover differences contribute more to the observed BRF variation than differences in soil reflectance.

## Reproduce the three tower-site runs

From the repository root with Python 3.12, install the simulator and the optional Earth Engine sampler:

    python -m pip install -r rt/requirements.txt
    python -m pip install earthengine-api
    earthengine authenticate

Replace PROJECT with an Earth Engine project you can use:

    python rt/sample_modis.py PROJECT 32.5849 -106.6032 2024-06-01 Jornada.csv
    python rt/sample_modis.py PROJECT 38.4309 -120.9660 2024-06-01 Tonzi.csv
    python rt/sample_modis.py PROJECT 60.08649722 17.47950278 2024-06-01 Norunda.csv
    python rt/acclim_rt.py Jornada.csv results/Jornada
    python rt/acclim_rt.py Tonzi.csv results/Tonzi
    python rt/acclim_rt.py Norunda.csv results/Norunda

The sampler uses a 10 km square (--km changes its side), native MODIS Aqua MYD09A1 red and NIR pixels, land and clear-sky screening, good band QA, and the actual pixel observation day and angles. The positional date is the composite start. The fitter selects the actual day with the most valid pixels, breaking ties by earliest date. Set --date YYYY-MM-DD to select another actual day. Inputs from another source can use the same CSV header with unitless reflectance, grid indices and angles in degrees, one area per file.

## Forward model and spatial test

The library evaluates 270 combinations of soil brightness, dry/wet soil mix, within-patch LAI, chlorophyll and mean leaf angle. The same explicit Lambertian soil spectrum is passed to PROSAIL beneath the patch and used as the exposed-soil endpoint. Directional BRF uses median solar zenith, view zenith and relative azimuth for the selected pixels; --angles overrides these. Red and NIR use equal-weight averages in the 620–670 and 841–876 nm windows. Sensor response weighting and within-box angular variability are approximations.

For each band, pixel BRF = soil BRF + cover × (vegetated patch BRF − soil BRF). Cover is inferred from red alone and bounded to [0,1]. The program selects the common soil and patch that minimize total squared red and NIR residuals on training pixels. Including the red residual penalizes segments that cannot reach the observed red range. It withholds each quadrant of the sampled grid in turn. Test NIR is used only for scoring. The benchmark is the training mean NIR. Lower held-out RMSE than this benchmark indicates predictive skill on that date. Variable substrate, litter, canopy properties and cover can produce similar two-band signals. Independent cover or soil observations are needed for attribution.

Outputs include the simulated spectral library, full-data endmembers and scatterplot, per-quadrant errors and fitted endmembers, and held-out cover and NIR predictions. Spatial testing is independent of tower GPP and ET.

[PROSAIL implementation](https://github.com/jgomezdans/prosail) · [MODIS Aqua reflectance product](https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MYD09A1)

## Norunda cover sensitivity

    python rt/cover_sensitivity.py Norunda.csv results/Norunda-cover --date 2024-06-02

Fits both observed bands with equal weights. PROSPECT-D leaf parameters are fixed (N 1.5, Cab 40, Car 8, brown pigment 0, water 0.015 cm, dry mass 0.009 g/cm²). Each independent run fixes patch LAI at 3, 5, or 7, uses a spherical leaf-angle distribution and hotspot parameter 0.05, and evaluates directional FourSAIL at each pixel's MODIS geometry. The same common Lambertian soil pair is used under the patch and in the exposed fraction. Soil BRF is constant within each band. Per-pixel cover minimizes the two-band squared residual by projection onto the soil-to-patch segment, constrained to [0,1]. A second optimization constrains every pixel to [0.8,1].

The broad soil grid spans red 0.01–0.40 and NIR from red to 0.50, spaced by 0.01. A restricted subset lies within Euclidean distance 0.01 of a PROSAIL wet/dry soil mixture with brightness 0.5–1.5. The reported ranges include soil pairs within 0.001, 0.002 or 0.005 absolute BRF of the best joint RMSE **within each family and LAI**. They are sensitivity intervals. The fits measure spectral reconstruction; independent canopy cover is needed to assess retrieval accuracy.
