/**
 * ============================================================================
 * LANDFLUX / ACCLIM v1.0 — Google Earth Engine port
 * ============================================================================
 *
 * Faithful GEE translation of gpp_acclimate_v3a_GRD_K2.Rmd
 * (Assimilation Capacity Coordinates with Light Interception Model, K. Tu).
 *
 * Architecture goals (intentionally over-modular for a v1 — see SECTION ORDER):
 *   - Pure model functions, separated from data I/O.
 *   - A single "forcing FeatureCollection" schema that any data source must
 *     conform to. Adding a new data source = writing one adapter that returns
 *     that schema.
 *   - Parameters live in a single dictionary that can be overridden by callers
 *     (the entry point for MCMC calibration and Monte Carlo uncertainty).
 *   - Uncertainty / calibration / ML hooks are stubbed with documented
 *     signatures so future modules plug in without rewiring the core.
 *
 * SECTION ORDER
 *   0. CONFIG       - constants, dataset IDs, defaults
 *   1. UTILS        - vapor pressure, daylength, smoothing, pressure helpers
 *   2. CO2          - MVP: constant 400 ppm (replace with lookup later)
 *   3. ADAPTERS     - one function per gridded data source
 *   4. ASSEMBLY     - merge adapters into the canonical forcing FC
 *   5. MODEL        - ACCLIM core, sequential running means via ee.List.iterate
 *   6. RUNNER       - top-level runAtPoint + uncertainty/ML stubs
 *   7. UI           - panel, inputs, results, charts
 *   8. EXPORTS      - to Drive
 *   9. ENTRY POINT  - assemble UI, attach to root
 *
 * RAW INPUTS (same data sources as PT-JPL — alignment was an explicit goal)
 *   TerraClimate (IDAHO_EPSCOR/TERRACLIMATE, monthly):
 *     srad (W/m2, scale 0.1), tmmn (degC, scale 0.1), tmmx (degC, scale 0.1),
 *     vap (kPa, scale 0.001) = actual vapor pressure ea; pr (mm/month, scale 1)
 *   MODIS (MODIS/061/MOD13Q1, 16-day, 250m):
 *     sur_refl_b01 (red), sur_refl_b02 (NIR), NDVI
 *
 * MVP DESIGN DECISIONS (each labeled so they're easy to find and revisit):
 *   [MVP-CO2]      ca_daily = 400 ppm (constant). Per user; replace with
 *                  NOAA/OCO-2 lookup when needed.
 *   [MVP-MONTHLY]  ACCLIM iterates at TerraClimate's native monthly cadence.
 *                  Daily carbon products interpolate monthly results afterward.
 *                  MODIS 16-day observations within each
 *                  month are mean-aggregated to one monthly value. With tc=7
 *                  days and dt≈30 days, w = exp(-30/7) ≈ 0.013, so each
 *                  running mean ≈ current month's value (the seed only
 *                  matters for month 1). The model is cadence-agnostic —
 *                  each forcing feature carries its own `dt` (days since
 *                  previous), so swapping to a finer-cadence source later
 *                  requires no model changes, only a new adapter.
 *   [MVP-PRESSURE] Pdaily = 101.325 kPa (sea-level standard, 1013.25 mbar)
 *                  for all sites. Per user; replace with elevation-derived
 *                  barometric formula or a gridded pressure product later.
 *   [MVP-Rn]       Rndaily derived internally using the PT-JPL formula
 *                  (isothermal net radiation), since neither TerraClimate
 *                  nor MODIS provide Rn directly and the user requested
 *                  "same inputs as PT-JPL".
 *   [MVP-SMOOTH]   At monthly cadence, MODIS observations are mean-aggregated
 *                  within each month — smoothing is implicit in the average.
 *                  The gap-filler smoother (smoothMovingAverage in Section 1,
 *                  port of vegetation_gap_filler.html) is kept available for
 *                  finer-cadence pipelines we'll add later.
 *
 * CANONICAL FORCING SCHEMA  (one feature per monthly timestep in v1; every
 *                             adapter must contribute toward this and the
 *                             assembled FC must have no nulls)
 *   date         string  'YYYY-MM-DD' (typically month-start)
 *   millis       number  date in millis since epoch (sort key)
 *   doy          number  day of year of the month's anchor date (1-366)
 *   dt           number  days since previous timestep (model reads this)
 *   Tmin         number  monthly mean of daily minimum air temperature, deg C
 *   Tmax         number  monthly mean of daily maximum air temperature, deg C
 *   Tdaily       number  monthly mean air temperature, deg C       [(Tmin+Tmax)/2]
 *   SWRdaily     number  monthly mean shortwave radiation, W/m2
 *   ea           number  monthly mean actual water vapor pressure, kPa
 *   VPDdaily     number  vapor pressure deficit, kPa               [es(Tdaily) - ea]
 *   ca_daily     number  atmospheric CO2 concentration, umol/mol (ppm) [MVP: 400]
 *   Pdaily       number  air pressure, kPa                              [MVP: 101.325]
 *   Rndaily      number  monthly mean net radiation, W/m2          [PT-JPL formula]
 *   EVI2         number  monthly mean enhanced vegetation index 2
 *   NDVI         number  monthly mean normalized difference vegetation index
 *   red          number  monthly mean red surface reflectance (0-1)
 *   nir          number  monthly mean NIR surface reflectance (0-1)
 *   daylength    number  astronomical daylength at the month's anchor doy, hours
 *
 *   Note on R-script deviations:
 *   - The R script uses FLUXNET's Tnight as input and derives Tmin = 2*Tnight
 *     - Tdaily. With gridded data we use Tmin/Tmax directly from TerraClimate
 *     and compute Tday consistently — only Tmin/Tmax/Tday enter the model
 *     after that point.
 *   - The R script's loop runs daily over FLUXNET data. We iterate monthly
 *     because the climate forcing is monthly. The model code is unchanged;
 *     the running-mean weight w = exp(-dt/tc) now uses dt ≈ 30 days instead
 *     of 1, so x_mean[i] ≈ current month's value (correct behavior per the
 *     EWMA, given that tc < dt).
 * ============================================================================
 */


// ============================================================================
// SECTION 0 — CONFIG
// ============================================================================

var ACCLIM_PARAMS = {
  faparMethod: 'NDVI',     // NDVI squared or EVI2 calibrated to NDVI squared
  k: 0.5,                  // canopy light extinction coefficient for LAI retrieval (must be > 0)
  CUEplant:       0.5,      // total plant NPP / GPP
  carbonToBiomass: 2.1276, // exact RAP conversion supplied by user (~47% C)
  fANPP_intercept: 0.171,   // fANPP = intercept + slope * MAT (deg C)
  fANPP_slope:    0.0129,
  tc:             7,        // running-mean time constant (days)
  GPPmean2GPPmax: 0.67,      // ratio of mean GPP to GPPmax asymptote
  Oi:             20900,    // leaf internal O2 partial pressure (Pa)
  QRC3:           12,       // C3 quantum requirement (mol photons / mol CO2)
  QRC4:           16.67,    // C4 quantum requirement
  EVI2min:        0.1,     // EVI2 of bare ground
  EVI2max:        0.8,     // EVI2 of full green canopy
  fresp:          0.1,      // leaf respiration as fraction of GPP capacity
  PT_coef:        1.26,     // Priestley-Taylor coefficient
  rho_air:        1.292,    // air density (kg/m3)
  Cp_air:         1010,     // specific heat of air (J/kg/K)
  gamma_psy:      0.066,    // psychrometric constant (kPa/K)
  cicaC4:         0.4,      // C4 ci/ca ratio (fixed)
  // C4frac logistic coefficients (from R: 1/(1+exp(-4.04*GPPmax_C4/GPPmax_C3 + 5.99)))
  C4frac_a:       -4.04,
  C4frac_b:       5.99,
  // fcold logistic coefficients: 1/(1+exp(0.1*(12-Tday_mean)))
  fcold_slope:    0.1,
  fcold_T0:       12,
  // tau temperature scaling: 2600 * 0.57^((T-25)/10)
  tau_base:       2600,
  tau_q10:        0.57,
  // Magnus vapor pressure
  magnus_a:       17.502,
  magnus_b:       240.97,
  magnus_es0:     0.61121,
  // PT-JPL net radiation derivation (used to fill ACCLIM's Rndaily input)
  sb:             5.67e-8,  // Stefan-Boltzmann (W/m2/K4)
  em_v:           0.97,     // emissivity of green vegetation
  em_nv:          0.92,     // emissivity of non-vegetation
  ndvi_min:       0.15,     // NDVI of bare soil
  ndvi_max:       0.95,     // NDVI of maximum green cover
  ndvi_pow:       1.25,     // exponent of fractional green cover ~ NDVI
  // Intercepted PAR fraction used for water and energy partitioning.
  fipar_ndvi_min: 0.15,
  fipar_ndvi_max: 0.97,
  fipar_ndvi_exp: 1.25
};

var DATASETS = {
  RAP_COVER:    'projects/rap-data-365417/assets/vegetation-cover-v3',
  TERRACLIMATE: 'IDAHO_EPSCOR/TERRACLIMATE',
  MODIS_VI:     'MODIS/061/MOD13Q1'
};

// MVP forcing constants — see the [MVP-*] notes in the header for rationale.
var MVP_CONSTANTS = {
  ca_daily: 400,            // CO2 (ppm), constant
  Pdaily:   101.325,        // air pressure (kPa) = 1013.25 mbar, sea-level standard
  // MODIS gap-filler smoothing (port of vegetation_gap_filler.html)
  viSmoothWindow:     5,    // centered window (days)
  viSmoothIterations: 20    // number of smoothing passes
};

var DEFAULTS = {
  lat:       38.41,         // Tonzi Ranch, CA
  lon:       -120.63,
  startDate: '2023-01-01',
  endDate:   '2024-12-31',
  // Spatial scale for sampling (m). TerraClimate native is ~4 km;
  // MODIS MOD13Q1 is 250 m. We sample both at this scale at the point.
  sampleScale: 1000
};


// ============================================================================
// SECTION 1 — UTILS
// ============================================================================
// Pure math helpers. All accept ee.Number (server-side) or plain JS numbers
// where indicated. Each function is testable in isolation.

/**
 * Magnus saturation vapor pressure (kPa) at temperature T (deg C).
 * Matches the R script's es formula exactly.
 *   esT = 0.61121 * exp(17.502 * T / (T + 240.97))
 */
function magnusVP(T) {
  T = ee.Number(T);
  var p = ACCLIM_PARAMS;
  return ee.Number(p.magnus_es0).multiply(
    T.multiply(p.magnus_a).divide(T.add(p.magnus_b)).exp()
  );
}

/**
 * Astronomical daylength (hours) from latitude (degrees) and day of year (1-366).
 * Standard formula (Forsythe et al. 1995 simplified): based on solar declination.
 *   delta = 0.4093 * sin(2*pi*(doy - 81)/365)    [solar declination, radians]
 *   cos(h0) = -tan(lat) * tan(delta)
 *   daylength = (2/pi) * acos(cos(h0)) * (12/1)  hours, clamped to [0.5, 23.5]
 *
 * R script used FLUXNET's DAY_D field; we compute from astronomy since DAY_D
 * is not available gridded. This is a necessary substitution, not a shortcut.
 */
function astronomicalDaylength(latDeg, doy) {
  var lat_rad = ee.Number(latDeg).multiply(Math.PI / 180);
  var doy_n = ee.Number(doy);
  // Declination (radians)
  var delta = doy_n.subtract(81).multiply(2 * Math.PI / 365).sin().multiply(0.4093);
  // Cosine of half-day-angle
  var cos_h0 = lat_rad.tan().multiply(-1).multiply(delta.tan());
  cos_h0 = cos_h0.max(-1).min(1);
  // Convert to hours: daylength = (2/pi) * acos(cos_h0) * 12
  var h_hours = cos_h0.acos().multiply(2).multiply(12).divide(Math.PI);
  return h_hours.max(0.5).min(23.5);
}

/**
 * Day-of-year (1-366) from a date string or ee.Date.
 */
function doyFromDate(date) {
  var d = ee.Date(date);
  var jan1 = ee.Date.fromYMD(d.get('year'), 1, 1);
  return d.difference(jan1, 'day').add(1);
}

/**
 * Build a daily date list (ee.List of ee.Date) spanning [startDate, endDate] inclusive.
 */
function dailyDateList(startDate, endDate) {
  var start = ee.Date(startDate);
  var end = ee.Date(endDate);
  var nDays = end.difference(start, 'day').add(1);
  return ee.List.sequence(0, nDays.subtract(1)).map(function(i) {
    return start.advance(ee.Number(i), 'day');
  });
}

/**
 * Iterative centered moving average smoother (port of
 * vegetation_gap_filler.html). Each pass replaces every value with the mean
 * of itself and ±halfWindow neighbors; repeated `iterations` times.
 *
 * Implemented server-side with ee.List.iterate so the smoothed series remains
 * in the compute graph and can feed straight into ACCLIM's iterate without a
 * client-side roundtrip.
 *
 * @param {ee.List} values         list of numbers (one per day)
 * @param {number} windowDays      window size; halfWindow = floor(windowDays/2)
 * @param {number} iterations      number of smoothing passes
 * @return {ee.List}
 */
function smoothMovingAverage(values, windowDays, iterations) {
  var halfWindow = Math.floor(windowDays / 2);

  function listMean(list) {
    list = ee.List(list);
    var sum = ee.Number(list.iterate(function(v, acc) {
      return ee.Number(acc).add(ee.Number(v));
    }, ee.Number(0)));
    return sum.divide(list.size());
  }

  function onePass(prevList) {
    prevList = ee.List(prevList);
    var n = prevList.size();
    return ee.List.sequence(0, n.subtract(1)).map(function(i) {
      i = ee.Number(i);
      var start = i.subtract(halfWindow).max(0);
      var end = i.add(halfWindow).min(n.subtract(1)).add(1);
      var window = prevList.slice(start, end);
      return listMean(window);
    });
  }

  return ee.List(
    ee.List.sequence(1, iterations).iterate(function(_, prev) {
      return onePass(prev);
    }, ee.List(values))
  );
}


// ============================================================================
// SECTION 2 — CO2
// ============================================================================
// MVP: constant 400 ppm (MVP_CONSTANTS.ca_daily). The signature of co2ForDate
// is preserved as the swap-in point for a later lookup. To enable time-varying
// CO2, replace this function body with a NOAA Mauna Loa anchor table or an
// OCO-2 adapter — no other code needs to change.

/**
 * CO2 (ppm) for a given date. MVP returns a constant; signature is the swap
 * point for a future time-varying lookup.
 *
 * @param {ee.Date|string} date  (unused in MVP; kept for future-compat)
 * @return {ee.Number}
 */
function co2ForDate(date) {
  return ee.Number(MVP_CONSTANTS.ca_daily);
}


// ============================================================================
// SECTION 3 — DATA ADAPTERS
// ============================================================================
// Each adapter:
//   - Takes (point, startDate, endDate) and a config object.
//   - Returns an ee.FeatureCollection keyed on date.
//   - Produces a *subset* of the canonical forcing properties.
// The assembly layer (Section 4) joins multiple adapters into a complete forcing FC.
//
// To add a new gridded source: copy one of these, populate the standardized
// property names, register it in assembleForcing.

/**
 * TerraClimate monthly aggregates -> Tmin, Tmax, SWRdaily, ea (vapor pressure).
 *
 * Band scale factors (from IDAHO_EPSCOR/TERRACLIMATE catalog):
 *   tmmn (min temp):   stored as degC * 10  -> scale 0.1
 *   tmmx (max temp):   stored as degC * 10  -> scale 0.1
 *   srad (shortwave):  stored as W/m2 * 10  -> scale 0.1
 *   vap (vapor press): stored as kPa * 1000 -> scale 0.001
 *
 * Returns one feature per MONTH (TerraClimate's native cadence). ACCLIM
 * iterates at this cadence; carbon products are interpolated after the run.
 */
function fetchTerraClimate(geom, startDate, endDate, scale) {
  // Small pad so a date right at the boundary still finds a bracketing month.
  var pad = 31;
  var startPadded = ee.Date(startDate).advance(-pad, 'day');
  var endPadded   = ee.Date(endDate).advance(pad, 'day');

  var ic = ee.ImageCollection(DATASETS.TERRACLIMATE)
    .filterDate(startPadded, endPadded)
    .filterBounds(geom);

  return ic.map(function(img) {
    var d = ee.Date(img.get('system:time_start'));
    // ee.Reducer.mean works for both points (single pixel = identity) and
    // polygons (spatial average over enclosed pixels), so the same fetcher
    // serves drawn points, drawn polygons, and uploaded asset geometries.
    var s = img.select(['tmmn', 'tmmx', 'srad', 'vap', 'pr']).reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: geom,
      scale: scale,
      bestEffort: true,
      maxPixels: 1e10
    });
    var Tmin     = ee.Number(s.get('tmmn')).multiply(0.1);
    var Tmax     = ee.Number(s.get('tmmx')).multiply(0.1);
    var SWRdaily = ee.Number(s.get('srad')).multiply(0.1);
    var ea       = ee.Number(s.get('vap')).multiply(0.001);  // kPa

    return ee.Feature(geom, {
      date:         d.format('YYYY-MM-dd'),
      millis:       d.millis(),
      year:         d.get('year'),
      month:        d.get('month'),
      Tmin:         Tmin,
      Tmax:         Tmax,
      SWRdaily:     SWRdaily,
      ea:           ea,
      precip:       s.get('pr'), // TerraClimate monthly accumulation, mm (scale 1)
      period_days:  d.advance(1, 'month').difference(d, 'day')
    });
  });
}

/**
 * MODIS MOD13Q1 (16-day, 250m) -> NDVI, EVI2, red, nir at the point.
 * Reflectance bands sur_refl_b01/b02 use scale 0.0001 (stored as int * 10000).
 *
 * EVI2 is computed from red/NIR (2-band form) so the value matches what the
 * R script's EVI2 forcing represents:
 *   EVI2 = 2.5 * (NIR - RED) / (NIR + 2.4 * RED + 1)
 *
 * red/nir are kept in the output so the assembly layer can derive albedo for
 * the PT-JPL net-radiation calculation (albedo = (red + nir) / 2).
 */
function fetchModisVI(geom, startDate, endDate, scale) {
  // Pad the window so we have anchors for interpolation at start/end.
  var pad = 32; // days
  var startPadded = ee.Date(startDate).advance(-pad, 'day');
  var endPadded   = ee.Date(endDate).advance(pad, 'day');

  var ic = ee.ImageCollection(DATASETS.MODIS_VI)
    .filterDate(startPadded, endPadded)
    .filterBounds(geom);

  var sampled = ic.map(function(img) {
    var d = ee.Date(img.get('system:time_start'));
    var s = img.select(['NDVI', 'sur_refl_b01', 'sur_refl_b02']).reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: geom,
      scale: scale,
      bestEffort: true,
      maxPixels: 1e10
    });
    // Masked pixels reduce to null. Filter raw samples BEFORE any arithmetic.
    return ee.Feature(geom, s).set({date: d.format('YYYY-MM-dd'), millis: d.millis()});
  }).filter(ee.Filter.notNull(['NDVI', 'sur_refl_b01', 'sur_refl_b02']));

  return sampled.map(function(sample) {
    var s = ee.Feature(sample).toDictionary();
    var ndvi = ee.Number(s.get('NDVI')).divide(10000);
    var red  = ee.Number(s.get('sur_refl_b01')).divide(10000).max(0.0001);
    var nir  = ee.Number(s.get('sur_refl_b02')).divide(10000).max(0.0001);
    var evi2 = nir.subtract(red).multiply(2.5)
      .divide(nir.add(red.multiply(2.4)).add(1))
      .max(0).min(1);

    return ee.Feature(geom, {
      date:   s.get('date'),
      millis: s.get('millis'),
      NDVI:   ndvi,
      EVI2:   evi2,
      red:    red,
      nir:    nir
    });
  });
}


// ============================================================================
// SECTION 4 — FORCING ASSEMBLY
// ============================================================================

/**
 * Linearly interpolate a sparse property (e.g., EVI2 from 16-day MODIS) to
 * the daily date list. For each daily date we find the bracketing observations
 * by millis and linearly interpolate. Edge dates are extrapolated by holding
 * the nearest anchor (flat extrapolation).
 *
 * NOTE — Not used in the v1 monthly pipeline (MODIS is aggregated within each
 * TerraClimate month instead). Kept available for finer-cadence pipelines
 * we'll add later (e.g., FLUXNET daily, or a daily climate product).
 *
 * @param {ee.FeatureCollection} sparseFC  with 'millis' and the named property
 * @param {ee.List} dailyDates             list of ee.Date
 * @param {string} propName
 * @return {ee.List} list of {date, millis, value}-shaped dictionaries
 */
function interpolateSparseToDaily(sparseFC, dailyDates, propName) {
  var sorted = sparseFC.sort('millis');
  return dailyDates.map(function(d) {
    d = ee.Date(d);
    var t = d.millis();
    var beforeFC = sorted.filter(ee.Filter.lte('millis', t)).sort('millis', false);
    var afterFC  = sorted.filter(ee.Filter.gt('millis', t)).sort('millis', true);
    var hasBefore = beforeFC.size().gt(0);
    var hasAfter  = afterFC.size().gt(0);

    // Bracketing anchors. If only one side exists (start/end of MODIS
    // collection beyond the padded window), fall back to the existing side.
    var prev = ee.Feature(ee.Algorithms.If(hasBefore, beforeFC.first(), afterFC.first()));
    var next = ee.Feature(ee.Algorithms.If(hasAfter,  afterFC.first(),  beforeFC.first()));

    var t0 = ee.Number(prev.get('millis'));
    var t1 = ee.Number(next.get('millis'));
    var v0 = ee.Number(prev.get(propName));
    var v1 = ee.Number(next.get(propName));
    var dtSpan = t1.subtract(t0);
    // When prev == next (one-sided edge case), fall back to v0 (no extrap).
    var frac = ee.Number(ee.Algorithms.If(dtSpan.gt(0),
      ee.Number(t).subtract(t0).divide(dtSpan),
      0));
    var v = v0.add(v1.subtract(v0).multiply(frac));

    return ee.Dictionary({
      date:   d.format('YYYY-MM-dd'),
      millis: t,
      value:  v
    });
  });
}

/**
 * Mean-aggregate MODIS observations that fall within a given month.
 * Returns a feature with EVI2/NDVI/red/nir mean values, or nulls if no
 * MODIS observation lies within the month (caller must handle).
 *
 * @param {ee.FeatureCollection} modisFC  fetched MODIS features
 * @param {ee.Date} monthStart            month start date
 * @param {ee.Date} monthEnd              first day of next month
 * @return {ee.Feature}
 */
function aggregateModisToMonth(modisFC, monthStart, monthEnd) {
  var inMonth = modisFC.filter(ee.Filter.and(
    ee.Filter.gte('millis', monthStart.millis()),
    ee.Filter.lt('millis', monthEnd.millis())
  ));
  // If no obs land in the month (rare at MODIS's 16-day cadence over a full
  // month), fall back to the single closest observation to month midpoint.
  var midMillis = monthStart.millis().add(monthEnd.millis()).divide(2);
  var fallback = modisFC.map(function(f) {
    return f.set('distFromMid', f.getNumber('millis').subtract(midMillis).abs());
  }).sort('distFromMid').first();

  var hasObs = inMonth.size().gt(0);
  var EVI2 = ee.Number(ee.Algorithms.If(hasObs,
    inMonth.aggregate_mean('EVI2'),
    ee.Feature(fallback).getNumber('EVI2')));
  var NDVI = ee.Number(ee.Algorithms.If(hasObs,
    inMonth.aggregate_mean('NDVI'),
    ee.Feature(fallback).getNumber('NDVI')));
  var red  = ee.Number(ee.Algorithms.If(hasObs,
    inMonth.aggregate_mean('red'),
    ee.Feature(fallback).getNumber('red')));
  var nir  = ee.Number(ee.Algorithms.If(hasObs,
    inMonth.aggregate_mean('nir'),
    ee.Feature(fallback).getNumber('nir')));

  return ee.Feature(null, {
    EVI2: EVI2, NDVI: NDVI, red: red, nir: nir
  });
}

/**
 * Build the canonical forcing FeatureCollection over a geometry for [start, end].
 * One feature per TerraClimate month, with MODIS aggregated to that month
 * and derived quantities (Tdaily, VPD, albedo, Rn, daylength, dt) filled in.
 *
 * Accepts any ee.Geometry — a point (drawn or programmatic), a drawn polygon,
 * or a polygon from an uploaded GEE asset. The data adapters mean-aggregate
 * across the geometry's enclosed pixels; for a point this is identity, for a
 * polygon it's a spatial average. Latitude (for astronomical daylength) is
 * derived from the geometry's centroid.
 *
 * @param {ee.Geometry} geom    point, polygon, or asset geometry
 * @param {string} startDate
 * @param {string} endDate
 * @param {Object} config       { scale }
 * @return {ee.FeatureCollection}  sorted by millis
 */
function assembleForcing(geom, startDate, endDate, config) {
  config = config || {};
  var scale = config.scale || DEFAULTS.sampleScale;
  var p = ACCLIM_PARAMS;

  // Centroid latitude for daylength. coordinates() returns [lon, lat] so we
  // pull index 1. maxError keeps the centroid computation tractable on big
  // polygons.
  var centroidCoords = geom.centroid({maxError: 1}).coordinates();
  var latDeg = ee.Number(centroidCoords.get(1));

  // Filter TerraClimate to the requested window (TC features land on the
  // first of each month, so this gives one feature per month).
  var towerRows = config.towerWeatherRows;
  var terra = (towerRows ? ee.FeatureCollection(towerRows.map(function(row){return towerWeatherFeature(row,geom);})) : fetchTerraClimate(geom, startDate, endDate, scale))
    .filter(ee.Filter.and(
      ee.Filter.gte('millis', ee.Date(startDate).millis()),
      ee.Filter.lte('millis', ee.Date(endDate).millis())
    ))
    .sort('millis');

  var modis = config.modis || fetchModisVI(geom, startDate, endDate, scale);

  // For each TerraClimate month, aggregate MODIS and compute derived fields.
  var enriched = terra.map(function(f) {
    var monthStart = ee.Date(f.getString('date'));
    var monthEnd   = monthStart.advance(1, 'month');

    var Tmin     = f.getNumber('Tmin');
    var Tmax     = f.getNumber('Tmax');
    var SWRdaily = f.getNumber('SWRdaily');
    var ea       = f.getNumber('ea');

    // Day of year of the month's anchor (mid-month for daylength fidelity)
    var midMonth = monthStart.advance(15, 'day');
    var doy = doyFromDate(midMonth);
    var dl  = astronomicalDaylength(latDeg, doy);

    // Temperatures (per R script)
    var Tdaily = towerRows ? f.getNumber('Tdaily') : Tmin.add(Tmax).divide(2);
    var Tday   = Tmax.add(Tmax.add(Tmin).multiply(0.5)).divide(2);

    // Saturation VPs and VPD (using observed ea, not derived from dewpoint)
    var esTday = magnusVP(Tday);
    var VPDdaily = esTday.subtract(ea).max(0);

    // MODIS monthly aggregation
    var vi = aggregateModisToMonth(modis, monthStart, monthEnd);
    var EVI2 = vi.getNumber('EVI2');
    var NDVI = vi.getNumber('NDVI');
    var red  = vi.getNumber('red');
    var nir  = vi.getNumber('nir');

    // PT-JPL net-radiation derivation (Rndaily, W/m2, 24-hour mean)
    //   albedo     = (red + nir) / 2
    //   greencover = ((NDVI - ndvi_min) / (ndvi_max - ndvi_min))^ndvi_pow, clamped
    //   em_sur     = greencover * em_v + (1 - greencover) * em_nv
    //   em_air     = 1.31 * (10 * ea / (Tdaily + 273.15))^(1/7)
    //   Rndaily    = (1 - albedo) * SWRdaily + sb * (Tdaily + 273.15)^4 * (em_air - em_sur)
    var albedo = red.add(nir).divide(2).max(0).min(1);
    var greencover = NDVI.subtract(p.ndvi_min).max(0)
      .divide(p.ndvi_max - p.ndvi_min)
      .pow(p.ndvi_pow).min(1);
    var em_sur = greencover.multiply(p.em_v)
      .add(ee.Number(1).subtract(greencover).multiply(p.em_nv));
    var em_air = ee.Number(10).multiply(ea).divide(Tdaily.add(273.15))
      .pow(1 / 7).multiply(1.31);
    var T_K4 = Tdaily.add(273.15).pow(4);
    var Rndaily = ee.Number(1).subtract(albedo).multiply(SWRdaily)
      .add(T_K4.multiply(p.sb).multiply(em_air.subtract(em_sur)));

    if (towerRows) Rndaily = f.getNumber('Rndaily');
    return f.set({
      weather_source: towerRows ? 'tower' : 'gridded',
      doy:       doy,
      daylength: dl,
      Tdaily:    Tdaily,
      Tday:      Tday,
      esTday:    esTday,
      VPDdaily:  VPDdaily,
      EVI2:      EVI2,
      NDVI:      NDVI,
      red:       red,
      nir:       nir,
      albedo:    albedo,
      greencover: greencover,
      Rndaily:   Rndaily,
      // MVP constants
      Pdaily:    towerRows ? f.getNumber('Pdaily') : MVP_CONSTANTS.Pdaily,
      ca_daily:  co2ForDate(monthStart)
    });
  });

  // Attach `dt` (days since previous timestep) to each feature so the model
  // is cadence-agnostic. First feature gets dt = 30 (typical month).
  var sorted = enriched.sort('millis');
  var sortedList = sorted.toList(sorted.size());
  var dayMs = 86400000;
  var withDt = ee.FeatureCollection(
    ee.List.sequence(0, sortedList.size().subtract(1)).map(function(i) {
      i = ee.Number(i);
      var f = ee.Feature(sortedList.get(i));
      var prevMillis = ee.Number(ee.Algorithms.If(
        i.gt(0),
        ee.Feature(sortedList.get(i.subtract(1).max(0))).getNumber('millis'),
        f.getNumber('millis').subtract(dayMs * 30)  // synthetic prev for step 0
      ));
      var dt = f.getNumber('millis').subtract(prevMillis).divide(dayMs);
      return f.set('dt', dt);
    })
  );

  // Do not fabricate vegetation values when the entire sampling window is masked.
  return ee.FeatureCollection(ee.Algorithms.If(modis.size().gt(0), withDt, ee.FeatureCollection([])));
}


// ============================================================================
// SECTION 5 — MODEL CORE
// ============================================================================
// The ACCLIM model in three pieces:
//   computeInitialState  - mean of first 5 timesteps, used to seed running means
//   acclimStep           - one timestep's update; appends output to history
//   runACCLIM            - drives ee.List.iterate, returns output FC
//
// Running-mean rule (from R script):
//     x_mean[i] = w[i] * x_mean[i-1] + (1 - w[i]) * x[i]      if PAR >= 10
//     x_mean[i] = x_mean[i-1]                                  if PAR < 10
//     x_mean[1] = mean(x[1:5])                                 (initial seed)
// where  dt[i] = days since previous timestep (read from feature property)
//        w[i]  = exp(-dt[i] / tc)
//
// At monthly cadence (TerraClimate), dt ≈ 30, so w ≈ exp(-30/7) ≈ 0.013 and
// each running mean ≈ current month's value — exactly the asymptotic limit
// expected of an EWMA when the sampling interval far exceeds the time
// constant. The model code is cadence-agnostic: feed it a daily forcing and
// dt drops to 1; feed it 16-day forcing and dt becomes 16.

/**
 * Seed the running means from the first 5 timesteps of forcing.
 * The R script does mean(x[1:5]) for each running-mean variable.
 *
 * @param {ee.List} forcingList  list of ee.Feature (one per timestep, sorted)
 * @param {Object} params
 * @return {ee.Dictionary} initial state for ee.List.iterate
 */
function computeInitialState(forcingList, params) {
  var first5 = ee.FeatureCollection(forcingList.slice(0, 5));

  // PAR is derived from SWR each day; we need to derive it for the first 5 too
  // so we can seed PARday_mean. We'll compute on the fly here.
  var first5_with_par = first5.map(function(f) {
    var swr = f.getNumber('SWRdaily');
    var dl = f.getNumber('daylength');
    var par = swr.multiply(2.2).multiply(ee.Number(24).divide(dl));
    var Tmin = f.getNumber('Tmin');
    var Tmax = f.getNumber('Tmax');
    var Tday = Tmax.add(Tmax.add(Tmin).multiply(0.5)).divide(2);
    var esTmin = magnusVP(Tmin);
    var esTmax = magnusVP(Tmax);
    var esTday = magnusVP(Tday);
    var useTowerHumidity = ee.Algorithms.IsEqual(f.get('weather_source'), 'tower');
    var humidityPressure = ee.Number(ee.Algorithms.If(useTowerHumidity, f.get('ea'), esTmin));
    var VPDmax = ee.Number(ee.Algorithms.If(useTowerHumidity, esTmax.subtract(humidityPressure).max(0), esTmax.subtract(esTmin)));
    var VPDday = ee.Number(ee.Algorithms.If(useTowerHumidity, esTday.subtract(humidityPressure).max(0), esTday.subtract(esTmin)));
    var RHmin = humidityPressure.divide(esTmax).max(0).min(1);
    var Rnday = f.getNumber('Rndaily').multiply(ee.Number(24).divide(dl));
    return f.set({
      Tday: Tday,
      esTmin: esTmin, esTmax: esTmax, esTday: esTday,
      VPDmax: VPDmax, VPDday: VPDday,
      PARday: par, RHmin: RHmin, Rnday: Rnday
    });
  });

  // Pmean in R is stored in Pa (1000 * Pdaily where Pdaily is kPa)
  var Pmean0      = first5_with_par.aggregate_mean('Pdaily').multiply(1000);
  var ca_mean0    = first5_with_par.aggregate_mean('ca_daily');
  var Tday_mean0  = first5_with_par.aggregate_mean('Tday');
  var Tmin_mean0  = first5_with_par.aggregate_mean('Tmin');
  var Tmax_mean0  = first5_with_par.aggregate_mean('Tmax');
  var esTday_mean0 = first5_with_par.aggregate_mean('esTday');
  var esTmax_mean0 = first5_with_par.aggregate_mean('esTmax');
  var VPDday_mean0 = first5_with_par.aggregate_mean('VPDday');
  var VPDmax_mean0 = first5_with_par.aggregate_mean('VPDmax');
  var PARday_mean0 = first5_with_par.aggregate_mean('PARday');
  var Rnday_mean0  = first5_with_par.aggregate_mean('Rnday');
  var RHmin_mean0  = first5_with_par.aggregate_mean('RHmin');

  return ee.Dictionary({
    // running means (carry across iterations)
    Pmean:        Pmean0,
    ca_mean:      ca_mean0,
    Tday_mean:    Tday_mean0,
    Tmin_mean:    Tmin_mean0,
    Tmax_mean:    Tmax_mean0,
    esTday_mean:  esTday_mean0,
    esTmax_mean:  esTmax_mean0,
    VPDday_mean:  VPDday_mean0,
    VPDmax_mean:  VPDmax_mean0,
    PARday_mean:  PARday_mean0,
    Rnday_mean:   Rnday_mean0,
    RHmin_mean:   RHmin_mean0,
    // dt/w state
    dt:           1,
    // step counter (so the first call knows it's the seed day)
    step:         0,
    // accumulating output history (one feature per timestep)
    history:      ee.List([])
  });
}

// FIPAR is independent of the selected photosynthetic FAPAR method.
function computeFIPAR(ndvi, p) {
  return ee.Number(ndvi).subtract(p.fipar_ndvi_min).max(0)
    .divide(p.fipar_ndvi_max - p.fipar_ndvi_min).pow(p.fipar_ndvi_exp);
}

/**
 * Compute mean-state-dependent diagnostics from a state dictionary +
 * today's instantaneous values. Pure function of inputs; no side effects.
 * Returns a dictionary of all derived quantities to attach to the output feature.
 */
function computeDiagnostics(running, instant, params) {
  var p = params;

  var Tday_mean   = ee.Number(running.get('Tday_mean'));
  var Tmax_mean   = ee.Number(running.get('Tmax_mean'));
  var esTday_mean = ee.Number(running.get('esTday_mean'));
  var esTmax_mean = ee.Number(running.get('esTmax_mean'));
  var VPDday_mean = ee.Number(running.get('VPDday_mean'));
  var VPDmax_mean = ee.Number(running.get('VPDmax_mean'));
  var PARday_mean = ee.Number(running.get('PARday_mean'));
  var Rnday_mean  = ee.Number(running.get('Rnday_mean'));
  var RHmin_mean  = ee.Number(running.get('RHmin_mean'));
  var Pmean       = ee.Number(running.get('Pmean'));
  var ca_mean     = ee.Number(running.get('ca_mean'));

  var FAPAR_in = ee.Number(instant.get('EVI2'));

  // fcold: 1 / (1 + exp(0.1 * (12 - Tday_mean)))
  var fcold_mean = ee.Number(1).divide(
    ee.Number(1).add(
      ee.Number(p.fcold_T0).subtract(Tday_mean).multiply(p.fcold_slope).exp()
    )
  );

  // tau: 2600 * 0.57^((Tday_mean - 25) / 10)
  var tau_mean = ee.Number(p.tau_base).multiply(
    ee.Number(p.tau_q10).pow(Tday_mean.subtract(25).divide(10))
  );

  // gamma_star (Pa): Oi / (2 * tau)
  var gamma_star_mean = ee.Number(p.Oi).divide(tau_mean.multiply(2));

  // Daytime saturation-vapor-pressure slope (kPa/K), using running means.
  var s_mean = esTday_mean.multiply(p.magnus_a).multiply(p.magnus_b)
    .divide(Tday_mean.add(p.magnus_b).pow(2));

  // zC3 (R formula, units arranged for the quadratic):
  //   zC3 = ((1.5*GPPmean2GPPmax*2500*18)/(QRC3*1.26*(1+fresp))) *
  //         ((s + gamma_psy) * VPDday_mean * PARday_mean * 1000) /
  //         (s * Rnday_mean * Pmean)
  var zCoef = ee.Number(1.5).multiply(p.GPPmean2GPPmax).multiply(2500).multiply(18)
    .divide(ee.Number(p.QRC3).multiply(p.PT_coef).multiply(ee.Number(1).add(p.fresp)));
  var zC3 = zCoef.multiply(
    s_mean.add(p.gamma_psy).multiply(VPDday_mean).multiply(PARday_mean).multiply(1000)
  ).divide(s_mean.multiply(Rnday_mean).multiply(Pmean));

  // Quadratic solves ciC3 in ppm, using Gamma* converted from Pa to ppm.
  // a=1, b=2*gamma_ppm+zC3-ca_mean, c=-gamma_ppm*(zC3+2*ca_mean).
  // Divide the positive root by ca_mean to obtain dimensionless ci/ca.
  var gamma_ppm = gamma_star_mean.multiply(1e6).divide(Pmean);
  var bC3 = gamma_ppm.multiply(2).add(zC3).subtract(ca_mean);
  var cC3 = gamma_ppm.multiply(zC3.multiply(-1).subtract(ca_mean.multiply(2)));
  // Guard against negative discriminant (numerical noise) by clamping to 0.
  var disc = bC3.pow(2).subtract(cC3.multiply(4)).max(0);
  var ciC3_ppm = bC3.multiply(-1).add(disc.sqrt()).divide(2);
  var cicaC3_mean = ciC3_ppm.divide(ca_mean);

  var cicaC4_mean = ee.Number(p.cicaC4);

  // pi (Pa): cica * ca (umol/mol) * Pmean (Pa) / 1e6
  var piC3_mean = cicaC3_mean.multiply(ca_mean).multiply(Pmean).divide(1e6);
  var piC4_mean = cicaC4_mean.multiply(ca_mean).multiply(Pmean).divide(1e6);

  // Quantum yields
  var QYC3_mean = fcold_mean.multiply(ee.Number(1).divide(p.QRC3))
    .multiply(piC3_mean.subtract(gamma_star_mean))
    .divide(piC3_mean.add(gamma_star_mean.multiply(2)));
  var QYC4_mean = fcold_mean.multiply(ee.Number(1).divide(p.QRC4));

  // Optional FAPAR calibration (enabled through the run controls).
  // Fit paired, quality-screened observations: NDVI^2 = a + b * EVI2 (b > 0).
  // Derive effective endpoints EVI2min = -a/b and EVI2max = (1-a)/b, then
  // use FAPAR = clamp((EVI2-EVI2min)/(EVI2max-EVI2min), 0, 1).
  // Intent: anchor the scale to NDVI^2 while retaining EVI2 temporal variation
  // and its potential sensitivity to green foliage amid non-green material.
  // These are calibration endpoints, not necessarily observed soil/full-canopy values.
  // If revisited, fit fixed endpoints across multiple years with balanced VI coverage;
  // check slope/endpoints and validate on withheld years against independent data.
  // Fitting to NDVI^2 can inherit its saturation/bias; improvement is not guaranteed.
  // Fit once across the selected record; coefficients stay fixed within each run.
  // Previous FAPAR from scaled EVI2 (retained for comparison).
  // var FAPAR = FAPAR_in.subtract(p.EVI2min)
  //   .divide(p.EVI2max - p.EVI2min)
  //   .max(0).min(1);

  // FAPAR from squared NDVI, using the current monthly NDVI input.
  var FIPAR = computeFIPAR(instant.get('NDVI'), p);
  var FAPAR = ee.Number(instant.get('NDVI')).pow(2);
  if (p.faparMethod === 'EVI2') {
    var scaledEVI2 = FAPAR_in.subtract(p.faparFit.min)
      .divide(p.faparFit.max.subtract(p.faparFit.min)).max(0).min(1);
    FAPAR = ee.Number(ee.Algorithms.If(p.faparFit.valid, scaledEVI2, FAPAR));
  }

  // GPPmax (umolCO2/m2s) = QY * FAPAR * PARday_mean
  var GPPmax_C3 = QYC3_mean.multiply(FAPAR).multiply(PARday_mean);
  var GPPmax_C4 = QYC4_mean.multiply(FAPAR).multiply(PARday_mean);

  // Mean daytime GPP
  var GPPday_mean_C3 = GPPmax_C3.multiply(p.GPPmean2GPPmax);
  var GPPday_mean_C4 = GPPmax_C4.multiply(p.GPPmean2GPPmax);

  // C4 fraction: 1 / (1 + exp(-4.04 * GPPmax_C4/GPPmax_C3 + 5.99))
  //   guard against division by zero on GPPmax_C3
  var ratio = GPPmax_C4.divide(GPPmax_C3.max(1e-9));
  var C4frac = ee.Number(1).divide(
    ee.Number(1).add(
      ratio.multiply(p.C4frac_a).add(p.C4frac_b).exp()
    )
  );

  // Maximum daytime GPP, with the same C3/C4 mixture as mean daytime GPP.
  var GPPmax = ee.Number(1).subtract(C4frac).multiply(GPPmax_C3)
    .add(C4frac.multiply(GPPmax_C4));

  // LAI is leaf area / ground area; Amax is per m2 leaf area.
  // Zero cover gives LAI=0 and undefined Amax; FAPAR>=1 has no finite LAI.
  var laiValid = FAPAR.gte(0).and(FAPAR.lt(1));
  var safeFAPAR = ee.Number(ee.Algorithms.If(laiValid, FAPAR, 0));
  var laiValue = ee.Number(1).subtract(safeFAPAR).log().multiply(-1).divide(p.k);
  var LAI = ee.Algorithms.If(laiValid, laiValue, null);
  var Amax = ee.Algorithms.If(laiValid.and(FAPAR.gt(0)),
    GPPmax.multiply(0.9).divide(ee.Number(ee.Algorithms.If(laiValue.gt(0), laiValue, 1))), null);

  // PredGPPmean (umolCO2/m2s)
  var PredGPPmean = ee.Number(1).subtract(C4frac).multiply(GPPday_mean_C3)
    .add(C4frac.multiply(GPPday_mean_C4));

  // Arithmetic C3/C4 fraction-weighted ci/ca diagnostic.
  var cica_mean = ee.Number(1).subtract(C4frac).multiply(cicaC3_mean)
    .add(C4frac.multiply(cicaC4_mean));

  // Water-vapor canopy conductance, mol H2O/m2 ground/s.
  // GPP is umol CO2/m2/s and ca is umol/mol, so the micro factors cancel.
  // Compute each pathway before fraction weighting; do not invert mean ci/ca.
  function pathwayConductance(gpp, cica) {
    var gradient = ca_mean.multiply(ee.Number(1).subtract(cica));
    var valid = ca_mean.gt(0).and(cica.gte(0)).and(cica.lt(1)).and(gpp.gte(0));
    var value = gpp.multiply(1.6 * 0.9)
      .divide(ee.Number(ee.Algorithms.If(valid, gradient, 1)));
    return {valid: valid, value: value};
  }
  var gs3 = pathwayConductance(GPPday_mean_C3, cicaC3_mean);
  var gs4 = pathwayConductance(GPPday_mean_C4, cicaC4_mean);
  var gsValid = gs3.valid.or(C4frac.eq(1)).and(gs4.valid.or(C4frac.eq(0)));
  var gsValue = gs3.value.multiply(ee.Number(1).subtract(C4frac)).add(gs4.value.multiply(C4frac));
  var gs = ee.Algorithms.If(gsValid, gsValue, null);
  // Convert molar conductance to m/s before inversion: rs = molar air density / gs.
  // Zero/invalid conductance has no finite resistance and is displayed as a gap.
  var temperatureK = Tday_mean.add(273.15);
  var rsValid = gsValid.and(gsValue.gt(0)).and(Pmean.gt(0)).and(temperatureK.gt(0));
  var rs = ee.Algorithms.If(rsValid, Pmean.divide(
    ee.Number(ee.Algorithms.If(rsValid, temperatureK.multiply(8.314462618).multiply(gsValue), 1))), null);

  // ET partitioning
  var s_ratio_mean = s_mean.divide(s_mean.add(p.gamma_psy));
  var fwet_mean = RHmin_mean.pow(10);
  var fSM_mean = RHmin_mean.pow(VPDmax_mean);
  var PT_mean = ee.Number(p.PT_coef).multiply(s_ratio_mean).multiply(Rnday_mean);

  var PredLEcday_mean = PT_mean.multiply(FIPAR).multiply(fcold_mean)
    .multiply(ee.Number(1).subtract(fwet_mean));
  var PredLEiday_mean = fwet_mean.multiply(PT_mean).multiply(FIPAR);
  var PredLEsday_mean = fwet_mean.add(ee.Number(1).subtract(fwet_mean).multiply(fSM_mean))
    .multiply(PT_mean).multiply(ee.Number(1).subtract(FIPAR));
  var PredLEday_mean = PredLEcday_mean.add(PredLEsday_mean).add(PredLEiday_mean);

  return ee.Dictionary({
    fcold_mean:       fcold_mean,
    tau_mean:         tau_mean,
    gamma_star_mean:  gamma_star_mean,
    s_mean:           s_mean,
    zC3:              zC3,
    bC3:              bC3,
    cC3:              cC3,
    cicaC3_mean:      cicaC3_mean,
    cicaC4_mean:      cicaC4_mean,
    cica_mean:        cica_mean,
    gs_C3: ee.Algorithms.If(gs3.valid, gs3.value, null),
    gs_C4: ee.Algorithms.If(gs4.valid, gs4.value, null),
    gs: gs, gs_units: 'mol H2O/m2 ground/s',
    rs: rs, rs_units: 's/m',
    piC3_mean:        piC3_mean,
    piC4_mean:        piC4_mean,
    QYC3_mean:        QYC3_mean,
    QYC4_mean:        QYC4_mean,
    FAPAR:            FAPAR,
    FIPAR:            FIPAR,
    GPPmax_C3:        GPPmax_C3,
    GPPmax_C4:        GPPmax_C4,
    GPPday_mean_C3:   GPPday_mean_C3,
    GPPday_mean_C4:   GPPday_mean_C4,
    C4frac:           C4frac,
    PredGPPmean:      PredGPPmean,
    GPPmax:           GPPmax,
    LAI: LAI, Amax: Amax, k: p.k,
    LAI_units: 'm2 leaf/m2 ground', Amax_units: 'umol CO2/m2 leaf/s',
    s_ratio_mean:     s_ratio_mean,
    fwet_mean:        fwet_mean,
    fSM_mean:         fSM_mean,
    PT_mean:          PT_mean,
    PredLEcday_mean:  PredLEcday_mean,
    PredLEiday_mean:  PredLEiday_mean,
    PredLEsday_mean:  PredLEsday_mean,
    PredLEday_mean:   PredLEday_mean
  });
}

/**
 * One iteration step. Called by ee.List.iterate.
 * Updates running means conditionally on PARday >= 10, computes diagnostics,
 * appends a feature to history.
 */
function acclimStep(currentFeat, prevState, params) {
  var state = ee.Dictionary(prevState);
  var f = ee.Feature(currentFeat);
  var p = params;

  // Today's instantaneous forcings (derive PAR, Tday, vapor pressures, RH, Rn)
  var Tmin = f.getNumber('Tmin');
  var Tmax = f.getNumber('Tmax');
  var Tdaily = f.getNumber('Tdaily');
  var SWRdaily = f.getNumber('SWRdaily');
  var Pdaily = f.getNumber('Pdaily');
  var Rndaily = f.getNumber('Rndaily');
  var ca_daily = f.getNumber('ca_daily');
  var EVI2 = f.getNumber('EVI2');
  var NDVI = f.getNumber('NDVI');
  var dl = f.getNumber('daylength');
  var doy = f.getNumber('doy');

  // Tday = (Tmax + 0.5*(Tmax+Tmin)) / 2   (from R script)
  var Tday = Tmax.add(Tmax.add(Tmin).multiply(0.5)).divide(2);

  var esTmin = magnusVP(Tmin);
  var esTmax = magnusVP(Tmax);
  var esTday = magnusVP(Tday);
  var useTowerHumidity = ee.Algorithms.IsEqual(f.get('weather_source'), 'tower');
  var humidityPressure = ee.Number(ee.Algorithms.If(useTowerHumidity, f.get('ea'), esTmin));
  var VPDmax = ee.Number(ee.Algorithms.If(useTowerHumidity, esTmax.subtract(humidityPressure).max(0), esTmax.subtract(esTmin)));
  var VPDday = ee.Number(ee.Algorithms.If(useTowerHumidity, esTday.subtract(humidityPressure).max(0), esTday.subtract(esTmin)));

  // PARday (umol/m2/s) = 2.2 * SWR(W/m2) * (24/dl)
  var PARday = SWRdaily.multiply(2.2).multiply(ee.Number(24).divide(dl));

  // Rnday (W/m2 over daylight) = Rndaily * (24/dl)
  var Rnday = Rndaily.multiply(ee.Number(24).divide(dl));

  // RHmin (fraction)
  var RHmin = humidityPressure.divide(esTmax).max(0).min(1);

  // dt and w. dt is supplied by the forcing FC (days since previous timestep),
  // so the model is cadence-agnostic — daily, 16-day, monthly all work.
  // When PARday < 10 (low-light timestep) we carry forward the previous dt
  // to keep w consistent with the running-mean carry-forward rule.
  var step = ee.Number(state.get('step'));
  var isFirst = step.eq(0);
  var lowLight = PARday.lt(10);
  var dt_supplied = f.getNumber('dt');
  var dt_prev = ee.Number(state.get('dt'));
  var dt = ee.Number(ee.Algorithms.If(
    isFirst, dt_supplied,
    ee.Algorithms.If(lowLight, dt_prev, dt_supplied)
  ));
  var w = dt.multiply(-1).divide(p.tc).exp();

  // Running-mean update closure
  function updateMean(prevKey, todayVal) {
    var prev = ee.Number(state.get(prevKey));
    return ee.Number(ee.Algorithms.If(
      isFirst, prev,                      // seed value retained on first step
      ee.Algorithms.If(
        lowLight, prev,                   // carry forward
        w.multiply(prev).add(ee.Number(1).subtract(w).multiply(todayVal))
      )
    ));
  }

  // R script multiplies Pdaily by 1000 (kPa -> Pa) when feeding Pmean
  var Pmean       = updateMean('Pmean',       Pdaily.multiply(1000));
  var ca_mean     = updateMean('ca_mean',     ca_daily);
  var Tday_mean   = updateMean('Tday_mean',   Tday);
  var Tmin_mean   = updateMean('Tmin_mean',   Tmin);
  var Tmax_mean   = updateMean('Tmax_mean',   Tmax);
  var esTday_mean = updateMean('esTday_mean', esTday);
  var esTmax_mean = updateMean('esTmax_mean', esTmax);
  var VPDday_mean = updateMean('VPDday_mean', VPDday);
  var VPDmax_mean = updateMean('VPDmax_mean', VPDmax);
  var PARday_mean = updateMean('PARday_mean', PARday);
  var Rnday_mean  = updateMean('Rnday_mean',  Rnday);
  var RHmin_mean  = updateMean('RHmin_mean',  RHmin);

  var runningNow = ee.Dictionary({
    Pmean: Pmean, ca_mean: ca_mean, Tday_mean: Tday_mean,
    Tmin_mean: Tmin_mean, Tmax_mean: Tmax_mean, esTmax_mean: esTmax_mean, esTday_mean: esTday_mean,
    VPDday_mean: VPDday_mean, VPDmax_mean: VPDmax_mean,
    PARday_mean: PARday_mean, Rnday_mean: Rnday_mean,
    RHmin_mean: RHmin_mean
  });
  var instant = ee.Dictionary({EVI2: EVI2, NDVI: NDVI});

  var diag = computeDiagnostics(runningNow, instant, p);

  // Output feature: forcings + instantaneous + running means + diagnostics
  var outProps = ee.Dictionary({
    date: f.getString('date'),
    millis: f.getNumber('millis'),
    doy: doy,
    // Forcings (carried through for reference)
    weather_source: f.get('weather_source'), weather_status: f.get('weather_status'),
    weather_energy_basis: f.get('weather_energy_basis'), weather_expected_intervals: f.get('weather_expected_intervals'), ea: f.get('ea'),
    Tmin: Tmin, Tmax: Tmax, Tdaily: Tdaily, Tday: Tday,
    SWRdaily: SWRdaily, Pdaily: Pdaily, Rndaily: Rndaily,
    precip: f.get('precip'), period_days: f.get('period_days'),
    ca_daily: ca_daily, EVI2: EVI2, NDVI: NDVI, daylength: dl,
    // Instantaneous derived
    esTmin: esTmin, esTmax: esTmax, esTday: esTday,
    VPDmax: VPDmax, VPDday: VPDday, PARday: PARday, Rnday: Rnday, RHmin: RHmin,
    dt: dt, w: w
  }).combine(runningNow).combine(diag);

  // Daytime mean flux integrated over daylight, zero GPP at night.
  // 12 g C/mol (not 44 g CO2/mol); no extra mean-to-maximum adjustment.
  var GPPdaily = ee.Number(diag.get('PredGPPmean')).multiply(dl).multiply(3600).multiply(12e-6);
  outProps = outProps.combine(ee.Dictionary({
    GPPdaily: GPPdaily, NPPdaily: GPPdaily.multiply(p.CUEplant), CUEplant: p.CUEplant
  }));
  var outFeat = ee.Feature(f.geometry(), outProps);

  // Return updated state
  return ee.Dictionary({
    Pmean: Pmean, ca_mean: ca_mean, Tday_mean: Tday_mean,
    Tmin_mean: Tmin_mean, Tmax_mean: Tmax_mean, esTmax_mean: esTmax_mean, esTday_mean: esTday_mean,
    VPDday_mean: VPDday_mean, VPDmax_mean: VPDmax_mean,
    PARday_mean: PARday_mean, Rnday_mean: Rnday_mean,
    RHmin_mean: RHmin_mean,
    dt: dt,
    step: step.add(1),
    history: ee.List(state.get('history')).add(outFeat)
  });
}

/**
 * Run the full ACCLIM model over a forcing FeatureCollection.
 *
 * @param {ee.FeatureCollection} forcingFC  sorted by millis, with canonical schema
 * @param {Object} params                   parameter overrides merged into ACCLIM_PARAMS
 * @return {ee.FeatureCollection}           one feature per timestep with all outputs
 */
// Retrospective effective evaporative storage calibration. Monthly totals in mm.
// Start each candidate full; no spin-up. Surface/litter evaporation (fwet) is
// external to this bucket. All precipitation enters it: this is an empirical
// moisture index, not a closed whole-system water balance or volumetric SM.
// Search a bounded log grid including 100 mm, then refine three times around the best.
// The original RH/VPD prediction remains the calibration target and fallback.
var SOIL_BUCKET = {initial: 100, minimum: 1, maximum: 2000, gridSize: 25};

// Exact period integration with constant mean forcing, using period totals only.
// dSM/du = precip - (potentialMM/capacity)*SM, u in [0,1], capped at capacity.
// Precipitation is uniformly available within the period; storm timing is unknown.
// Integrate evaporation as storage changes, rather than holding start-SM for a month.
function soilBucketPeriod(sm, capacity, precip, potentialMM) {
  sm = ee.Number(sm); capacity = ee.Number(capacity);
  precip = ee.Number(precip); potentialMM = ee.Number(potentialMM);
  var q = potentialMM.divide(capacity);
  var surplus = precip.subtract(potentialMM);
  var gapOverSurplus = capacity.subtract(sm).max(0).divide(surplus.max(1e-12));
  var x = q.multiply(gapOverSurplus);
  // Stable log(1+x)/x and exponential integrals near zero.
  var logRatio = ee.Number(ee.Algorithms.If(x.lt(1e-5),
    ee.Number(1).subtract(x.divide(2)).add(x.pow(2).divide(3)),
    x.add(1).log().divide(x.max(1e-30))));
  var t = ee.Number(ee.Algorithms.If(surplus.gt(0), gapOverSurplus.multiply(logRatio).min(1), 1));
  var z = q.multiply(t);
  var decay = z.multiply(-1).exp();
  var phi = ee.Number(ee.Algorithms.If(z.lt(1e-4),
    ee.Number(1).subtract(z.divide(2)).add(z.pow(2).divide(6)).subtract(z.pow(3).divide(24)),
    ee.Number(1).subtract(decay).divide(z.max(1e-30))));
  var psi = ee.Number(ee.Algorithms.If(z.lt(1e-4),
    ee.Number(0.5).subtract(z.divide(6)).add(z.pow(2).divide(24)).subtract(z.pow(3).divide(120)),
    ee.Number(1).subtract(phi).divide(z.max(1e-30))));
  var meanSM = sm.multiply(t).multiply(phi).add(precip.multiply(t.pow(2)).multiply(psi))
    .add(capacity.multiply(ee.Number(1).subtract(t))).max(0).min(capacity);
  var fraction = meanSM.divide(capacity);
  var evaporation = potentialMM.multiply(fraction);
  var overflow = surplus.max(0).multiply(ee.Number(1).subtract(t));
  var endSM = sm.add(precip).subtract(evaporation).subtract(overflow).max(0).min(capacity);
  return ee.Dictionary({meanSM: meanSM, fraction: fraction,
    evaporation: evaporation, overflow: overflow, endSM: endSM});
}

function soilBucketTrajectory(rows, capacity, keepHistory) {
  capacity = ee.Number(capacity);
  return ee.Dictionary(ee.List(rows).iterate(function(item, previous) {
    var r = ee.Feature(item), st = ee.Dictionary(previous);
    var sm = ee.Number(st.get('SM'));
    var pr = r.getNumber('precip'); // mm over this entire calendar month
    var days = r.getNumber('period_days');
    var conversion = r.getNumber('daylength').multiply(3600).multiply(days).divide(2450000);
    var wet = r.getNumber('fwet_mean');
    var potential = r.getNumber('PT_mean').multiply(ee.Number(1).subtract(r.getNumber('FIPAR')));
    var surfaceLE = wet.multiply(potential); // unchanged original surface-water term
    var potentialMM = ee.Number(1).subtract(wet).multiply(potential.max(0)).multiply(conversion);
    var period = soilBucketPeriod(sm, capacity, pr, potentialMM);
    var fraction = ee.Number(period.get('fraction'));
    var soilMM = ee.Number(period.get('evaporation'));
    var soilLE = soilMM.divide(conversion);
    var totalLE = surfaceLE.add(soilLE);
    var overflow = ee.Number(period.get('overflow'));
    var nextSM = ee.Number(period.get('endSM'));
    var error = totalLE.subtract(r.getNumber('PredLEsday_mean'));
    // Weight monthly mean flux errors by represented days; no fitting on
    // condensation/zero-energy periods, though their precipitation updates SM.
    var weight = ee.Number(ee.Algorithms.If(potential.gt(0), days, 0));
    var out = r.set({SM: sm, SM_mean: period.get('meanSM'), SM_end: nextSM, SMmax: capacity,
      fSM_start: sm.divide(capacity),
      SM_units: 'mm effective evaporative storage', fSM_mean: fraction,
      fSM_RHVPD: r.get('fSM_mean'), PredLEsday_RHVPD: r.get('PredLEsday_mean'),
      PredLEday_RHVPD: r.get('PredLEday_mean'),
      PredLEs_surface_mean: surfaceLE, PredLEs_soil_mean: soilLE,
      Esoil_period: soilMM, Esurface_period: surfaceLE.multiply(conversion),
      bucket_overflow: overflow, bucket_water_limited: potentialMM.gt(0).and(fraction.lt(1)),
      SM_integration: 'analytic period mean; uniform precipitation within period',
      PredLEsday_mean: totalLE,
      PredLEday_mean: r.getNumber('PredLEcday_mean').add(r.getNumber('PredLEiday_mean')).add(totalLE)});
    return ee.Dictionary({SM: nextSM,
      sse: ee.Number(st.get('sse')).add(error.pow(2).multiply(weight)),
      weight: ee.Number(st.get('weight')).add(weight),
      history: keepHistory ? ee.List(st.get('history')).add(out) : ee.List([])});
  }, ee.Dictionary({SM: capacity, sse: 0, weight: 0, history: ee.List([])})));
}

function fitSoilBucket(baseline, options) {
  var cfg = mergeParams(SOIL_BUCKET, options || {});
  if (!(cfg.minimum > 0 && cfg.maximum > cfg.minimum && cfg.initial >= cfg.minimum &&
        cfg.initial <= cfg.maximum && cfg.gridSize >= 3 && cfg.gridSize % 1 === 0))
    throw new Error('Invalid soil bucket capacity search settings.');
  baseline = baseline.sort('millis');
  var checked = baseline.map(function(f) {
    var required = ['precip', 'period_days', 'daylength', 'fwet_mean', 'PT_mean',
      'FIPAR', 'fSM_mean', 'PredLEsday_mean', 'PredLEday_mean', 'PredLEcday_mean', 'PredLEiday_mean'];
    var complete = f.toDictionary(required).size().eq(required.length);
    // Defaults only protect validation expressions; missing values are never filled.
    var d = ee.Dictionary.fromLists(required, ee.List.repeat(-9999, required.length))
      .combine(f.toDictionary(required));
    var valid = complete.and(ee.Number(d.get('precip')).gte(0))
      .and(ee.Number(d.get('period_days')).gte(28)).and(ee.Number(d.get('period_days')).lte(31))
      .and(ee.Number(d.get('daylength')).gt(0)).and(ee.Number(d.get('daylength')).lte(24))
      .and(ee.Number(d.get('fwet_mean')).gte(0)).and(ee.Number(d.get('fwet_mean')).lte(1))
      .and(ee.Number(d.get('FIPAR')).gte(0)).and(ee.Number(d.get('FIPAR')).lte(1))
      .and(ee.Number(d.get('fSM_mean')).gte(0)).and(ee.Number(d.get('fSM_mean')).lte(1));
    return f.set('_bucket_valid', valid);
  });
  var validRows = checked.filter(ee.Filter.eq('_bucket_valid', 1));
  var n = baseline.size();
  var first = ee.Date(baseline.aggregate_min('millis'));
  var last = ee.Date(baseline.aggregate_max('millis')).advance(1, 'month');
  var contiguous = last.difference(first, 'month').round().eq(n);
  var eligible = n.gte(12).and(validRows.size().eq(n)).and(contiguous);
  var rows = validRows.toList(validRows.size());
  function score(cap) {
    cap = ee.Number(cap);
    var trajectory = soilBucketTrajectory(rows, cap, false);
    return ee.Feature(null, {capacity: cap,
      mse: ee.Number(trajectory.get('sse')).divide(ee.Number(trajectory.get('weight')).max(1))});
  }
  var logMin = Math.log(cfg.minimum), logMax = Math.log(cfg.maximum);
  var step = (logMax - logMin) / (cfg.gridSize - 1);
  var candidates = ee.List.sequence(0, cfg.gridSize - 1).map(function(i) {
    return ee.Number(i).multiply(step).add(logMin).exp();
  });
  // Initial guess wins exact ties on an uninformative objective.
  var initialScore = score(cfg.initial);
  var coarse = ee.FeatureCollection(candidates.map(score));
  var coarseBest = ee.Feature(coarse.sort('mse').first());
  var best = ee.Feature(ee.Algorithms.If(initialScore.getNumber('mse').lte(coarseBest.getNumber('mse')),
    initialScore, coarseBest));
  var spread = ee.Number(coarse.aggregate_max('mse')).subtract(coarse.aggregate_min('mse'));
  for (var pass = 0; pass < 3; pass++) {
    var center = best.getNumber('capacity').log();
    var refined = ee.List.sequence(-4, 4).map(function(i) {
      return ee.Number(i).multiply(step / 4).add(center).max(logMin).min(logMax).exp();
    });
    var challenger = ee.Feature(ee.FeatureCollection(refined.map(score)).sort('mse').first());
    best = ee.Feature(ee.Algorithms.If(challenger.getNumber('mse').lt(best.getNumber('mse')), challenger, best));
    step /= 4;
  }
  var fitted = soilBucketTrajectory(rows, best.getNumber('capacity'), true);
  var informative = spread.gt(1e-8).and(ee.Number(fitted.get('weight')).gt(0));
  // Full initialization + precipitation >= full-bucket loss in EVERY period
  // implies fSM=1 for every capacity, including zero-demand records.
  var fullRows = validRows.map(function(f) {
    var loss = f.getNumber('PT_mean').multiply(ee.Number(1).subtract(f.getNumber('FIPAR'))).max(0)
      .multiply(ee.Number(1).subtract(f.getNumber('fwet_mean')))
      .multiply(f.getNumber('daylength')).multiply(3600).multiply(f.getNumber('period_days')).divide(2450000);
    return f.set('_always_full', f.getNumber('precip').gte(loss));
  }).filter(ee.Filter.eq('_always_full', 1));
  var allFull = eligible.and(fullRows.size().eq(n));
  var apply = eligible.and(informative).and(allFull.eq(0));
  var boundary = best.getNumber('capacity').lte(cfg.minimum * 1.001)
    .or(best.getNumber('capacity').gte(cfg.maximum / 1.001));
  var status = ee.String(ee.Algorithms.If(allFull, 'unconstrained_full_bucket', ee.Algorithms.If(eligible,
    ee.Algorithms.If(informative, ee.Algorithms.If(boundary, 'fitted_at_search_boundary', 'fitted'),
      'baseline_unidentifiable_capacity'), 'baseline_requires_complete_contiguous_year')));
  var fallback = baseline.map(function(f) {
    var potential = f.getNumber('PT_mean').multiply(ee.Number(1).subtract(f.getNumber('FIPAR')));
    var surface = f.getNumber('fwet_mean').multiply(potential);
    var soil = ee.Number(1).subtract(f.getNumber('fwet_mean')).multiply(potential.max(0));
    var fullLE = surface.add(soil);
    return f.set({fSM_RHVPD: f.get('fSM_mean'), PredLEsday_RHVPD: f.get('PredLEsday_mean'),
      PredLEday_RHVPD: f.get('PredLEday_mean'),
      fSM_mean: ee.Algorithms.If(allFull, 1, f.get('fSM_mean')),
      fSM_start: ee.Algorithms.If(allFull, 1, null),
      PredLEsday_mean: ee.Algorithms.If(allFull, fullLE, f.get('PredLEsday_mean')),
      PredLEday_mean: ee.Algorithms.If(allFull, fullLE.add(f.getNumber('PredLEcday_mean'))
        .add(f.getNumber('PredLEiday_mean')), f.get('PredLEday_mean'))});
  });
  return ee.FeatureCollection(ee.Algorithms.If(apply, ee.FeatureCollection(ee.List(fitted.get('history'))), fallback))
    .map(function(f) { return f.set({SM_fit_status: status,
      SM_fit_RMSE: ee.Algorithms.If(apply, best.getNumber('mse').sqrt(), null),
      SM_initial_guess: cfg.initial, SM_search_min: cfg.minimum, SM_search_max: cfg.maximum,
      SM_initialization: 'full; no spin-up', precip_units: 'mm/month'}); });
}


// Separate ET-effective diagnostic. Fit the FIPAR-free moisture equation to
// ORIGINAL full-canopy total LE, not the shallow-bucket-adjusted total LE.
// Only the isolated fitting proxy sets FIPAR=0. Copy diagnostics back explicitly,
// preserving actual FIPAR, FAPAR, GPP, transpiration, and all near-surface outputs.
function fitETBucket(nearSurface, options) {
  nearSurface = nearSurface.sort('millis');
  var proxy = nearSurface.map(function(f) {
    var properties = f.toDictionary(['date', 'millis', 'precip', 'period_days',
      'daylength', 'fwet_mean', 'PT_mean', 'fSM_mean']);
    return ee.Feature(null, properties.combine(ee.Dictionary({FIPAR: 0,
      PredLEsday_mean: f.get('PredLEday_RHVPD'),
      PredLEday_mean: f.get('PredLEday_RHVPD'),
      PredLEcday_mean: 0, PredLEiday_mean: 0})));
  });
  var fitted = fitSoilBucket(proxy, options);
  var pairs = nearSurface.toList(nearSurface.size()).zip(fitted.toList(fitted.size()));
  return ee.FeatureCollection(pairs.map(function(pair) {
    pair = ee.List(pair);
    var original = ee.Feature(pair.get(0)), et = ee.Feature(pair.get(1));
    var hasFit = ee.Algorithms.If(ee.Algorithms.IsEqual(et.get('SMmax'), null), false, true);
    function fittedValue(name) { return ee.Algorithms.If(hasFit, et.get(name), null); }
    var full = ee.Algorithms.IsEqual(et.get('SM_fit_status'), 'unconstrained_full_bucket');
    function moistureValue(name) { return ee.Algorithms.If(full, et.get(name), fittedValue(name)); }
    return original.set({SM_et: fittedValue('SM'), SM_et_mean: fittedValue('SM_mean'),
      SM_et_end: fittedValue('SM_end'), SMmax_et: fittedValue('SMmax'),
      fSM_et_mean: moistureValue('fSM_mean'), fSM_et_start: moistureValue('fSM_start'),
      PredLE_ETeffective: moistureValue('PredLEsday_mean'),
      PredLE_ETreference: original.get('PredLEday_RHVPD'),
      ETbucket_period: fittedValue('Esoil_period'), ETsurface_period: fittedValue('Esurface_period'),
      ETbucket_overflow: fittedValue('bucket_overflow'),
      SM_et_fit_status: et.get('SM_fit_status'), SM_et_fit_RMSE: et.get('SM_fit_RMSE'),
      SM_et_initial_guess: et.get('SM_initial_guess'), SM_et_search_min: et.get('SM_search_min'),
      SM_et_search_max: et.get('SM_search_max'), SM_et_initialization: et.get('SM_initialization'),
      SM_et_integration: fittedValue('SM_integration'),
      SM_et_units: 'mm effective ET storage', ETbucket_period_units: 'mm/period'});
  }));
}

// Absolute AWC regression supplied by the user (Saxton-model texture data).
// Sand is percent by mass; AWC is percent by volume. No bulk-density conversion.
// Extrapolates 0-30 cm texture over the effective profile; diagnostic depth only.
function effectiveDepthFields(row, sandPercent) {
  var valid = carbonFinite(sandPercent) && sandPercent >= 0 && sandPercent <= 100;
  var awc = valid ? (18.57 - 0.14 * sandPercent) / 100 : null;
  var result = {sand_0_30cm_pct: valid ? sandPercent : null, AWC_fraction: awc,
    AWC_pct_vol: valid ? awc * 100 : null};
  [['SMmax', 'effective_soil_evaporation_depth_cm'],
   ['SMmax_et', 'effective_root_depth_et_cm']].forEach(function(pair) {
    result[pair[1]] = valid && carbonFinite(row[pair[0]]) && row[pair[0]] > 0 ?
      row[pair[0]] / awc / 10 : null;
  });
  return result;
}
// Source-depth mixing uses the same ORIGINAL flux partition as the ET fit.
// Interception remains in total ET with zero soil-source depth.
// Integrate daytime LE over each period; do not average monthly flux ratios.
function transpirationSourceDepth(rows) {
  var tSum = 0, etSum = 0, eSum = 0, depthMoment = 0;
  var valid = rows.length > 0;
  var excluded = 0, included = 0;
  rows.forEach(function(r) {
    var t = r.PredLEcday_mean, e = r.PredLEsday_RHVPD;
    var et = r.PredLE_ETreference;
    var zET = r.effective_root_depth_et_cm;
    var zE = r.effective_soil_evaporation_depth_cm;
    // No evaporative source depth is defined when all component fluxes are nonpositive.
    // Exclude these periods from the diagnostic only; retain original model fluxes.
    if ([t, e, et].every(carbonFinite) && t <= 0 && e <= 0 && et <= 0) {
      excluded++; return;
    }
    var values = [t, e, et, r.daylength, r.period_days, zET, zE];
    if (!values.every(carbonFinite) || t < 0 || e < 0 || et < t + e - 1e-8 ||
        r.daylength <= 0 || r.period_days <= 0 || zET < 0 || zE < 0) {
      valid = false; return;
    }
    included++;
    var factor = r.daylength * 3600 * r.period_days / 2450000;
    tSum += t * factor; eSum += e * factor; etSum += et * factor;
    depthMoment += (zET * et - zE * e) * factor;
  });
  var depth = valid && tSum > 0 && etSum > 0 ? depthMoment / tSum : null;
  return {
    root_depth_periods_included: included,
    root_depth_periods_non_evaporating: excluded,
    effective_root_depth_t_cm: carbonFinite(depth) && depth >= 0 ? depth : null,
    root_depth_T_fraction: valid && etSum > 0 ? tSum / etSum : null,
    root_depth_E_fraction: valid && etSum > 0 ? eSum / etSum : null
  };
}
function depthChartLabel(info, field, label) {
  return label + ': ' + (carbonFinite(info[field]) ? (info[field] / 100).toFixed(2) + ' m (' + (info[field] / 30.48).toFixed(2) + ' ft)' : 'unavailable');
}

// Squared Pearson correlation over paired finite monthly fluxes (equal weight).
// Undefined for fewer than three pairs or a constant series; not 1-SSE/SST.
function bucketR2(rows, predicted, reference) {
  var pairs = rows.filter(function(r) { return carbonFinite(r[predicted]) && carbonFinite(r[reference]); });
  if (pairs.length < 3) return null;
  var mx = 0, my = 0;
  pairs.forEach(function(r) { mx += r[predicted] / pairs.length; my += r[reference] / pairs.length; });
  var xx = 0, yy = 0, xy = 0;
  pairs.forEach(function(r) { var x = r[predicted] - mx, y = r[reference] - my; xx += x*x; yy += y*y; xy += x*y; });
  return xx > 0 && yy > 0 ? Math.max(0, Math.min(1, (xy / Math.sqrt(xx) / Math.sqrt(yy)) * (xy / Math.sqrt(xx) / Math.sqrt(yy)))) : null;
}
function soilChartFooter(info, field) {
  var texture = carbonFinite(info.sand_0_30cm_pct) && carbonFinite(info.AWC_pct_vol) ?
    info.sand_0_30cm_pct.toFixed(1) + '%; AWC: ' + info.AWC_pct_vol.toFixed(2) + '%vol' : 'unavailable';
  return ui.Label({value: 'SoilGrids sand (0–30 cm): ' + texture + '\n' +
    'r² (bucket vs ACCLIM): ' + (carbonFinite(info[field]) ? info[field].toFixed(2) : 'unavailable'),
    style: {fontSize: '11px', color: '#666', whiteSpace: 'pre-wrap'}});
}

// Soil moisture here denotes effective extractable soil water, not total water content.
function addETMoistureCharts(output, info, regionName) {
//  var message = 'ET-effective bucket: ' + (info.SM_et_fit_status || 'unavailable');
  if (carbonFinite(info.SMmax_et)) {
//    message += '; SMmax,et = ' + info.SMmax_et.toFixed(2) + ' mm';
    addResultChart(ui.Chart.feature.byFeature(output, 'date', ['SM_et_mean', 'SMmax_et', 'precip'])
      .setChartType('ComboChart')
      .setSeriesNames(['SM,ET mean', 'SMmax,ET', 'precip'])
      .setOptions({title: regionChartTitle('Soil Moisture - Root Zone', regionName) + '\n' + depthChartLabel(info, 'effective_root_depth_t_cm', 'effective_root_depth_t_cm'),
        chartArea: {top: 110}, hAxis: {title: 'Date'},
        vAxes: {0: {title: 'Extractable soil moisture (mm)', viewWindow: {min: 0}},
                1: {title: 'Precipitation (mm)', viewWindow: {min: 0},
                    gridlines: {color: 'transparent'}, minorGridlines: {color: 'transparent'}}},
        seriesType: 'line',
        series: {0: {targetAxisIndex: 0, color: '#7B1FA2'},
                 1: {targetAxisIndex: 0, color: '#757575', lineDashStyle: [4, 4]},
                 2: {targetAxisIndex: 1, type: 'bars', color: '#64B5F6'}},
        lineWidth: 1.5, pointSize: 0}));
  }
  if (carbonFinite(info.SMmax_et)) addChartNote(soilChartFooter(info, 'SM_et_r2'));
//  UI.chartHolder.add(ui.Label(message));
}

// Retired transpiration-only bucket: see archive/transpiration_bucket_2026-09-18/README.md.
// Potential fluxes and diagnostic Omega use the final bucket-updated LE (W/m2).
function addPotentialFluxDiagnostics(f, p) {
  var pt = f.getNumber('PT_mean'), fipar = f.getNumber('FIPAR');
  var le = f.getNumber('PredLEday_mean'), t = f.getNumber('PredLEcday_mean');
  var e = f.getNumber('PredLEsday_mean');
  var pet = pt.multiply(2).subtract(le);
  var tpot = pt.multiply(fipar).multiply(2).subtract(t);
  var epot = pt.multiply(ee.Number(1).subtract(fipar)).multiply(2).subtract(e);
  function ratio(actual, potential) {
    // Nonpositive potential does not define a daytime evaporative coefficient.
    var valid = potential.gt(0);
    return ee.Algorithms.If(valid, actual.divide(ee.Number(ee.Algorithms.If(valid, potential, 1))), null);
  }
  var slope = f.getNumber('s_mean'), rn = f.getNumber('Rnday');
  var numerator = f.getNumber('VPDday_mean').multiply(p.rho_air).multiply(p.Cp_air);
  function aerodynamicResistance(potential, fraction) {
    // VPD and both vapor-pressure slopes use kPa consistently; result is s/m.
    var denominator = potential.multiply(slope.add(p.gamma_psy))
      .subtract(slope.multiply(rn).multiply(fraction));
    // Retain signed equation results for diagnosis; zero denominator is undefined.
    var valid = denominator.neq(0);
    return ee.Algorithms.If(valid, numerator.divide(ee.Number(ee.Algorithms.If(valid, denominator, 1))), null);
  }
  return f.set({PET_mean: pet, Tpot_mean: tpot, Epot_mean: epot,
    ra: aerodynamicResistance(pet, ee.Number(1)),
    rac: aerodynamicResistance(tpot, fipar),
    ras: aerodynamicResistance(epot, ee.Number(1).subtract(fipar)),
    aerodynamic_resistance_units: 's/m',
    Omega_surface: ratio(le, pet), Omega_canopy: ratio(t, tpot), Omega_soil: ratio(e, epot),
    potential_LE_units: 'W/m2 (mean daytime)', Omega_units: 'dimensionless'});
}

// Energy-balance temperature diagnostics (degrees C), after final LE updates.
function addTemperatureDiagnostics(f, p) {
  var air = f.getNumber('Tday_mean'), rn = f.getNumber('Rnday_mean');
  var fipar = f.getNumber('FIPAR');
  function temperature(resistanceKey, fraction, leKey) {
    var resistance = f.get(resistanceKey);
    var missing = ee.Algorithms.IsEqual(resistance, null);
    var safeResistance = ee.Number(ee.Algorithms.If(missing, 0, resistance));
    var value = air.add(safeResistance.multiply(rn.multiply(fraction).subtract(f.getNumber(leKey)))
      .divide(p.rho_air * p.Cp_air));
    return ee.Algorithms.If(missing, null, value);
  }
  var surface = temperature('ra', ee.Number(1), 'PredLEday_mean');
  var missingSurface = ee.Algorithms.IsEqual(surface, null);
  // User-specified maximum-temperature relationship; Tmin uses the same running mean.
  var maximum = ee.Number(ee.Algorithms.If(missingSurface, 0, surface)).multiply(1.33)
    .subtract(f.getNumber('Tmin_mean').multiply(0.33));
  return f.set({Tsurface_mean: surface,
    Tcanopy_mean: temperature('rac', fipar, 'PredLEcday_mean'),
    Tsoil_mean: temperature('ras', ee.Number(1).subtract(fipar), 'PredLEsday_mean'),
    Tsurface_max: ee.Algorithms.If(missingSurface, null, maximum),
    surface_temperature_units: 'degrees C'});
}

// Fit the paired geometry-mean monthly VI values used by this model run.
function fitEVI2FAPAR(forcingFC) {
  var pairs = forcingFC.filter(ee.Filter.notNull(['NDVI', 'EVI2']))
    .filter(ee.Filter.gte('NDVI', -1)).filter(ee.Filter.lte('NDVI', 1))
    .map(function(f) { return f.set('NDVI_squared', f.getNumber('NDVI').pow(2)); });
  var fit = ee.Dictionary(pairs.reduceColumns(ee.Reducer.linearFit(), ['EVI2', 'NDVI_squared']));
  var a = ee.Number(ee.Algorithms.If(ee.Algorithms.IsEqual(fit.get('offset'), null), 0, fit.get('offset')));
  var b = ee.Number(ee.Algorithms.If(ee.Algorithms.IsEqual(fit.get('scale'), null), 0, fit.get('scale')));
  var valid = pairs.size().gte(3).and(b.gt(1e-8));
  var safeB = ee.Number(ee.Algorithms.If(valid, b, 1));
  return {a: a, b: b, valid: valid, count: pairs.size(),
    min: a.multiply(-1).divide(safeB), max: ee.Number(1).subtract(a).divide(safeB)};
}

function runACCLIM(forcingFC, params) {
  // Merge param overrides (ES5-safe; GEE Code Editor lacks Object.assign).
  var p = mergeParams(ACCLIM_PARAMS, params || {});
  if (p.faparMethod !== 'NDVI' && p.faparMethod !== 'EVI2') throw new Error('Unknown FAPAR method.');
  if (p.faparMethod === 'EVI2' && !p.faparFit) p.faparFit = fitEVI2FAPAR(forcingFC);
  if (typeof p.k !== 'number' || !isFinite(p.k) || p.k <= 0) throw new Error('k must be a positive finite number.');
  var forcingList = forcingFC.sort('millis').toList(forcingFC.size());
  var initial = computeInitialState(forcingList, p);
  var finalState = ee.Dictionary(forcingList.iterate(function(curr, prev) {
    return acclimStep(curr, prev, p);
  }, initial));
  var baseline = ee.FeatureCollection(ee.List(finalState.get('history')));
  var nearSurface = fitSoilBucket(baseline, p.soilBucket);
  // Convert final partitioned fluxes after the soil-bucket update. Monthly rows
  // report equivalent mean daily water depths using their computed daylength;
  // these are not monthly totals. Nighttime flux is assumed zero.
  return fitETBucket(nearSurface, p.etBucket).map(function(f) {
    f = f.set('FAPAR_method_requested', p.faparMethod);
    if (p.faparMethod === 'EVI2') {
      f = f.set({FAPAR_method: ee.Algorithms.If(p.faparFit.valid, 'EVI2', 'NDVI'),
        FAPAR_fit_status: ee.Algorithms.If(p.faparFit.valid, 'fitted', 'NDVI fallback: insufficient pairs or nonpositive/degenerate slope'),
        FAPAR_fit_n: p.faparFit.count,
        FAPAR_fit_a: p.faparFit.a, FAPAR_fit_b: p.faparFit.b,
        FAPAR_EVI2min: ee.Algorithms.If(p.faparFit.valid, p.faparFit.min, null),
        FAPAR_EVI2max: ee.Algorithms.If(p.faparFit.valid, p.faparFit.max, null)});
    } else { f = f.set({FAPAR_method: 'NDVI', FAPAR_fit_status: 'not requested'}); }
    f = addPotentialFluxDiagnostics(f, p);
    f = addTemperatureDiagnostics(f, p);
    var factor = f.getNumber('daylength').multiply(3600).divide(2450000);
    // LE / (J/kg * kg/mol) gives mol H2O/m2/s; express carbon uptake in mmol.
    var transpirationMol = f.getNumber('PredLEcday_mean').divide(2450000 * 0.018);
    var netCarbonMmol = f.getNumber('PredGPPmean').multiply(0.9).divide(1000);
    var validWUE = transpirationMol.gt(0).and(netCarbonMmol.gte(0));
    var wue = ee.Algorithms.If(validWUE, netCarbonMmol.divide(
      ee.Number(ee.Algorithms.If(validWUE, transpirationMol, 1))), null);
    return f.set({
      Tair: f.get('Tdaily'), Tair_F: f.getNumber('Tdaily').multiply(1.8).add(32),
      VPD_kPa: f.get('VPDday_mean'), WUE: wue,
      transpiration_mol_m2_s: transpirationMol, WUE_units: 'mmol CO2/mol H2O',
      ET_mm_day: f.getNumber('PredLEday_mean').multiply(factor),
      T_mm_day: f.getNumber('PredLEcday_mean').multiply(factor),
      E_mm_day: f.getNumber('PredLEsday_mean').multiply(factor),
      I_mm_day: f.getNumber('PredLEiday_mean').multiply(factor)
    });
  });
}

/**
 * ES5-safe shallow merge: returns {...base, ...overrides}.
 * GEE's Code Editor runtime is an ES5 subset; Object.assign is unavailable.
 */
function mergeParams(base, overrides) {
  var out = {};
  var k;
  for (k in base) {
    if (base.hasOwnProperty(k)) out[k] = base[k];
  }
  for (k in overrides) {
    if (overrides.hasOwnProperty(k)) out[k] = overrides[k];
  }
  return out;
}


// ============================================================================
// SECTION 6 — RUNNER + UNCERTAINTY/ML HOOKS (stubs)
// ============================================================================

/**
 * Top-level runner: any geometry -> output FC.
 *
 * The forcing layer mean-aggregates across the geometry, so the same call
 * handles a single drawn point, a drawn polygon, or an uploaded asset.
 *
 * @param {ee.Geometry} geom
 * @param {string} startDate
 * @param {string} endDate
 * @param {Object} [paramOverrides]
 * @param {Object} [config]   { scale }
 * @return {ee.FeatureCollection}
 */
// Observational comparison only: MOD21C3 does not alter model forcing or fluxes.
function attachMonthlyLST(output, geom, start, end) {
  var source = ee.ImageCollection('MODIS/061/MOD11A1');
  var projection = ee.Image(source.first()).select('LST_Day_1km').projection();
  var collection = source.filterDate(start, end).filterBounds(geom).map(function(img) {
    var raw = img.select('LST_Day_1km');
    var good = img.select('QC_Day').bitwiseAnd(3).eq(0).and(raw.gt(0));
    return raw.multiply(0.02).subtract(273.15).updateMask(good).rename('LST')
      .copyProperties(img, ['system:time_start']);
  });
  var area = geom.area(1);
  var scale = ee.Number(ee.Algorithms.If(area.gt(0),
    projection.nominalScale().min(area.sqrt().divide(20).max(1)), projection.nominalScale()));
  return output.map(function(f) {
    var month = ee.Date(f.getString('date'));
    var images = collection.filterDate(month, month.advance(1, 'month'));
    // A masked placeholder guarantees typed bands even outside dataset coverage.
    var empty = ee.Image.constant(0).rename('LST').updateMask(ee.Image.constant(0));
    var safe = images.merge(ee.ImageCollection.fromImages([empty]));
    // First average valid days per pixel, then average over the selected geometry.
    // Only one spatial reduction per month, rather than one per daily image.
    var monthly = safe.mean().rename('LST').addBands(safe.count().rename('LST_valid_days'));
    var values = monthly.reduceRegion({reducer: ee.Reducer.mean(), geometry: geom,
      crs: projection, scale: scale, maxPixels: 1e8, tileScale: 4});
    return f.set({LST: values.get('LST'), LST_valid_days: values.get('LST_valid_days'),
      LST_image_count: images.size(), LST_units: 'degrees C', LST_source: 'MODIS/061/MOD11A1 LST_Day_1km',
      LST_sampling: 'Monthly mean of good-quality daily daytime LST; valid days averaged over AOI'});
  });
}

function attachMonthlyAquaLST(output, geom, start, end) {
  var source = ee.ImageCollection('MODIS/061/MYD11A1');
  var projection = ee.Image(source.first()).select('LST_Day_1km').projection();
  var collection = source.filterDate(start, end).filterBounds(geom).map(function(img) {
    var raw = img.select('LST_Day_1km');
    var good = img.select('QC_Day').bitwiseAnd(3).eq(0).and(raw.gt(0));
    return raw.multiply(0.02).subtract(273.15).updateMask(good).rename('Aqua_LST')
      .copyProperties(img, ['system:time_start']);
  });
  var area = geom.area(1);
  var scale = ee.Number(ee.Algorithms.If(area.gt(0),
    projection.nominalScale().min(area.sqrt().divide(20).max(1)), projection.nominalScale()));
  return output.map(function(f) {
    var month = ee.Date(f.getString('date'));
    var images = collection.filterDate(month, month.advance(1, 'month'));
    // A masked placeholder guarantees typed bands even outside dataset coverage.
    var empty = ee.Image.constant(0).rename('Aqua_LST').updateMask(ee.Image.constant(0));
    var safe = images.merge(ee.ImageCollection.fromImages([empty]));
    // First average valid days per pixel, then average over the selected geometry.
    // Only one spatial reduction per month, rather than one per daily image.
    var monthly = safe.mean().rename('Aqua_LST').addBands(safe.count().rename('Aqua_LST_valid_days'));
    var values = monthly.reduceRegion({reducer: ee.Reducer.mean(), geometry: geom,
      crs: projection, scale: scale, maxPixels: 1e8, tileScale: 4});
    return f.set({Aqua_LST: values.get('Aqua_LST'), Aqua_LST_valid_days: values.get('Aqua_LST_valid_days'),
      Aqua_LST_image_count: images.size(), Aqua_LST_units: 'degrees C', Aqua_LST_source: 'MODIS/061/MYD11A1 LST_Day_1km',
      Aqua_LST_sampling: 'Monthly mean of good-quality daily daytime LST; valid days averaged over AOI'});
  });
}

// Monthly mean of valid scene-level AOI means; counts are valid scenes, not pixels.
// Comparison only: Landsat observations never alter ACCLIM calculations.
function attachLandsatLST(output, geom, start, end) {
  var scenes = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'))
    .filterBounds(geom).filterDate(start, end).filter(ee.Filter.eq('PROCESSING_LEVEL', 'L2SP'));
  var area = geom.area(1);
  var samplingScale = ee.Number(ee.Algorithms.If(area.gt(0), area.sqrt().divide(20).max(1).min(30), 30));
  var sampled = scenes.map(function(image) {
    // Exclude fill, dilated cloud, cirrus, cloud, shadow, snow and terrain occlusion.
    var clear = image.select('QA_PIXEL').bitwiseAnd(63).eq(0)
      .and(image.select('QA_RADSAT').bitwiseAnd(1 << 11).eq(0));
    var raw = image.select('ST_B10');
    var lst = raw.multiply(0.00341802).add(149.0).subtract(273.15)
      .updateMask(clear.and(raw.gt(0))).rename('Landsat_LST');
    var mean = lst.reduceRegion({reducer: ee.Reducer.mean(), geometry: geom,
      crs: raw.projection(), scale: samplingScale, maxPixels: 1e8, tileScale: 4}).get('Landsat_LST');
    return image.set('_AOI_LST', mean);
  });
  return output.map(function(f) {
    var month = ee.Date(f.getString('date'));
    var candidates = sampled.filterDate(month, month.advance(1, 'month'));
    var valid = candidates.filter(ee.Filter.notNull(['_AOI_LST']));
    return f.set({Landsat_LST: ee.Algorithms.If(valid.size().gt(0), valid.aggregate_mean('_AOI_LST'), null),
      Landsat_LST_scene_count: valid.size(), Landsat_LST_candidate_count: candidates.size(),
      Landsat_LST_units: 'degrees C', Landsat_LST_source: 'Landsat 8/9 Collection 2 Tier 1 Level 2 ST_B10',
      Landsat_LST_sampling: 'Monthly mean of quality-screened scene AOI means; 30 m grid, native thermal 100 m'});
  });
}

function runForGeometry(geom, startDate, endDate, paramOverrides, config) {
  config = config || {};
  if (config.weatherPlan) {
    // Use the same full-window vegetation samples and FAPAR fit in both weather modes.
    var firstTowerMonth=ee.Date(startDate.slice(0,7)+'-01');
    var towerEnd=ee.Date(endDate.slice(0,7)+'-01').advance(1,'month');
    var fullModis=fetchModisVI(geom,firstTowerMonth,towerEnd.advance(-1,'day'),config.scale||DEFAULTS.sampleScale);
    var allVI=ee.FeatureCollection(ee.Algorithms.If(fullModis.size().gt(0),
      ee.FeatureCollection(ee.List.sequence(0,towerEnd.difference(firstTowerMonth,'month').subtract(1)).map(function(i){
        var month=firstTowerMonth.advance(i,'month');return aggregateModisToMonth(fullModis,month,month.advance(1,'month'));
      })),ee.FeatureCollection([])));
    var groups=[];
    config.weatherPlan.segments.forEach(function(rows){
      var forcing=assembleForcing(geom,rows[0].month+'-01',ee.Date(rows[rows.length-1].month+'-01').advance(1,'month').advance(-1,'day'),{towerWeatherRows:rows,modis:fullModis,scale:config.scale});
      groups.push(forcing);
    });
    var overrides=mergeParams(paramOverrides||{},{});
    if(overrides.faparMethod==='EVI2')overrides.faparFit=fitEVI2FAPAR(allVI);
    var result=ee.FeatureCollection([]);
    groups.forEach(function(forcing){result=result.merge(ee.FeatureCollection(ee.Algorithms.If(forcing.size().gt(0),runACCLIM(forcing,overrides),ee.FeatureCollection([]))));});
    var missing=ee.FeatureCollection(config.weatherPlan.missing.map(function(month){return ee.Feature(geom,{date:month+'-01',millis:ee.Date(month+'-01').millis(),weather_source:'tower',weather_status:'missing_weather'});}));
    return result.merge(missing).sort('millis').map(function(f){return f.set('weather_fields',JSON.stringify(config.weatherPlan.source));});
  }

  // Include every month intersecting the inclusive requested date range.
  // This prevents mid-month requests from losing their first month's forcing.
  var firstMonth = ee.Date(ee.Date(startDate).format('YYYY-MM-01'));
  var lastDay = ee.Date(ee.Date(endDate).format('YYYY-MM-01')).advance(1, 'month').advance(-1, 'day');
  var forcing = assembleForcing(geom, firstMonth, lastDay, config);
  return ee.FeatureCollection(ee.Algorithms.If(forcing.size().gt(0),
    runACCLIM(forcing, paramOverrides), ee.FeatureCollection([])));
}

/**
 * Convenience wrapper for the common lat/lon case.
 *
 * @param {number} lat
 * @param {number} lon
 * @param {string} startDate
 * @param {string} endDate
 * @param {Object} [paramOverrides]
 * @param {Object} [config]
 * @return {ee.FeatureCollection}
 */
function runAtPoint(lat, lon, startDate, endDate, paramOverrides, config) {
  var point = ee.Geometry.Point([lon, lat]);
  return runForGeometry(point, startDate, endDate, paramOverrides, config);
}

/**
 * STUB: Run an ensemble of ACCLIM realizations under sampled parameter sets.
 * Used for Monte Carlo uncertainty propagation. Each sample produces one
 * output FC; the ensemble can be summarized (mean, quantiles) externally.
 *
 * Future implementation:
 *   - paramSamples: array of partial-param overrides, one per ensemble member
 *   - returns: ee.List of output FeatureCollections, indexed by sample id
 *
 * @param {number} lat
 * @param {number} lon
 * @param {string} startDate
 * @param {string} endDate
 * @param {Array<Object>} paramSamples
 * @return {ee.List}
 */
function runEnsemble(lat, lon, startDate, endDate, paramSamples) {
  throw new Error(
    'runEnsemble is a stub. Provide paramSamples sampled from prior, then ' +
    'call runAtPoint for each. Aggregate client-side after .getInfo().'
  );
}

/**
 * STUB: MCMC calibration of ACCLIM parameters against observed daily GPP/LE.
 * Future implementation will likely run client-side (after .getInfo() to pull
 * forcing + obs locally) since GEE does not support stochastic sampling
 * efficiently server-side. The interface freezes here so that downstream UI
 * code can be written in advance.
 *
 * @param {ee.FeatureCollection} forcingFC
 * @param {ee.FeatureCollection} observationsFC  with GPP_obs, LE_obs, date
 * @param {Object} priors         { paramName: {dist: 'normal'|'uniform', ...} }
 * @param {Object} mcmcConfig     { iterations, burnin, thin, ... }
 * @return {Object}               { samples: [...], posterior: {...} }
 */
function mcmcCalibrate(forcingFC, observationsFC, priors, mcmcConfig) {
  throw new Error('mcmcCalibrate is a stub. See JSDoc for the planned interface.');
}

/**
 * STUB: Knowledge-guided ML correction of ACCLIM predictions.
 * Trains/applies a model that learns the residual between ACCLIM and
 * observations, with physics-derived features (FAPAR, fcold, etc.) as
 * inputs alongside raw forcings.
 *
 * @param {ee.FeatureCollection} acclimOutputFC
 * @param {ee.FeatureCollection} observationsFC
 * @param {Object} mlConfig    { features: [...], model: 'rf'|'gbm'|... }
 * @return {ee.FeatureCollection}  acclimOutputFC + corrected_GPP, corrected_LE
 */
function kgmlPredict(acclimOutputFC, observationsFC, mlConfig) {
  throw new Error('kgmlPredict is a stub. See JSDoc for the planned interface.');
}


// ============================================================================
// SECTION 7 — UI
// ============================================================================

// ============================================================================
// CARBON PRODUCTS — client-side processing of evaluated monthly model output.
// No daily meteorology is inferred: these are interpolated monthly estimates.
// ============================================================================
var CARBON_UNITS = { dayMs: 86400000, gramsCPerMol: 12,
  squareMetresPerAcre: 4046.86, poundsPerKilogram: 2.20462 };

function carbonDateMillis(s) {
  if (!/^\d{4}-\d{2}-\d{2}$/.test(s)) throw new Error('Use YYYY-MM-DD dates.');
  var parts = s.split('-');
  var ms = Date.UTC(+parts[0], +parts[1] - 1, +parts[2]);
  if (carbonDateString(ms) !== s) throw new Error('Invalid calendar date: ' + s);
  return ms;
}
function carbonDateString(ms) { return new Date(ms).toISOString().slice(0, 10); }
function carbonFinite(x) { return typeof x === 'number' && isFinite(x); }
function carbonMonthKey(ms) {
  var d = new Date(ms);
  return d.getUTCFullYear() * 12 + d.getUTCMonth();
}
function carbonClamp(x) { return Math.max(0, Math.min(1, x)); }

/** Annual cover is a proxy for the NPP share, NOT a PFT productivity model.
 * The ratio of spatial mean covers allocates already ground-area-based NPP.
 * C4frac is deliberately not used: C3/C4 and woody/herbaceous are different.
 */
function carbonAllocation(context, p) {
  context = context || {};
  var bands = ['AFG', 'PFG', 'SHR', 'TRE'];
  var coversOK = bands.every(function(k) {
    return carbonFinite(context[k]) && context[k] >= 0 && context[k] <= 100;
  });
  var herb = coversOK ? context.AFG + context.PFG : null;
  var total = coversOK ? herb + context.SHR + context.TRE : null;
  var share = total > 0 ? herb / total : null;
  var matOK = context.temperatureMonths === 12 && carbonFinite(context.MAT);
  // Hui & Jackson grassland allocation, also used in the local Ecoexplorer:
  // fANPP = 0.171 + 0.0129 * MAT (degrees C). Bound the empirical fraction.
  var above = matOK ? carbonClamp(p.fANPP_intercept + p.fANPP_slope * context.MAT) : null;
  return {
    RAP_AFG_pct: coversOK ? context.AFG : null,
    RAP_PFG_pct: coversOK ? context.PFG : null,
    RAP_SHR_pct: coversOK ? context.SHR : null,
    RAP_TRE_pct: coversOK ? context.TRE : null,
    RAP_herbaceousCover: coversOK ? herb / 100 : null,
    herbaceousNPPFraction: share,
    shrubNPPFraction: total > 0 ? context.SHR / total : null,
    treeNPPFraction: total > 0 ? context.TRE / total : null, MAT: matOK ? context.MAT : null,
    fANPP: above, fBNPP: above === null ? null : 1 - above,
    rootShootRatio: above > 0 ? (1 - above) / above : null,
    carbonToBiomass: p.carbonToBiomass,
    carbonFractionDryMatter: 1 / p.carbonToBiomass,
    allocationStatus: !coversOK ? 'missing_RAP_cover' :
      !(total > 0) ? 'zero_vegetation_cover' : !matOK ? 'missing_full_year_temperature' : 'cover_proxy'
  };
}

/** Inputs: plain property dictionaries returned by Earth Engine evaluate().
 * Dates are inclusive. Monthly daily-carbon estimates are anchored on day 16,
 * matching the daylength anchor used by assembleForcing. Interpolation is on
 * daily carbon totals, not accumulated stocks. It is not constrained to conserve
 * the original monthly means. Edges are held only within an available month;
 * missing months stay missing and interpolation never bridges them.
 */
function buildCarbonProducts(monthlyRows, annualRows, startDate, endDate, overrides) {
  var p = mergeParams(ACCLIM_PARAMS, overrides || {});
  if (!carbonFinite(p.CUEplant) || p.CUEplant < 0 || p.CUEplant > 1)
    throw new Error('CUEplant must be between 0 and 1.');
  if (!carbonFinite(p.carbonToBiomass) || p.carbonToBiomass < 1)
    throw new Error('carbonToBiomass must be finite and at least 1.');
  var start = carbonDateMillis(startDate), end = carbonDateMillis(endDate);
  if (end < start) throw new Error('End date must be on or after start date.');
  var dayMs = CARBON_UNITS.dayMs;
  // Unlike RAP's scaled integer images, these values are already g C/m2.
  // g C -> kg C -> lb C, using the same conversion constants as RAP.
  var conversion = 0.001 * CARBON_UNITS.poundsPerKilogram * CARBON_UNITS.squareMetresPerAcre;
  var monthly = {}, annual = {}, allocations = {};
  annualRows.forEach(function(r) { annual[r.year] = r; });
  monthlyRows.forEach(function(r) {
    var ms = carbonDateMillis(r.date), key = carbonMonthKey(ms);
    if (monthly.hasOwnProperty(key)) throw new Error('Duplicate model month: ' + r.date);
    var d = new Date(ms);
    var valid = carbonFinite(r.PredGPPmean) && r.PredGPPmean >= 0 &&
      carbonFinite(r.daylength) && r.daylength >= 0 && r.daylength <= 24;
    monthly[key] = { weather_source:r.weather_source || 'gridded', anchor: Date.UTC(d.getUTCFullYear(), d.getUTCMonth(), 16),
      gpp: valid ? r.PredGPPmean * r.daylength * 3600 * CARBON_UNITS.gramsCPerMol / 1e6 : null };
  });
  var daily = [], windows = {}, periods = [];
  for (var ms = start; ms <= end; ms += dayMs) {
    var date = new Date(ms), year = date.getUTCFullYear(), key = carbonMonthKey(ms);
    var current = monthly[key], gpp = null, mode = 'missing_month';
    if (current && carbonFinite(current.gpp)) {
      var left = ms < current.anchor ? monthly[key - 1] : current;
      var right = ms < current.anchor ? current : monthly[key + 1];
      if (left && right && carbonFinite(left.gpp) && carbonFinite(right.gpp)) {
        var weight = (ms - left.anchor) / (right.anchor - left.anchor);
        gpp = left.gpp + weight * (right.gpp - left.gpp);
        mode = 'linear_monthly';
      } else {
        gpp = current.gpp;
        mode = 'held_within_month';
      }
    }
    if (!allocations[year]) allocations[year] = carbonAllocation(annual[year], p);
    var a = allocations[year];
    var npp = carbonFinite(gpp) ? gpp * p.CUEplant : null;
    var biomass = carbonFinite(npp) && a.allocationStatus === 'cover_proxy' ?
      npp * a.herbaceousNPPFraction * a.fANPP * p.carbonToBiomass * conversion : null;
    var row = mergeParams(a, {
      date: carbonDateString(ms), millis: ms, 'system:time_start': ms,
      year: year, weather_source:current ? current.weather_source : null, GPPdaily: gpp, NPPdaily: npp, CUEplant: p.CUEplant,
      // Total (above + belowground) PFT NPP; cover allocation needs no MAT.
      NPPgrass: carbonFinite(npp) && carbonFinite(a.herbaceousNPPFraction) ? npp * a.herbaceousNPPFraction : null,
      NPPshrub: carbonFinite(npp) && carbonFinite(a.shrubNPPFraction) ? npp * a.shrubNPPFraction : null,
      NPPtree: carbonFinite(npp) && carbonFinite(a.treeNPPFraction) ? npp * a.treeNPPFraction : null,
      NPPgrass_units: 'gC/m2/day', NPPshrub_units: 'gC/m2/day', NPPtree_units: 'gC/m2/day',
      HerbaceousAGBdaily: biomass, interpolation: mode,
      GPPdaily_units: 'gC/m2/day', NPPdaily_units: 'gC/m2/day',
      HerbaceousAGBdaily_units: 'lb dry matter/acre/day'
    });
    daily.push(row);
    // RAP API labels the first two windows Jan 16 and Feb 1 (DOY 16,32).
    // Non-overlapping Jan-1 anchored windows, final window ends Dec 31.
    var jan1 = Date.UTC(year, 0, 1), nextJan1 = Date.UTC(year + 1, 0, 1);
    var index = Math.floor((ms - jan1) / dayMs / 16);
    var windowStart = jan1 + index * 16 * dayMs;
    var windowEnd = Math.min(windowStart + 16 * dayMs, nextJan1);
    var windowKey = String(windowStart);
    if (!windows[windowKey]) {
      windows[windowKey] = { start: windowStart, end: windowEnd, year: year,
        index: index, requested: 0, valid: 0, biomassValid: 0, held: 0,
        gpp: 0, npp: 0, biomass: 0 };
      periods.push(windows[windowKey]);
    }
    var w = windows[windowKey];
    w.requested++;
    if (mode === 'held_within_month') w.held++;
    if (carbonFinite(npp)) { w.valid++; w.gpp += gpp; w.npp += npp; }
    if (carbonFinite(biomass)) { w.biomassValid++; w.biomass += biomass; }
  }
  var npp16day = periods.map(function(w) {
    var expected = (w.end - w.start) / dayMs;
    var complete = w.valid === expected;
    var biomassComplete = w.biomassValid === expected;
    return mergeParams(allocations[w.year], {
      date: carbonDateString(w.end - dayMs), year: w.year,
      millis: w.end - dayMs, 'system:time_start': w.start,
      period: w.index + 1, windowStart: carbonDateString(w.start),
      windowEnd: carbonDateString(w.end - dayMs),
      windowMidpoint: carbonDateString(w.start + (w.end - w.start) / 2),
      windowEndExclusive: carbonDateString(w.end),
      windowDays: expected, requestedDays: w.requested, validDays: w.valid,
      biomassValidDays: w.biomassValid, heldDays: w.held,
      completeWindow: complete ? 1 : 0, biomassCompleteWindow: biomassComplete ? 1 : 0,
      CUEplant: p.CUEplant,
      GPP16day_gC_m2: complete ? w.gpp : null,
      NPP16day_gC_m2: complete ? w.npp : null,
      NPP16day: complete ? w.npp * conversion : null,
      HerbaceousAGB16day: biomassComplete ? w.biomass : null,
      // Annual RAP cover and MAT allocation are constant within each calendar window.
      HerbaceousBGB16day: biomassComplete ? w.npp * conversion *
        allocations[w.year].herbaceousNPPFraction * allocations[w.year].fBNPP * p.carbonToBiomass : null,
      HerbaceousBGB16day_units: 'lb dry matter/acre/window',
      // Partial sums are explicitly separate; never present them as full totals.
      NPPpartial_lbC_acre: w.valid ? w.npp * conversion : null,
      HerbaceousAGBpartial_lb_acre: w.biomassValid ? w.biomass : null,
      NPP16day_units: 'lb C/acre/window',
      HerbaceousAGB16day_units: 'lb dry matter/acre/window'
    });
  });
  // Annual cumulative aboveground herbaceous carbon; reset before Jan 1 production.
  // Undo dry-matter and lb/acre conversions. This is production, not standing stock.
  // A partial first year starts at the selected start date. Missing days make
  // subsequent totals unknown only through the end of that calendar year.
  var cumulative = 0, cumulativeComplete = true, cumulativeByDate = {}, cumulativeYear = null;
  daily.forEach(function(r) {
    if (r.year !== cumulativeYear) {
      cumulativeYear = r.year;
      cumulative = 0;
      cumulativeComplete = true;
    }
    if (!carbonFinite(r.HerbaceousAGBdaily)) cumulativeComplete = false;
    if (cumulativeComplete) cumulative += r.HerbaceousAGBdaily / (r.carbonToBiomass * conversion);
    r.HerbaceousAGBcumulative_gC_m2 = cumulativeComplete ? cumulative : null;
    r.HerbaceousAGBcumulative_units = 'gC/m2';
    r.cumulativeStart = startDate > r.year + '-01-01' ? startDate : r.year + '-01-01';
    cumulativeByDate[r.date] = r.HerbaceousAGBcumulative_gC_m2;
  });
  npp16day.forEach(function(r) {
    r.HerbaceousAGBcumulative_gC_m2 = cumulativeByDate.hasOwnProperty(r.date) ? cumulativeByDate[r.date] : null;
    // Accumulation through the midpoint calendar date, for midpoint charting.
    r.HerbaceousAGBcumulative_mid_gC_m2 = cumulativeByDate.hasOwnProperty(r.windowMidpoint) ?
      cumulativeByDate[r.windowMidpoint] : null;
    r.HerbaceousAGBcumulative_units = 'gC/m2';
    r.cumulativeStart = startDate > r.year + '-01-01' ? startDate : r.year + '-01-01';
  });
  return { daily: daily, npp16day: npp16day };
}

/** Annual RAP v3 cover + complete-year TerraClimate MAT for biomass allocation.
 * Use matching years only. Missing cover or incomplete temperature => null
 * biomass, while GPP/NPP carbon outputs remain available. RAP spatial means
 * are sampled at 30m; the model's climate sampling scale is unchanged.
 */
function fetchCarbonAnnualContext(geom, startDate, endDate, config) {
  config = config || {};
  var startYear = Number(startDate.slice(0, 4)), endYear = Number(endDate.slice(0, 4));
  var rows = [];
  for (var y = startYear; y <= endYear; y++) {
    var yearStart = ee.Date.fromYMD(y, 1, 1), yearEnd = yearStart.advance(1, 'year');
    var coverIC = ee.ImageCollection(DATASETS.RAP_COVER).filterDate(yearStart, yearEnd);
    var blankCover = ee.Image.constant([0, 0, 0, 0]).rename(['AFG', 'PFG', 'SHR', 'TRE'])
      .updateMask(ee.Image.constant(0));
    var cover = ee.Image(ee.Algorithms.If(coverIC.size().gt(0), coverIC.first(), blankCover))
      .select(['AFG', 'PFG', 'SHR', 'TRE']).reduceRegion({
        reducer: ee.Reducer.mean(), geometry: geom, scale: 30, maxPixels: 1e10, tileScale: 4
      });
    if (config.weatherPlan) {
      rows.push(ee.Feature(null,cover.combine(ee.Dictionary(towerAnnualTemperature(config.weatherPlan,y))).combine(ee.Dictionary({year:y,coverImages:coverIC.size()}))));
      continue;
    }
    var temperature = ee.ImageCollection(DATASETS.TERRACLIMATE).filterDate(yearStart, yearEnd)
      .map(function(img) {
        return img.select('tmmn').add(img.select('tmmx')).multiply(0.05).rename('MAT');
      });
    var blankMAT = ee.Image.constant(0).rename('MAT').updateMask(ee.Image.constant(0));
    var matImage = ee.Image(ee.Algorithms.If(temperature.size().eq(12),
      temperature.mean().updateMask(temperature.count().eq(12)), blankMAT));
    var mat = matImage.reduceRegion({reducer: ee.Reducer.mean(), geometry: geom,
      scale: config.scale || DEFAULTS.sampleScale, bestEffort: true, maxPixels: 1e10});
    rows.push(ee.Feature(null, cover.combine(mat).combine(ee.Dictionary({
      year: y, temperatureMonths: temperature.size(), coverImages: coverIC.size()
    }))));
  }
  return ee.FeatureCollection(rows);
}

function carbonRowsToFC(rows) {
  return ee.FeatureCollection(rows.map(function(r) { return ee.Feature(null, r); }));
}

// SoilGrids 2.0 sand: g/kg fine earth, 250 m; divide by 10 for percent.
// Thickness-weight the 0-5, 5-15 and 15-30 cm mean predictions (5:10:15).
// Arithmetic intersects depth masks, so incomplete profiles are not filled.
// https://docs.isric.org/globaldata/soilgrids/access_on_gee.html
function fetchSoilGridsSand(geom) {
  var sand = ee.Image('projects/soilgrids-isric/sand_mean')
    .select(['sand_0-5cm_mean', 'sand_5-15cm_mean', 'sand_15-30cm_mean']);
  sand = sand.updateMask(sand.gte(0).and(sand.lte(1000)));
  var sand030 = sand.select('sand_0-5cm_mean').multiply(5)
    .add(sand.select('sand_5-15cm_mean').multiply(10))
    .add(sand.select('sand_15-30cm_mean').multiply(15))
    .divide(30).divide(10).rename('sand_0_30cm_pct');
  // Point: local pixel. Polygon: spatial mean of the depth-weighted pixels.
  return sand030.reduceRegion({reducer: ee.Reducer.mean(), geometry: geom,
    crs: sand.select('sand_0-5cm_mean').projection(), scale: 250,
    maxPixels: 1e10, tileScale: 4});
}

// User-supplied Relative AWC versus sand (% by mass) plot.
// The plotted line is 1 at 0% sand and 0.25 at 100% sand:
// RAWHC = 1 - 0.0075*sandPercent = 1 - 0.75*sandFraction.
// The figure legend's additional /100 is inconsistent with the plotted line.
// Dimensionless relative AWC, not absolute theta_FC-theta_WP or water depth.
// Display-only context; this does not modify the fitted evaporative bucket.
function relativeAWCFromSand(sandPercent) {
  return carbonFinite(sandPercent) && sandPercent >= 0 && sandPercent <= 100 ?
    1 - 0.0075 * sandPercent : null;
}

function addSoilMoistureChart(output, geom, requestId, regionName, info) {
  var baseTitle = regionChartTitle('Soil Moisture - Near Surface', regionName);
  var options = {title: baseTitle, chartArea: {top: 100},
    hAxis: {title: 'Date'},
    vAxes: {0: {title: 'Extractable soil moisture (mm)', viewWindow: {min: 0}},
                1: {title: 'Precipitation (mm)', viewWindow: {min: 0},
                    gridlines: {color: 'transparent'}, minorGridlines: {color: 'transparent'}}},
    seriesType: 'line',
    series: {0: {targetAxisIndex: 0, color: '#8D6E63'},
             1: {targetAxisIndex: 0, color: '#757575', lineDashStyle: [4, 4]},
             2: {targetAxisIndex: 1, type: 'bars', color: '#64B5F6'}},
    lineWidth: 1.5, pointSize: 0};
  var chart = ui.Chart.feature.byFeature(output, 'date', ['SM_mean', 'SMmax', 'precip'])
    .setChartType('ComboChart').setOptions(options);
  addResultChart(chart);
  options.title = baseTitle + '\n' + depthChartLabel(info, 'effective_soil_evaporation_depth_cm', 'Effective source depth');
  options.chartArea.top = 110;
  chart.setOptions(options);
  addChartNote(soilChartFooter(info, 'SM_r2'));
}

// Aligned Celsius/Fahrenheit axes: identical physical bounds and tick positions.
// Annual climate summaries use 24-hour air temperature, not daytime Tair.
// Sum daily carbon fluxes without annualizing partial years or filling missing days.
function annualCarbonNote(rows, fields, labels) {
  var years = {};
  rows.forEach(function(r) {
    var year = r.date.slice(0, 4);
    if (!years[year]) years[year] = {};
    years[year][r.date] = r;
  });
  return fields.map(function(field, index) {
    var totals = Object.keys(years).sort().map(function(year) {
      var sum = 0, complete = true;
      var end = Date.UTC(Number(year) + 1, 0, 1);
      for (var ms = Date.UTC(Number(year), 0, 1); ms < end; ms += 86400000) {
        var key = new Date(ms).toISOString().slice(0, 10);
        var row = years[year][key];
        if (!row || !carbonFinite(row[field])) complete = false;
        else sum += row[field];
      }
      return year + ': ' + (complete ? sum.toFixed(1) : 'unavailable (incomplete year/data)');
    });
    return 'Annual ' + labels[index] + ' (gC/m²/yr): ' + totals.join('; ');
  }).join('\n');
}

function annualClimateNote(rows, start, end, variable) {
  var years = {};
  rows.forEach(function(row) {
    var year = row.date.slice(0, 4);
    if (!years[year]) years[year] = {};
    years[year][row.date.slice(5, 7)] = row;
  });
  var values = Object.keys(years).sort().map(function(year) {
    var months = years[year], sum = 0, days = 0;
    var complete = start <= year + '-01-01' && end >= year + '-12-31';
    for (var m = 1; m <= 12; m++) {
      var row = months[('0' + m).slice(-2)];
      var n = new Date(Date.UTC(Number(year), m, 0)).getUTCDate();
      if (!row || !carbonFinite(row[variable])) { complete = false; continue; }
      sum += row[variable] * (variable === 'Tdaily' ? n : 1);
      days += n;
    }
    return year + ': ' + (complete ? (variable === 'Tdaily' ? sum / days : sum).toFixed(1) +
      (variable === 'Tdaily' ? ' °C' : ' mm/year') : 'unavailable (incomplete year/data)');
  });
  return (variable === 'Tdaily' ? 'Mean annual air temperature (24-hour): ' : 'Annual precipitation: ') + values.join('; ');
}

function addWeatherCharts(output, rows, regionName) {
  addResultChart(ui.Chart.feature.byFeature(output, 'date', ['VPD_kPa', 'WUE'])
    .setSeriesNames(['VPD', 'WUE'])
    .setOptions({title: regionChartTitle('VPD & WUE', regionName), chartArea: {top: 80},
      hAxis: {title: 'Date'},
      vAxes: {0: {title: 'VPD (kPa)'}, 1: {title: 'WUE (mmol/mol)', viewWindow: {min: 0},
        gridlines: {color: 'transparent'}, minorGridlines: {color: 'transparent'}}},
      series: {0: {targetAxisIndex: 0, color: '#EF6C00'}, 1: {targetAxisIndex: 1, color: '#1976D2', lineDashStyle: [6, 4]}},
      lineWidth: 1.5, pointSize: 0, interpolateNulls: false}));
  // Paired monthly values; omit missing pairs rather than manufacturing zeros.
  function scatter(metric) {
    var label = metric === 'WUE' ? 'WUE (mmol/mol)' : 'gs (mol/m2/s)';
    var pairs = rows.filter(function(r) { return carbonFinite(r.VPD_kPa) && carbonFinite(r[metric]); });
    var table = {cols: [{id: 'vpd', label: 'VPD (kPa)', type: 'number'},
      {id: metric, label: label, type: 'number'}],
      rows: pairs.map(function(r) { return {c: [{v: r.VPD_kPa}, {v: r[metric]}]}; })};
    return ui.Chart(table).setChartType('ScatterChart').setOptions({
      title: regionChartTitle(metric + ' vs VPD', regionName), chartArea: {top: 80},
      hAxis: {title: 'VPD (kPa)'}, vAxis: {title: label, viewWindow: {min: 0}},
      legend: {position: 'none'}, colors: [metric === 'WUE' ? '#1976D2' : '#2E7D32'],
      pointSize: 4, pointShape: 'circle', lineWidth: 0
    });
  }
  addResultChart(scatter('WUE'));
  var replaceScatter = chartReplacer();
  UI._chartCard.add(ui.Select({items: ['WUE', 'gs'], value: 'WUE',
    onChange: function(metric) { replaceScatter(scatter(metric)); }, style: {margin: '4px 8px'}}));
}

// BEGIN CARBON ASYNC UI
function loadCarbonProducts(output, geom, start, end, requestId, regionName, c4Chart) {
  function ready(texture, error) {
    if (requestId !== UI._requestId) return;
    loadCarbonProductsWithTexture(output, geom, start, end, requestId, regionName,
      error || !texture ? null : texture.sand_0_30cm_pct, c4Chart);
  }
  // Texture failure leaves depths missing but does not suppress other outputs.
  try { fetchSoilGridsSand(geom).evaluate(ready); } catch (error) { ready(null, error); }
}
function loadCarbonProductsWithTexture(output, geom, start, end, requestId, regionName, sandPercent, c4Chart) {
  var context = fetchCarbonAnnualContext(geom, start, end, {weatherPlan:UI._runWeatherPlan});
  var monthlyProperties = output.toList(output.size()).map(function(f) { return ee.Feature(f).toDictionary(); });
  var annualProperties = context.toList(context.size()).map(function(f) { return ee.Feature(f).toDictionary(); });
  ee.Dictionary({monthly: monthlyProperties, annual: annualProperties}).evaluate(function(result, error) {
    if (requestId !== UI._requestId) return;
    UI.runButton.setDisabled(false);
    if (error) { setStatus('Carbon products could not be calculated: ' + error); return; }
    if (!result.monthly || !result.monthly.length) {
      setStatus('No usable forcing data for this location and date range. MODIS vegetation samples may all be masked. Check the coordinates or choose a different date range.');
      return;
    }
    try {
      var soilR2 = bucketR2(result.monthly, 'PredLEsday_mean', 'PredLEsday_RHVPD');
      var etR2 = bucketR2(result.monthly, 'PredLE_ETeffective', 'PredLE_ETreference');
      result.monthly.forEach(function(row) {
        row.SM_r2 = carbonFinite(row.SMmax) ? soilR2 : null;
        row.SM_et_r2 = carbonFinite(row.SMmax_et) ? etR2 : null;
        var depths = effectiveDepthFields(row, sandPercent);
        Object.keys(depths).forEach(function(key) { row[key] = depths[key]; });
        var tDepth = transpirationSourceDepth([row]);
        Object.keys(tDepth).forEach(function(key) { row[key] = tDepth[key]; });
      });
      output = carbonRowsToFC(result.monthly);
      UI._lastOutput = output; // monthly CSV includes AWC and ET, E, and T source depths
      var bucketInfo = {};
      var bucketRow=result.monthly.filter(function(r){return carbonFinite(r.PredGPPmean);})[0];
      if (bucketRow) Object.keys(bucketRow).forEach(function(key) {
        bucketInfo[key] = bucketRow[key];
      });
      var selectedPeriodDepth = transpirationSourceDepth(result.monthly);
      Object.keys(selectedPeriodDepth).forEach(function(key) {
        bucketInfo[key] = selectedPeriodDepth[key];
      });
      var products = buildCarbonProducts(result.monthly, result.annual, start, end);
      UI._lastDaily = carbonRowsToFC(products.daily);
      if (UI._annualGPPNote) UI._annualGPPNote.setValue(annualCarbonNote(products.daily, ['GPPdaily'], ['GPP']));
      UI._last16day = carbonRowsToFC(products.npp16day);
      addOptionalSeriesChart({rows: products.daily},
        ['NPPdaily', 'NPPgrass', 'NPPshrub', 'NPPtree'],
        ['NPPeco', 'NPPgrass', 'NPPshrub', 'NPPtree'],
        ['#000000', '#2E7D32', '#EF6C00', '#1976D2'], [1, 2, 3], {
          title: regionChartTitle('Daily NPP', regionName), chartArea: {top: 80},
          hAxis: {title: 'Date'}, vAxis: {title: 'g C/m2/day'}, lineWidth: 1.5, pointSize: 0
        }, annualCarbonNote);
      var hasBiomass = products.npp16day.some(function(r) { return carbonFinite(r.HerbaceousAGB16day) || carbonFinite(r.HerbaceousBGB16day); });
      if (hasBiomass) {
        // Chart-only Jan 1 anchors show zero before that day's production.
        // Daily and 16-day export totals retain their end-of-period meaning.
        var annualChartRows = products.npp16day.map(function(r) {
          // Above- and belowground herbaceous production share dry-matter units.
          return {date: r.windowMidpoint, HerbaceousAGB16day: r.HerbaceousAGB16day,
            HerbaceousBGB16day: r.HerbaceousBGB16day,
            HerbaceousAGBcumulative_gC_m2: r.HerbaceousAGBcumulative_mid_gC_m2};
        });
        products.daily.forEach(function(r) {
          if (r.date.slice(5) === '01-01') annualChartRows.push({date: r.date,
            HerbaceousAGB16day: null, HerbaceousAGBcumulative_gC_m2: 0});
        });
        annualChartRows.sort(function(a, b) { return carbonDateMillis(a.date) - carbonDateMillis(b.date); });
        addOptionalSeriesChart({rows: annualChartRows},
          ['HerbaceousAGB16day', 'HerbaceousAGBcumulative_gC_m2', 'HerbaceousBGB16day'],
          ['ANPP', 'Aboveground Biomass', 'BNPP'], ['#2E7D32', '#2E7D32', '#000000'], [2], {
            title: regionChartTitle('Herbaceous Biomass Production', regionName), chartArea: {top: 80},
            hAxis: {title: 'Window midpoint date'},
            vAxes: {0: {title: 'Dry Matter Production (lbs/acre/16days)'},
              1: {title: 'Aboveground Biomass (gC/m2)', gridlines: {color: 'transparent'},
                minorGridlines: {color: 'transparent'}}},
            series: {0: {targetAxisIndex: 0, color: '#2E7D32'},
              1: {targetAxisIndex: 1, color: '#2E7D32', lineWidth: 0, pointSize: 4},
              2: {targetAxisIndex: 0, color: '#000000'}},
            lineWidth: 1.5, pointSize: 2
          });
        var dryBiomassLbAcrePerGCm2 = 0.001 * CARBON_UNITS.poundsPerKilogram *
          CARBON_UNITS.squareMetresPerAcre * ACCLIM_PARAMS.carbonToBiomass;
        addChartNote(ui.Label('1 gC/m² = ' +
          dryBiomassLbAcrePerGCm2.toFixed(2) + ' lb DM/acre (DM=47%C).'));
//        UI.chartHolder.add(ui.Label('Missing data make totals unknown ' +
//          'for the rest of that year. This is accumulated production, not standing biomass.'));
      }
      addETMoistureCharts(output, bucketInfo, regionName);
      if (carbonFinite(bucketInfo.SMmax_et)) addChartNote(ui.Label(annualClimateNote(result.monthly, start, end, 'precip')));
      if (carbonFinite(bucketInfo.SMmax)) {
        addSoilMoistureChart(output, geom, requestId, regionName, bucketInfo);
      }
      c4Chart(result.monthly);
      if (bucketInfo.FAPAR_method_requested === 'EVI2') {
        addChartNote(ui.Label(bucketInfo.FAPAR_method === 'EVI2'
          ? 'FAPAR = EVI2 scaled between EVI2min=' + bucketInfo.FAPAR_EVI2min.toFixed(3) +
            ' and EVI2max=' + bucketInfo.FAPAR_EVI2max.toFixed(3)
          : 'FAPAR: NDVI² fallback. EVI2 fit needs at least 3 paired months and a positive slope.'));
      }
      addWeatherCharts(output, result.monthly, regionName);
      var valid = products.daily.filter(function(r) { return carbonFinite(r.NPPdaily); }).length;
      var biomassDays = products.daily.filter(function(r) { return carbonFinite(r.HerbaceousAGBdaily); }).length;
      setStatus('Ready: ' + valid + '/' + products.daily.length + ' days have GPP/NPP; ' +
        biomassDays + ' have biomass estimates. Incomplete windows have blank full totals in CSV.');
    } catch (e) { setStatus('Carbon products error: ' + e.message); }
  });
}
// END CARBON ASYNC UI


var UI = {};  // namespace

UI.mainPanel = ui.Panel({
  layout: ui.Panel.Layout.flow('vertical'),
  style: { stretch: 'horizontal', padding: '10px', backgroundColor: '#f5f5f5' }
});

UI.title = ui.Label({
  value: 'ACCLIM (v2)',
  style: { fontSize: '22px', fontWeight: 'bold', color: '#2E7D32',
           margin: '6px 0px', textAlign: 'center' }
});

UI.subtitle = ui.Label({
  value: '"Assimilation Capacity Coordinates with Light Interception Model"\n' +
         'A dynamic physiology model of photosynthesis and evapotranspiration driven with global satellite observations',
  style: { fontSize: '12px', color: '#666', textAlign: 'center', whiteSpace: 'pre-wrap',
           margin: '0px 0px 12px 0px' }
});

UI.inputPanel = ui.Panel({
  layout: ui.Panel.Layout.flow('vertical'),
  style: { border: '1px solid #ddd', padding: '12px',
           backgroundColor: 'white', margin: '8px 0px' }
});

function makeLabeledInput(labelText, defaultVal) {
  var label = ui.Label({
    value: labelText,
    style: { fontSize: '13px', fontWeight: 'bold', margin: '8px 0px 2px 0px' }
  });
  var input = ui.Textbox({
    value: String(defaultVal),
    style: { width: '100%', margin: '2px 0px' }
  });
  return { label: label, input: input };
}

// --- Geometry inputs: asset dropdown + drawing tools (added below to the map) ---
// The asset dropdown is populated synchronously from the ecoexplorer asset
// folder. ee.data.listAssets is a blocking client call, so we can build the
// items list at script-load time.
UI.ASSET_FOLDER = 'projects/ee-landflux/assets/public';
UI.assetList = [];   // [{name, assetId}, ...]
try {
  var rawAssets = ee.data.listAssets(UI.ASSET_FOLDER).assets || [];
  for (var ai = 0; ai < rawAssets.length; ai++) {
    var aid = rawAssets[ai].id;
    UI.assetList.push({ name: aid.split('/').pop(), assetId: aid });
  }
} catch (e) {
  // If the folder is unreadable from the running account, just leave the
  // dropdown empty — drawn geometries still work.
  print('Warning: could not list assets in ' + UI.ASSET_FOLDER + ': ' + e.message);
}

UI.assetLabel = ui.Label({
  value: 'Pre-loaded geometry asset (optional)',
  style: { fontSize: '13px', fontWeight: 'bold', margin: '8px 0px 2px 0px' }
});
UI.assetSelector = ui.Select({
  items: UI.assetList.map(function(a) { return a.name; }),
  placeholder: UI.assetList.length ? 'Choose asset…' : '(no assets found)',
  style: { width: '100%', margin: '2px 0px' }
});

UI.uploadLink = ui.Label({
  value: '↗ Upload a polygon (KML / shapefile / GeoJSON)',
  style: { fontSize: '12px', color: '#1976D2', margin: '4px 0px 8px 0px' },
  targetUrl: 'https://ecoexplorer-uploader-public.streamlit.app/'
});

UI.geomHelp = ui.Label({
  value: 'Or use the map: draw a point, draw a polygon, or pick an asset above.',
  style: { fontSize: '11px', color: '#666', margin: '0px 0px 6px 0px' }
});

// --- Date range + optional site name ---
UI.startI     = makeLabeledInput('Start date (YYYY-MM-DD)', DEFAULTS.startDate);
UI.endI       = makeLabeledInput('End date (YYYY-MM-DD)',   DEFAULTS.endDate);

UI.inputPanel.add(UI.assetLabel);
UI.inputPanel.add(UI.assetSelector);
UI.inputPanel.add(UI.uploadLink);
UI.inputPanel.add(UI.geomHelp);
UI.useCoordinates = ui.Checkbox({label: 'Use latitude / longitude instead of map or asset', value: false});
UI.latitude = ui.Textbox({placeholder: 'Latitude (−90 to 90)', style: {stretch: 'horizontal'}});
UI.longitude = ui.Textbox({placeholder: 'Longitude (−180 to 180)', style: {stretch: 'horizontal'}});
UI.inputPanel.add(UI.useCoordinates);
UI.inputPanel.add(ui.Panel({widgets: [UI.latitude, UI.longitude],
  layout: ui.Panel.Layout.flow('horizontal'), style: {stretch: 'horizontal'}}));

[UI.startI, UI.endI].forEach(function(pair) {
  UI.inputPanel.add(pair.label);
  UI.inputPanel.add(pair.input);
});

// NOTE: ui.Label has no onClick in GEE. Buttons must be ui.Button. The
// styling options on ui.Button are limited to: stretch, shown, position,
// padding, margin, height, width — no backgroundColor / color / fontSize /
// textAlign. We get GEE's default button chrome in exchange for working
// click handlers. (onRun/onReset are function declarations defined below,
// so hoisting makes the forward reference safe.)
UI.runButton = ui.Button({
  label: 'Run ACCLIM',
  onClick: function() { onRun(); },
  style: { width: '140px', height: '30px', padding: '0px',
           margin: '6px 0px 0px 0px' }
});

UI.faparSelector = ui.Select({items: [
  {label: 'FAPAR = NDVI²', value: 'NDVI'},
  {label: 'FAPAR = scaled EVI2', value: 'EVI2'}], value: 'EVI2',
  style: {stretch: 'horizontal', margin: '4px 0px'}});
UI.inputPanel.add(ui.Label('FAPAR method (applied on Run ACCLIM)'));
UI.inputPanel.add(UI.faparSelector);
setupTowerImport();
UI.inputPanel.add(UI.runButton);

UI.statusPanel = ui.Panel({
  layout: ui.Panel.Layout.flow('vertical'),
  style: { border: '1px solid #ddd', padding: '10px',
           backgroundColor: '#f0f8ff', margin: '8px 0px', shown: false }
});
UI.statusLabel = ui.Label({
  value: 'Ready.',
  style: { fontSize: '12px', color: '#333' }
});
UI.statusPanel.add(UI.statusLabel);

UI.resultsPanel = ui.Panel({
  layout: ui.Panel.Layout.flow('vertical'),
  style: { border: '1px solid #ddd', padding: '10px',
           backgroundColor: 'white', margin: '8px 0px', shown: false }
});

UI.resultsTitle = ui.Label({
  value: 'Results',
  style: { fontSize: '16px', fontWeight: 'bold', margin: '0px 0px 8px 0px' }
});
// Results title is provided by the always-visible floating-panel header.

UI.chartHolder = ui.Panel({
  layout: ui.Panel.Layout.flow('vertical'),
  style: { margin: '8px 0px' }
});
UI.resultsPanel.add(UI.chartHolder);

// Keep chart notes inside their card; fill rows left-to-right as results arrive.
UI._resultCharts = [];
UI._chartCount = 0;
UI._chartRow = null;
UI._chartCard = null;
function addResultChart(chart) {
  UI._resultCharts.push(chart);
  if (UI._chartCount % 2 === 0) {
    UI._chartRow = ui.Panel({layout: ui.Panel.Layout.flow('horizontal'),
      style: {stretch: 'horizontal', margin: '0px', padding: '0px'}});
    UI.chartHolder.add(UI._chartRow);
  }
  chart.style().set({width: '100%', height: '380px', margin: '0px'});
  UI._chartCard = ui.Panel({widgets: [chart], layout: ui.Panel.Layout.flow('vertical'),
    style: {width: '49%', margin: '4px', padding: '0px'}});
  UI._chartRow.add(UI._chartCard);
  UI._chartCount++;
}
// Keep replacement charts registered for redraw after Results expands.
function chartReplacer() {
  var card = UI._chartCard, index = UI._resultCharts.length - 1;
  var requestId = UI._requestId;
  return function(chart) {
    if (requestId !== UI._requestId) return;
    chart.style().set({width: '100%', height: '380px', margin: '0px'});
    card.widgets().set(0, chart);
    UI._resultCharts[index] = chart;
  };
}
// Optional chart series use the monthly values already fetched for carbon products.
// No server-side collection or evaluate call is used when the selection changes.
function addOptionalSeriesChart(cache, properties, names, colors, optional, options, annualFooter) {
  var choices = [{label: 'Optional', value: '0'}];
  for (var mask = 1; mask < Math.pow(2, optional.length); mask++) {
    var labels = [];
    optional.forEach(function(index, bit) { if (mask & (1 << bit)) labels.push(names[index]); });
    choices.push({label: 'Add ' + labels.join(' + '), value: String(mask)});
  }
  var selection = 0;
  var annualNote = annualFooter ? ui.Label({value: '', style: {whiteSpace: 'pre-wrap'}}) : null;
  function build() {
    var selected = [], palette = [], seriesOptions = {};
    var cols = [{id: 'date', label: 'Date', type: 'date'}];
    properties.forEach(function(property, index) {
      var bit = optional.indexOf(index);
      if (bit < 0 || (selection & (1 << bit))) {
        if (options.series && options.series[index]) seriesOptions[selected.length] = options.series[index];
        selected.push(property); palette.push(colors[index]);
        cols.push({id: property, label: names[index], type: 'number'});
      }
    });
    var rows = (cache.rows || []).map(function(row) {
      var parts = row.date.split('-');
      var cells = [{v: new Date(Number(parts[0]), Number(parts[1]) - 1, Number(parts[2]))}];
      selected.forEach(function(property) {
        cells.push({v: carbonFinite(row[property]) ? row[property] : null});
      });
      return {c: cells};
    });
    if (annualNote) annualNote.setValue(annualFooter(cache.rows || [], selected, selected.map(function(p) { return names[properties.indexOf(p)]; })));
    var chartOptions = {};
    Object.keys(options).forEach(function(key) { chartOptions[key] = options[key]; });
    chartOptions.colors = palette;
    if (options.series) chartOptions.series = seriesOptions;
    chartOptions.interpolateNulls = false;
    return ui.Chart({cols: cols, rows: rows}).setChartType('LineChart').setOptions(chartOptions);
  }
  addResultChart(build());
  var redraw = chartReplacer();
  var selector = ui.Select({items: choices, value: '0', disabled: !cache.rows,
    onChange: function(value) { selection = Number(value); redraw(build()); },
    style: {margin: '4px 8px', stretch: 'horizontal'}});
  UI._chartCard.add(selector);
  if (annualNote) addChartNote(annualNote);
  if (!cache.rows) cache.ready.push(function() { redraw(build()); selector.setDisabled(false); });
  if (cache.refresh) cache.refresh.push(function() { redraw(build()); });
}

// Tower monthly values are prepared locally with tower-csv-importer.html.
// No raw tower records or metadata are uploaded to Earth Engine assets.
function parseTowerMonthly(text) {
  var data = JSON.parse(text), seen = {};
  if ((data.version !== 2 && data.version !== 3) || data.units !== 'micromol/m2/s' ||
      data.filter !== 'INCOMING_SHORTWAVE_RADIATION > 0' ||
      !Array.isArray(data.daily) || !data.daily.length || data.daily.length > 40000) {
    throw new Error('Use the updated tower-csv-importer.html output.');
  }
  data.daily.forEach(function(row) {
    if (!/^\d{4}-\d{2}-\d{2}$/.test(row.date) ||
        !isFinite(Date.parse(row.date)) || new Date(row.date).toISOString().slice(0,10) !== row.date || seen[row.date] ||
        !carbonFinite(row.sum) || !carbonFinite(row.n) || row.n < 0 || row.n % 1 !== 0) {
      throw new Error('Invalid or duplicate daily tower summary.');
    }
    if(data.version===3) {
      ['daylight','total','unknownLight','absent','secondN'].forEach(function(k){if(!carbonFinite(row[k])||row[k]<0||row[k]%1!==0)throw new Error('Invalid coverage count.');});
      if(!carbonFinite(row.secondSum)||row.n>row.daylight||row.secondN>row.daylight||row.daylight+row.unknownLight>row.total||(!row.n&&row.sum!==0)||(!row.secondN&&row.secondSum!==0))throw new Error('Inconsistent tower coverage.');
    }
    if(data.leVariable==='LE_CORR'&&(!carbonFinite(row.leSum)||!carbonFinite(row.leN)||row.leN<0||row.leN%1!==0||row.leN>row.total||(!row.leN&&row.leSum!==0)))throw new Error('Invalid LE_CORR summary.');
    seen[row.date] = true;
  });
  data.daily.sort(function(a,b) { return a.date < b.date ? -1 : 1; });
  var valid=data.daily.filter(function(row){return row.n>0 || row.secondN>0;});
  if(!valid.length)throw new Error('No valid daytime GPP.');
  data.start=valid[0].date;data.end=valid[valid.length-1].date;
  data.site = String(data.site || 'Tower CSV');
  if (data.latitude !== null && (!carbonFinite(data.latitude) || Math.abs(data.latitude) > 90)) throw new Error('Invalid latitude.');
  if (data.longitude !== null && (!carbonFinite(data.longitude) || Math.abs(data.longitude) > 180)) throw new Error('Invalid longitude.');
  return data;
}
function towerMonthsInRange(tower, start, end) {
  var months = {};
  tower.daily.forEach(function(day) {
    if (day.date < start || day.date > end) return;
    var key = day.date.slice(0,7), m = months[key] || (months[key] = {leSum:0,leN:0,sum:0,n:0,days:0,secondSum:0,secondN:0,daylight:0,total:0,unknownLight:0,absent:0});
    ['sum','n','secondSum','secondN','daylight','total','unknownLight','absent','leSum','leN'].forEach(function(k){m[k] += day[k] || 0;});
    if(tower.version === 2 || day.total > 0) m.days++;
  });
  Object.keys(months).forEach(function(key) {
    var m=months[key]; m.gpp=m.n ? m.sum/m.n : null; m.secondGpp=m.secondN ? m.secondSum/m.secondN : null;
    m.leET=tower.leVariable==='LE_CORR'&&m.leN ? m.leSum/m.leN*86400/2450000 : null;
    m.leMissingPct=tower.leVariable==='LE_CORR'&&m.total ? 100*(m.total-m.leN)/m.total : null;
    m.missingPct=tower.version===3 && m.daylight ? 100*(m.daylight-m.n)/m.daylight : null;
    m.secondMissingPct=tower.version===3 && m.daylight ? 100*(m.daylight-m.secondN)/m.daylight : null;
  });
  return months;
}
// Validate imported monthly weather independently of the summary's complete flag.
function towerWeatherPlan(tower, start, end) {
  if (!tower || !tower.weather || tower.weather.version !== 2)
    throw new Error('Reprocess the original tower CSV/ZIP with the updated importer to include P_F, PA_F and LE_CORR + H_CORR weather inputs.');
  var seen = {}, byMonth = {}, fields = ['Tmin','Tmax','Tdaily','SWRdaily','ea','Pdaily','precip','Rndaily'];
  var counts = ['ta','sw','ea','pa','rain','rn'];
  if (!Array.isArray(tower.weather.months)) throw new Error('Invalid tower weather summary.');
  tower.weather.months.forEach(function(m) {
    if (!/^\d{4}-\d{2}$/.test(m.month) || seen[m.month]) throw new Error('Invalid or duplicate weather month.');
    var date = new Date(m.month + '-01T00:00:00Z');
    if (!isFinite(date.getTime()) || date.toISOString().slice(0,7) !== m.month) throw new Error('Invalid weather month.');
    seen[m.month] = true;
    var days = new Date(Date.UTC(date.getUTCFullYear(),date.getUTCMonth()+1,0)).getUTCDate();
    var interval = tower.intervalMinutes;
    if (interval !== 30 && interval !== 60) throw new Error('Unsupported weather interval.');
    var expected = days*1440/interval;
    var valid = m.complete === true && m.expected === expected && fields.every(function(k){return carbonFinite(m[k]);}) &&
      counts.every(function(k){return m.counts && m.counts[k] === expected;}) &&
      m.Tmin <= m.Tdaily && m.Tdaily <= m.Tmax && m.ea >= 0 && m.Pdaily > 0 && m.precip >= 0;
    byMonth[m.month] = {data:m,valid:valid};
  });
  var segments=[], missing=[], rejected=[], current=null;
  var first = new Date(start.slice(0,7)+'-01T00:00:00Z'), last=end.slice(0,7);
  for(var d=first;d.toISOString().slice(0,7)<=last;d=new Date(Date.UTC(d.getUTCFullYear(),d.getUTCMonth()+1,1))){
    var key=d.toISOString().slice(0,7), entry=byMonth[key];
    if(entry && entry.valid){if(!current){current=[];segments.push(current);}current.push(entry.data);}
    else {
      missing.push(key);current=null;
      var reasons=[],labels={ta:'temperature',sw:'shortwave',ea:'humidity',pa:'pressure',rain:'precipitation',rn:'LE_CORR + H_CORR'};
      if(!entry)reasons.push('no weather records');
      else {
        counts.forEach(function(k){var n=entry.data.counts&&entry.data.counts[k];if(n!==entry.data.expected)reasons.push(labels[k]+': '+(carbonFinite(n)?n:0)+'/'+entry.data.expected+' valid intervals');});
        if(!reasons.length)reasons.push('invalid monthly values or incomplete daily temperature extrema');
      }
      rejected.push(key+' — '+reasons.join('; '));
    }
  }
  return {segments:segments,missing:missing,rejected:rejected,byMonth:byMonth,source:tower.weather.source||{}};
}
function towerWeatherFeature(row, geom) {
  var d=ee.Date(row.month+'-01');
  return ee.Feature(geom,{date:row.month+'-01',millis:d.millis(),year:d.get('year'),month:d.get('month'),
    Tmin:row.Tmin,Tmax:row.Tmax,Tdaily:row.Tdaily,SWRdaily:row.SWRdaily,ea:row.ea,
    Pdaily:row.Pdaily,precip:row.precip,Rndaily:row.Rndaily,period_days:d.advance(1,'month').difference(d,'day'),
    weather_source:'tower',weather_energy_basis:'LE_CORR + H_CORR',weather_status:'complete',weather_expected_intervals:row.expected});
}
function towerAnnualTemperature(plan, year) {
  var sum=0,days=0,count=0;
  for(var m=1;m<=12;m++){
    var key=year+'-'+(m<10?'0':'')+m, entry=plan.byMonth[key], row=entry&&entry.data;
    if(row && carbonFinite(row.Tdaily) && row.counts && row.counts.ta===row.expected){
      var n=new Date(Date.UTC(year,m,0)).getUTCDate();sum+=row.Tdaily*n;days+=n;count++;
    }
  }
  return {MAT:count===12?sum/days:null,temperatureMonths:count,MAT_source:'tower'};
}

function setupTowerImport() {
  UI.inputPanel.add(ui.Label({value: 'Open flux tower CSV importer',
    targetUrl: 'https://kevinptu.github.io/ACCLIM/tower-csv-importer.html', style: {color: '#1565C0'}}));
  UI.towerInput = ui.Textbox({placeholder: 'Tower CSV data: paste importer output', style: {stretch:'horizontal'}});
  UI.towerNote = ui.Label({value:'Prepare your CSV with tower-csv-importer.html. Import before running ACCLIM.', style:{whiteSpace:'pre-wrap'}});
  function load(defaults) {
    try {
      var data = parseTowerMonthly(UI.towerInput.getValue());
      if (defaults && (data.latitude === null || data.longitude === null)) throw new Error('CSV coordinates missing. Enter coordinates manually and use Load with current settings.');
      UI._towerData = data;
      if (defaults) {
        UI.latitude.setValue(String(data.latitude)); UI.longitude.setValue(String(data.longitude));
        UI.useCoordinates.setValue(true);
        UI.startI.input.setValue(data.start); UI.endI.input.setValue(data.end);
      }
      UI.towerNote.setValue(data.site + ': ' + data.start + ' to ' + data.end +
        '. Review/edit coordinates and dates, then Run ACCLIM. Existing results are unchanged until you run again.');
    } catch(error) { UI.towerNote.setValue('Tower import error: ' + error.message); }
  }
  UI.inputPanel.add(UI.towerInput);
  UI.inputPanel.add(ui.Button({label:'Use CSV coordinates and dates', onClick:function(){load(true);}}));
  UI.inputPanel.add(ui.Button({label:'Load with current settings', onClick:function(){load(false);}}));
  UI.inputPanel.add(ui.Button({label:'Clear tower import',onClick:function(){UI._towerData=null;UI.towerInput.setValue('');UI.towerNote.setValue('Tower import cleared. Run again to update the chart.');}}));
  UI.inputPanel.add(UI.towerNote);
  UI.useTowerWeather=ui.Checkbox({label:'Use Tower Weather Data',value:false});
  UI.inputPanel.add(UI.useTowerWeather);
  UI.weatherNote=ui.Label({value:'Optional: tower temperature, humidity, shortwave, pressure, precipitation and LE_CORR + H_CORR. Complete months only; no gridded substitutions. CO2 and satellite vegetation remain unchanged.',style:{whiteSpace:'pre-wrap'}});
  UI.inputPanel.add(UI.weatherNote);
}
function addTowerControls(cache, tower, start, end, lat, lon, coordinateMode) {
  var note=ui.Label({value:'Import tower data in the inputs before running ACCLIM.',style:{whiteSpace:'pre-wrap'}});
  UI._chartCard.add(note); if(!tower)return;
  var months=towerMonthsInRange(tower,start,end);
  var coverage=ui.Label({value:'',style:{whiteSpace:'pre-wrap'}});
  var monthChoice=ui.Select({items:[],placeholder:'Monthly tower coverage',onChange:function(key){
    var m=months[key]; if(!m)return;
    function pct(x){return x===null?'unavailable':x.toFixed(1)+'%';}
    coverage.setValue(key+': missing GPP '+pct(m.missingPct)+(tower.secondVariable?'; missing '+tower.secondVariable+' '+pct(m.secondMissingPct):'')+
      '; LE_CORR missing '+pct(m.leMissingPct)+'. Valid primary/daytime records: '+m.n+'/'+m.daylight+'. Unknown radiation: '+(m.total?(100*m.unknownLight/m.total).toFixed(1)+'%':'unavailable')+
      '; absent intervals: '+m.absent+'. No coverage cutoff. Missing percentages exclude records whose daytime status is unknown.');
  }});
  UI._chartCard.add(monthChoice).add(coverage);
  function update(){
    var count=0,partial=[],keys=[];
    cache.rows.forEach(function(row){
      var key=row.date.slice(0,7),m=months[key];
      row.TowerLE_CORR_mm_day=m?m.leET:null;row.TowerLE_CORR_count=m?m.leN:null;row.TowerLE_CORR_missing_pct=m?m.leMissingPct:null;
      row.TowerGPP=m?m.gpp:null;row.TowerGPP2=m?m.secondGpp:null;row.TowerGPP_count=m?m.n:null;
      row.TowerGPP_missing_pct=m?m.missingPct:null;row.TowerGPP2_missing_pct=m?m.secondMissingPct:null;
      if(m){keys.push(key);if(m.n||m.secondN)count++;var p=key.split('-'),days=new Date(Number(p[0]),Number(p[1]),0).getDate();if(m.days<days)partial.push(key);}
    });
    monthChoice.items().reset(keys);if(keys.length)monthChoice.setValue(keys[0]);
    var mismatch=!coordinateMode||tower.latitude===null||tower.longitude===null||Math.abs(lat-tower.latitude)>0.00001||Math.abs(lon-tower.longitude)>0.00001;
    var provenance=tower.provenance||{};
    note.setValue(tower.site+': '+count+' overlapping months with GPP. Primary: '+(tower.gppVariable||'GPP')+(tower.secondVariable?'; additional: '+tower.secondVariable:'')+
      '. '+(mismatch?'Verify selected location against tower coordinates. ':'')+'Whole-footprint observations; dates '+start+' through '+end+'. '+
      (partial.length?'Partial-month coverage: '+partial.join(', ')+'. Model points remain full-month estimates. ':'')+
      (provenance.citation?'Citation: '+provenance.citation+' ':'')+(provenance.productId?'Product: '+provenance.productId:''));
    UI._lastOutput=carbonRowsToFC(cache.rows);
    cache.refresh.forEach(function(refresh){refresh();});
  }
  if(cache.rows)update();else cache.ready.push(update);
}

// Independent observation requests: never include model graphs in LST reductions.
function loadThermalComparisons(cache, geom, start, end, requestId) {
  var dates = ee.FeatureCollection(cache.rows.map(function(row) { return ee.Feature(null, {date: row.date}); }));
  var first = ee.Date(ee.Date(start).format('YYYY-MM-01'));
  var last = ee.Date(ee.Date(end).format('YYYY-MM-01')).advance(1, 'month');
  function load(name, field, fetch) {
    var note = ui.Label(name + ': loading observations… Model results are ready.');
    addChartNote(note);
    function finish(result, error) {
      if (requestId !== UI._requestId) return;
      if (error || !result) { note.setValue(name + ': unavailable — ' + (error || 'empty response')); return; }
      var byDate = {};
      result.forEach(function(row) { byDate[row.date] = row; });
      cache.rows.forEach(function(row) {
        var observed = byDate[row.date];
        if (observed) Object.keys(observed).forEach(function(key) { if (key !== 'date') row[key] = observed[key]; });
      });
      var count = cache.rows.filter(function(row) { return carbonFinite(row[field]); }).length;
      var details = field === 'Landsat_LST' ? ' Valid scenes/month: ' + cache.rows.map(function(row) {
        return row.date.slice(0, 7) + ': ' + (row.Landsat_LST_scene_count || 0);
      }).join(', ') : ' Quality-screened observations.';
      note.setValue(name + ': ' + count + '/' + cache.rows.length + ' months available.' + details);
      UI._lastOutput = carbonRowsToFC(cache.rows); // include completed observations in monthly export
      cache.refresh.forEach(function(refresh) { refresh(); });
    }
    try {
      var observations = fetch(dates, geom, first, last);
      observations.toList(observations.size()).map(function(f) { return ee.Feature(f).toDictionary(); }).evaluate(finish);
    } catch (error) { finish(null, error.message); }
  }
  load('Terra LST (MOD11A1)', 'LST', attachMonthlyLST);
  load('Aqua LST (MYD11A1)', 'Aqua_LST', attachMonthlyAquaLST);
  // Landsat comparison retained above but disabled for this MOD11A1 trial.
}

function addChartNote(note) {
  if (UI._chartCard) UI._chartCard.add(note);
}


UI.exportPanel = ui.Panel({
  layout: ui.Panel.Layout.flow('vertical'),
  style: { margin: '8px 0px 0px 0px' }
});
UI.exportCsv = ui.Button({
  label: 'Export CSV to Drive',
  onClick: function() { onExportCsv(); },
  style: { stretch: 'horizontal', padding: '6px',
           margin: '0px 4px 0px 0px' }
});
UI.exportPanel.add(UI.exportCsv);
UI.exportDaily = ui.Button({label: 'Export daily CSV',
  onClick: function() { onExportCarbon('daily'); }, style: {stretch: 'horizontal'}});
UI.export16day = ui.Button({label: 'Export 16-day CSV',
  onClick: function() { onExportCarbon('16day'); }, style: {stretch: 'horizontal'}});
UI.exportPanel.add(UI.exportDaily);
UI.exportPanel.add(UI.export16day);
UI.resultsPanel.add(UI.exportPanel);

// Title is attached to the floating control header below.
UI.mainPanel.add(UI.subtitle);
UI.mainPanel.add(UI.inputPanel);
UI.mainPanel.add(UI.statusPanel);
// Results live in a separate map overlay.


// Event handlers ------------------------------------------------------------

// Module-scope reference to the latest output FC (for export, chart re-render).
UI._lastOutput = null;
UI._lastDaily = null;
UI._last16day = null;
UI._lastRun = null;
UI._requestId = 0;

function setStatus(msg) {
  UI.statusPanel.style().set({shown: true});
  UI.statusLabel.setValue(msg);
}

function clearResults() {
  UI._requestId++;
  UI._lastOutput = null;
  UI._lastDaily = null;
  UI._last16day = null;
  UI._lastRun = null;
  UI.runButton.setDisabled(false);
  UI.chartHolder.clear();
  UI._resultCharts = [];
  UI._chartCount = 0;
  UI._chartRow = null;
  UI._chartCard = null;
  UI.resultsPanel.style().set({shown: false});
  if (UI.resultsOverlay) UI.resultsOverlay.setExpanded(false);
}

/**
 * Resolve the user's geometry choice into an ee.Geometry.
 *
 * Precedence (explicit latitude/longitude mode overrides the following):
 *   1. A geometry drawn on the map (last drawn layer wins).
 *   2. A pre-loaded asset selected in the dropdown.
 *   3. null (caller surfaces an error).
 *
 * Drawn points are passed through unbuffered — the data fetchers use
 * mean-reduction with bestEffort, so a bare point samples a single pixel.
 * (If you want a point to represent a buffered footprint instead, .buffer()
 * the result here.)
 */
function readEnteredCoordinates() {
  if (!UI.useCoordinates.getValue()) return null;
  var latText = String(UI.latitude.getValue() || '').trim();
  var lonText = String(UI.longitude.getValue() || '').trim();
  if (!latText || !lonText) throw new Error('Enter both latitude and longitude in decimal degrees.');
  var lat = Number(latText), lon = Number(lonText);
  if (!isFinite(lat) || lat < -90 || lat > 90) throw new Error('Latitude must be a number from −90 to 90.');
  if (!isFinite(lon) || lon < -180 || lon > 180) throw new Error('Longitude must be a number from −180 to 180.');
  return [lon, lat]; // Earth Engine point coordinate order
}

function resolveGeometry() {
  // Explicit coordinate mode takes precedence without deleting existing drawings.
  var coordinates = readEnteredCoordinates();
  if (coordinates) return ee.Geometry.Point(coordinates);
  // 1. Drawing tools
  var drawn = UI.drawingTools && UI.drawingTools.layers().length()
    ? UI.drawingTools.layers().get(0).getEeObject()
    : null;
  if (drawn) return drawn;

  // 2. Asset dropdown
  var assetName = UI.assetSelector.getValue();
  if (assetName) {
    for (var i = 0; i < UI.assetList.length; i++) {
      if (UI.assetList[i].name === assetName) {
        return ee.FeatureCollection(UI.assetList[i].assetId).geometry();
      }
    }
  }
  return null;
}

// Match resolveGeometry's precedence and capture the label when the run starts.
// A later dropdown change must not relabel an in-flight result.
function resolveGeometryName(geom, done) {
  var coordinates = readEnteredCoordinates();
  if (coordinates) {
    done('Center: ' + coordinates[1].toFixed(5) + '° latitude, ' + coordinates[0].toFixed(5) + '° longitude');
    return;
  }
  if (UI.drawingTools && UI.drawingTools.layers().length() &&
      UI.drawingTools.layers().get(0).getEeObject()) {
    // Resolve the captured run geometry, not a later edited drawing.
    geom.centroid(1, ee.Projection('EPSG:4326')).coordinates().evaluate(function(coords, error) {
      if (error || !coords || !carbonFinite(coords[0]) || !carbonFinite(coords[1])) {
        done(null, error || 'Could not determine the drawn geometry center.');
        return;
      }
      done('Center: ' + coords[1].toFixed(5) + '° latitude, ' + coords[0].toFixed(5) + '° longitude');
    });
    return;
  }
  done(UI.assetSelector.getValue() || 'Selected geometry');
}

function regionChartTitle(title, regionName) {
  return title + '\n' + regionName;
}

function onRun() {
  clearResults();
  setStatus('Resolving geometry…');

  var startDate = UI.startI.input.getValue();
  var endDate = UI.endI.input.getValue();
  var runTower = UI._towerData;
  var towerLat = Number(UI.latitude.getValue()), towerLon = Number(UI.longitude.getValue());
  var towerCoordinateMode = UI.useCoordinates.getValue();
  var geom;
  try { geom = resolveGeometry(); }
  catch (coordinateError) { setStatus(coordinateError.message); return; }

  if (!geom) {
    setStatus('Please enter coordinates and enable coordinate mode, draw a point/polygon, or choose a pre-loaded asset.');
    return;
  }
  if (!startDate || !endDate) {
    setStatus('Please fill both start and end dates (YYYY-MM-DD).');
    return;
  }

  setStatus('Assembling forcing data…');

  try {
    if (carbonDateMillis(endDate) < carbonDateMillis(startDate))
      throw new Error('End date must be on or after start date.');
    UI._runWeatherPlan=null;
    if(UI.useTowerWeather.getValue()) {
      if(!runTower || !towerCoordinateMode || !carbonFinite(runTower.latitude) || !carbonFinite(runTower.longitude) || Math.abs(towerLat-runTower.latitude)>0.00001 || Math.abs(towerLon-runTower.longitude)>0.00001)
        throw new Error('Use the imported tower coordinates for a tower-weather run.');
      UI._runWeatherPlan=towerWeatherPlan(runTower,startDate,endDate);
      if(!UI._runWeatherPlan.segments.length){
        UI.weatherNote.setValue('Tower weather unavailable:\n'+UI._runWeatherPlan.rejected.join('\n'));
        throw new Error('No complete tower-weather months. '+UI._runWeatherPlan.rejected.slice(0,2).join(' | ')+'. See the weather note for all months.');
      }
      UI.weatherNote.setValue('Tower weather selected. Unavailable months: '+(UI._runWeatherPlan.missing.join(', ')||'none')+'. Acclimation and soil-water calculations restart after gaps. CO2 and satellite vegetation remain unchanged. Supplied gap-filled fields are retained.');
    } else {UI.weatherNote.setValue('Gridded weather selected; tower weather is not used.');}
    UI.runButton.setDisabled(true);
    var runRequestId = UI._requestId;
    resolveGeometryName(geom, function(regionName, locationError) {
      if (runRequestId !== UI._requestId) return;
      if (locationError) {
        UI.runButton.setDisabled(false);
        setStatus('Location error: ' + locationError);
        return;
      }
      try {
    var output = runForGeometry(geom, startDate, endDate, {faparMethod: UI.faparSelector.getValue()}, {weatherPlan:UI._runWeatherPlan});
    UI._lastOutput = output;
    UI._lastRun = {start: startDate, end: endDate, site: regionName, regionName: regionName, weatherSource:UI._runWeatherPlan?'tower':'gridded'};

    regionName += UI._runWeatherPlan ? '\nTower weather' : '\nGridded weather';
    // Center the map on the active geometry so the user sees where they ran.
    if (UI.map) UI.map.centerObject(geom);

    setStatus('Building chart (server-side)...');

    var conductanceChart;
    var chartCache = {rows: null, ready: [], refresh: []};
    var c4Chart = function(rows) {
      chartCache.rows = rows;
      addRemainingDiagnosticCharts();
      var resistance = rows.map(function(row) { return row.rs; }).filter(carbonFinite);
      var conductanceOptions = conductanceChart.getOptions();
      conductanceOptions.vAxes[1].viewWindow.max = resistance.length ? 10 * Math.min.apply(null, resistance) : 1000;
      conductanceChart.setOptions(conductanceOptions);
      chartCache.ready.forEach(function(ready) { ready(); });
      chartCache.ready = [];
      addOptionalSeriesChart(chartCache, ['FAPAR', 'C4frac', 'NDVI', 'EVI2', 'LAI'],
        ['FAPAR', 'C4 fraction', 'NDVI', 'EVI2', 'Green LAI'],
        ['#000000', '#1976D2', '#2E7D32', '#8E24AA', '#EF6C00'], [2, 3, 4], {
          series: {0: {color: '#000000', lineDashStyle: [6, 4]},
            1: {color: '#1976D2', lineDashStyle: [6, 4]},
            2: {color: '#2E7D32'}, 3: {color: '#8E24AA'},
            4: {color: '#EF6C00', targetAxisIndex: 1}},
          title: regionChartTitle('FAPAR & C4 Fraction', regionName), chartArea: {top: 80},
          hAxis: {title: 'Date'},
          vAxes: {0: {title: 'Fraction / vegetation index'},
            1: {title: 'Green LAI (m² leaf/m² ground)', viewWindow: {min: 0},
              gridlines: {color: 'transparent'}, minorGridlines: {color: 'transparent'}}},
          lineWidth: 1.2, pointSize: 0
        });
    };
    var extraTower=runTower && runTower.secondVariable;
    addOptionalSeriesChart(chartCache, ['PredGPPmean', 'GPPmax', 'Amax', 'TowerGPP'].concat(extraTower?['TowerGPP2']:[]),
      ['GPPmean', 'GPPmax', 'Amax', 'Tower '+(runTower && runTower.gppVariable || 'GPP')].concat(extraTower?['Tower '+extraTower]:[]),
      ['#2E7D32', '#1976D2', '#EF6C00', '#111111'].concat(extraTower?['#8E24AA']:[]), extraTower?[1,2,3,4]:[1,2,3], {
        series: {3: {lineWidth: 0, pointSize: 5},4:{lineWidth:0,pointSize:5}},
        title: regionChartTitle('GPP', regionName), chartArea: {top: 80},
        hAxis: {title: 'Date'}, vAxis: {title: 'umol/m2/s (daytime)'}, lineWidth: 1.5, pointSize: 0
      });
    addTowerControls(chartCache, runTower, startDate, endDate, towerLat, towerLon, towerCoordinateMode);
    UI._annualGPPNote = ui.Label({value: 'Annual GPP: calculating...', style: {whiteSpace: 'pre-wrap'}});
    addChartNote(UI._annualGPPNote);
    addOptionalSeriesChart(chartCache,['ET_mm_day','T_mm_day','E_mm_day','I_mm_day','TowerLE_CORR_mm_day'],
      ['ET','T','E','I','Tower LE_CORR'],['#000000','#2E7D32','#1976D2','#808080','#D81B60'],[4],{
        title:regionChartTitle('Evapotranspiration',regionName),chartArea:{top:80},
        hAxis:{title:'Date'},vAxis:{title:'mm/day'},lineWidth:1.2,pointSize:0,
        series:{4:{lineWidth:0,pointSize:5}}
      });
    addChartNote(ui.Label(runTower&&runTower.leVariable==='LE_CORR'?
      'Tower LE_CORR: available full-day flux mean × 86400 / 2450000 (mm/day). Includes night and negative fluxes; missing percentages are in Monthly tower coverage.':
      'To add tower LE_CORR, reprocess a CSV/ZIP containing LE_CORR and import the new summary.'));

    function addRemainingDiagnosticCharts() {
    conductanceChart = ui.Chart.feature.byFeature(output, 'date', ['gs', 'rs'])
      .setSeriesNames(['gs', 'rs'])
      .setOptions({title: regionChartTitle('Canopy Stomatal Conductance', regionName), chartArea: {top: 80},
        hAxis: {title: 'Date'},
        vAxes: {0: {title: 'gs (mol/m2/s)', viewWindow: {min: 0}},
          1: {title: 'rs (s/m)', viewWindow: {min: 0, max: 1000},
            gridlines: {color: 'transparent'}, minorGridlines: {color: 'transparent'}}},
        //series: {0: {targetAxisIndex: 0, color: '#2E7D32'}, 1: {targetAxisIndex: 1, color: '#EF6C00'}},
		series: {0: {targetAxisIndex: 0, color: '#2E7D32'}, 1: {targetAxisIndex: 1, color: '#EF6C00', lineDashStyle: [6, 4]}},		
        lineWidth: 1.5, pointSize: 0, interpolateNulls: false});
    addResultChart(conductanceChart);
    addOptionalSeriesChart(chartCache, ['cica_mean', 'cicaC3_mean', 'cicaC4_mean'],
      ['ci/ca', 'C3', 'C4'], ['#000000', '#2E7D32', '#1976D2'], [1, 2], {
        title: regionChartTitle('ci/ca', regionName), chartArea: {top: 80},
        hAxis: {title: 'Date'}, vAxis: {title: 'ci/ca', viewWindow: {min: 0, max: 1}}, lineWidth: 0, pointSize: 4, pointShape: 'circle'
      });
    addOptionalSeriesChart(chartCache, ['Omega_surface', 'Omega_canopy', 'Omega_soil'],
      ['Omega_surface', 'Omega_canopy', 'Omega_soil'], ['#000000', '#2E7D32', '#1976D2'], [2], {
        title: regionChartTitle('Decoupling Coefficient', regionName), chartArea: {top: 80},
        hAxis: {title: 'Date'}, vAxis: {title: 'Omega (dimensionless)'}, lineWidth: 1.5, pointSize: 0
      });
    addOptionalSeriesChart(chartCache,
      ['Tsurface_mean', 'Tcanopy_mean', 'Tsoil_mean', 'Tday_mean', 'Tsurface_max', 'Tmax_mean', 'LST', 'Aqua_LST'],
      ['Tsurface_mean', 'Tcanopy_mean', 'Tsoil_mean', 'Tair', 'Tsurface_max', 'Tmax', 'Terra LST', 'Aqua LST'],
      ['#000000', '#2E7D32', '#8D6E63', '#1976D2', '#EF6C00', '#8E24AA', '#000000', '#D81B60'], [0, 1, 2, 3, 5, 6], {
        series: {6: {color: '#000000', lineWidth: 0, pointSize: 3, pointShape: 'circle'},
          7: {color: '#D81B60', lineWidth: 0, pointSize: 3, pointShape: 'circle'}},
        title: regionChartTitle('Surface, Canopy & Soil Temperatures', regionName), chartArea: {top: 80},
        hAxis: {title: 'Date'}, vAxis: {title: 'Temperature (°C)'}, lineWidth: 1.5, pointSize: 0
      });
    addChartNote(ui.Label(annualClimateNote(chartCache.rows, startDate, endDate, 'Tdaily')));
    loadThermalComparisons(chartCache, geom, startDate, endDate, runRequestId);

    }
    UI.resultsPanel.style().set({shown: true});
    UI.resultsOverlay.setExpanded(true);
    setStatus('Calculating daily carbon totals and loading annual RAP cover...');
    loadCarbonProducts(output, geom, startDate, endDate, runRequestId, regionName, c4Chart);
      } catch (runError) {
        UI.runButton.setDisabled(false);
        setStatus('Error: ' + runError.message);
      }
    });
  } catch (e) {
    UI.runButton.setDisabled(false);
    setStatus('Error: ' + e.message);
  }
}

function onReset() {
  UI.startI.input.setValue(DEFAULTS.startDate);
  UI.endI.input.setValue(DEFAULTS.endDate);
  // Clear drawn geometry (asset selection is left alone — ui.Select has no
  // reliable way to reset to "no selection" without re-creating the widget).
  if (UI.drawingTools) {
    while (UI.drawingTools.layers().length() > 0) {
      UI.drawingTools.layers().remove(UI.drawingTools.layers().get(0));
    }
  }
  UI.statusPanel.style().set({shown: false});
  clearResults();
}

// (No trailing .onClick() calls — click handlers are wired into the ui.Button
//  constructors above. ui.Label has no onClick method.)


// ============================================================================
// SECTION 8 — EXPORTS
// ============================================================================

function onExportCarbon(cadence) {
  var collection = cadence === 'daily' ? UI._lastDaily : UI._last16day;
  if (!collection || !UI._lastRun) { setStatus('Run the model and wait for carbon products before exporting.'); return; }
  var name = ('acclim_' + UI._lastRun.site + '_' + cadence + '_' + UI._lastRun.start + '_to_' + UI._lastRun.end)
    .replace(/[^A-Za-z0-9_-]/g, '_');
  Export.table.toDrive({collection: collection, description: name, fileNamePrefix: name, fileFormat: 'CSV'});
  setStatus('Export task queued: "' + name + '". Open the Tasks tab to start it.');
}

function onExportCsv() {
  if (!UI._lastOutput) {
    setStatus('Run the model before exporting.');
    return;
  }
  var site = UI._lastRun.site;
  var start = UI._lastRun.start;
  var end = UI._lastRun.end;
  var name = 'acclim_' + site + '_' + start + '_to_' + end;

  // Drop the geometry to keep CSV tabular.
  var tabular = UI._lastOutput.map(function(f) { return f.setGeometry(null); });

  Export.table.toDrive({
    collection: tabular,
    description: name,
    fileNamePrefix: name,
    fileFormat: 'CSV'
  });
  setStatus('Export task queued: "' + name + '". Open the Tasks tab to start it.');
}


// ============================================================================
// SECTION 9 — ENTRY POINT
// ============================================================================

// ----------------------------------------------------------------------------
// Build the map + drawing tools (pattern adapted from ecoexplorer_v1.txt).
// ----------------------------------------------------------------------------
UI.map = ui.Map();
UI.map.setOptions('SATELLITE');
UI.map.setControlVisibility({ layerList: false, mapTypeControl: false, zoomControl: false });
UI.map.style().set('cursor', 'crosshair');
UI.map.setCenter(DEFAULTS.lon, DEFAULTS.lat, 6);

UI.drawingTools = UI.map.drawingTools();
UI.drawingTools.setShown(false);            // hide the default toolbar
UI.drawingTools.setDrawModes(['point', 'polygon']);
UI.drawingTools.setLinked(false);

function _clearDrawnGeometry() {
  while (UI.drawingTools.layers().length() > 0) {
    UI.drawingTools.layers().remove(UI.drawingTools.layers().get(0));
  }
}

function _startDrawPoint() {
  _clearDrawnGeometry();
  UI.drawingTools.setShape('point');
  UI.drawingTools.draw();
}

function _startDrawPolygon() {
  _clearDrawnGeometry();
  UI.drawingTools.setShape('polygon');
  UI.drawingTools.draw();
}

// Floating control panel docked to the map. The buttons mirror ecoexplorer's
// pattern: draw point, draw polygon, clear. ui.Button is required for onClick.
UI.drawPanel = ui.Panel({
  widgets: [
    ui.Label({value: 'Geometry', style: {fontWeight: 'bold', fontSize: '13px', margin: '0px 0px 4px 0px'}}),
    ui.Panel({widgets: [
      ui.Button({label: 'Draw point', onClick: _startDrawPoint,
        style: {width: '100px', height: '30px', padding: '0px', margin: '0px 4px 0px 0px'}}),
      ui.Button({label: 'Draw polygon', onClick: _startDrawPolygon,
        style: {width: '110px', height: '30px', padding: '0px', margin: '0px 4px 0px 0px'}}),
      ui.Button({label: 'Clear drawn', onClick: _clearDrawnGeometry,
        style: {width: '100px', height: '30px', padding: '0px', margin: '0px'}})
    ], layout: ui.Panel.Layout.flow('horizontal'), style: {margin: '0px', padding: '0px'}})
  ], style: {padding: '4px', backgroundColor: 'white'}
});
// Geometry controls belong to the collapsible input panel.
UI.drawPanel.style().set({position: 'top-left', width: '100%'});
UI.mainPanel.add(UI.drawPanel);

// ----------------------------------------------------------------------------
// Attach panels to the UI root.
// ----------------------------------------------------------------------------
ui.root.clear();
// Floating overlays preserve a full-size map; no sidebar or split panel.
function makeCollapsibleMapPanel(titleWidget, content, position, width, expanded, fullHeight) {
  var body = ui.Panel({widgets: [content], layout: ui.Panel.Layout.flow('vertical'),
    style: {stretch: 'horizontal', height: '480px', maxHeight: '60vh', shown: expanded}});
  var toggle = ui.Button({label: expanded ? '▼' : '▶',
    style: {width: '36px', margin: '2px 6px 2px 0px'}});
  var header = ui.Panel({widgets: [toggle, titleWidget], layout: ui.Panel.Layout.flow('horizontal'),
    style: {stretch: 'horizontal', padding: '4px', backgroundColor: '#f5f5f5', borderRadius: '8px'}});
  var panel = ui.Panel({widgets: [header, body], layout: ui.Panel.Layout.flow('vertical'),
    style: {position: position, width: width, maxHeight: '95%', margin: '8px',
      padding: '4px', backgroundColor: 'white', border: '1px solid #bbb', borderRadius: '12px'}});
  function setExpanded(value) {
    expanded = value;
    body.style().set('shown', expanded);
    toggle.setLabel(expanded ? '▼' : '▶');
    panel.style().set('width', expanded ? width : '190px');
    if (fullHeight) panel.style().set('height', expanded ? '95%' : '46px');
  }
  toggle.onClick(function() { setExpanded(!expanded); });
  setExpanded(expanded);
  return {panel: panel, setExpanded: setExpanded};
}
UI.map.style().set('stretch', 'both');
UI.controlsOverlay = makeCollapsibleMapPanel(UI.title, UI.mainPanel, 'top-left', '420px', true, false);
UI.controlsOverlay.panel.style().set({maxHeight: '100%', margin: '0px 0px 0px 8px'});
// Keep the chart-bearing panel at its expanded size. Only visibility changes.
// A separate launcher provides the compact collapsed state.
function makeResultsOverlay() {
  var timer = null, expanded = false, revision = 0;
  var close = ui.Button({label: '▼', onClick: function() { setExpanded(false); },
    style: {width: '36px', margin: '2px 6px 2px 0px'}});
  var header = ui.Panel({widgets: [close, UI.resultsTitle], layout: ui.Panel.Layout.flow('horizontal'),
    style: {stretch: 'horizontal', padding: '4px', backgroundColor: '#f5f5f5', borderRadius: '8px'}});
  var panel = ui.Panel({widgets: [header, UI.resultsPanel], layout: ui.Panel.Layout.flow('vertical'),
    style: {position: 'top-right', width: '70%', height: '95%', maxHeight: '95%', shown: false,
      margin: '8px', padding: '4px', backgroundColor: 'white', border: '1px solid #bbb', borderRadius: '12px'}});
  var launcher = ui.Panel({widgets: [ui.Button({label: '▶ Results', onClick: function() { setExpanded(true); },
    style: {stretch: 'horizontal', margin: '0px'}})],
    style: {position: 'top-right', width: '190px', margin: '8px', padding: '4px',
      backgroundColor: 'white', border: '1px solid #bbb', borderRadius: '12px'}});
  function setExpanded(value) {
    expanded = value;
    revision++;
    if (timer !== null) { ui.util.clearTimeout(timer); timer = null; }
    panel.style().set('shown', expanded);
    launcher.style().set('shown', !expanded);
    if (expanded) {
      var currentRevision = revision, request = UI._requestId;
      timer = ui.util.setTimeout(function() {
        timer = null;
        if (!expanded || revision !== currentRevision || request !== UI._requestId) return;
        // Redraw after layout settles, preserving each chart's options and data.
        UI._resultCharts.forEach(function(chart) { chart.setOptions(chart.getOptions()); });
      }, 250);
    }
  }
  return {panel: panel, launcher: launcher, setExpanded: setExpanded};
}
UI.resultsOverlay = makeResultsOverlay();
// Custom navigation controls: map types above zoom, inputs beside zoom.
function selectBaseMap(type) {
  UI.map.setOptions(type);
  UI.mapButton.setDisabled(type === 'ROADMAP');
  UI.satelliteButton.setDisabled(type === 'SATELLITE');
}
UI.mapButton = ui.Button({label: 'Map', onClick: function() { selectBaseMap('ROADMAP'); },
  style: {margin: '0px', width: '70px', height: '36px'}});
UI.satelliteButton = ui.Button({label: 'Satellite', onClick: function() { selectBaseMap('SATELLITE'); },
  style: {margin: '0px', width: '90px', height: '36px'}});
UI.baseMapPanel = ui.Panel({widgets: [UI.mapButton, UI.satelliteButton],
  layout: ui.Panel.Layout.flow('horizontal'),
  style: {position: 'top-left', margin: '8px', padding: '4px', backgroundColor: 'white', borderRadius: '8px'}});
function changeMapZoom(delta) {
  UI.map.setZoom(Math.max(0, Math.min(24, UI.map.getZoom() + delta)));
}
UI.zoomPanel = ui.Panel({widgets: [
  ui.Button({label: '+', onClick: function() { changeMapZoom(1); }, style: {width: '36px', height: '36px', margin: '0px'}}),
  ui.Button({label: '−', onClick: function() { changeMapZoom(-1); }, style: {width: '36px', height: '36px', margin: '0px'}})
], layout: ui.Panel.Layout.flow('vertical'),
  style: {position: 'top-left', width: '44px', minWidth: '44px', height: '80px', maxHeight: '80px', margin: '0px', padding: '4px', backgroundColor: 'white', borderRadius: '8px'}});
selectBaseMap('SATELLITE');
UI.map.add(UI.baseMapPanel);
// A shared horizontal parent prevents the map from stacking these controls.
UI.leftControlsRow = ui.Panel({widgets: [UI.zoomPanel, UI.controlsOverlay.panel],
  layout: ui.Panel.Layout.flow('horizontal'),
  style: {position: 'top-left', margin: '8px', padding: '0px',
    width: '95%', maxWidth: '480px', maxHeight: '85%', backgroundColor: 'rgba(0,0,0,0)'}});
UI.map.add(UI.leftControlsRow);
UI.map.add(UI.resultsOverlay.panel);
UI.map.add(UI.resultsOverlay.launcher);
ui.root.add(UI.map);

// Expose the API on the global scope so the Code Editor console can poke at it.
exports = exports || {};
exports.ACCLIM_PARAMS = ACCLIM_PARAMS;
exports.runAtPoint = runAtPoint;
exports.runForGeometry = runForGeometry;
exports.runACCLIM = runACCLIM;
exports.assembleForcing = assembleForcing;
exports.magnusVP = magnusVP;
exports.astronomicalDaylength = astronomicalDaylength;
exports.co2ForDate = co2ForDate;
exports.runEnsemble = runEnsemble;
exports.mcmcCalibrate = mcmcCalibrate;
exports.kgmlPredict = kgmlPredict;

exports.buildCarbonProducts = buildCarbonProducts;
exports.fetchCarbonAnnualContext = fetchCarbonAnnualContext;
