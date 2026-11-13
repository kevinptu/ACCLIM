/**
 * ACCLIM research web app — selectable weather and satellite reflectance.
 * One app Timestep selection aggregates all inputs before shared physiology.
 * Charts and CSV use the selected daily, weekly or monthly intervals.
 * The importer preserves daily tower weather and flux sums/counts.
 * Satellite red/NIR: quality screening, actual dates where available, otherwise
 * composite midpoints; paired linear interpolation with bounded gaps.
 * Bucket-only interpolation preserves storage; missing carbon stays unavailable.
 * All weather sources supply actual vapor pressure to RH/VPD calculations.
 * ERA5-Land supplies pressure; other gridded sources retain constant pressure.
 * CO2 remains the configured constant. Thermal diagnostics are unchanged.
 * See docs/data-sources.md for units, source bands, coverage flags and exports.
 * Legacy monthly adapter helpers are retained for research comparisons;
 * the app uses assembleSelectedForcing through runForGeometryAsync.
 */

// ============================================================================
// SECTION 0 — CONFIG
// ============================================================================

var ACCLIM_PARAMS = {
  kndvi_sigma: 0.15, // RBF bandwidth in fractional reflectance units (IPL-UV GEE example)
  GE_scale: 0.07,        // GE reference used by fge = clamp01(GE / GE_scale)^fge_exp
  blend_exp: 1.6,          // NDVI exponent controlling the weight in both blend options
  fiparMethod: 'NDVI',    // NDVI scaling or NDVI-weighted GE blend
  fexp: 1,                // Exponent applied once to every selected FAPAR method
  fge_exp: 1,             // Raw GE / GE_scale exponent; zero disables it
  fitVIMinP05: false,     // Optional ten-year observed fifth-percentile baseline
  fitNDVIminGE: false,     // Optional FAPAR menu experiment; configured NDVImin by default
  fMgpp_exp: 0,  // GPP seasonal FAPAR response
  FAPARmult: 0.8,          // Multiplier for the direct FIPAR option
  faparMethod: 'EVI2',     // Fixed bounds by default; optional SoilGrids minima resolve per location
  fitEVI2maxFixed: false, // Experiment: fit only the lower endpoint; fix upper at 0.6
  useSharedFPAR: false,   // Use photosynthetic FAPAR for water and energy partitioning too
  fitEVI2: false,         // Fit EVI2 bounds against FIPAR; otherwise use fixed bounds
  k: 0.5,                  // canopy light extinction coefficient for LAI retrieval (must be > 0)
  CUEplant:       0.5,      // total plant NPP / GPP
  carbonToBiomass: 2.1276, // exact RAP conversion supplied by user (~47% C)
  fANPP_intercept: 0.171,   // fANPP = intercept + slope * MAT (deg C)
  fANPP_slope:    0.0129,
  tc_sm:          30,       // humidity-derived soil-moisture running-mean time constant (days)
  tc:             7,        // running-mean time constant (days)
  fDtrans_exp: 0,  // transpiration VPD response
  fMtrans_exp: 0,  // transpiration seasonal FAPAR response
  fDgpp_exp: 0,  // GPP VPD response
  GPPmean2GPPmax: 0.67,      // ratio of mean GPP to GPPmax asymptote
  Oi:             20900,    // leaf internal O2 partial pressure (Pa)
  QRC3:           12,       // C3 quantum requirement (mol photons / mol CO2)
  QRC4:           16.67,    // C4 quantum requirement
  EVI2min:        0,     // EVI2 of bare ground
  powerFixedEVI2: true,    // Scaled EVI2 uses EVI2exp for both interactive and batch runs
  EVI2exp: 1.2,            // Exponent for the powered fixed-bounds EVI2 option
  EVI2max:        0.6,     // EVI2 of full green canopy
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
  // Intercepted PAR fraction used for water and energy partitioning.
  ndvi_min: 0, // Used unless the GE′=0 NDVImin option is selected.
  ndvi_max: 1,
  ndvi_exp: 1.6
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
// Daytime radiation: only shortwave is converted from a 24-hour mean.
// Temperatures are Celsius on input; vapor pressure ea is kPa.
function griddedDaytimeRn(swrDaily, daylength, tday, ea, albedo, emSurface, p, periodSeconds, daylightSeconds) {
  var kelvin = tday.add(273.15);
  var seconds=periodSeconds===undefined?ee.Number(86400):ee.Number(periodSeconds);
  var light=daylightSeconds===undefined?daylength.multiply(3600):ee.Number(daylightSeconds);
  var swDay=swrDaily.multiply(seconds).divide(light);
  var emAirDay = ee.Number(10).multiply(ea).divide(kelvin).pow(1 / 7).multiply(1.31);
  return ee.Number(1).subtract(albedo).multiply(swDay)
    .add(kelvin.pow(4).multiply(p.sb).multiply(emAirDay.subtract(emSurface)));
}

// Dataset adapters. Weather is sampled at its native grid; reflectance at 30/250 m.
// These settings are research choices, not claims about observed revisit intervals.
// Product-wide catalog year ranges checked 2026-10-27; final years may be partial.
// Combined satellites use the union: HLS Landsat from 2013, Sentinel from 2015; MODIS Terra from 2000, Aqua from 2002.
var WEATHER_SOURCES = {
  TERRACLIMATE: {label:'TerraClimate (global/4.6km/monthly/1958-2024)', id:'IDAHO_EPSCOR/TERRACLIMATE', daily:false, scale:4638},
  ERA5: {label:'ERA5-Land (global/11km/daily/1950-2026)', id:'ECMWF/ERA5_LAND/DAILY_AGGR', daily:true, scale:11132},
  GRIDMET: {label:'GRIDMET (contiguous US/4km/daily/1979-2026)', id:'IDAHO_EPSCOR/GRIDMET', daily:true, scale:4638}
};
var SATELLITE_SOURCES = {
  LANDSAT8: {label:'Landsat (global/30m/8-day/1984-2026)', ids:['LANDSAT/COMPOSITES/C02/T1_L2_8DAY'], days:8, scale:30, gap:32, kind:'landsat'},
  LANDSAT32: {label:'Landsat (global/30m/32-day/1984-2026)', ids:['LANDSAT/COMPOSITES/C02/T1_L2_32DAY'], days:32, scale:30, gap:64, kind:'landsat'},
  HLS: {label:'HLS (global/30m/daily/2013-2026)', ids:['NASA/HLS/HLSL30/v002','NASA/HLS/HLSS30/v002'], days:1, scale:30, gap:32, kind:'hls'},
  MODIS1: {label:'MODIS (global/250m/daily/2000-2026)', ids:['MODIS/061/MOD09GQ','MODIS/061/MYD09GQ'], days:1, scale:250, gap:32, kind:'gq'},
  MODIS8: {label:'MODIS (global/250m/8-day/2000-2026)', ids:['MODIS/061/MOD09Q1','MODIS/061/MYD09Q1'], days:8, scale:250, gap:32, kind:'q1'},
  MODIS16: {label:'MODIS (global/250m/16-day/2000-2026)', ids:['MODIS/061/MOD13Q1','MODIS/061/MYD13Q1'], days:16, scale:250, gap:48, kind:'vi'}
};
// One temporal contract for every weather adapter. Flux density is kept separate
// from mass integrated over the represented weather interval.
function withPeriodTiming(f) {
  var seconds=ee.Number(ee.Algorithms.If(ee.Algorithms.IsEqual(f.get('period_seconds'),null),
    f.getNumber('period_days').multiply(86400),f.get('period_seconds')));
  var light=ee.Number(ee.Algorithms.If(ee.Algorithms.IsEqual(f.get('daylight_seconds'),null),
    seconds.multiply(f.getNumber('daylength')).divide(24),f.get('daylight_seconds')));
  return f.set({period_seconds:seconds,daylight_seconds:light,period_days:seconds.divide(86400),
    period_start_millis:f.get('millis'),period_end_millis:f.getNumber('millis').add(seconds.multiply(1000)),
    period_end_exclusive:ee.Date(f.getNumber('millis').add(seconds.multiply(1000))).format('YYYY-MM-dd'),
    flux_time_basis:'daytime mean',period_time_units:'seconds'});
}
function nativePeriodSeconds(r) {
  if(carbonFinite(r.period_seconds))return r.period_seconds;
  if(carbonFinite(r.period_days))return r.period_days*86400;
  // Compatibility only for older exported rows without an interval contract.
  var d=new Date(carbonDateMillis(r.date));
  return r.model_cadence==='daily'?86400:(Date.UTC(d.getUTCFullYear(),d.getUTCMonth()+1,1)-d.getTime())/1000;
}
function nativeDaylightSeconds(r) {
  return carbonFinite(r.daylight_seconds)?r.daylight_seconds:nativePeriodSeconds(r)*r.daylength/24;
}
function periodCarbonByDay(rows,start,end) {
  var byDay={},previousEnd=-Infinity;
  rows.slice().sort(function(a,b){return (carbonFinite(a.period_start_millis)?a.period_start_millis:carbonDateMillis(a.date))-
    (carbonFinite(b.period_start_millis)?b.period_start_millis:carbonDateMillis(b.date));}).forEach(function(r) {
    var a=carbonFinite(r.period_start_millis)?r.period_start_millis:carbonDateMillis(r.date);
    var seconds=nativePeriodSeconds(r),light=nativeDaylightSeconds(r),b=a+seconds*1000;
    if(!(seconds>0))throw new Error('Weather period must have a positive duration: '+r.date);
    if(a<previousEnd)throw new Error('Overlapping weather periods: '+r.date);
    previousEnd=b;
    var valid=carbonFinite(r.PredGPPmean)&&r.PredGPPmean>=0&&carbonFinite(light)&&light>=0&&light<=seconds;
    var mass=valid?r.PredGPPmean*light*12e-6:null;
    var first=Math.max(start,Math.floor(a/86400000)*86400000),last=Math.min(end,b-1);
    for(var t=first;t<=last;t+=86400000) {
      var overlap=Math.max(0,Math.min(b,t+86400000)-Math.max(a,t))/1000;
      if(!overlap)continue;
      var key=carbonDateString(t);
      if(!byDay[key])byDay[key]={gpp:0,validSeconds:0,coveredSeconds:0,weather_source:r.weather_source||'gridded',mode:'native_daily'};
      var day=byDay[key];day.coveredSeconds+=overlap;
      if(valid){day.gpp+=mass*overlap/seconds;day.validSeconds+=overlap;}
      if(seconds!==86400||a!==t)day.mode='allocated_from_weather_period';
      if(day.weather_source!==(r.weather_source||'gridded'))day.weather_source='mixed';
    }
  });
  Object.keys(byDay).forEach(function(k){if(Math.abs(byDay[k].validSeconds-86400)>1e-6)byDay[k].gpp=null;});
  return byDay;
}

// Shared, clipped intervals for every model input, chart and export.
function selectedTimestepPeriods(start,end,cadence) {
  if(['daily','weekly','monthly'].indexOf(cadence)<0)throw new Error('Choose Daily, Weekly or Monthly.');
  var a=carbonDateMillis(start),stop=carbonDateMillis(end)+86400000,out=[];
  if(!isFinite(a)||!isFinite(stop)||stop<=a||(stop-a)/86400000>40000)throw new Error('Invalid timestep date range.');
  while(a<stop) {
    var d=new Date(a),y=d.getUTCFullYear(),b;
    if(cadence==='monthly')b=Date.UTC(y,d.getUTCMonth()+1,1);
    else if(cadence==='weekly')b=Math.min(Date.UTC(y,0,1)+(Math.floor((a-Date.UTC(y,0,1))/604800000)+1)*604800000,Date.UTC(y+1,0,1));
    else b=a+86400000;
    b=Math.min(b,stop);
    out.push({date:carbonDateString(a),millis:a,year:y,month:d.getUTCMonth()+1,
      period_days:(b-a)/86400000,period_seconds:(b-a)/1000,period_start_millis:a,period_end_millis:b,
      period_end_exclusive:carbonDateString(b),model_cadence:cadence});
    a=b;
  }
  return out;
}
// Integrate amounts, duration-weight means, and retain missingness and source resolution.
function aggregateWeatherPeriods(rows,periods,sourceCadence,intervalMinutes) {
  var source=rows.map(function(r){var a=carbonDateMillis(r.date||r.month+'-01');
    var seconds=carbonFinite(r.period_seconds)?r.period_seconds:carbonFinite(r.period_days)?r.period_days*86400:
      sourceCadence==='monthly'?(Date.UTC(new Date(a).getUTCFullYear(),new Date(a).getUTCMonth()+1,1)-a)/1000:86400;
    return {r:r,a:a,b:a+seconds*1000,seconds:seconds};}).sort(function(a,b){return a.a-b.a;});
  source.forEach(function(s,i){if(!(s.seconds>0)||!isFinite(s.a)||(i&&s.a<source[i-1].b))throw new Error('Invalid or overlapping source weather periods.');});
  return periods.map(function(period){
    var a=period.millis,b=period.period_end_millis,seconds=period.period_seconds;
    var parts=source.filter(function(s){return s.a<b&&s.b>a;}).map(function(s){return {s:s,w:(Math.min(b,s.b)-Math.max(a,s.a))/1000};});
    var r=mergeParams(parts.length?parts[0].s.r:{},period),covered=parts.reduce(function(n,p){return n+p.w;},0);
    ['Tmin','Tmax','Tdaily','SWRdaily','ea','Pdaily','precip','Rndaily'].forEach(function(k){
      r[k]=Math.abs(covered-seconds)<0.01&&parts.every(function(p){return carbonFinite(p.s.r[k]);})?
        parts.reduce(function(n,p){return n+p.s.r[k]*p.w/(k==='precip'?p.s.seconds:seconds);},0):null;
    });
    r.weather_source_cadence=sourceCadence;
    r.weather_temporally_allocated_pct=100*parts.reduce(function(n,p){return n+(p.s.seconds>seconds?p.w:0);},0)/seconds;
    r.weather_temporal_method=r.weather_temporally_allocated_pct?'coarser source means held constant; precipitation allocated by duration':'duration-weighted means; precipitation summed';
    r.weather_missing_pct=100*(1-parts.reduce(function(n,p){return n+(p.s.r.complete===false||!['Tmin','Tmax','Tdaily','SWRdaily','ea','Pdaily','precip'].every(function(k){return carbonFinite(p.s.r[k]);})?0:p.w);},0)/seconds);
    r.weather_missing_day_pct=r.weather_missing_pct;
    if(intervalMinutes){
      r.expected=seconds/(60*intervalMinutes);r.counts={};
      ['ta','sw','ea','pa','rain','rn'].forEach(function(k){r.counts[k]=parts.reduce(function(n,p){return n+((p.s.r.counts&&p.s.r.counts[k])||0)*p.w/p.s.seconds;},0);});
      r.complete=Math.abs(covered-seconds)<0.01&&parts.every(function(p){return p.s.r.complete===true;})&&
        Object.keys(r.counts).every(function(k){return Math.abs(r.counts[k]-r.expected)<0.01;})&&
        ['Tmin','Tmax','Tdaily','SWRdaily','ea','Pdaily','precip','Rndaily'].every(function(k){return carbonFinite(r[k]);});
    }
    return r;
  });
}
function selectedTowerWeatherPlan(tower,start,end,cadence) {
  if(!tower||!tower.weather||tower.weather.version!==2)throw new Error('Reprocess the source CSV/ZIP with the current importer for tower weather.');
  var w=tower.weather,source=w.daily||w.periods||w.months,sourceCadence=w.daily?'daily':w.periods?(w.cadence||'daily'):'monthly';
  if(!Array.isArray(source)||!source.length)throw new Error('No tower weather available.');
  if(tower.intervalMinutes!==30&&tower.intervalMinutes!==60)throw new Error('Unsupported tower interval.');
  var rows=aggregateWeatherPeriods(source,selectedTimestepPeriods(start,end,cadence),sourceCadence,tower.intervalMinutes);
  var plan={cadence:cadence,records:[],segments:[],missing:[],rejected:[],byMonth:{},source:w.source||{}},segment=null;
  rows.forEach(function(r){var valid=r.complete&&r.Tmin<=r.Tdaily&&r.Tdaily<=r.Tmax&&r.ea>=0&&r.Pdaily>0&&r.precip>=0;
    plan.records.push({data:r,valid:valid});
    if(valid){if(!segment){segment=[];plan.segments.push(segment);}segment.push(r);}
    else{segment=null;plan.missing.push(r.date);plan.rejected.push(r.date+': incomplete '+cadence+' tower weather');}
  });
  (w.months||[]).forEach(function(r){plan.byMonth[r.month]={data:r,valid:r.complete===true};});
  return plan;
}

function sourceCalendar(start, end, daily) {
  start = ee.Date(start); end = ee.Date(end);
  if (!daily) start = ee.Date(start.format('YYYY-MM-01'));
  var unit = daily ? 'day' : 'month';
  return ee.FeatureCollection(ee.List.sequence(0, end.difference(start, unit).floor()).map(function(i) {
    var d = start.advance(i, unit), next = d.advance(1, unit);
    return ee.Feature(null, {date:d.format('YYYY-MM-dd'), millis:d.millis(), year:d.get('year'), month:d.get('month'),
      period_days:next.difference(d,'day'),period_seconds:next.difference(d,'second'),
      period_start_millis:d.millis(),period_end_millis:next.millis(), model_cadence:daily?'daily':'monthly'});
  }));
}
function sampleWeather(geom, start, end, key) {
  var spec = WEATHER_SOURCES[key];
  if (!spec) throw new Error('Unknown weather source: '+key);
  var names = key==='TERRACLIMATE' ? ['tmmn','tmmx','srad','vap','pr'] : key==='GRIDMET' ?
    ['tmmn','tmmx','srad','rmax','pr'] : ['temperature_2m_min','temperature_2m_max','temperature_2m',
      'dewpoint_temperature_2m','surface_solar_radiation_downwards_sum','total_precipitation_sum','surface_pressure'];
  var images = ee.ImageCollection(spec.id).filterBounds(geom).filterDate(start,ee.Date(end).advance(1,'day'));
  // Retain nulls in the calendar. Never turn masked pixels into valid zero weather.
  return sourceCalendar(start,end,spec.daily).map(function(row) {
    var date=ee.Date(row.get('millis'));
    var matches=images.filterDate(date,date.advance(1,'day'));
    var blank=ee.Image.constant(ee.List.repeat(0,names.length)).rename(names).updateMask(ee.Image.constant(0));
    var img=ee.Image(ee.Algorithms.If(matches.size(),matches.first(),blank));
    var s=ee.Feature(null,img.select(names).reduceRegion({reducer:ee.Reducer.mean(),geometry:geom,
      scale:spec.scale,maxPixels:1e10,tileScale:4}));
    function number(name, factor, offset) {
      return ee.Algorithms.If(ee.Algorithms.IsEqual(s.get(name),null),null,
        ee.Number(s.get(name)).multiply(factor===undefined?1:factor).add(offset||0));
    }
    var lo=key==='ERA5'?'temperature_2m_min':'tmmn', hi=key==='ERA5'?'temperature_2m_max':'tmmx';
    var tmin=number(lo,key==='TERRACLIMATE'?0.1:1,key==='TERRACLIMATE'?0:-273.15);
    var tmax=number(hi,key==='TERRACLIMATE'?0.1:1,key==='TERRACLIMATE'?0:-273.15);
    var ea=key==='TERRACLIMATE'?number('vap',0.001): key==='GRIDMET'?
      ee.Algorithms.If(s.toDictionary(['tmmn','rmax']).size().eq(2),magnusVP(ee.Number(tmin)).multiply(ee.Number(s.get('rmax'))).divide(100),null):
      ee.Algorithms.If(ee.Algorithms.IsEqual(s.get('dewpoint_temperature_2m'),null),null,magnusVP(ee.Number(s.get('dewpoint_temperature_2m')).subtract(273.15)));
    var tmean=key==='ERA5'?number('temperature_2m',1,-273.15):
      ee.Algorithms.If(s.toDictionary([lo,hi]).size().eq(2),ee.Number(tmin).add(ee.Number(tmax)).divide(2),null);
    var rain=key==='ERA5'?number('total_precipitation_sum',1000):number('pr');
    rain=ee.Algorithms.If(ee.Algorithms.IsEqual(rain,null),null,ee.Algorithms.If(ee.Number(rain).gte(0),rain,null));
    return row.setGeometry(geom).set({Tmin:tmin,Tmax:tmax,Tdaily:tmean,ea:ea,
      SWRdaily:key==='ERA5'?number('surface_solar_radiation_downwards_sum',1/86400):number('srad',key==='TERRACLIMATE'?0.1:1),
      precip:rain,
      Pdaily:key==='ERA5'?number('surface_pressure',0.001):MVP_CONSTANTS.Pdaily,
      weather_source:key,weather_dataset:spec.id,weather_scale_m:spec.scale,
      weather_gapfill_provenance:'supplied product values retained; original observation gaps not supplied by catalog',
      weather_source_missing_pct:null});
  });
}
function clearModisState(state) {
  // Cloud-state 3 (not set/assumed clear) is explicitly excluded.
  return state.bitwiseAnd(3).eq(0).and(state.bitwiseAnd(4).eq(0))
    .and(state.bitwiseAnd(3<<8).eq(0)).and(state.bitwiseAnd((1<<10)|(1<<12)|(1<<13)|(1<<15)).eq(0));
}
function satelliteImages(geom,start,end,key,landOnly) {
  var spec=SATELLITE_SOURCES[key], merged=ee.ImageCollection([]);
  if (!spec) throw new Error('Unknown satellite source: '+key);
  spec.ids.forEach(function(id) {
    var images=ee.ImageCollection(id).filterBounds(geom).filterDate(start,end);
    if(spec.kind==='gq') {
      var qa=ee.ImageCollection(id.replace('09GQ','09GA')).filterDate(start,end);
      images=ee.ImageCollection(ee.Join.saveFirst('cloudQA').apply(images,qa,
        ee.Filter.equals({leftField:'system:time_start',rightField:'system:time_start'})))
        .filter(ee.Filter.notNull(['cloudQA']));
    }
    images=images.map(function(img) {
      img=img.clip(geom.buffer(spec.scale*2));
      var d=ee.Date(ee.Date(img.get('system:time_start')).format('YYYY-MM-dd'));
      var dateDay=ee.Image.constant(d.millis().divide(86400000));
      var sampleDay=d;
      var red,nir,green,blue,mask=ee.Image.constant(1);
      if(spec.kind==='landsat') {
        red=img.select('red');nir=img.select('nir');green=img.select('green');blue=img.select('blue');
        dateDay=dateDay.add(spec.days/2); // Catalog composites have no selected-pixel date band.
        sampleDay=d.advance(spec.days/2,'day');
      } else if(spec.kind==='hls') {
        red=img.select('B4');nir=img.select(id.indexOf('HLSS30')>=0?'B8A':'B5');green=img.select('B3');blue=img.select('B2');
        var fm=img.select('Fmask');mask=fm.bitwiseAnd(30).eq(0).and(fm.rightShift(6).bitwiseAnd(3).lt(3));
      } else {
        red=img.select('sur_refl_b01').multiply(0.0001);nir=img.select('sur_refl_b02').multiply(0.0001);
        if(spec.kind==='vi') {
          var qa=img.select('DetailedQA');
          mask=img.select('SummaryQA').eq(0).and(qa.bitwiseAnd((1<<8)|(1<<10)|(1<<14)|(1<<15)).eq(0));
          var jan=ee.Date.fromYMD(d.get('year'),1,1);
          dateDay=img.select('DayOfYear').subtract(1).add(jan.millis().divide(86400000));
        } else {
          var qc=img.select(spec.kind==='gq'?'QC_250m':'QA');
          mask=qc.bitwiseAnd(3).eq(0).and(qc.rightShift(4).bitwiseAnd(255).eq(0));
          var state=spec.kind==='gq'?ee.Image(img.get('cloudQA')).select('state_1km'):img.select('State');
          mask=mask.and(clearModisState(state));
          if(spec.kind==='q1') {
            // MODIS final composite stops at year end, unlike overlapping Landsat composites.
            var nextYear=ee.Date.fromYMD(d.get('year').add(1),1,1);
            dateDay=dateDay.add(nextYear.difference(d,'day').min(8).divide(2).floor());
            sampleDay=d.advance(nextYear.difference(d,'day').min(8).divide(2).floor(),'day');
          }
        }
      }
      if(landOnly){
        if(spec.kind==='hls')mask=mask.and(fm.bitwiseAnd(32).eq(0));
        else if(spec.kind==='vi')mask=mask.and(qa.rightShift(11).bitwiseAnd(7).eq(1));
        else if(spec.kind!=='landsat')mask=mask.and(state.rightShift(3).bitwiseAnd(7).eq(1));
        else mask=mask.and(ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select('occurrence').unmask(0).lt(90));
      }
      mask=mask.and(red.gte(0)).and(nir.gte(0)).and(red.lte(1)).and(nir.lte(1)).and(red.add(nir).gt(0));
      // Keep the RGB mask separate from model red/NIR; fit the GE baseline before the run.
      var ge=ee.Image.constant(0).rename('GE').updateMask(ee.Image.constant(0));
      if(green&&blue)ge=green.multiply(2).subtract(red).subtract(blue).rename('GE')
        .updateMask(green.gte(0).and(green.lte(1)).and(blue.gte(0)).and(blue.lte(1)));
      return red.rename('red').addBands(nir.rename('nir')).addBands(ge).addBands(dateDay.rename('observation_day'))
        .updateMask(mask).set({'system:time_start':img.get('system:time_start'),sample_millis:sampleDay.millis()});
    });
    merged=merged.merge(images);
  });
  return merged;
}
// Wang et al. (2023), p. 410: sigma = sum_i |NIR_i - red_i| / N for each image.
function addSpatialKNDVI(img,region,scale) {
  var red=img.select('red'), nir=img.select('nir');
  var distance=nir.subtract(red).abs().rename('sigma');
  var result=distance.reduceRegion({reducer:ee.Reducer.mean().unweighted(),geometry:region,
    crs:red.projection(),scale:scale,maxPixels:1e9,tileScale:4});
  var sigma=ee.Number(ee.Algorithms.If(result.get('sigma'),result.get('sigma'),0));
  var sigmaBand=red.multiply(0).add(sigma).rename('kndvi_sigma').updateMask(ee.Image.constant(sigma.gt(0)));
  var index=nir.subtract(red).divide(sigmaBand.multiply(2)).pow(2).tanh().rename('kNDVI');
  return img.addBands(index).addBands(sigmaBand);
}
function sampleSatellite(geom,start,end,key,unpadded,sigmaRegion) {
  var spec=SATELLITE_SOURCES[key], pad=spec.gap+32;
  var a=ee.Date(start).advance(unpadded?0:-pad,'day'), b=ee.Date(end).advance(unpadded?0:pad,'day');
  var images=satelliteImages(sigmaRegion?sigmaRegion.union(geom,1):geom,a,b.advance(1,'day'),key,!!sigmaRegion);
  if(sigmaRegion)images=images.map(function(img){return addSpatialKNDVI(img,sigmaRegion,spec.scale);});
  var sampleBands=sigmaRegion?['red','nir','GE','kNDVI','kndvi_sigma']:['red','nir','GE'];
  // Group pixels by their actual observation day BEFORE reducing spatially.
  // Same-day Terra/Aqua or HLS values are averaged with a common red/NIR mask.
  // Sample existing images only, rather than rebuilding every composite on every calendar day.
  var sampled=ee.FeatureCollection(images.toList(images.size().max(1)).map(function(item) {
    var img=ee.Image(item),startDay=ee.Date(img.get('sample_millis'));
    var count=spec.kind==='vi'?ee.Date.fromYMD(startDay.get('year').add(1),1,1).difference(startDay,'day').min(16):ee.Number(1);
    return ee.List.sequence(0,count.subtract(1)).map(function(i) {
      var d=startDay.advance(i,'day');
      var pair=img.select(sampleBands).updateMask(img.select('observation_day').eq(d.millis().divide(86400000)));
      var values=pair.reduceRegion({reducer:ee.Reducer.mean(),geometry:geom,scale:spec.scale,maxPixels:1e10,tileScale:4});
      return ee.Feature(null,values).set({date:d.format('YYYY-MM-dd'),millis:d.millis()});
    });
  }).flatten()).filter(ee.Filter.notNull(['red','nir']));
  var observed=sampled.distinct(['millis']).map(function(f) {
    var same=sampled.filter(ee.Filter.eq('millis',f.get('millis')));
    if(sigmaRegion)f=f.set({kNDVI:same.aggregate_mean('kNDVI'),kndvi_sigma:same.aggregate_mean('kndvi_sigma')});
    return f.set({red:same.aggregate_mean('red'),nir:same.aggregate_mean('nir'),GE:same.aggregate_mean('GE'),GE_observation_count:same.filter(ee.Filter.notNull(['GE'])).size(),satellite_observation_count:same.size()});
  });
  return observed;
}
function alignedSatellite(geom,start,end,key,observations) {
  var spec=SATELLITE_SOURCES[key],observed=observations||sampleSatellite(geom,start,end,key);
  var safe=ee.Feature(null,{millis:0,red:0,nir:0});
  return sourceCalendar(start,end,true).map(function(f) {
    var t=f.getNumber('millis');
    var left=observed.filter(ee.Filter.lte('millis',t)).sort('millis',false);
    var right=observed.filter(ee.Filter.gte('millis',t)).sort('millis');
    var l=ee.Feature(ee.Algorithms.If(left.size(),left.first(),safe));
    var r=ee.Feature(ee.Algorithms.If(right.size(),right.first(),safe));
    var span=r.getNumber('millis').subtract(l.getNumber('millis')).divide(86400000);
    var valid=left.size().gt(0).and(right.size().gt(0)).and(span.lte(spec.gap));
    var w=t.subtract(l.getNumber('millis')).divide(r.getNumber('millis').subtract(l.getNumber('millis')).max(1));
    var red=l.getNumber('red').add(r.getNumber('red').subtract(l.getNumber('red')).multiply(w));
    var nir=l.getNumber('nir').add(r.getNumber('nir').subtract(l.getNumber('nir')).multiply(w));
    var same=l.getNumber('millis').eq(t).and(valid);
    var geValid=valid.and(l.toDictionary(['GE']).size().eq(1)).and(r.toDictionary(['GE']).size().eq(1));
    var geLeft=ee.Number(ee.Algorithms.If(geValid,l.get('GE'),0)),geRight=ee.Number(ee.Algorithms.If(geValid,r.get('GE'),0));
    return f.set({GE:ee.Algorithms.If(geValid,geLeft.add(geRight.subtract(geLeft).multiply(w)),null),GE_available_pct:ee.Algorithms.If(geValid,100,0),red:ee.Algorithms.If(valid,red,null),nir:ee.Algorithms.If(valid,nir,null),
      NDVI:ee.Algorithms.If(valid,nir.subtract(red).divide(nir.add(red).max(1e-12)),null),
      EVI2:ee.Algorithms.If(valid,nir.subtract(red).multiply(2.5).divide(nir.add(red.multiply(2.4)).add(1)).max(0).min(1),null),
      satellite_status:ee.Algorithms.If(valid,ee.Algorithms.If(same,'observed','interpolated'),'unavailable'),
      satellite_missing_pct:ee.Algorithms.If(same,0,100),satellite_unavailable_pct:ee.Algorithms.If(valid,0,100),
      satellite_gap_days:ee.Algorithms.If(valid,span,null),
      satellite_left_date:ee.Algorithms.If(valid,ee.Date(l.get('millis')).format('YYYY-MM-dd'),null),
      satellite_right_date:ee.Algorithms.If(valid,ee.Date(r.get('millis')).format('YYYY-MM-dd'),null),
      satellite_date_basis:spec.kind==='vi'?'pixel observation date':spec.days===1?'observation date':'composite midpoint',
      satellite_source:key,satellite_scale_m:spec.scale,satellite_max_gap_days:spec.gap});
  });
}
function assembleSelectedForcing(geom,start,end,config) {
  config=config||{};
  var key=config.weatherSource||'TERRACLIMATE', sat=config.satelliteSource||'MODIS16';
  var cadence=config.timestep||(config.weatherPlan?config.weatherPlan.cadence:(WEATHER_SOURCES[key].daily?'daily':'monthly'));
  var daily=cadence==='daily';
  var cal=config.weatherPlan?ee.FeatureCollection(config.weatherPlan.records.map(function(e){return towerWeatherFeature(e.data,geom);})):config.periodCalendar||sourceCalendar(start,end,daily);
  var first=ee.Date(cal.aggregate_min('millis'));
  var last=ee.Feature(cal.sort('millis',false).first());
  var lastDate=ee.Date(last.get('millis')).advance(last.getNumber('period_days').subtract(1),'day');
  var weather;
  if(config.weatherPlan) {
    weather=ee.FeatureCollection(config.weatherPlan.records.map(function(entry) {
      var r=entry.data;
      var counts=['ta','sw','ea','pa','rain','rn'];
      var available=counts.reduce(function(sum,k){return sum+(r.counts&&carbonFinite(r.counts[k])?r.counts[k]:0);},0);
      return towerWeatherFeature(r,geom).set({tower_weather_valid:entry.valid?1:0,
        weather_source_missing_pct:r.expected?100*(1-available/(counts.length*r.expected)):100,
        weather_gapfill_provenance:'supplied tower gap-filled fields retained; coverage refers to supplied records'});
    }));
  } else weather=config.sampledWeather||sampleWeather(geom,first,lastDate,key);
  var satellite=config.alignedSatellite||alignedSatellite(geom,first,lastDate,sat,config.sampledSatellite), p=mergeParams(ACCLIM_PARAMS,config.modelParams||{});
  var lat=ee.Number(geom.centroid({maxError:1}).coordinates().get(1));
  var raw=cal.map(function(f) {
    var t=f.getNumber('millis'), days=f.getNumber('period_days');
    var matches=weather.filter(ee.Filter.eq('millis',t));
    var wf=ee.Feature(ee.Algorithms.If(matches.size(),matches.first(),ee.Feature(null,{})));
    var vi=satellite.filter(ee.Filter.gte('millis',t)).filter(ee.Filter.lt('millis',t.add(days.multiply(86400000))));
    var valid=vi.filter(ee.Filter.notNull(['red','nir','NDVI','EVI2']));
    var monthly = cadence === 'monthly';
    var usableSatellite = monthly ? valid.size().gt(0) : valid.size().eq(days);
    var props=ee.Dictionary({});
    // Unavailable satellite periods have null means; Dictionary.set rejects null.
    // Omit those properties so the existing completeness check flags the period.
    var fields=['red','nir','NDVI','EVI2','GE','GE_available_pct','satellite_missing_pct','satellite_unavailable_pct'];
    if(p.faparMethod==='KNDVI_SPATIAL')fields=fields.concat(['kNDVI','kndvi_sigma']);
    fields.forEach(function(k){
      var isReflectanceOrIndex = ['red','nir','NDVI','EVI2','kNDVI','kndvi_sigma'].indexOf(k) >= 0;
      var mean = monthly && isReflectanceOrIndex
        ? valid.reduceColumns(ee.Reducer.median(), [k]).get('median')
        : vi.aggregate_mean(k);
      props=ee.Dictionary(ee.Algorithms.If(ee.Algorithms.IsEqual(mean,null),props,props.set(k,mean)));
    });
    if (config.gePeriods) {
      var geMatches = config.gePeriods.filter(ee.Filter.eq('date', f.get('date')));
      var geValue = ee.Feature(ee.Algorithms.If(geMatches.size(), geMatches.first(), ee.Feature(null, {}))).get('GE');
      var gePrime = ee.Feature(ee.Algorithms.If(geMatches.size(), geMatches.first(), ee.Feature(null, {}))).get('GE_normalized');
      props = ee.Dictionary(ee.Algorithms.If(ee.Algorithms.IsEqual(gePrime, null), props, props.set('GE_normalized', gePrime)));
      props = props.select(props.keys().remove('GE'));
      props = ee.Dictionary(ee.Algorithms.If(ee.Algorithms.IsEqual(geValue, null), props, props.set('GE', geValue)));
    }
    if(config.rtPeriods) props=props.combine(ee.Dictionary(config.rtPeriods.get(f.getString('date'),ee.Dictionary({}))),true);
    var dl=ee.Number(ee.List.sequence(0,days.subtract(1)).map(function(i){
      return astronomicalDaylength(lat,doyFromDate(ee.Date(t).advance(i,'day')));
    }).reduce(ee.Reducer.mean()));
    var complete=wf.toDictionary(['Tmin','Tmax','Tdaily','SWRdaily','ea','Pdaily']).size().eq(6).and(usableSatellite).and(dl.gt(0));
    if(p.faparMethod==='KNDVI_SPATIAL'){
      var kernelDays=vi.filter(ee.Filter.notNull(['kNDVI'])).size();
      complete=complete.and(monthly?kernelDays.gt(0):kernelDays.eq(days));
    }
    if(p.faparMethod==='RT_CHL'||p.faparMethod==='RT_LEAF') complete=complete.and(ee.Number(ee.Algorithms.If(props.contains('FAPAR_leaf'),1,0)));
    if(p.faparMethod==='PROSAIL' || p.fiparMethod==='PROSAIL') complete=complete.and(ee.Number(ee.Algorithms.If(props.contains('rt_fipar_daylight'),1,0)));
    if(p.fiparMethod==='NDVI_BLEND' && p.fge_exp>0)complete=complete.and(ee.Number(ee.Algorithms.If(props.contains('GE'),1,0)));
    if(config.weatherPlan)complete=complete.and(ee.Number(ee.Algorithms.If(ee.Algorithms.IsEqual(wf.get('tower_weather_valid'),null),0,wf.get('tower_weather_valid'))).eq(1));
    return f.setGeometry(geom).set(wf.toDictionary()).set(props).set({
      latitude:lat,daylength:dl,daylight_seconds:dl.multiply(3600).multiply(days),
      doy:doyFromDate(ee.Date(t)),dt:days,ca_daily:co2ForDate(ee.Date(t)),
      model_cadence:cadence,weather_source:config.weatherPlan?'tower':key,
      model_valid:complete,weather_missing_pct:ee.Number(7).subtract(wf.toDictionary(['Tmin','Tmax','Tdaily','SWRdaily','ea','Pdaily','precip']).size()).multiply(100/7),
      precipitation_missing_pct:ee.Algorithms.If(ee.Algorithms.IsEqual(wf.get('precip'),null),100,0),
      satellite_source:sat,satellite_scale_m:SATELLITE_SOURCES[sat].scale,
      satellite_valid_days:valid.size(),
      satellite_missing_pct:ee.Number(1).subtract(valid.size().divide(days)).multiply(100),
      satellite_unavailable_pct:ee.Number(1).subtract(valid.size().divide(days)).multiply(100),
      satellite_unobserved_pct:ee.Number(1).subtract(valid.filter(ee.Filter.eq('satellite_status','observed')).size().divide(days)).multiply(100),
      satellite_aggregation:monthly?'median of available daily values after gap filling':'mean with complete daily coverage',
      satellite_date_basis:ee.Feature(vi.first()).get('satellite_date_basis'),
      satellite_status:ee.Algorithms.If(usableSatellite,ee.Algorithms.If(valid.size().lt(days),'partial',ee.Algorithms.If(vi.aggregate_sum('satellite_missing_pct'),'interpolated','observed')),'unavailable'),
      satellite_gap_days:vi.aggregate_max('satellite_gap_days'),satellite_left_date:ee.Feature(vi.sort('millis').first()).get('satellite_left_date'),
      satellite_right_date:ee.Feature(vi.sort('millis',false).first()).get('satellite_right_date'),satellite_max_gap_days:SATELLITE_SOURCES[sat].gap});
  });
  var ready=raw.filter(ee.Filter.eq('model_valid',1)).map(function(f) {
    var tday=f.getNumber('Tmax').multiply(0.75).add(f.getNumber('Tmin').multiply(0.25));
    var ea=f.getNumber('ea'), albedo=f.getNumber('red').add(f.getNumber('nir')).divide(2).max(0).min(1);
    var fipar=p.fiparMethod==='PROSAIL' ? f.getNumber('rt_fipar_daylight') : computeFIPAR(f.getNumber('NDVI'),p,f.get('GE'));
    var em=fipar.multiply(p.em_v).add(ee.Number(1).subtract(fipar).multiply(p.em_nv));
    var kelvin=f.getNumber('Tdaily').add(273.15);
    var rn=ee.Number(1).subtract(albedo).multiply(f.getNumber('SWRdaily')).add(kelvin.pow(4).multiply(p.sb)
      .multiply(ea.multiply(10).divide(kelvin).pow(1/7).multiply(1.31).subtract(em)));
    if(config.weatherPlan)rn=f.getNumber('Rndaily');
    return f.set({Tday:tday,albedo:albedo,FIPAR:fipar,Rndaily:rn,
      Rnday:config.weatherPlan?rn.multiply(f.getNumber('period_seconds')).divide(f.getNumber('daylight_seconds')):
        griddedDaytimeRn(f.getNumber('SWRdaily'),f.getNumber('daylength'),tday,ea,albedo,em,p,f.getNumber('period_seconds'),f.getNumber('daylight_seconds')),
      esTday:magnusVP(tday),VPDdaily:magnusVP(tday).subtract(ea).max(0)});
  });
  return ready.merge(raw.filter(ee.Filter.eq('model_valid',0))).sort('millis');
}

// Fill only bucket inputs. Carbon and meteorological model outputs stay missing.
function bridgeBucketForcing(baseline) {
  var prepared=baseline.map(function(f) {
    var valid=f.toDictionary(['PT_mean','FIPAR','fwet_mean','daylength']).size().eq(4);
    var demand=ee.Number(f.get('PT_mean')).multiply(ee.Number(1).subtract(ee.Number(f.get('FIPAR'))));
    return f.set({_bucket_rain:ee.Algorithms.If(ee.Algorithms.IsEqual(f.get('precip'),null),null,f.getNumber('precip').divide(f.getNumber('period_seconds'))),
      _bucket_demand:ee.Algorithms.If(valid,demand.max(0).multiply(ee.Number(1).subtract(f.getNumber('fwet_mean'))).multiply(f.getNumber('daylight_seconds')).divide(f.getNumber('period_seconds')).divide(2450000),null),
      _bucket_surface:ee.Algorithms.If(valid,demand.multiply(f.getNumber('fwet_mean')),null)});
  });
  var fields=['_bucket_rain','_bucket_demand','_bucket_surface','PT_mean','FIPAR','fwet_mean','fSM_mean',
    'PredLEsday_mean','PredLEday_mean','PredLEcday_mean','PredLEiday_mean'];
  return prepared.map(function(f) {
    var missing=f.toDictionary(fields).size().lt(fields.length).or(ee.Number(ee.Algorithms.If(ee.Algorithms.IsEqual(f.get('bucket_forcing_missing_pct'),null),0,f.get('bucket_forcing_missing_pct'))).gt(0)), filled=f;
    fields.forEach(function(k) {
      var valid=prepared.filter(ee.Filter.notNull([k]));
      var left=valid.filter(ee.Filter.lte('millis',f.get('millis'))).sort('millis',false);
      var right=valid.filter(ee.Filter.gte('millis',f.get('millis'))).sort('millis');
      var safe=ee.Feature(null,{millis:0}).set(k,0);
      var l=ee.Feature(ee.Algorithms.If(left.size(),left.first(),safe)),r=ee.Feature(ee.Algorithms.If(right.size(),right.first(),safe));
      var w=f.getNumber('millis').subtract(l.getNumber('millis')).divide(r.getNumber('millis').subtract(l.getNumber('millis')).max(1));
      var value=l.getNumber(k).add(r.getNumber(k).subtract(l.getNumber(k)).multiply(w));
      filled=filled.set(k,ee.Algorithms.If(ee.Algorithms.IsEqual(f.get(k),null),
        ee.Algorithms.If(left.size().gt(0).and(right.size().gt(0)),value,null),f.get(k)));
    });
    return filled.set({precip:ee.Algorithms.If(ee.Algorithms.IsEqual(filled.get('_bucket_rain'),null),null,
      filled.getNumber('_bucket_rain').multiply(f.getNumber('period_seconds'))),
      bucket_forcing_missing_pct:ee.Algorithms.If(missing,100,0),
      bucket_forcing_status:ee.Algorithms.If(missing,
        ee.Algorithms.If(filled.toDictionary(fields).size().eq(fields.length),'interpolated','unavailable'),'supplied'),
      bucket_fit_observed:ee.Algorithms.If(missing,0,1)});
  });
}

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
    //   FIPAR uses the shared NDVI interception function below.
    //   em_sur     = FIPAR * em_v + (1 - FIPAR) * em_nv
    //   em_air     = 1.31 * (10 * ea / (Tdaily + 273.15))^(1/7)
    //   Rndaily    = (1 - albedo) * SWRdaily + sb * (Tdaily + 273.15)^4 * (em_air - em_sur)
    var albedo = red.add(nir).divide(2).max(0).min(1);
    var FIPAR = computeFIPAR(NDVI, p,vi.get('GE'));
    var em_sur = FIPAR.multiply(p.em_v)
      .add(ee.Number(1).subtract(FIPAR).multiply(p.em_nv));
    var em_air = ee.Number(10).multiply(ea).divide(Tdaily.add(273.15))
      .pow(1 / 7).multiply(1.31);
    var T_K4 = Tdaily.add(273.15).pow(4);
    var Rndaily = ee.Number(1).subtract(albedo).multiply(SWRdaily)
      .add(T_K4.multiply(p.sb).multiply(em_air.subtract(em_sur)));

    if (towerRows) Rndaily = f.getNumber('Rndaily');
    var Rnday = towerRows ? Rndaily.multiply(ee.Number(24).divide(dl))
      : griddedDaytimeRn(SWRdaily, dl, Tday, ea, albedo, em_sur, p);
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
      FIPAR: FIPAR,
      Rndaily:   Rndaily,
      Rnday:     Rnday,
      // MVP constants
      Pdaily:    towerRows ? f.getNumber('Pdaily') : MVP_CONSTANTS.Pdaily,
      ca_daily:  co2ForDate(monthStart)
    });
  });

  // Attach `dt` (days since previous timestep) to each feature so the model
  // is cadence-agnostic. First feature gets dt = 30 (typical month).
  var sorted = enriched.sort('millis');
  var sortedList = sorted.toList(sorted.size().max(1));
  var dayMs = 86400000;
  var withDt = ee.FeatureCollection(
    ee.List.sequence(0, sortedList.size().subtract(1)).slice(0,sortedList.size()).map(function(i) {
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
//   computeInitialState  - time-weighted first five days, used to seed running means
//   acclimStep           - one timestep's update; appends output to history
//   runACCLIM            - drives ee.List.iterate, returns output FC
//
// Running-mean rule (from R script):
//     x_mean[i] = w[i] * x_mean[i-1] + (1 - w[i]) * x[i]      if PAR >= 10
//     x_mean[i] = x_mean[i-1]                                  if PAR < 10
//     x_mean[1] = duration-weighted mean over the initialization window                                 (initial seed)
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
  var initialSeconds=(params.initialization_days===undefined?5:params.initialization_days)*86400;
  if(!(initialSeconds>0&&isFinite(initialSeconds)))throw new Error('Initialization window must be positive.');
  var seedStart=ee.Feature(forcingList.get(0)).getNumber('period_start_millis');
  var seedEnd=seedStart.add(initialSeconds*1000);
  var first5=ee.FeatureCollection(forcingList).filter(ee.Filter.lt('period_start_millis',seedEnd)).map(function(f){
    return f.set('_seed_seconds',f.getNumber('period_end_millis').min(seedEnd)
      .subtract(f.getNumber('period_start_millis')).max(0).divide(1000));
  });

  // PAR is derived from SWR each day; we need to derive it for the first 5 too
  // so we can seed PARday_mean. We'll compute on the fly here.
  var first5_with_par = first5.map(function(f) {
    var swr = f.getNumber('SWRdaily');
    var dl = f.getNumber('daylength');
    var par=swr.multiply(2.2).multiply(f.getNumber('period_seconds')).divide(f.getNumber('daylight_seconds'));
    var Tmin = f.getNumber('Tmin');
    var Tmax = f.getNumber('Tmax');
    var Tday = Tmax.add(Tmax.add(Tmin).multiply(0.5)).divide(2);
    var esTmin = magnusVP(Tmin);
    var esTmax = magnusVP(Tmax);
    var esTday = magnusVP(Tday);
    var humidityPressure = f.getNumber('ea');
    var VPDmax = esTmax.subtract(humidityPressure).max(0);
    var VPDday = esTday.subtract(humidityPressure).max(0);
    var RHmin = humidityPressure.divide(esTmax).max(0).min(1);
    var Rnday = f.getNumber('Rnday');
    return f.set({
      Tday: Tday,
      esTmin: esTmin, esTmax: esTmax, esTday: esTday,
      VPDmax: VPDmax, VPDday: VPDday,
      PARday: par, RHmin: RHmin, Rnday: Rnday
    });
  });

  function seedMean(name) {
    var weighted=first5_with_par.map(function(f){return f.set('_seed_value',f.getNumber(name).multiply(f.getNumber('_seed_seconds')));});
    return ee.Number(weighted.aggregate_sum('_seed_value')).divide(ee.Number(weighted.aggregate_sum('_seed_seconds')));
  }
  // Pmean in R is stored in Pa (1000 * Pdaily where Pdaily is kPa)
  var Pmean0      = seedMean('Pdaily').multiply(1000);
  var ca_mean0    = seedMean('ca_daily');
  var Tday_mean0  = seedMean('Tday');
  var Tmin_mean0  = seedMean('Tmin');
  var Tmax_mean0  = seedMean('Tmax');
  var esTday_mean0 = seedMean('esTday');
  var esTmax_mean0 = seedMean('esTmax');
  var VPDday_mean0 = seedMean('VPDday');
  var VPDmax_mean0 = seedMean('VPDmax');
  var PARday_mean0 = seedMean('PARday');
  var Rnday_mean0  = seedMean('Rnday');
  var RHmin_mean0  = seedMean('RHmin');

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
    RHmin_sm_mean: RHmin_mean0,
    VPDmax_sm_mean: VPDmax_mean0,
    // dt/w state
    dt:           1,
    // step counter (so the first call knows it's the seed day)
    step:         0,
    // accumulating output history (one feature per timestep)
    history:      ee.List([])
  });
}

// FIPAR is independent of the selected photosynthetic FAPAR method.
function computeFIPAR(ndvi, p, ge) {
  // Site calibration: FIPAR = clamp01(m_cover * NDVI + b_cover).
  if (p.fiparMethod === 'RT_NDVI') return ee.Number(ndvi).multiply(p.m_cover).add(p.b_cover).max(0).min(1);
  if (p.fiparMethod === 'NDVI_BLEND') {
    var n = ee.Number(ndvi).max(0);
    var fge = greenExcessMultiplier(ge, p);
    var w = ndviBlendWeight(n, p);
    return w.multiply(n.pow(p.ndvi_exp)).add(ee.Number(1).subtract(w).multiply(fge));
  }
  return ee.Number(ndvi).subtract(p.ndvi_min).max(0)
    .divide(p.ndvi_max - p.ndvi_min).pow(p.ndvi_exp);
}

/**
 * Compute mean-state-dependent diagnostics from a state dictionary +
 * today's instantaneous values. Pure function of inputs; no side effects.
 * Returns a dictionary of all derived quantities to attach to the output feature.
 */
// FAPAR menu equations. Each option first calculates B, then FAPAR = B^fexp.
// clamp01(x) = min(1, max(0, x)); n = max(0, NDVI).
// GEprime = (GE - GEmin)/(GEmax - GEmin) when GEmin < 0, otherwise GEprime = GE.
// GE bounds use valid paired observations at the selected run timestep/date range.
// FIPAR: B = selected FIPAR * FAPARmult (before the optional shared-FAPAR override).
//   Scaled NDVI FIPAR = [max(0, NDVI - ndvi_min)/(ndvi_max - ndvi_min)]^ndvi_exp.
//   Blended FIPAR = w*n^ndvi_exp + (1-w)*fge; w = n^blend_exp.
// FIPAR_GE: B = selected FIPAR with ndvi_min = intercept of NDVI versus GEprime.
// FIPAR_P05: B = selected FIPAR with ndvi_min = P05(historical NDVI).
//   EVI2min = P05(historical EVI2), using observations in the ten years before start.
// NDVI: B = NDVI^2.
// FGE: B = fge = clamp01(GE / GE_scale)^fge_exp; see greenExcessMultiplier for GE_scale.
// NDVI_BLEND: B = w*n^ndvi_exp + (1-w)*fge; w = n^blend_exp.
// GE_NDVI_MIN: B = min(fge, n^ndvi_exp), using the same fge scaling as FGE.
// EVI2_FIXED_POWER: B = clamp01((EVI2 - EVI2min)/(EVI2max - EVI2min))^EVI2exp.
// EVI2_FITTED: fit FIPAR = a + b*EVI2; lo = -a/b; hi = (1-a)/b.
//   B = clamp01((EVI2 - lo)/(hi - lo)).
// EVI2_FITTED_MIN: same fitted lo, hi = 0.6; B = clamp01((EVI2-lo)/(0.6-lo)).
// Invalid EVI2 fits fall back to B = NDVI^2. fexp is applied after that fallback too.
function photosyntheticFAPAR(instant, p) {
  // IPL-UV general RBF kNDVI: k=exp(-(NIR-red)^2/(2*sigma^2)); B=(1-k)/(1+k).
  // FAPAR = B^fexp. Inputs are the selected timestep's fractional reflectances.
  if(p.faparMethod==='KNDVI_SPATIAL')return ee.Number(instant.get('kNDVI')).pow(p.fexp);
  if (p.faparMethod === 'KNDVI') {
    var distance = ee.Number(instant.get('nir')).subtract(instant.get('red')).divide(p.kndvi_sigma);
    var kernel = distance.pow(2).multiply(-0.5).exp();
    return ee.Number(1).subtract(kernel).divide(kernel.add(1)).pow(p.fexp);
  }
  // Explicit cover option: FAPAR equals retrieved cover at the selected model timestep.
  // Photon-weighted 400–700 nm absorptance from the fitted FourSAIL foliage state.
  // Fitted FAPAR = clamp01(m_abs_index * model index + b_abs_index).
  if (p.faparMethod === 'RT_CHL_NDVI') return ee.Number(instant.get('NDVI')).multiply(p.m_abs_ndvi).add(p.b_abs_ndvi).max(0).min(1);
  if (p.faparMethod === 'RT_CHL_EVI2') return ee.Number(instant.get('EVI2')).multiply(p.m_abs_evi2).add(p.b_abs_evi2).max(0).min(1);
  if (p.faparMethod === 'RT_CHL') return ee.Number(instant.get('FAPAR_chl'));
  if (p.faparMethod === 'RT_LEAF') return ee.Number(instant.get('FAPAR_leaf'));
  if (p.faparMethod === 'PROSAIL') return ee.Number(instant.get('rt_fipar_daylight'));
  // NIRV: NIRv = NIR * NDVI; B = max(0, NIRv); FAPAR = B^fexp.
  if (p.faparMethod === 'NIRV') {
    return ee.Number(instant.get('nir')).multiply(ee.Number(instant.get('NDVI'))).max(0).pow(p.fexp);
  }
  // NDVI_RED_EXP: B = max(0,NDVI)^(R*0.5 + (1-R)*2), with R = red reflectance.
  // FAPAR = B^fexp. Reflectance uses fractions (0-1).
  if (p.faparMethod === 'NDVI_RED_EXP') {
    var red = ee.Number(instant.get('red'));
    var exponent = red.multiply(0.5).add(ee.Number(1).subtract(red).multiply(2));
    return ee.Number(instant.get('NDVI')).max(0).pow(exponent).pow(p.fexp);
  }
  // NDVI_GE_OFFSET: B = max(0, max(0, NDVI)^ndvi_exp - ndvi_min).
  // ndvi_min is the NDVI-versus-GEprime intercept fitted before forcing assembly.
  if (p.faparMethod === 'NDVI_GE_OFFSET') {
    return ee.Number(instant.get('NDVI')).max(0).pow(p.ndvi_exp).subtract(p.ndvi_min).max(0).pow(p.fexp);
  }
  if (p.faparMethod === 'GE_NDVI_MIN') {
    var geForMinimum = p.fge_exp === 0 ? 0 : instant.get('GE');
    var fgeForMinimum = greenExcessMultiplier(geForMinimum, p);
    return fgeForMinimum.min(ee.Number(instant.get('NDVI')).max(0).pow(p.ndvi_exp)).pow(p.fexp);
  }
  if (p.faparMethod === 'FGE' || p.faparMethod === 'NDVI_BLEND') {
    var ge = p.fge_exp === 0 ? 0 : instant.get('GE');
    var fge = greenExcessMultiplier(ge, p);
    if (p.faparMethod === 'NDVI_BLEND') {
      var ndvi = ee.Number(instant.get('NDVI')).max(0);
      var w = ndviBlendWeight(ndvi, p);
      return w.multiply(ndvi.pow(p.ndvi_exp)).add(ee.Number(1).subtract(w).multiply(fge)).pow(p.fexp);
    }
    return fge.pow(p.fexp);
  }
  if (p.faparMethod === 'FIPAR') {
    var originalFIPAR = instant.get(p.useSharedFPAR ? 'FIPAR_NDVI' : 'FIPAR');
    return ee.Number(originalFIPAR).multiply(directFIPARMultiplier(p)).pow(p.fexp);
  }
  var FAPAR_in = ee.Number(instant.get('EVI2'));
  var FAPAR = ee.Number(instant.get('NDVI')).pow(2);
  if (p.faparMethod === 'EVI2' && !p.fitEVI2) {
    FAPAR = FAPAR_in.subtract(p.EVI2min).divide(p.EVI2max - p.EVI2min).max(0).min(1);
    if (p.powerFixedEVI2) FAPAR = FAPAR.pow(p.EVI2exp);
  } else if (p.faparMethod === 'EVI2') {
    var fitRange = ee.Number(ee.Algorithms.If(p.faparFit.valid, p.faparFit.max.subtract(p.faparFit.min), 1));
    var scaledEVI2 = FAPAR_in.subtract(p.faparFit.min)
      .divide(fitRange).max(0).min(1);
    FAPAR = ee.Number(ee.Algorithms.If(p.faparFit.valid, scaledEVI2, FAPAR));
  }

  return FAPAR.pow(p.fexp);
}

function ndviBlendWeight(ndvi, p) {
  if (typeof p.blend_exp !== 'number' || !isFinite(p.blend_exp) || p.blend_exp < 0) {
    throw new Error('blend_exp must be finite and nonnegative.');
  }
  return ee.Number(ndvi).max(0).pow(p.blend_exp);
}

function directFIPARMultiplier(p) {
  if (p.fitNDVIminGE || p.fitVIMinP05 || p.FAPARmult === undefined) return 1;
  return p.FAPARmult;
}

function greenExcessMultiplier(ge,p) {
  if(p.fge_exp===0)return ee.Number(1);
  // Bound the base before a fractional power to avoid negative-base NaNs.
  return ee.Number(ge).divide(p.GE_scale).max(0).min(1).pow(p.fge_exp);
}
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

  // Optional FAPAR calibration (enabled through the run controls).
  // Fit paired, quality-screened observations: FIPAR = a + b * EVI2 (b > 0).
  // Derive effective endpoints EVI2min = -a/b and EVI2max = (1-a)/b, then
  // use FAPAR = clamp((EVI2-EVI2min)/(EVI2max-EVI2min), 0, 1).
  // Intent: anchor the scale to FIPAR while retaining EVI2 temporal variation
  // and its potential sensitivity to green foliage amid non-green material.
  // These are calibration endpoints, not necessarily observed soil/full-canopy values.
  // If revisited, fit fixed endpoints across multiple years with balanced VI coverage;
  // check slope/endpoints and validate on withheld years against independent data.
  // Fitting to FIPAR can inherit its saturation/bias; improvement is not guaranteed.
  // Fit once across the selected record; coefficients stay fixed within each run.
  // Previous FAPAR from scaled EVI2 (retained for comparison).
  // var FAPAR = FAPAR_in.subtract(p.EVI2min)
  //   .divide(p.EVI2max - p.EVI2min)
  //   .max(0).min(1);

  // FAPAR from squared NDVI, using the current monthly NDVI input.
  var FIPAR = ee.Number(instant.get('FIPAR'));
  var FAPAR = photosyntheticFAPAR(instant, p);

  // zC3 (R formula, units arranged for the quadratic):
  //   zC3 = ((1.5*GPPmean2GPPmax*2500*18)/(QRC3*1.26*(1+fresp))) *
  //         ((s + gamma_psy) * VPDday_mean * PARday_mean * 1000) /
  //         (s * Rnday_mean * Pmean)
  // PAR is umol photons/m2/s; 2500 J/g * 18 g/mol gives J/mol H2O.
  // VPD * 1000 / Pmean is Pa/Pa, so zC3 is umol/mol (ppm).
  // The vegetation ratio no longer cancels when FAPAR differs from FIPAR.
  // A numerical floor avoids division by zero; flag these non-identifiable cases.
  // Each exponent acts directly on a raw environmental response, once per flux.
  var vpdResponse = ee.Number(1).divide(VPDday_mean.max(1)); // VPD in kPa
  var seasonalRatio = ee.Number(instant.get('seasonal_ratio', 1)); // FAPAR / running maximum
  var fDgpp = p.fDgpp_exp === 0 ? ee.Number(1) : vpdResponse.pow(p.fDgpp_exp);
  var fDtrans = p.fDtrans_exp === 0 ? ee.Number(1) : vpdResponse.pow(p.fDtrans_exp);
  var fMgpp = p.fMgpp_exp === 0 ? ee.Number(1) : seasonalRatio.pow(p.fMgpp_exp);
  var fMtrans = p.fMtrans_exp === 0 ? ee.Number(1) : seasonalRatio.pow(p.fMtrans_exp);
  var gppFactor = fDgpp.multiply(fMgpp);
  var transpirationFactor = fDtrans.multiply(fMtrans);
  var validTransRatio = transpirationFactor.gt(0);
  // At zero T multiplier, A/T is unidentified. Retain the flagged internal fallback.
  var carbonWaterFactor = ee.Number(ee.Algorithms.If(validTransRatio,
    gppFactor.divide(ee.Number(ee.Algorithms.If(validTransRatio, transpirationFactor, 1))), 1));
  var faparFiparValid = FIPAR.gt(1e-9);
  var faparFiparRatio = (p.useSharedFPAR || (p.faparMethod === 'FIPAR' && p.fexp === 1 && directFIPARMultiplier(p) === 1)) ? ee.Number(1) : FAPAR.divide(FIPAR.max(1e-9));
  var zCoef = ee.Number(1.5).multiply(p.GPPmean2GPPmax).multiply(2500).multiply(18)
    .divide(ee.Number(p.QRC3).multiply(p.PT_coef).multiply(ee.Number(1).add(p.fresp)));
  var zMet = zCoef.multiply(
    s_mean.add(p.gamma_psy).multiply(VPDday_mean).multiply(PARday_mean).multiply(1000)
  ).divide(s_mean.multiply(Rnday_mean).multiply(Pmean));
  // zC3 = meteorological term * FAPAR/FIPAR * GPP factors / transpiration factors.
  var zC3 = zMet.multiply(faparFiparRatio).multiply(carbonWaterFactor);

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

  // GPPmax (umolCO2/m2s) = QY * FAPAR * PARday_mean * fDgpp * fMgpp
  var FAPARmax=ee.Number(ee.Algorithms.If(instant.contains('FAPARmax'),instant.get('FAPARmax'),FAPAR));
  var fge=greenExcessMultiplier(instant.get('GE',0),p);
  var GPPmax_C3 = QYC3_mean.multiply(FAPAR).multiply(PARday_mean).multiply(gppFactor);
  var GPPmax_C4 = QYC4_mean.multiply(FAPAR).multiply(PARday_mean).multiply(gppFactor);

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
  C4frac = ee.Number(0); // C3-only test; comment out to restore modeled C4 fraction

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
  var fSM_mean = ee.Number(running.get('RHmin_sm_mean'))
    .pow(ee.Number(running.get('VPDmax_sm_mean')));
  var PT_mean = ee.Number(p.PT_coef).multiply(s_ratio_mean).multiply(Rnday_mean);

  var PredLEcday_mean = PT_mean.multiply(FIPAR).multiply(fcold_mean)
    .multiply(ee.Number(1).subtract(fwet_mean)).multiply(transpirationFactor);
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
    T_multiplier: transpirationFactor, GPP_multiplier: gppFactor,
    zMet: zMet, carbon_water_factor: carbonWaterFactor, seasonal_ratio: seasonalRatio,
    fDgpp: fDgpp, fDtrans: fDtrans, fMgpp: fMgpp, fMtrans: fMtrans,
    fDgpp_exp: p.fDgpp_exp, fDtrans_exp: p.fDtrans_exp,
    fMgpp_exp: p.fMgpp_exp, fMtrans_exp: p.fMtrans_exp,
    ci_T_multiplier_valid: validTransRatio,
    FAPAR_FIPAR_ratio: faparFiparRatio,
    FAPAR_FIPAR_valid: faparFiparValid,
    bC3:              bC3,
    cC3:              cC3,
    ciC3_ppm:        ciC3_ppm,
    gamma_ppm:       gamma_ppm,
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
    FAPAR:            FAPAR, fexp:p.fexp, blend_exp:p.blend_exp,
    FAPARmax: FAPARmax, fM: fMgpp, // historical output alias for the GPP seasonal factor
    fge:ee.Algorithms.If(instant.contains('GE'),fge,p.fge_exp===0?1:null), fge_exp:p.fge_exp, GE_for_GPP:ee.Algorithms.If(instant.contains('GE'),instant.get('GE',0),null),
    FIPAR:            FIPAR, FIPAR_method:p.fiparMethod,
    FPAR_mode: p.useSharedFPAR ? 'shared FAPAR' : 'separate FAPAR/FIPAR',
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
// Reset in the period containing Jan 1 (north/equator) or DOY 182 (south).
// Missing periods retain no invented peaks; first valid period seeds each season.
function seasonalFAPARState(fapar,latitude,periodEnd,state) {
  var end=ee.Date(periodEnd).advance(-1,'second'),year=ee.Number(end.get('year'));
  var doy=ee.Number(end.getRelative('day','year')).add(1);
  var season=year.subtract(ee.Number(latitude).lt(0).and(doy.lt(182)));
  var reset=ee.Number(state.get('step')).eq(0).or(season.neq(ee.Number(ee.Algorithms.If(state.contains('FAPAR_season'),state.get('FAPAR_season'),season))));
  var previous=ee.Number(ee.Algorithms.If(state.contains('FAPARmax'),state.get('FAPARmax'),fapar));
  var maximum=ee.Number(ee.Algorithms.If(reset,fapar,previous.max(fapar)));
  var ratio=ee.Number(ee.Algorithms.If(maximum.gt(0),ee.Number(fapar).divide(ee.Number(ee.Algorithms.If(maximum.gt(0),maximum,1))),0));
  return ee.Dictionary({FAPARmax:maximum,seasonal_ratio:ratio,FAPAR_season:season,FAPARmax_reset:reset});
}
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
  var humidityPressure = f.getNumber('ea');
    var VPDmax = esTmax.subtract(humidityPressure).max(0);
    var VPDday = esTday.subtract(humidityPressure).max(0);

  // PARday (umol/m2/s) = 2.2 * SWR(W/m2) * (24/dl)
  var PARday=SWRdaily.multiply(2.2).multiply(f.getNumber('period_seconds')).divide(f.getNumber('daylight_seconds'));

  // Use source-specific daytime energy prepared in assembleForcing.
  var Rnday = f.getNumber('Rnday');

  // RHmin (fraction)
  var RHmin = humidityPressure.divide(esTmax).max(0).min(1);

  // Acclimation always uses elapsed seconds. Low light holds the state, not time.
  var step=ee.Number(state.get('step')),isFirst=step.eq(0),lowLight=PARday.lt(10);
  var elapsedSeconds=f.getNumber('elapsed_seconds'),dt=elapsedSeconds.divide(86400);
  var w=elapsedSeconds.multiply(-1).divide(ee.Number(p.tc).multiply(86400)).exp();
  var w_sm=elapsedSeconds.multiply(-1).divide(ee.Number(p.tc_sm).multiply(86400)).exp();

  // Running-mean update closure
  function updateMean(prevKey, todayVal, weight) {
    weight = weight === undefined ? w : weight;
    var prev = ee.Number(state.get(prevKey));
    return ee.Number(ee.Algorithms.If(
      isFirst, prev,                      // seed value retained on first step
      ee.Algorithms.If(
        lowLight, prev,                   // carry forward
        weight.multiply(prev).add(ee.Number(1).subtract(weight).multiply(todayVal))
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
  var RHmin_sm_mean = updateMean('RHmin_sm_mean', RHmin, w_sm);
  var VPDmax_sm_mean = updateMean('VPDmax_sm_mean', VPDmax, w_sm);

  var runningNow = ee.Dictionary({
    Pmean: Pmean, ca_mean: ca_mean, Tday_mean: Tday_mean,
    Tmin_mean: Tmin_mean, Tmax_mean: Tmax_mean, esTmax_mean: esTmax_mean, esTday_mean: esTday_mean,
    VPDday_mean: VPDday_mean, VPDmax_mean: VPDmax_mean,
    PARday_mean: PARday_mean, Rnday_mean: Rnday_mean,
    RHmin_mean: RHmin_mean,
    RHmin_sm_mean: RHmin_sm_mean, VPDmax_sm_mean: VPDmax_sm_mean
  });
  var instant = ee.Dictionary({kNDVI:f.get('kNDVI'),FAPAR_leaf:f.get('FAPAR_leaf'),FAPAR_chl:f.get('FAPAR_chl'),rt_fipar_daylight:f.get('rt_fipar_daylight'), EVI2: EVI2, NDVI: NDVI, red:f.get('red'),nir:f.get('nir'),FIPAR: f.getNumber('FIPAR'),FIPAR_NDVI:p.useSharedFPAR?f.get('FIPAR_NDVI'):f.get('FIPAR'),GE:f.get('GE'),GE_normalized:f.get('GE_normalized')});

  var latitude=ee.Number(ee.Algorithms.If(f.propertyNames().contains('latitude'),f.get('latitude'),
    f.geometry().centroid(1).coordinates().get(1)));
  var seasonal=seasonalFAPARState(photosyntheticFAPAR(instant,p),latitude,f.getNumber('period_end_millis'),state);
  instant=instant.combine(seasonal);
  var diag = computeDiagnostics(runningNow, instant, p);

  // Output feature: forcings + instantaneous + running means + diagnostics
  var outProps = ee.Dictionary({
    date: f.getString('date'),
    millis: f.getNumber('millis'),
    doy: doy,
    latitude:latitude,FAPAR_season:seasonal.get('FAPAR_season'),FAPARmax_reset:seasonal.get('FAPARmax_reset'),
    // Forcings (carried through for reference)
    weather_source: f.get('weather_source'), weather_status: f.get('weather_status'),
    weather_energy_basis: f.get('weather_energy_basis'), weather_expected_intervals: f.get('weather_expected_intervals'), ea: f.get('ea'),
    Tmin: Tmin, Tmax: Tmax, Tdaily: Tdaily, Tday: Tday,
    SWRdaily: SWRdaily, Pdaily: Pdaily, Rndaily: Rndaily,
    precip: f.get('precip'), period_days: f.get('period_days'),
    ca_daily: ca_daily, EVI2: EVI2, NDVI: NDVI, daylength: dl,
    FIPAR_NDVI: p.useSharedFPAR ? f.get('FIPAR_NDVI') : f.get('FIPAR'),
    // Instantaneous derived
    esTmin: esTmin, esTmax: esTmax, esTday: esTday,
    VPDmax: VPDmax, VPDday: VPDday, PARday: PARday, Rnday: Rnday, RHmin: RHmin,
    RHmax: humidityPressure.divide(esTmin).max(0).min(1), VPDmin: esTmin.subtract(humidityPressure).max(0),
    dt: dt, w: w, w_sm: w_sm, tc: p.tc, tc_sm: p.tc_sm,
    initialization_seconds:ee.Number(p.initialization_days===undefined?5:p.initialization_days).multiply(86400)
  }).combine(runningNow).combine(diag);

  // Daytime mean flux integrated over daylight, zero GPP at night.
  // 12 g C/mol (not 44 g CO2/mol); no extra mean-to-maximum adjustment.
  var GPPperiod=ee.Number(diag.get('PredGPPmean')).multiply(f.getNumber('daylight_seconds')).multiply(12e-6);
  var GPPdaily=GPPperiod.multiply(86400).divide(f.getNumber('period_seconds'));
  outProps = outProps.combine(ee.Dictionary({
    GPPdaily: GPPdaily, NPPdaily: GPPdaily.multiply(p.CUEplant), CUEplant: p.CUEplant,
    GPP_period_gC_m2:GPPperiod,NPP_period_gC_m2:GPPperiod.multiply(p.CUEplant),
    carbon_integration:'daytime mean flux multiplied by daylight seconds'
  }));
  var outFeat = ee.Feature(f.geometry(), f.toDictionary().combine(outProps, true));

  // Return updated state
  return ee.Dictionary({
    Pmean: Pmean, ca_mean: ca_mean, Tday_mean: Tday_mean,
    Tmin_mean: Tmin_mean, Tmax_mean: Tmax_mean, esTmax_mean: esTmax_mean, esTday_mean: esTday_mean,
    VPDday_mean: VPDday_mean, VPDmax_mean: VPDmax_mean,
    PARday_mean: PARday_mean, Rnday_mean: Rnday_mean,
    RHmin_mean: RHmin_mean,
    RHmin_sm_mean: RHmin_sm_mean, VPDmax_sm_mean: VPDmax_sm_mean,
    dt: dt,
    FAPARmax:seasonal.get('FAPARmax'),FAPAR_season:seasonal.get('FAPAR_season'),
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
    var seconds=r.getNumber('period_seconds');
    var conversion=r.getNumber('daylight_seconds').divide(2450000);
    var wet = r.getNumber('fwet_mean');
    var potential = r.getNumber('PT_mean').multiply(ee.Number(1).subtract(r.getNumber('FIPAR')));
    var surfaceLE = r.getNumber('_bucket_surface');
    var potentialMM = r.getNumber('_bucket_demand').multiply(seconds);
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
    var weight = ee.Number(ee.Algorithms.If(potential.gt(0).and(r.getNumber('bucket_fit_observed').eq(1)), seconds, 0));
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
  baseline = bridgeBucketForcing(baseline.map(withPeriodTiming).sort('millis'));
  var checked = baseline.map(function(f) {
    var required = ['precip', 'period_days', 'daylength', 'fwet_mean', 'PT_mean',
      'FIPAR', 'fSM_mean', 'PredLEsday_mean', 'PredLEday_mean', 'PredLEcday_mean', 'PredLEiday_mean', '_bucket_demand', '_bucket_surface'];
    var complete = f.toDictionary(required).size().eq(required.length);
    // Defaults only protect validation expressions; missing values are never filled.
    var d = ee.Dictionary.fromLists(required, ee.List.repeat(-9999, required.length))
      .combine(f.toDictionary(required));
    var valid = complete.and(ee.Number(d.get('precip')).gte(0))
      .and(ee.Number(d.get('period_days')).gt(0))
      .and(ee.Number(d.get('daylength')).gt(0)).and(ee.Number(d.get('daylength')).lte(24))
      .and(ee.Number(d.get('fwet_mean')).gte(0)).and(ee.Number(d.get('fwet_mean')).lte(1))
      .and(ee.Number(d.get('FIPAR')).gte(0)).and(ee.Number(d.get('FIPAR')).lte(1))
      .and(ee.Number(d.get('fSM_mean')).gte(0)).and(ee.Number(d.get('fSM_mean')).lte(1));
    return f.set('_bucket_valid', valid);
  });
  var validRows = checked.filter(ee.Filter.eq('_bucket_valid', 1));
  var n = validRows.size();
  var first = ee.Date(ee.Algorithms.If(n,validRows.aggregate_min('millis'),0));
  var lastRow = ee.Feature(ee.Algorithms.If(n,validRows.sort('millis',false).first(),ee.Feature(null,{millis:0,period_days:0})));
  var last = ee.Date(lastRow.get('millis')).advance(lastRow.getNumber('period_days'),'day');
  var duration = ee.Number(validRows.aggregate_sum('period_days'));
  var contiguous=last.difference(first,'second').subtract(duration.multiply(86400)).abs().lt(1e-6);
  var eligible = duration.gte(365).and(contiguous);
  var missingDays = baseline.map(function(f){return f.set('_missing_days',f.getNumber('period_seconds').multiply(f.getNumber('bucket_forcing_missing_pct')).divide(100));}).aggregate_sum('_missing_days');
  var missingPct = ee.Number(missingDays).divide(ee.Number(baseline.aggregate_sum('period_seconds')).max(1e-9)).multiply(100);
  var rows = validRows.toList(validRows.size().max(1));
  var scoringRows=validRows.select(['precip','period_days','period_seconds','daylight_seconds','daylength','fwet_mean','PT_mean','FIPAR','PredLEsday_mean','_bucket_surface','_bucket_demand','bucket_fit_observed']).toList(validRows.size().max(1));
  function score(cap) {
    cap = ee.Number(cap);
    var trajectory = soilBucketTrajectory(scoringRows, cap, false);
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
  var fallback = validRows.map(function(f) {
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
  return ee.FeatureCollection(ee.Algorithms.If(apply, ee.FeatureCollection(ee.List(fitted.get('history'))), fallback)).merge(checked.filter(ee.Filter.eq('_bucket_valid',0))).sort('millis')
    .map(function(f) { return f.set({SM_fit_status: status,
      SM_fit_RMSE: ee.Algorithms.If(apply, best.getNumber('mse').sqrt(), null),
      SM_initial_guess: cfg.initial, SM_search_min: cfg.minimum, SM_search_max: cfg.maximum,
      SM_initialization: 'full at first bracketed forcing; no spin-up', precip_units: 'mm/period', bucket_record_missing_pct: missingPct}); });
}


// Separate ET-effective diagnostic. Fit the FIPAR-free moisture equation to
// ORIGINAL full-canopy total LE, not the shallow-bucket-adjusted total LE.
// Only the isolated fitting proxy sets FIPAR=0. Copy diagnostics back explicitly,
// preserving actual FIPAR, FAPAR, GPP, transpiration, and all near-surface outputs.
function fitETBucket(nearSurface, options) {
  nearSurface = nearSurface.sort('millis');
  var proxy = nearSurface.map(function(f) {
    var properties = f.toDictionary(['date', 'millis', 'precip', 'period_days','period_seconds','daylight_seconds',
      'daylength', 'fwet_mean', 'PT_mean', 'fSM_mean','model_valid','bucket_forcing_missing_pct','bucket_fit_observed']);
    return ee.Feature(null, properties.combine(ee.Dictionary({FIPAR: 0,
      PredLEsday_mean: f.get('PredLEday_RHVPD'),
      PredLEday_mean: f.get('PredLEday_RHVPD'),
      PredLEcday_mean: 0, PredLEiday_mean: 0})));
  });
  var fitted = fitSoilBucket(proxy, options);
  var pairs = nearSurface.toList(nearSurface.size().max(1)).zip(fitted.toList(fitted.size().max(1)));
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
  var texture = ['sand', 'silt', 'clay'].map(function(name) {
    var pct = info[name + '_0_30cm_pct'];
    return name + ' ' + (carbonFinite(pct) ? pct.toFixed(1) + '%' : 'unavailable');
  }).join(', ');
  return ui.Label({value: 'SoilGrids (0–30 cm): ' + texture + '\nAWC: ' +
    (carbonFinite(info.AWC_pct_vol) ? info.AWC_pct_vol.toFixed(2) + '%vol' : 'unavailable') + '\n' +
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
function fitEVI2FAPAR(forcingFC, params) {
  var p = mergeParams(ACCLIM_PARAMS, params || {});
  var pairs = forcingFC.filter(ee.Filter.notNull(['NDVI', 'EVI2']))
    .filter(ee.Filter.gte('NDVI', -1)).filter(ee.Filter.lte('NDVI', 1))
    .map(function(f) { return f.set('FIPAR_fit_target', p.fiparMethod==='PROSAIL' ? f.getNumber('rt_fipar_daylight') : computeFIPAR(f.getNumber('NDVI'), p,f.get('GE'))); });
  var fit = ee.Dictionary(pairs.reduceColumns(ee.Reducer.linearFit(), ['EVI2', 'FIPAR_fit_target']));
  var a = ee.Number(ee.Algorithms.If(ee.Algorithms.IsEqual(fit.get('offset'), null), 0, fit.get('offset')));
  var b = ee.Number(ee.Algorithms.If(ee.Algorithms.IsEqual(fit.get('scale'), null), 0, fit.get('scale')));
  var valid = pairs.size().gte(3).and(b.gt(1e-8));
  var safeB = ee.Number(ee.Algorithms.If(valid, b, 1));
  return {a: a, b: b, valid: valid, count: pairs.size(),
    min: a.multiply(-1).divide(safeB), max: ee.Number(1).subtract(a).divide(safeB)};
}

// Preserve the regression-derived lower endpoint while fixing the upper endpoint.
function fixedUpperEVI2Fit(fit) {
  var upper = ee.Number(0.6);
  return {a: fit.a, b: fit.b, count: fit.count, min: fit.min, max: upper,
    valid: ee.Number(fit.valid).and(upper.subtract(fit.min).gt(1e-9))};
}

// Called after fitting, before initialization: use the same fraction everywhere.
function applySharedFPAR(f, p) {
  var rawFIPAR = f.getNumber('FIPAR');
  var shared = photosyntheticFAPAR(f, mergeParams(p,{useSharedFPAR:false}));
  var emSurface = shared.multiply(p.em_v).add(ee.Number(1).subtract(shared).multiply(p.em_nv));
  var kelvin = f.getNumber('Tdaily').add(273.15);
  var emAir = ee.Number(10).multiply(f.getNumber('ea')).divide(kelvin).pow(1 / 7).multiply(1.31);
  var rnDaily = ee.Number(1).subtract(f.getNumber('albedo')).multiply(f.getNumber('SWRdaily'))
    .add(kelvin.pow(4).multiply(p.sb).multiply(emAir.subtract(emSurface)));
  var rnDay = griddedDaytimeRn(f.getNumber('SWRdaily'), f.getNumber('daylength'),
    f.getNumber('Tday'), f.getNumber('ea'), f.getNumber('albedo'), emSurface, p);
  var tower = ee.Algorithms.IsEqual(f.get('weather_source'), 'tower');
  return f.set({FIPAR_NDVI: rawFIPAR, FIPAR: shared,
    Rndaily: ee.Algorithms.If(tower, f.get('Rndaily'), rnDaily),
    Rnday: ee.Algorithms.If(tower, f.get('Rnday'), rnDay)});
}

function runACCLIM(forcingFC, params) {
  // Merge param overrides (ES5-safe; GEE Code Editor lacks Object.assign).
  var p = mergeParams(ACCLIM_PARAMS, params || {});
  if(typeof p.fexp!=='number'||!isFinite(p.fexp)||p.fexp<0)throw new Error('fexp must be finite and nonnegative.');
  if(typeof p.fge_exp!=='number'||!isFinite(p.fge_exp)||p.fge_exp<0)throw new Error('fge_exp must be finite and nonnegative.');
  forcingFC=forcingFC.map(withPeriodTiming);
  if((p.faparMethod==='FGE' || p.faparMethod==='GE_NDVI_MIN' || p.faparMethod==='NDVI_BLEND' || p.fiparMethod==='NDVI_BLEND') && p.fge_exp>0)forcingFC=forcingFC.map(function(f){
    var hasGE=f.toDictionary(['GE']).size().gt(0);
    return f.set({model_valid:ee.Algorithms.If(hasGE,f.get('model_valid'),0),GE_missing_for_GPP:hasGE.not()});
  });
  if (p.faparMethod === 'KNDVI' && !(typeof p.kndvi_sigma==='number' && isFinite(p.kndvi_sigma) && p.kndvi_sigma>0))
    throw new Error('kNDVI sigma must be finite and positive.');
  var bands = {NIRV:['nir'], NDVI_RED_EXP:['red'], KNDVI:['red','nir']}[p.faparMethod];
  if (bands) forcingFC = forcingFC.map(function(f) {
    var available = f.toDictionary(bands).size().eq(bands.length);
    return f.set('model_valid', ee.Algorithms.If(available, f.get('model_valid'), 0));
  });
  var missingForcing = forcingFC.filter(ee.Filter.eq('model_valid',0));
  if (bands) forcingFC = forcingFC.filter(ee.Filter.notNull(bands));
  forcingFC = forcingFC.filter(ee.Filter.notNull(['Tmin','Tmax','SWRdaily','ea','NDVI','EVI2','Rnday']));
  if((p.faparMethod==='FGE' || p.faparMethod==='GE_NDVI_MIN' || p.faparMethod==='NDVI_BLEND' || p.fiparMethod==='NDVI_BLEND') && p.fge_exp>0)forcingFC=forcingFC.filter(ee.Filter.notNull(['GE']));
  if (p.faparMethod !== 'NDVI' && p.faparMethod !== 'EVI2' && p.faparMethod !== 'FIPAR' && p.faparMethod !== 'FGE' && p.faparMethod !== 'NDVI_BLEND' && p.faparMethod !== 'GE_NDVI_MIN' && p.faparMethod !== 'NDVI_GE_OFFSET' && p.faparMethod !== 'NIRV' && p.faparMethod !== 'NDVI_RED_EXP' && p.faparMethod !== 'PROSAIL' && p.faparMethod !== 'KNDVI' && p.faparMethod !== 'KNDVI_SPATIAL' && !rtAbsorptionMethod(p.faparMethod)) throw new Error('Unknown FAPAR method.');
  if (p.faparMethod === 'EVI2' && p.fitEVI2 && !p.faparFit) p.faparFit = fitEVI2FAPAR(forcingFC, p);
  if (p.faparMethod === 'EVI2' && p.fitEVI2 && p.fitEVI2maxFixed) p.faparFit = fixedUpperEVI2Fit(p.faparFit);
  if (p.faparMethod === 'EVI2' && !p.fitEVI2 &&
      (!isFinite(p.EVI2min) || !isFinite(p.EVI2max) || p.EVI2max <= p.EVI2min)) {
    throw new Error('Fixed EVI2 bounds must be finite with EVI2max > EVI2min.');
  }
  ['fDgpp_exp','fDtrans_exp','fMgpp_exp','fMtrans_exp'].forEach(function(key) {
    if (typeof p[key] !== 'number' || !isFinite(p[key]) || p[key] < 0) throw new Error(key + ' must be finite and nonnegative.');
  });
  if (typeof p.k !== 'number' || !isFinite(p.k) || p.k <= 0) throw new Error('k must be a positive finite number.');
  if (typeof p.tc_sm !== 'number' || !isFinite(p.tc_sm) || p.tc_sm <= 0) throw new Error('tc_sm must be positive and finite.');
  if (p.useSharedFPAR) forcingFC = forcingFC.map(function(f) { return applySharedFPAR(f, p); });
  var ordered = forcingFC.sort('millis').toList(forcingFC.size().max(1));
  // Retain elapsed time across unavailable physiological periods; do not invent carbon.
  var forcingList = ee.List.sequence(0,ordered.size().subtract(1)).slice(0,ordered.size()).map(function(i) {
    i=ee.Number(i);var f=ee.Feature(ordered.get(i));
    var previous=ee.Feature(ordered.get(i.subtract(1).max(0)));
    var elapsed=ee.Number(ee.Algorithms.If(i.eq(0),f.get('period_seconds'),
      f.getNumber('period_end_millis').subtract(previous.getNumber('period_end_millis')).divide(1000)));
    return f.set({dt:elapsed.divide(86400),elapsed_seconds:elapsed});
  });
  var initial = computeInitialState(forcingList, p);
  var finalState = ee.Dictionary(ee.Algorithms.If(forcingList.size().gt(0), forcingList.iterate(function(curr, prev) {
    return acclimStep(curr, prev, p);
  }, initial), ee.Dictionary({history:ee.List([])})));
  var baseline = ee.FeatureCollection(ee.List(finalState.get('history'))).merge(missingForcing).sort('millis');
  if (p.deferBuckets) return baseline;
  return finishACCLIM(baseline,p);
}

function finishACCLIM(baseline, params) {
  var p=mergeParams(ACCLIM_PARAMS,params||{});
  if(p.faparMethod==='EVI2' && p.fitEVI2 && !p.faparFit) {
    p.faparFit=fitEVI2FAPAR(baseline.filter(ee.Filter.notNull(['PredGPPmean'])),p);
    if(p.fitEVI2maxFixed)p.faparFit=fixedUpperEVI2Fit(p.faparFit);
  }
  var nearSurface = p.fittedBuckets ? null : fitSoilBucket(baseline, p.soilBucket);
  // Convert final partitioned fluxes after the soil-bucket update. Monthly rows
  // report equivalent mean daily water depths using their computed daylength;
  // these are not monthly totals. Nighttime flux is assumed zero.
  var bucketOutput = p.fittedBuckets || fitETBucket(nearSurface, p.etBucket);
  var missingOutput = bucketOutput.filter(ee.Filter.eq('model_valid',0)).map(function(f) {
    // Bucket interpolation must never manufacture carbon or other physiology.
    var keep=f.propertyNames().filter(ee.Filter.stringStartsWith('item','SM'))
      .cat(f.propertyNames().filter(ee.Filter.stringStartsWith('item','bucket_')))
      .cat(['date','millis','period_days','daylength','period_seconds','daylight_seconds','period_start_millis','period_end_millis','period_end_exclusive','elapsed_seconds','model_cadence','model_valid','weather_source',
        'weather_missing_pct','weather_dataset','weather_scale_m','weather_source_missing_pct','weather_gapfill_provenance',
        'weather_source_cadence','weather_temporally_allocated_pct','weather_temporal_method','weather_missing_day_pct',
        'precipitation_missing_pct','satellite_source','satellite_scale_m','satellite_missing_pct','satellite_unavailable_pct','satellite_unobserved_pct','satellite_valid_days','satellite_aggregation','satellite_status',
        'satellite_date_basis','satellite_gap_days','satellite_left_date','satellite_right_date','satellite_max_gap_days',
        'Tmin','Tmax','Tdaily','SWRdaily','ea','Pdaily','NDVI','EVI2','red','nir','GE','GE_available_pct',
        'FAPAR_leaf','FAPAR_chl','rt_diffuse_par_fraction','rt_abs_observed_days','rt_abs_interpolated_days','rt_abs_missing_pct','rt_fipar_daylight','rt_fcover','rt_observed_days','rt_interpolated_days','rt_missing_pct','precip','precip_units','fSM_mean','fSM_et_mean']);
    return ee.Feature(f.geometry(),f.toDictionary(keep));
  });
  return bucketOutput.filter(ee.Filter.notNull(['PredGPPmean'])).map(function(f) {
    f = f.set('FAPAR_method_requested', p.faparMethod);
    if (p.faparMethod === 'EVI2' && !p.fitEVI2) {
      f = f.set({FAPAR_method: 'EVI2', FAPAR_fit_status: 'fixed bounds',
        FAPAR_EVI2min: p.EVI2min, FAPAR_EVI2max: p.EVI2max,
        FAPAR_EVI2exp_used: p.powerFixedEVI2 ? p.EVI2exp : 1});
    } else if (p.faparMethod === 'EVI2') {
      f = f.set({FAPAR_method: ee.Algorithms.If(p.faparFit.valid, 'EVI2', 'NDVI'),
        FAPAR_fit_status: ee.Algorithms.If(p.faparFit.valid, p.fitEVI2maxFixed ? 'fitted minimum; fixed maximum 0.6' : 'fitted', 'NDVI fallback: invalid fit or endpoint range'),
        FAPAR_fit_target: 'FIPAR',
        FAPAR_fit_n: p.faparFit.count,
        FAPAR_fit_a: p.faparFit.a, FAPAR_fit_b: p.faparFit.b,
        FAPAR_EVI2min: ee.Algorithms.If(p.faparFit.valid, p.faparFit.min, null),
        FAPAR_EVI2max: ee.Algorithms.If(p.faparFit.valid, p.faparFit.max, null)});
    } else { f = f.set({FAPAR_method: p.faparMethod, FAPAR_fit_status: 'not requested'}); }
    f = addPotentialFluxDiagnostics(f, p);
    f = addTemperatureDiagnostics(f, p);
    var periodFactor=f.getNumber('daylight_seconds').divide(2450000);
    var factor=periodFactor.multiply(86400).divide(f.getNumber('period_seconds'));
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
      ET_period_mm:f.getNumber('PredLEday_mean').multiply(periodFactor),
      T_period_mm:f.getNumber('PredLEcday_mean').multiply(periodFactor),
      E_period_mm:f.getNumber('PredLEsday_mean').multiply(periodFactor),
      I_period_mm:f.getNumber('PredLEiday_mean').multiply(periodFactor),
      ET_mm_day: f.getNumber('PredLEday_mean').multiply(factor),
      T_mm_day: f.getNumber('PredLEcday_mean').multiply(factor),
      E_mm_day: f.getNumber('PredLEsday_mean').multiply(factor),
      I_mm_day: f.getNumber('PredLEiday_mean').multiply(factor)
    });
  }).merge(missingOutput).sort('millis');
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
    var images = collection.filterDate(month, month.advance(f.getNumber('period_days'), 'day'));
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
      LST_sampling: 'Selected-period mean of good-quality daily daytime LST; valid days averaged over AOI'});
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
    var images = collection.filterDate(month, month.advance(f.getNumber('period_days'), 'day'));
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
      Aqua_LST_sampling: 'Selected-period mean of good-quality daily daytime LST; valid days averaged over AOI'});
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

// Historical minima use observed location-mean red/NIR, not interpolated timesteps.
function historicalVIWindow(start) {
  var end=carbonDateMillis(start),d=new Date(end),year=d.getUTCFullYear()-10;
  var first=Date.UTC(year,d.getUTCMonth(),d.getUTCDate());
  if(new Date(first).getUTCMonth()!==d.getUTCMonth())first=Date.UTC(year,d.getUTCMonth()+1,0);
  return {start:carbonDateString(first),end:carbonDateString(end-86400000),stop:end};
}
function historicalVIMinima(rows,start,params) {
  var window=historicalVIWindow(start),first=carbonDateMillis(window.start),groups={};
  rows.forEach(function(r){
    if(!carbonFinite(r.millis)||r.millis<first||r.millis>=window.stop||!carbonFinite(r.red)||!carbonFinite(r.nir)||r.red<0||r.nir<0||r.red>1||r.nir>1||r.red+r.nir<=0)return;
    var k=String(r.millis);if(!groups[k])groups[k]={millis:r.millis,red:0,nir:0,n:0};
    groups[k].red+=r.red;groups[k].nir+=r.nir;groups[k].n++;
  });
  var values=Object.keys(groups).map(function(k){var g=groups[k],r=g.red/g.n,n=g.nir/g.n;
    return {millis:g.millis,ndvi:(n-r)/(n+r),evi2:2.5*(n-r)/(n+2.4*r+1)};
  }).sort(function(a,b){return a.millis-b.millis;});
  if(values.length<3)throw new Error('Ten-year percentile baseline requires at least three valid historical observation dates. Choose a later start date or another satellite source.');
  function p05(field){var v=values.map(function(r){return r[field];}).sort(function(a,b){return a-b;}),i=.05*(v.length-1),j=Math.floor(i);return v[j]+(v[Math.ceil(i)]-v[j])*(i-j);}
  var ndvi=p05('ndvi'),evi2=p05('evi2');
  if(ndvi>=params.ndvi_max)throw new Error('Historical NDVI fifth percentile must be below NDVImax.');
  return {ndvi:ndvi,evi2:evi2,n:values.length,start:window.start,end:window.end,
    first:carbonDateString(values[0].millis),last:carbonDateString(values[values.length-1].millis)};
}
function fetchHistoricalVIMinima(geom,start,key,params,request,done) {
  var window=historicalVIWindow(start),cursor=carbonDateMillis(window.start)-32*86400000,all=[];
  function batch(){
    if(request!==UI._requestId)return;
    try{
      if(cursor>=window.stop){done(historicalVIMinima(all,start,params));return;}
      var next=Math.min(cursor+90*86400000,window.stop),a=carbonDateString(cursor),b=carbonDateString(next-86400000);
      setStatus('Sampling ten-year vegetation-index baseline: '+a+' to '+b+'...');
      var fc=sampleSatellite(geom,a,b,key,true);
      fc.toList(fc.size().max(1)).map(function(f){return ee.Feature(f).toDictionary();}).evaluate(function(rows,error){
        if(request!==UI._requestId)return;
        if(error){done(null,error);return;}all=all.concat(rows||[]);cursor=next;batch();
      });
    }catch(error){done(null,error.message);}
  }
  batch();
}

// Keep up to five input bundles in this app session. The batch runner supplies a disk store.
var recentInputBundles = [];
var sessionInputStore = {
  get:function(key) {
    for (var i=0;i<recentInputBundles.length;i++) {
      if (recentInputBundles[i].key===key) return recentInputBundles[i].bundle;
    }
    return null;
  },
  put:function(key,bundle) {
    recentInputBundles=recentInputBundles.filter(function(item){return item.key!==key;});
    recentInputBundles.push({key:key,bundle:bundle});
    if (recentInputBundles.length>5) recentInputBundles.shift();
  }
};
// BEGIN ACCLIM-RT RETRIEVAL
// ACCLIM-RT Earth Engine retrieval library.
// Generated by rt/gee/build.py. Edit core.js and app.js, then rebuild.
// Fixed PROSPECT-D: N=1.5, Cab=40, Car=8, Cbrown=0, Cw=.015, Cm=.009.
// FourSAIL: spherical Verhoef LIDF (-.35,-.15), hotspot=.05.
// Spectral means over 620-670 and 841-876 nm, sampled at 1 nm.
// Equations follow Verhoef et al. 2007, doi:10.1109/TGRS.2007.895844.
// Numerical reference: jgomezdans/prosail 2.0.5, FourSAIL.py.
var RT_LEAF = {"rho":[0.0625896963911335,0.06202388870234529,0.06153597796090371,0.0611113551360732,0.06074170463782174,0.06043424261372378,0.06017215571707808,0.059977556992232214,0.05979697518626102,0.05964384057629166,0.059464169266544895,0.05927388999524043,0.05902263173689722,0.058704407156019156,0.0582947065985025,0.05779627119596954,0.05720186413486978,0.056486668666981016,0.055700191351999474,0.054823037104860826,0.0539039530606307,0.052939127459364646,0.05195835094703228,0.05098354709222833,0.050022093976397446,0.04908603145425414,0.04820812132472298,0.04741468702119819,0.046678413192517317,0.04605843895493955,0.045495770730433725,0.0450274333622101,0.044588090683059654,0.04418506838154737,0.04376721723402183,0.04333428966211632,0.04286219621618144,0.04234431045504196,0.04178700870839548,0.04119506662037898,0.04058487682761451,0.03996103721911245,0.03936903516208812,0.03879441316205869,0.03828219121447477,0.03781826883843458,0.03739731733013242,0.03706794115372812,0.03677434100016923,0.0365523748634502,0.036352019518955815,0.4419675291509636,0.4419599918128174,0.4419492246041353,0.4419446486890788,0.4419398036991777,0.4419290374107764,0.4419258076199554,0.44188930365644075,0.44187450362504066,0.44187584904411686,0.4418693911024509,0.44186104985581254,0.4418562066865009,0.44185405419868634,0.4418470587487745,0.4418325303994726,0.44183010909481163,0.4417968352082262,0.44178231026565906,0.4417766619181345,0.4417734343516311,0.4417629450651276,0.44172644053182236,0.44171971806135313,0.44170546705682073,0.4416912169122524,0.44168234461510714,0.44167858071096294,0.4416764299355298,0.4416659451858067,0.44164417065451594,0.4416350313260488,0.4415971909466057,0.44157945413661015,0.4415670929053522,0.44155526972441816],"tau":[0.05123347951585097,0.05047239418041808,0.04982068247863875,0.04925834502642028,0.04877417065064444,0.048379512074061874,0.04805035937328354,0.047791716617632755,0.047581325467070314,0.04738300238943259,0.04717307437890459,0.046918871963291926,0.04660144932268637,0.046155153259127515,0.04559834130719383,0.0448775350656202,0.04400774917996745,0.04294839740895938,0.04177131197242107,0.04044317189869365,0.03903572105565749,0.03753988652519439,0.03599983729109728,0.03442150753587642,0.032870729055173964,0.03133942642446945,0.029883386312195195,0.028523205694879122,0.027270535590735886,0.026178183946982847,0.025202106871682027,0.02435724254834073,0.023584076739956262,0.02284194151139271,0.02209237926123047,0.02127781765690792,0.02040641058352095,0.019405150432256615,0.018337664668817293,0.017150498636146868,0.015927884114876778,0.014615042183395759,0.013333208337050333,0.012080620551148282,0.01089888123218498,0.009791511549737876,0.008783422568998339,0.00794053641526767,0.007198751321097137,0.006593586763541768,0.006067960860071246,0.47382915359944644,0.47382145608477455,0.4738104600509623,0.4738057868846835,0.47380083892254754,0.47378984380566,0.47378654536594406,0.47382203152133745,0.4738069158576984,0.4738082899707533,0.4738016942977375,0.4737931751471601,0.4737882286782981,0.47378603027945715,0.47377888561821857,0.4737640473660736,0.4737615744106243,0.47380037491039145,0.47378553909042825,0.47377976984527265,0.4737764731941909,0.4737657593817154,0.47380127725037435,0.47379441040229187,0.47377985331496747,0.473765297084547,0.4737562342040988,0.4737523894460736,0.47375019246832656,0.47373948248132675,0.47371724014424244,0.4737079044507195,0.47374207054686834,0.4737239513661494,0.4737113236313029,0.47369924553175496],"lidf":[0.018624621096667782,0.01926721219784622,0.020582745805722245,0.02263420000333073,0.025521564202268546,0.02938720659490064,0.03441893557919243,0.04084070738727483,0.04886538547152547,0.0585532050175378,0.06949362034336437,0.08034136286927757,0.08874789320807291,0.09261736167203649,0.09196705668043492,0.0888580555091335,0.08560546372702504,0.0836734026343885]};
var RT_PAR = {"mu":[0.0001,0.01,0.025,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0],"direct":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,87718,71606,87141,71135,75922,61977,55564,45358,42811,34948,34647,28283,25003,20411,19539,15950,16033,13088,13589,11093,11789,9624,10405,8494,9306,7596,8411,6866,7667,6258,7036,5744,6496,5303,6027,4920,5617,4585,5255,4289,4935,4028,4653,3798,4404,3595,4192,3422,87764,71643,87193,71177,76075,62101,55900,45632,43262,35316,35172,28711,25614,20909,20199,16489,16724,13652,14302,11675,12518,10219,11147,9099,10057,8210,9171,7486,8433,6884,7808,6374,7273,5937,6808,5558,6401,5226,6042,4933,5726,4674,5446,4446,5200,4245,4990,4073,87809,71680,87244,71219,76228,62226,56237,45907,43713,35684,35696,29139,26225,21408,20859,17028,17416,14217,15016,12258,13249,10815,11889,9706,10810,8824,9931,8107,9200,7510,8581,7005,8051,6572,7590,6196,7187,5867,6831,5577,6517,5320,6240,5094,5996,4895,5788,4725,87855,71717,87296,71261,76382,62352,56574,46182,44165,36053,36221,29568,26837,21908,21521,17568,18109,14782,15731,12841,13979,11412,12633,10312,11563,9439,10692,8728,9968,8137,9355,7636,8829,7207,8373,6835,7973,6509,7621,6221,7310,5967,7035,5743,6793,5546,6587,5377,87901,71755,87348,71303,76535,62477,56911,46457,44618,36422,36747,29997,27450,22408,22183,18108,18802,15349,16447,13426,14711,12009,13377,10920,12317,10055,11455,9351,10737,8765,10129,8269,9609,7844,9157,7475,8761,7152,8412,6867,8103,6615,7831,6393,7592,6197,7387,6030,88177,71979,88188,71988,86605,70695,76350,62325,65085,53129,55916,45644,43115,35195,34917,28503,29295,23914,25211,20580,22115,18052,19684,16068,17721,14466,16104,13145,14744,12035,13582,11087,12580,10269,11705,9555,10933,8925,10250,8367,9645,7873,9109,7435,8634,7048,8229,6717,88257,72044,88269,72054,86716,70786,76639,62560,65568,53523,56556,46167,43976,35897,35918,29320,30393,24809,26379,21533,23335,19049,20946,17098,19017,15523,17427,14226,16091,13135,14949,12203,13964,11399,13104,10697,12346,10078,11674,9530,11079,9044,10552,8614,10086,8233,9687,7908,88338,72109,88351,72120,86827,70876,76929,62796,66052,53918,57198,46691,44838,36601,36921,30138,31493,25707,27549,22488,24559,20047,22211,18131,20316,16584,18754,15309,17440,14237,16319,13321,15351,12531,14506,11841,13761,11233,13102,10695,12516,10217,11999,9795,11541,9421,11149,9101,88418,72174,88433,72187,86939,70967,77219,63033,66537,54313,57842,47216,45702,37306,37927,30959,32595,26607,28721,23445,25785,21048,23479,19165,21617,17646,20083,16394,18793,15341,17692,14442,16741,13666,15911,12988,15180,12391,14532,11862,13957,11393,13449,10978,12999,10611,12614,10297,88499,72240,88516,72254,87051,71058,77510,63270,67023,54710,58487,47742,46568,38013,38934,31781,33699,27508,29896,24404,27013,22050,24749,20202,22921,18710,21415,17481,20149,16447,19067,15564,18134,14803,17319,14137,16601,13551,15965,13032,15401,12571,14901,12164,14460,11803,14082,11495,88600,72322,88622,72340,88428,72181,84288,68803,76762,62659,69024,56343,56257,45921,47029,38389,40270,32871,35151,28693,31158,25433,27956,22819,25329,20675,23136,18885,21273,17365,19669,16055,18275,14917,17049,13917,15963,13031,14997,12242,14137,11540,13373,10916,12694,10362,12113,9888,88705,72407,88730,72427,88543,72275,84510,68983,77172,62993,69625,56833,57173,46669,48172,39321,41579,33940,36587,29865,32691,26685,29568,24136,27006,22044,24867,20298,23050,18815,21485,17538,20125,16428,18930,15452,17871,14587,16928,13818,16089,13133,15344,12525,14682,11985,14115,11522,88811,72493,88837,72515,88660,72370,84733,69165,77582,63328,70228,57324,58092,47418,49319,40257,42893,35012,38026,31040,34230,27941,31185,25456,28688,23417,26603,21715,24832,20270,23307,19025,21981,17943,20816,16991,19784,16149,18865,15399,18047,14731,17321,14138,16676,13612,16123,13161,88917,72579,88946,72603,88776,72464,84957,69346,77994,63663,70832,57818,59013,48170,50469,41196,44210,36087,39471,32218,35773,29200,32808,26780,30375,24794,28345,23137,26620,21729,25134,20516,23843,19462,22708,18536,21702,17715,20808,16985,20011,16334,19303,15757,18675,15244,18137,14805,89024,72666,89054,72691,88894,72560,85181,69529,78408,64000,71439,58312,59937,48924,51623,42137,45532,37166,40919,33400,37320,30463,34435,28107,32067,26175,30091,24562,28412,23192,26967,22012,25710,20986,24605,20084,23627,19286,22756,18575,21981,17942,21292,17380,20680,16881,20157,16453,88988,72637,89018,72662,89022,72664,87469,71397,82967,67722,77153,62977,65813,53720,56578,46182,49355,40286,43653,35632,39076,31896,35329,28837,32204,26287,29563,24131,27296,22281,25327,20673,23602,19265,22077,18021,20719,16912,19505,15921,18419,15035,17451,14245,16589,13541,15849,12937,89112,72737,89144,72763,89152,72770,87655,71548,83298,67992,77669,63397,66687,54433,57743,47132,50747,41422,45225,36914,40791,33296,37161,30333,34135,27862,31577,25774,29381,23982,27473,22425,25803,21061,24325,19856,23010,18782,21834,17822,20782,16963,19844,16198,19009,15516,18292,14931,89235,72838,89270,72866,89283,72876,87841,71699,83631,68263,78187,63819,67565,55149,58912,48087,52144,42562,46802,38201,42512,34700,39000,31834,36072,29444,33597,27424,31473,25689,29627,24183,28011,22864,26582,21697,25309,20658,24171,19729,23153,18899,22246,18158,21438,17499,20744,16932,89360,72939,89397,72969,89415,72983,88028,71851,83965,68534,78707,64243,68446,55868,60086,49044,53547,43707,48384,39493,44240,36110,40846,33340,38017,31031,35625,29079,33572,27403,31789,25947,30227,24673,28846,23545,27616,22542,26517,21644,25533,20841,24657,20126,23876,19489,23205,18941,89485,73040,89525,73072,89547,73090,88217,72004,84300,68807,79230,64669,69330,56589,61264,50005,54955,44855,49973,40789,45974,37525,42699,34852,39969,32624,37661,30740,35680,29123,33959,27718,32452,26488,31119,25400,29932,24432,28871,23565,27922,22791,27076,22100,26323,21485,25676,20957,89345,72927,89382,72957,89422,72989,88876,72544,86340,70474,82241,67127,72782,59407,64115,52333,56880,46427,50928,41569,46011,37555,41899,34200,38417,31357,35435,28923,32849,26812,30583,24963,28584,23331,26806,21880,25214,20580,23783,19413,22499,18364,21350,17426,20323,16588,19439,15867,89481,73036,89520,73068,89564,73104,89047,72682,86613,70695,82671,67478,73572,60051,65233,53245,58271,47562,52542,42886,47810,39024,43853,35794,40501,33058,37631,30716,35143,28684,32962,26904,31038,25334,29327,23937,27794,22686,26417,21562,25181,20553,24075,19651,23087,18844,22236,18149,89617,73147,89659,73181,89707,73220,89219,72822,86887,70918,83104,67831,74366,60699,66356,54161,59667,48701,54164,44209,49617,40498,45816,37395,42595,34767,39838,32516,37446,30564,35350,28854,33502,27345,31858,26003,30385,24801,29062,23721,27874,22752,26812,21884,25862,21109,25044,20442,89754,73258,89798,73294,89851,73337,89392,72963,87162,71142,83539,68185,75164,61349,67484,55080,61070,49846,55792,45538,51432,41979,47786,39004,44698,36482,42053,34324,39759,32452,37749,30811,35977,29365,34400,28077,32987,26925,31718,25889,30579,24959,29560,24127,28649,23384,27865,22744,89892,73369,89939,73408,89996,73454,89567,73104,87439,71367,83976,68541,75965,62002,68616,56004,62479,50995,57428,46872,53255,43467,49766,40619,46809,38206,44278,36140,42083,34348,40159,32778,38462,31393,36952,30161,35600,29057,34386,28066,33295,27176,32320,26380,31448,25668,30697,25055,89673,73193,89716,73228,89769,73271,89611,73141,88242,72024,85467,69759,77882,63568,70074,57195,63117,51517,57153,46648,52083,42511,47755,38978,44028,35936,40797,33299,37966,30988,35464,28946,33241,27132,31251,25507,29459,24045,27842,22725,26384,21535,25074,20466,23901,19508,22887,18681,89816,73308,89862,73345,89919,73392,89778,73276,88476,72214,85827,70052,78576,64134,71108,58039,64454,52607,58748,47950,53898,43991,49756,40611,46191,37701,43100,35178,40391,32967,37996,31013,35869,29277,33965,27722,32250,26323,30703,25060,29308,23921,28055,22899,26932,21982,25962,21190,89960,73425,90008,73464,90070,73514,89945,73413,88711,72405,86188,70346,79274,64702,72148,58886,65797,53703,60351,49258,55721,45479,51768,42252,48364,39474,45413,37066,42827,34955,40541,33089,38510,31432,36692,29948,35055,28612,33578,27406,32245,26319,31049,25342,29977,24467,29051,23711,90104,73542,90156,73583,90222,73637,90114,73549,88948,72597,86552,70642,79975,65274,73193,59739,67147,54804,61961,50572,57553,46974,53789,43901,50547,41256,47737,38962,45274,36952,43097,35175,41163,33597,39432,32184,37873,30911,36466,29763,35197,28727,34058,27798,33037,26964,32155,26244,90250,73660,90304,73704,90375,73761,90284,73687,89186,72791,86918,70939,80680,65848,74243,60595,68503,55910,63580,51892,59394,48476,55819,45558,52741,43046,50072,40868,47733,38959,45666,37271,43829,35772,42185,34430,40704,33222,39368,32131,38163,31148,37081,30265,36111,29473,35273,28789,89974,73437,90023,73476,90086,73527,90080,73523,89373,72945,87551,71458,81633,66628,74793,61045,68291,55738,62482,50997,57401,46850,52972,43235,49097,40072,45694,37295,42681,34836,39995,32644,37592,30682,35427,28915,33467,27315,31690,25865,30081,24552,28631,23368,27327,22304,26198,21382,90121,73555,90172,73597,90239,73652,90246,73657,89582,73115,87856,71706,82233,67117,75728,61808,69543,56760,64016,52249,59182,48303,54967,44863,51279,41853,48041,39210,45174,36870,42618,34784,40330,32917,38270,31235,36404,29713,34713,28332,33182,27082,31802,25956,30561,24943,29486,24066,90268,73674,90323,73719,90394,73777,90412,73792,89793,73286,88163,71955,82837,67608,76669,62575,70803,57787,65559,53507,60972,49764,56973,46499,53474,43643,50400,41135,47679,38914,45254,36935,43083,35163,41127,33567,39357,32122,37752,30812,36298,29626,34988,28557,33811,27595,32790,26763,90417,73795,90474,73841,90550,73903,90580,73927,90005,73458,88471,72206,83444,68103,77615,63346,72069,58819,67110,54772,62772,51232,58989,48144,55679,45443,52772,43070,50198,40969,47903,39097,45850,37421,43999,35911,42325,34544,40806,33304,39431,32182,38192,31171,37077,30261,36112,29474,90566,73916,90627,73965,90707,74030,90749,74064,90219,73631,88782,72459,84054,68600,78567,64121,73342,59857,68669,56044,64581,52708,61016,49798,57896,47252,55156,45015,52730,43035,50567,41270,48631,39690,46887,38266,45308,36978,43876,35809,42580,34752,41412,33798,40361,32941,39451,32198,90251,73660,90305,73704,90376,73762,90436,73811,90096,73533,88932,72583,84406,68890,78538,64100,72586,59242,67045,54720,62059,50650,57620,47027,53673,43806,50163,40942,47023,38379,44200,36075,41656,33998,39349,32115,37250,30402,35337,28841,33598,27422,32026,26138,30607,24980,29375,23975,90398,73780,90455,73826,90530,73887,90600,73944,90289,73690,89195,72797,84922,69310,79371,64779,73739,60183,68494,55902,63773,52049,59570,48618,55833,45568,52509,42856,49536,40429,46862,38247,44452,36280,42267,34497,40279,32874,38468,31396,36821,30051,35331,28836,33988,27739,32820,26787,90546,73900,90606,73948,90686,74013,90765,74078,90483,73847,89461,73013,85441,69732,80210,65463,74899,61128,69951,57090,65497,53455,61531,50218,58004,47340,54868,44780,52062,42490,49538,40430,47264,38574,45202,36891,43325,35360,41615,33964,40061,32696,38655,31548,37387,30513,36285,29614,90696,74021,90758,74072,90843,74140,90932,74212,90679,74006,89728,73230,85963,70157,81054,66150,76065,62079,71416,58285,67231,54869,63503,51827,60188,49121,57240,46715,54602,44562,52229,42626,50091,40881,48152,39299,46388,37859,44781,36547,43319,35354,41997,34275,40805,33302,39769,32457,90847,74143,90912,74196,91001,74268,91100,74348,90877,74166,89998,73449,86489,70585,81902,66842,77238,63035,72890,59486,68974,56291,65486,53444,62384,50913,59625,48661,57156,46646,54936,44834,52934,43201,51120,41720,49468,40372,47963,39144,46595,38027,45358,37018,44241,36107,43272,35315,90504,73866,90563,73914,90642,73978,90734,74053,90596,73940,89874,73351,86472,70574,81517,66530,76157,62155,70955,57910,66139,53979,61761,50406,57805,47177,54242,44270,51022,41641,48101,39258,45450,37094,43031,35119,40818,33313,38792,31660,36943,30150,35264,28780,33745,27541,32422,26461,90650,73984,90712,74034,90795,74102,90896,74183,90778,74087,90108,73540,86915,70934,82253,67129,77204,63009,72304,59009,67765,55305,63638,51937,59910,48895,56551,46153,53515,43676,50762,41429,48262,39388,45982,37527,43895,35824,41985,34265,40241,32842,38658,31550,37226,30381,35978,29363,90797,74103,90862,74155,90949,74226,91059,74315,90961,74236,90343,73731,87361,71297,82993,67732,78258,63868,73660,60115,69401,56639,65527,53478,62027,50622,58874,48048,56023,45722,53438,43612,51091,41696,48949,39949,46990,38350,45196,36886,43559,35549,42072,34336,40727,33238,39555,32282,90945,74222,91013,74277,91105,74352,91223,74448,91146,74385,90580,73923,87810,71662,83738,68339,79319,64732,75025,61228,71046,57981,67427,55028,64157,52359,61210,49954,58546,47780,56130,45808,53936,44018,51934,42384,50103,40890,48426,39521,46895,38272,45506,37138,44248,36112,43153,35218,91095,74343,91165,74401,91262,74479,91388,74582,91332,74536,90819,74117,88262,72030,84487,68949,80386,65602,76398,62348,72702,59331,69339,56587,66299,54106,63560,51871,61084,49850,58838,48017,56798,46353,54937,44834,53234,43444,51675,42172,50252,41010,48960,39956,47790,39002,46772,38171,90737,74054,90800,74106,90886,74175,90998,74266,90973,74245,90543,73895,88023,71838,83894,68468,79127,64578,74307,60644,69714,56895,65450,53415,61535,50220,57964,47306,54703,44644,51719,42209,48991,39983,46486,37939,44183,36058,42064,34329,40122,32744,38352,31300,36746,29989,35343,28844,90880,74169,90946,74223,91036,74296,91156,74394,91145,74385,90754,74065,88404,72148,84539,68993,80070,65346,75549,61656,71239,58139,67237,54873,63563,51874,60210,49138,57148,46640,54347,44354,51786,42263,49434,40344,47271,38579,45282,36955,43458,35467,41796,34111,40288,32879,38970,31804,91023,74286,91092,74342,91187,74419,91315,74522,91319,74526,90966,74237,88789,72460,85189,69522,81019,66120,76798,62675,72774,59391,69036,56340,65603,53538,62471,50982,59609,48647,56992,46511,54598,44557,52400,42764,50378,41114,48519,39596,46814,38205,45261,36938,43851,35787,42619,34782,91169,74403,91240,74461,91339,74542,91475,74652,91495,74668,91180,74410,89176,72775,85842,70054,81975,66898,78056,63700,74318,60650,70845,57816,67656,55213,64745,52837,62085,50667,59652,48681,57427,46865,55384,45198,53504,43664,51776,42254,50191,40960,48747,39782,47436,38712,46291,37777,91315,74521,91390,74582,91493,74666,91637,74783,91672,74811,91395,74585,89566,73092,86500,70590,82936,67682,79322,64732,75872,61917,72666,59301,69721,56897,67033,54703,64577,52699,62329,50865,60274,49187,58386,47647,56650,46230,55052,44927,53588,43732,52254,42643,51043,41655,49984,40791,90950,74227,91018,74282,91110,74357,91236,74459,91275,74491,91038,74297,89199,72796,85797,70020,81603,66596,77182,62988,72847,59451,68737,56096,64902,52967,61359,50076,58090,47407,55073,44946,52295,42678,49728,40583,47355,38646,45161,36856,43142,35208,41296,33702,39615,32330,38142,31128,91089,74338,91159,74395,91255,74474,91389,74582,91440,74624,91230,74452,89529,73064,86361,70478,82446,67283,78316,63913,74265,60607,70422,57471,66837,54545,63524,51841,60466,49346,57645,47044,55046,44923,52645,42964,50425,41152,48373,39477,46484,37936,44757,36526,43184,35243,41807,34118,91228,74451,91301,74510,91402,74592,91543,74707,91605,74758,91424,74610,89863,73335,86928,70940,83295,67975,79458,64843,75692,61770,72119,58855,68785,56134,65703,53618,62858,51297,60233,49155,57815,47182,55581,45359,53515,43672,51605,42114,49848,40680,48240,39368,46776,38173,45494,37127,91369,74565,91444,74626,91549,74712,91698,74833,91773,74893,91620,74768,90199,73608,87500,71405,84149,68671,80607,65780,77128,62941,73827,60247,70745,57732,67895,55407,65265,53260,62838,51280,60601,49455,58535,47769,56624,46209,54858,44768,53232,43441,51745,42228,50391,41122,49205,40154,91511,74680,91589,74743,91698,74832,91855,74959,91942,75030,91817,74928,90537,73883,88076,71874,85010,69372,81764,66723,78574,64120,75545,61648,72717,59340,70102,57206,67687,55236,65459,53418,63405,51742,61508,50194,59753,48762,58131,47438,56638,46219,55272,45105,54028,44090,52939,43201,91146,74385,91217,74443,91315,74523,91453,74634,91531,74697,91419,74606,90100,73530,87326,71265,83669,68281,79649,65000,75593,61690,71665,58485,67941,55446,64456,52602,61207,49950,58183,47482,55377,45192,52769,43064,50344,41085,48092,39247,46011,37549,44101,35990,42356,34566,40824,33316,91279,74492,91353,74552,91455,74635,91599,74753,91687,74824,91597,74750,90389,73764,87818,71665,84418,68891,80678,65839,76902,62757,73243,59772,69774,56940,66526,54290,63498,51819,60680,49519,58064,47385,55633,45401,53373,43556,51274,41843,49333,40260,47553,38807,45926,37479,44497,36313,91413,74601,91489,74663,91596,74749,91748,74873,91845,74952,91776,74895,90680,74000,88313,72069,85174,69507,81715,66684,78219,63831,74832,61067,71619,58445,68610,55990,65805,53700,63193,51569,60769,49591,58516,47753,56421,46043,54476,44455,52677,42988,51027,41641,49519,40410,48194,39330,91549,74710,91628,74774,91738,74864,91897,74993,92005,75081,91957,75041,90974,74238,88812,72475,85934,70126,82758,67534,79546,64913,76431,62371,73476,59959,70708,57701,68126,55594,65723,53633,63492,51812,61418,50120,59489,48546,57698,47085,56043,45733,54523,44494,53135,43360,51915,42365,91685,74821,91767,74887,91882,74980,92048,75115,92166,75211,92139,75189,91270,74479,89316,72884,86700,70750,83809,68390,80881,66001,78041,63684,75345,61483,72820,59423,70463,57500,68269,55709,66233,54048,64339,52502,62578,51065,60942,49731,59430,48497,58042,47364,56774,46329,55660,45420,91325,74530,91400,74590,91504,74675,91651,74795,91753,74878,91725,74854,90801,74100,88559,72270,85395,69688,81768,66728,78001,63653,74274,60613,70684,57683,67280,54905,64074,52288,61064,49832,58250,47536,55619,45389,53159,43381,50865,41509,48735,39771,46774,38170,44976,36703,43392,35411,91452,74632,91529,74695,91637,74782,91791,74908,91902,74998,91890,74988,91054,74305,88988,72619,86060,70229,82697,67485,79202,64633,75743,61810,72409,59089,69248,56510,66269,54079,63473,51797,60859,49664,58414,47669,56128,45804,53996,44063,52017,42448,50194,40961,48522,39597,47051,38396,91580,74736,91660,74801,91772,74891,91933,75022,92052,75119,92056,75122,91311,74513,89421,72971,86729,70774,83633,68248,80411,65618,77221,63015,74146,60505,71228,58125,68479,55882,65898,53775,63485,51806,61227,49964,59117,48242,57147,46635,55320,45143,53636,43769,52093,42510,50733,41400,91710,74840,91792,74907,91907,75001,92076,75138,92204,75242,92225,75258,91569,74723,89857,73325,87404,71324,84576,69016,81629,66611,78710,64229,75894,61931,73223,59751,70704,57696,68339,55766,66128,53962,64059,52274,62125,50695,60320,49223,58645,47856,57101,46596,55686,45442,54440,44425,91841,74946,91926,75015,92045,75112,92220,75254,92357,75366,92395,75396,91830,74934,90296,73683,88084,71877,85525,69789,82856,67611,80209,65451,77655,63367,75231,61388,72945,59523,70797,57771,68789,56132,66911,54599,65154,53166,63514,51828,61992,50586,60590,49442,59304,48393,58172,47469,91489,74662,91567,74726,91676,74814,91832,74941,91952,75038,91978,75060,91352,74548,89558,73084,86841,70866,83590,68213,80112,65375,76599,62508,73159,59701,69855,57005,66710,54439,63733,52009,60929,49721,58291,47568,55811,45544,53487,43648,51321,41880,49318,40246,47477,38743,45851,37416,91610,74760,91690,74825,91803,74916,91966,75049,92093,75152,92132,75184,91577,74730,89933,73389,87429,71345,84425,68894,81210,66270,77959,63617,74774,61018,71715,58522,68803,56145,66045,53895,63448,51775,61003,49781,58705,47906,56552,46148,54544,44510,52689,42996,50982,41603,49475,40373,91732,74858,91815,74926,91931,75020,92100,75158,92235,75267,92288,75310,91804,74914,90311,73696,88021,71827,85267,69579,82315,67170,79328,64733,76401,62344,73588,60049,70909,57863,68372,55793,65983,53843,63734,52008,61619,50282,59637,48665,57789,47157,56081,45764,54511,44482,53123,43350,91855,74958,91940,75027,92060,75124,92236,75267,92379,75383,92445,75437,92033,75100,90693,74006,88618,72312,86115,70270,83428,68077,80707,65857,78038,63679,75474,61587,73030,59593,70715,57704,68535,55925,66483,54250,64553,52675,62743,51199,61057,49823,59497,48550,58063,47380,56796,46346,91980,75058,92067,75129,92191,75230,92374,75378,92524,75501,92604,75565,92264,75287,91077,74318,89219,72802,86969,70965,84549,68990,82095,66989,79688,65024,77373,63135,75166,61334,73075,59628,71105,58021,69250,56507,67506,55084,65871,53750,64346,52506,62937,51355,61640,50297,60495,49363,91640,74784,91721,74850,91834,74942,91999,75075,92131,75183,92195,75235,91792,74905,90371,73746,88054,71855,85156,69490,81965,66886,78671,64198,75392,61522,72203,58920,69135,56416,66205,54025,63426,51758,60795,49610,58308,47581,55966,45670,53774,43881,51741,42222,49865,40691,48203,39335,91754,74876,91837,74944,91954,75039,92125,75177,92264,75290,92339,75351,91993,75068,90700,74013,88573,72277,85905,70100,82963,67699,79925,65220,76898,62750,73953,60347,71120,58035,68413,55826,65846,53731,63415,51747,61116,49872,58952,48106,56926,46453,55047,44919,53313,43504,51777,42251,91870,74969,91955,75039,92075,75136,92252,75280,92398,75399,92484,75468,92195,75232,91032,74282,89096,72702,86660,70715,83969,68519,81187,66249,78415,63986,75716,61784,73119,59665,70637,57640,68283,55719,66052,53899,63944,52179,61959,50558,60100,49042,58376,47635,56784,46336,55375,45186,91986,75063,92074,75135,92197,75235,92381,75384,92534,75509,92631,75587,92400,75398,91366,74553,89623,73131,87421,71334,84982,69344,82459,67285,79942,65231,77491,63231,75131,61306,72876,59466,70736,57719,68708,56065,66792,54501,64986,53028,63296,51648,61727,50369,60280,49188,58998,48142,92104,75158,92194,75231,92321,75335,92511,75489,92672,75619,92779,75707,92607,75565,91703,74827,90154,73563,88187,71958,86003,70176,83740,68329,81480,66485,79279,64689,77158,62958,75131,61304,73207,59734,71383,58246,69659,56839,68034,55514,66513,54273,65102,53121,63800,52058,62646,51117,91778,74895,91861,74963,91979,75059,92151,75198,92293,75314,92383,75387,92149,75195,91037,74288,89073,72685,86505,70589,83591,68211,80518,65703,77408,63165,74343,60665,71365,58234,68495,55892,65753,53655,63141,51523,60658,49498,58309,47581,56102,45780,54046,44102,52144,42550,50455,41171,91886,74982,91971,75052,92092,75150,92270,75294,92418,75415,92518,75496,92329,75341,91326,74522,89531,73057,87175,71135,84497,68949,81669,66642,78806,64306,75984,62003,73240,59764,70596,57606,68069,55544,65661,53579,63373,51712,61207,49945,59172,48284,57277,46738,55523,45306,53965,44035,91994,75070,92082,75141,92206,75242,92390,75391,92545,75517,92654,75605,92511,75489,91617,74759,89993,73433,87850,71684,85409,69692,82829,67587,80215,65454,77637,63350,75129,61304,72712,59332,70402,57446,68199,55649,66106,53942,64125,52325,62264,50806,60530,49391,58925,48082,57499,46918,92104,75159,92194,75232,92321,75335,92511,75489,92673,75621,92791,75716,92695,75637,91911,74997,90458,73811,88531,72238,86329,70441,83998,68539,81634,66610,79301,64707,77031,62855,74843,61069,72750,59362,70755,57734,68859,56187,67064,54722,65377,53345,63805,52063,62350,50876,61058,49822,92215,75248,92307,75323,92438,75429,92634,75588,92803,75725,92930,75828,92881,75788,92207,75237,90927,74192,89216,72796,87255,71195,85175,69498,83063,67775,80978,66074,78948,64417,76989,62819,75116,61291,73330,59833,71631,58447,70024,57136,68512,55902,67104,54753,65800,53690,64642,52745,91904,74997,91990,75067,92111,75166,92290,75311,92441,75433,92550,75522,92442,75433,91586,74734,89932,73385,87668,71537,85019,69375,82163,67045,79226,64648,76294,62256,73414,59905,70615,57621,67921,55423,65338,53315,62870,51301,60524,49387,58309,47580,56240,45891,54319,44324,52608,42927,92005,75079,92093,75150,92217,75251,92402,75401,92559,75528,92675,75623,92605,75565,91840,74940,90336,73713,88266,72023,85838,70042,83218,67904,80521,65704,77828,63506,75180,61346,72607,59246,70130,57225,67755,55286,65485,53434,63327,51673,61290,50011,59386,48458,57618,47015,56044,45731,92108,75161,92197,75234,92325,75338,92515,75492,92678,75624,92802,75726,92770,75698,92097,75148,90743,74044,88868,72514,86664,70715,84281,68770,81826,66767,79372,64765,76960,62796,74614,60882,72355,59039,70189,57272,68118,55582,66149,53976,64291,52459,62554,51042,60941,49726,59504,48553,92211,75245,92303,75319,92433,75425,92629,75584,92798,75721,92931,75829,92936,75833,92355,75358,91154,74378,89476,73008,87496,71392,85352,69642,83140,67838,80928,66033,78752,64257,76635,62530,74596,60867,72640,59271,70771,57745,68992,56294,67313,54924,65744,53644,64287,52455,62988,51395,92315,75329,92409,75405,92543,75513,92744,75677,92919,75819,93060,75934,93104,75968,92617,75570,91568,74714,90088,73506,88334,72075,86430,70521,84464,68916,82495,67310,80557,65729,78671,64190,76854,62707,75109,61284,73442,59923,71855,58629,70358,57407,68957,56264,67656,55203,66498,54258,92020,75091,92108,75162,92233,75264,92418,75414,92576,75542,92699,75642,92687,75632,92040,75103,90658,73975,88671,72354,86273,70397,83631,68241,80867,65986,78072,63705,75297,61441,72577,59221,69940,57069,67395,54993,64950,52998,62615,51093,60402,49287,58327,47593,56393,46016,54667,44607,92115,75167,92205,75240,92332,75344,92523,75498,92686,75631,92816,75737,92835,75752,92265,75285,91014,74265,89204,72787,87013,70999,84594,69026,82063,66960,79501,64870,76956,62794,74461,60758,72042,58783,69706,56877,67462,55047,65319,53298,63287,51640,61381,50085,59606,48636,58020,47343,92210,75244,92302,75319,92433,75425,92628,75584,92798,75721,92935,75832,92985,75872,92491,75469,91374,74557,89741,73224,87758,71606,85565,69817,83267,67942,80940,66043,78628,64156,76359,62305,74158,60509,72033,58776,69992,57110,68041,55518,66192,54009,64457,52594,62841,51275,61398,50098,92307,75322,92401,75398,92534,75507,92735,75670,92910,75812,93054,75929,93136,75994,92720,75655,91736,74851,90283,73665,88509,72218,86544,70614,84481,68931,82390,67225,80311,65528,78270,63863,76290,62247,74378,60687,72540,59187,70783,57754,69118,56396,67555,55121,66100,53933,64799,52872,92405,75401,92500,75479,92637,75589,92843,75757,93024,75904,93176,76027,93288,76117,92952,75842,92102,75148,90829,74109,89266,72834,87530,71417,85704,69927,83851,68415,82007,66910,80196,65433,78438,63998,76739,62612,75106,61280,73545,60006,72065,58799,70676,57665,69381,56609,68225,55666,92126,75176,92216,75249,92344,75353,92535,75508,92699,75642,92833,75751,92895,75800,92418,75411,91272,74475,89537,73059,87375,71295,84939,69307,82348,67193,79693,65026,77028,62852,74394,60702,71820,58602,69321,56563,66908,54594,64591,52704,62387,50905,60311,49212,58372,47629,56635,46212,92214,75247,92306,75322,92437,75428,92633,75587,92802,75725,92943,75839,93029,75909,92618,75572,91587,74731,90012,73445,88042,71838,85818,70023,83449,68090,81021,66109,78583,64119,76171,62151,73814,60228,71525,58361,69314,56557,67191,54825,65171,53176,63270,51625,61492,50174,59900,48875,92304,75319,92397,75395,92531,75504,92731,75667,92907,75809,93053,75928,93165,76019,92819,75735,91904,74988,90490,73834,88713,72384,86703,70744,84559,68994,82359,67199,80148,65395,77960,63610,75822,61866,73745,60171,71738,58533,69810,56961,67976,55464,66249,54055,64634,52737,63188,51557,92394,75392,92489,75470,92626,75580,92831,75747,93012,75894,93165,76018,93303,76130,93022,75900,92224,75248,90972,74226,89390,72935,87595,71470,85677,69905,83706,68297,81725,66681,79763,65080,77845,63515,75981,61994,74179,60524,72448,59112,70801,57768,69250,56502,67799,55319,66500,54259,92486,75466,92583,75545,92722,75657,92932,75828,93118,75979,93278,76109,93442,76242,93227,76066,92547,75510,91459,74622,90073,73491,88494,72202,86804,70823,85064,69404,83314,67976,81580,66561,79883,65176,78233,63830,76638,62529,75106,61279,73647,60088,72272,58967,70987,57918,69836,56979,92223,75254,92315,75329,92446,75435,92642,75595,92812,75733,92955,75849,93072,75944,92736,75668,91794,74899,90286,73669,88344,72085,86106,70258,83684,68282,81169,66230,78619,64149,76075,62073,73571,60030,71124,58033,68748,56095,66457,54225,64267,52438,62198,50750,60259,49168,58517,47746,92305,75321,92399,75397,92532,75505,92733,75668,92908,75811,93057,75931,93195,76043,92913,75812,92072,75125,90708,74012,88944,72573,86905,70909,84696,69107,82401,67234,80072,65333,77747,63436,75458,61569,73221,59744,71048,57971,68953,56261,66950,54626,65057,53082,63282,51634,61688,50333,92389,75388,92484,75465,92620,75576,92825,75742,93005,75889,93160,76014,93320,76143,93092,75957,92352,75353,91134,74358,89548,73064,87711,71565,85717,69938,83642,68245,81535,66526,79431,64809,77359,63119,75333,61466,73365,59860,71466,58311,69651,56830,67936,55430,66327,54117,64882,52939,92473,75456,92570,75534,92708,75647,92918,75817,93104,75968,93263,76098,93445,76245,93273,76104,92635,75582,91563,74707,90158,73560,88523,72226,86745,70775,84892,69264,83010,67728,81128,66193,79274,64680,77461,63201,75699,61763,73998,60375,72373,59049,70836,57795,69394,56619,68099,55562,92558,75524,92657,75604,92798,75719,93012,75893,93203,76048,93369,76182,93573,76347,93456,76251,92921,75814,91996,75059,90772,74060,89341,72892,87780,71619,86152,70290,84495,68938,82838,67586,81204,66253,79605,64948,78050,63680,76549,62455,75114,61284,73757,60177,72484,59138,71340,58205,92311,75326,92405,75402,92539,75510,92740,75674,92916,75816,93066,75938,93226,76068,93004,75886,92237,75260,90934,74197,89197,72779,87147,71106,84890,69264,82515,67327,80080,65340,77630,63341,75201,61359,72812,59409,70479,57506,68218,55661,66049,53891,63992,52213,62057,50634,60315,49212,92388,75388,92483,75465,92619,75575,92825,75742,93005,75888,93160,76014,93339,76159,93162,76014,92484,75460,91310,74502,89736,73217,87873,71697,85819,70021,83655,68256,81436,66445,79200,64621,76984,62812,74803,61033,72673,59295,70609,57611,68628,55995,66749,54461,64981,53019,63389,51720,92466,75450,92563,75529,92701,75641,92910,75811,93095,75961,93256,76091,93453,76251,93322,76144,92732,75661,91688,74810,90279,73659,88605,72293,86755,70784,84804,69192,82801,67558,80783,65911,78780,64277,76809,62669,74883,61098,73017,59575,71225,58113,69525,56726,67926,55421,66485,54246,92544,75513,92643,75593,92783,75707,92997,75881,93187,76035,93352,76169,93568,76344,93484,76274,92982,75865,92070,75120,90826,74105,89342,72894,87698,71552,85962,70135,84177,68679,82376,67210,80589,65752,78830,64316,77110,62913,75442,61552,73841,60246,72322,59007,70892,57840,69605,56790,92624,75577,92723,75658,92866,75774,93085,75951,93279,76109,93450,76248,93684,76438,93647,76406,93235,76070,92456,75433,91379,74554,90086,73499,88649,72326,87128,71085,85563,69808,83982,68519,82412,67238,80866,65976,79353,64742,77886,63545,76477,62395,75139,61304,73881,60277,72747,59352,92393,75391,92488,75469,92624,75579,92830,75746,93010,75893,93166,76019,93361,76177,93231,76071,92616,75568,91496,74654,89947,73390,88075,71862,85978,70151,83741,68326,81423,66434,79069,64514,76718,62595,74391,60696,72106,58832,69881,57017,67737,55268,65697,53603,63771,52032,62033,50613,92464,75449,92561,75527,92699,75639,92908,75809,93093,75959,93254,76090,93464,76260,93373,76185,92834,75745,91830,74925,90431,73783,88734,72399,86829,70844,84795,69185,82685,67463,80541,65714,78399,63966,76278,62236,74195,60536,72167,58881,70212,57286,68351,55767,66594,54334,65008,53040,92536,75507,92634,75586,92775,75700,92988,75873,93177,76027,93343,76162,93569,76345,93517,76301,93054,75924,92167,75199,90918,74180,89398,72939,87687,71543,85857,70050,83957,68499,82025,66923,80093,65347,78180,63786,76300,62253,74469,60759,72704,59319,71024,57947,69437,56653,68004,55484,92609,75565,92708,75646,92851,75762,93068,75938,93262,76095,93433,76234,93674,76430,93661,76418,93277,76104,92506,75475,91410,74580,90067,73484,88551,72247,86927,70921,85238,69544,83519,68141,81800,66738,80096,65348,78421,63981,76789,62650,75215,61366,73716,60143,72301,58989,71023,57946,92683,75624,92784,75706,92928,75824,93150,76004,93348,76164,93523,76307,93781,76516,93808,76536,93501,76285,92849,75753,91906,74983,90742,74033,89422,72956,88005,71800,86529,70595,85025,69368,83519,68140,82026,66921,80557,65723,79126,64555,77745,63428,76429,62355,75186,61341,74064,60425,92467,75451,92564,75530,92702,75642,92912,75812,93097,75962,93258,76094,93480,76273,93426,76228,92941,75832,91983,75050,90608,73928,88903,72537,86960,70951,84858,69236,82656,67440,80400,65599,78129,63746,75868,61901,73635,60079,71451,58297,69337,56572,67317,54924,65405,53364,63673,51951,92533,75504,92631,75584,92771,75698,92985,75871,93174,76024,93340,76159,93574,76349,93553,76331,93134,75989,92280,75292,91042,74281,89500,73023,87739,71585,85831,70029,83830,68396,81778,66722,79713,65037,77654,63357,75622,61699,73633,60076,71707,58505,69868,57004,68126,55583,66548,54296,92600,75558,92699,75639,92842,75754,93058,75930,93251,76087,93422,76225,93670,76427,93682,76435,93330,76147,92580,75535,91479,74636,90102,73513,88523,72224,86811,70827,85012,69360,83167,67854,81308,66337,79454,64825,77623,63331,75831,61868,74095,60453,72437,59099,70866,57818,69444,56657,92668,75612,92768,75694,92913,75811,93133,75990,93330,76150,93505,76292,93767,76505,93812,76541,93527,76307,92882,75780,91920,74995,90709,74006,89314,72868,87798,71631,86204,70331,84566,68994,82915,67647,81268,66303,79640,64975,78045,63674,76501,62414,75025,61210,73627,60069,72360,59036,92736,75667,92838,75750,92984,75869,93209,76051,93410,76214,93589,76360,93865,76584,93944,76647,93726,76468,93187,76028,92365,75357,91320,74504,90112,73517,88793,72442,87405,71309,85976,70143,84534,68967,83095,67792,81671,66631,80276,65493,78925,64390,77633,63336,76408,62337,75299,61432,92535,75506,92633,75586,92774,75699,92987,75873,93176,76026,93342,76161,93585,76358,93593,76364,93219,76058,92407,75395,91190,74402,89643,73139,87846,71673,85876,70066,83789,68363,81631,66602,79443,64817,77250,63027,75072,61250,72931,59504,70851,57806,68856,56179,66961,54633,65241,53229,92596,75555,92696,75636,92838,75751,93054,75927,93247,76084,93418,76222,93672,76428,93708,76457,93391,76198,92671,75609,91579,74718,90183,73579,88558,72252,86772,70796,84879,69251,82920,67652,80932,66030,78938,64404,76959,62789,75012,61200,73119,59656,71304,58175,69580,56769,68014,55491,92658,75605,92759,75686,92903,75803,93123,75982,93319,76141,93494,76283,93760,76499,93824,76551,93565,76338,92937,75826,91971,75037,90728,74022,89275,72836,87675,71531,85977,70146,84218,68710,82432,67253,80640,65791,78859,64338,77107,62909,75404,61520,73771,60187,72218,58920,70808,57770,92721,75655,92822,75738,92969,75856,93192,76038,93392,76200,93571,76345,93849,76571,93942,76645,93740,76480,93206,76044,92367,75358,91277,74469,89998,73425,88586,72273,87084,71047,85526,69776,83943,68485,82354,67188,80774,65899,79219,64630,77706,63396,76255,62212,74876,61087,73623,60065,92784,75706,92887,75790,93035,75910,93262,76094,93466,76259,93649,76408,93940,76644,94060,76741,93917,76623,93478,76264,92766,75683,91831,74919,90726,74018,89503,73020,88199,71956,86844,70850,85466,69726,84081,68596,82703,67471,81346,66364,80026,65287,78758,64253,77553,63270,76458,62377,92597,75556,92697,75637,92839,75752,93056,75928,93248,76085,93419,76223,93679,76434,93739,76482,93459,76253,92775,75694,91703,74819,90303,73677,88646,72325,86803,70821,84829,69210,82770,67530,80665,65813,78542,64081,76422,62351,74329,60643,72285,58976,70318,57371,68444,55842,66737,54449,92654,75602,92754,75683,92898,75800,93118,75978,93314,76138,93489,76279,93759,76499,93842,76565,93612,76377,93010,75885,92052,75103,90792,74074,89295,72853,87628,71493,85840,70034,83973,68510,82063,66953,80136,65380,78211,63810,76309,62258,74453,60743,72665,59285,70962,57895,69410,56629,92711,75648,92813,75730,92958,75848,93181,76029,93381,76191,93559,76336,93840,76564,93947,76650,93767,76502,93247,76078,92403,75388,91284,74475,89950,73386,88460,72170,86859,70864,85185,69498,83472,68101,81742,66689,80014,65279,78305,63885,76636,62524,75030,61213,73498,59963,72103,58825,92769,75694,92872,75777,93019,75897,93245,76081,93448,76245,93630,76393,93922,76630,94053,76735,93923,76629,93487,76272,92758,75676,91781,74879,90610,73924,89298,72853,87885,71700,86406,70494,84892,69258,83360,68008,81829,66759,80315,65524,78836,64317,77412,63155,76053,62047,74815,61037,92828,75741,92931,75825,93081,75947,93310,76132,93517,76300,93702,76451,94005,76696,94160,76821,94081,76756,93729,76468,93116,75967,92282,75287,91275,74465,90143,73541,88920,72543,87637,71496,86322,70423,84991,69337,83659,68251,82341,67175,81053,66124,79811,65111,78627,64146,77548,63265,92654,75602,92755,75683,92899,75800,93118,75979,93315,76138,93489,76279,93764,76502,93865,76584,93666,76421,93096,75955,92156,75188,90892,74156,89368,72913,87648,71509,85785,69989,83822,68388,81802,66740,79751,65066,77691,63385,75647,61717,73643,60082,71707,58503,69856,56993,68166,55614,92707,75644,92808,75726,92954,75844,93176,76025,93376,76187,93553,76331,93838,76562,93959,76659,93803,76532,93305,76125,92468,75441,91334,74516,89960,73395,88407,72127,86721,70752,84944,69302,83114,67809,81254,66291,79386,64767,77531,63253,75712,61770,73954,60336,72274,58965,70739,57713,92760,75686,92862,75770,93009,75889,93235,76072,93437,76236,93619,76384,93912,76622,94053,76736,93941,76643,93517,76296,92783,75697,91779,74878,90557,73881,89172,72750,87665,71521,86075,70224,84435,68886,82768,67526,81093,66159,79428,64801,77796,63469,76219,62182,74710,60951,73331,59827,92813,75729,92917,75813,93065,75934,93294,76120,93499,76286,93684,76436,93987,76682,94149,76813,94081,76756,93730,76469,93101,75955,92228,75243,91159,74371,89942,73378,88616,72296,87214,71152,85767,69971,84294,68769,82813,67561,81341,66360,79896,65181,78500,64042,77163,62952,75943,61956,92867,75773,92972,75858,93122,75980,93354,76168,93562,76337,93751,76490,94064,76743,94246,76891,94221,76870,93945,76644,93421,76216,92681,75611,91766,74865,90719,74010,89575,73076,88363,72087,87109,71064,85832,70022,84546,68973,83267,67930,82012,66906,80798,65915,79636,64967,78573,64101,92706,75644,92808,75726,92953,75844,93176,76025,93375,76187,93553,76331,93840,76564,93975,76673,93846,76567,93376,76183,92555,75513,91419,74585,90020,73443,88418,72136,86662,70704,84796,69181,82860,67602,80881,65987,78883,64357,76890,62731,74928,61130,73026,59578,71202,58090,69531,56727,92755,75683,92857,75766,93004,75885,93229,76068,93432,76232,93613,76379,93908,76618,94060,76742,93968,76666,93562,76334,92835,75740,91818,74909,90559,73883,89115,72704,87529,71410,85841,70033,84089,68603,82297,67141,80486,65664,78680,64190,76901,62739,75176,61331,73521,59981,72005,58745,92804,75722,92907,75806,93055,75926,93283,76111,93488,76278,93673,76427,93977,76673,94146,76811,94092,76765,93751,76487,93117,75969,92220,75237,91103,74326,89817,73276,88403,72122,86895,70891,85328,69613,83723,68304,82102,66981,80483,65660,78889,64360,77342,63098,75858,61887,74498,60777,92853,75761,92957,75846,93107,75968,93338,76155,93546,76324,93734,76476,94046,76729,94233,76881,94216,76866,93941,76641,93402,76200,92626,75566,91652,74772,90526,73853,89283,72839,87956,71756,86576,70630,85161,69475,83730,68308,82300,67141,80891,65992,79524,64877,78212,63806,77009,62825,92903,75802,93008,75887,93160,76010,93393,76199,93604,76370,93795,76525,94116,76785,94321,76952,94342,76968,94133,76796,93689,76433,93035,75899,92205,75222,91240,74434,90171,73562,89026,72628,87833,71654,86609,70656,85370,69645,84131,68634,82909,67637,81723,66669,80583,65740,79539,64888,92754,75682,92856,75765,93003,75885,93229,76067,93431,76231,93612,76378,93909,76619,94073,76752,94002,76693,93621,76381,92908,75800,91889,74968,90608,73922,89119,72707,87468,71360,85696,69914,83843,68403,81937,66848,80002,65269,78062,63686,76145,62122,74278,60599,72482,59134,70833,57788,92799,75718,92902,75802,93050,75922,93278,76107,93483,76273,93667,76422,93971,76670,94150,76814,94112,76782,93787,76516,93159,76003,92249,75261,91099,74322,89759,73228,88269,72013,86669,70707,84994,69341,83270,67934,81518,66505,79761,65072,78024,63655,76333,62275,74705,60947,73210,59727,92844,75754,92948,75839,93098,75960,93328,76147,93535,76315,93722,76467,94034,76720,94228,76877,94222,76871,93955,76652,93411,76208,92612,75556,91595,74726,90404,73754,89077,72671,87650,71506,86154,70286,84612,69028,83046,67750,81473,66467,79918,65198,78403,63963,76945,62773,75605,61680,92889,75791,92994,75876,93145,75998,93378,76187,93588,76358,93778,76512,94098,76771,94307,76940,94334,76961,94124,76790,93666,76415,92979,75854,92094,75132,91054,74283,89892,73334,88638,72312,87323,71238,85965,70131,84585,69004,83198,67873,81826,66754,80489,65663,79202,64613,78018,63647,92935,75828,93041,75913,93194,76037,93429,76228,93642,76401,93835,76557,94163,76823,94387,77004,94446,77052,94296,76928,93924,76624,93348,76154,92598,75542,91710,74816,90713,74003,89634,73123,88500,72198,87329,71242,86136,70268,84936,69290,83749,68321,82591,67376,81475,66466,80449,65629,92798,75717,92901,75801,93049,75921,93277,76106,93482,76272,93666,76422,93972,76670,94159,76821,94138,76804,93835,76555,93220,76053,92309,75310,91139,74355,89757,73227,88207,71963,86528,70592,84758,69148,82925,67653,81053,66126,79168,64587,77296,63060,75467,61568,73702,60128,72076,58801,92839,75750,92943,75835,93092,75956,93322,76142,93529,76311,93716,76462,94029,76716,94229,76878,94236,76883,93983,76675,93444,76235,92634,75574,91586,74719,90345,73705,88948,72566,87433,71330,85835,70026,84178,68674,82485,67293,80780,65901,79085,64519,77430,63169,75831,61864,74358,60663,92880,75783,92985,75868,93136,75991,93368,76179,93578,76349,93767,76503,94087,76762,94300,76935,94335,76963,94133,76797,93671,76419,92962,75840,92037,75086,90936,74187,89695,73174,88345,72073,86920,70910,85441,69703,83928,68469,82403,67225,80888,65989,79407,64781,77976,63613,76657,62537,92922,75817,93028,75903,93180,76026,93414,76216,93627,76389,93819,76545,94145,76809,94371,76992,94436,77044,94284,76919,93899,76605,93292,76109,92492,75455,91533,74673,90447,73787,89265,72822,88013,71801,86713,70740,85382,69654,84040,68559,82705,67470,81399,66404,80137,65375,78974,64426,92965,75851,93071,75938,93225,76062,93461,76254,93676,76428,93871,76586,94205,76857,94444,77051,94537,77125,94437,77043,94130,76792,93626,76380,92950,75828,92134,75162,91206,74405,90191,73577,89115,72698,87995,71785,86848,70848,85688,69903,84535,68962,83406,68041,82315,67150,81308,66329,92838,75749,92942,75834,93091,75955,93321,76141,93528,76310,93715,76461,94028,76715,94235,76883,94258,76900,94022,76707,93496,76278,92685,75615,91618,74745,90339,73701,88886,72516,87297,71219,85609,69841,83848,68405,82040,66930,80210,65437,78385,63948,76596,62488,74863,61074,73262,59768,92875,75780,92980,75865,93131,75987,93363,76175,93572,76345,93761,76498,94081,76758,94299,76934,94346,76971,94155,76815,93697,76441,92978,75854,92025,75076,90877,74139,89570,73073,88138,71904,86615,70662,85026,69365,83392,68032,81738,66682,80088,65336,78469,64015,76900,62736,75451,61554,92914,75810,93019,75896,93171,76019,93405,76209,93617,76381,93808,76536,94134,76800,94363,76986,94435,77043,94289,76923,93900,76605,93274,76094,92435,75410,91420,74581,90260,73634,88986,72595,87630,71488,86212,70332,84754,69142,83277,67937,81803,66734,80356,65554,78953,64410,77657,63352,92952,75841,93059,75927,93212,76052,93448,76243,93662,76417,93856,76574,94188,76844,94429,77039,94524,77115,94424,77032,94104,76771,93572,76336,92849,75746,91967,75026,90955,74200,89841,73291,88652,72321,87408,71306,86127,70261,84828,69201,83531,68143,82257,67103,81021,66096,79880,65164,92992,75873,93098,75959,93253,76085,93491,76278,93707,76453,93904,76613,94243,76887,94495,77092,94615,77189,94560,77143,94311,76939,93873,76581,93266,76085,92518,75475,91656,74771,90702,73993,89682,73160,88613,72288,87510,71388,86391,70475,85272,69562,84172,68665,83105,67795,82119,66990,92874,75779,92979,75864,93130,75986,93361,76174,93571,76344,93760,76497,94080,76757,94303,76937,94362,76985,94187,76841,93740,76476,93021,75889,92051,75098,90869,74133,89510,73023,88008,71798,86399,70486,84711,69108,82967,67685,81193,66238,79416,64788,77667,63362,75968,61975,74393,60691,92909,75807,93014,75892,93166,76015,93400,76205,93611,76376,93802,76531,94128,76796,94361,76984,94441,77049,94305,76937,93920,76622,93285,76103,92421,75398,91362,74534,90141,73538,88789,72435,87340,71252,85817,70009,84242,68725,82640,67418,81034,66108,79453,64818,77917,63564,76493,62402,92944,75835,93050,75921,93203,76045,93439,76236,93652,76409,93846,76566,94177,76835,94420,77032,94521,77113,94425,77034,94102,76769,93552,76320,92794,75702,91860,74939,90777,74056,89577,73076,88287,72024,86931,70917,85527,69772,84098,68606,82665,67437,81253,66285,79880,65165,78607,64127,92980,75863,93086,75950,93240,76075,93478,76267,93694,76442,93890,76601,94227,76875,94479,77079,94602,77178,94546,77131,94285,76918,93821,76539,93170,76007,92361,75347,91419,74578,90371,73723,89243,72802,88054,71832,86822,70827,85567,69803,84307,68776,83066,67763,81858,66778,80738,65864,93016,75892,93123,75979,93278,76105,93518,76299,93735,76476,93934,76637,94277,76915,94540,77128,94684,77244,94668,77230,94470,77068,94092,76759,93549,76316,92866,75758,92066,75104,91171,74374,90205,73586,89185,72754,88127,71891,87047,71010,85963,70125,84893,69252,83851,68402,82885,67614,92908,75805,93013,75891,93165,76014,93398,76203,93609,76375,93801,76530,94127,76795,94364,76986,94454,77059,94332,76958,93956,76652,93321,76133,92442,75416,91352,74526,90082,73489,88665,72333,87135,71085,85517,69765,83837,68394,82119,66993,80392,65583,78684,64190,77020,62832,75473,61571,92940,75831,93046,75917,93198,76041,93434,76232,93647,76405,93840,76562,94171,76830,94417,77029,94525,77117,94438,77044,94118,76782,93559,76327,92779,75689,91804,74893,90664,73964,89390,72924,88013,71800,86555,70611,85039,69374,83489,68109,81929,66837,80387,65578,78882,64351,77485,63211,92972,75857,93079,75943,93232,76068,93469,76261,93684,76435,93880,76594,94216,76866,94470,77072,94597,77175,94545,77131,94280,76914,93800,76522,93117,75965,92259,75264,91251,74442,90121,73520,88897,72521,87601,71463,86251,70362,84869,69234,83478,68099,82102,66977,80759,65881,79511,64863,93005,75883,93112,75970,93267,76096,93505,76290,93723,76466,93920,76626,94262,76903,94525,77116,94670,77233,94653,77218,94445,77048,94042,76719,93459,76242,92718,75637,91843,74923,90858,74120,89789,73247,88654,72321,87472,71357,86260,70368,85038,69371,83830,68385,82650,67423,81553,66528,93038,75910,93146,75997,93301,76124,93542,76319,93761,76497,93961,76659,94308,76940,94580,77160,94744,77293,94762,77306,94611,77182,94287,76917,93804,76522,93181,76014,92439,75409,91601,74725,90688,73979,89716,73187,88702,72359,87661,71510,86611,70653,85571,69805,84555,68976,83610,68204,92938,75830,93044,75916,93197,76040,93432,76231,93645,76404,93838,76560,94170,76829,94418,77030,94535,77125,94459,77062,94148,76807,93590,76351,92796,75703,91792,74884,90607,73917,89272,72827,87818,71641,86271,70379,84654,69060,82992,67704,81314,66335,79649,64977,78021,63648,76503,62410,92968,75854,93074,75940,93228,76065,93464,76257,93679,76431,93875,76589,94211,76862,94466,77069,94599,77176,94554,77138,94292,76924,93805,76526,93101,75951,92205,75220,91144,74354,89945,73376,88637,72309,87244,71172,85787,69983,84288,68761,82774,67525,81271,66299,79800,65099,78429,63981,92997,75877,93104,75964,93259,76090,93497,76283,93714,76459,93911,76619,94252,76895,94515,77109,94664,77229,94650,77216,94438,77043,94021,76702,93409,76202,92622,75559,91685,74795,90623,73928,89463,72982,88225,71971,86928,70914,85594,69825,84245,68724,82905,67631,81593,66561,80370,65563,93028,75902,93135,75989,93291,76115,93531,76310,93749,76487,93948,76648,94294,76929,94565,77148,94730,77282,94747,77294,94586,77162,94240,76880,93719,76454,93042,75901,92231,75239,91307,74485,90295,73660,89213,72777,88079,71851,86910,70897,85727,69932,84551,68973,83399,68033,82325,67157,93058,75926,93166,76014,93323,76141,93564,76337,93785,76515,93986,76678,94336,76962,94615,77189,94797,77335,94844,77373,94734,77282,94461,77059,94032,76708,93466,76246,92781,75686,91996,75046,91134,74343,90209,73588,89238,72796,88236,71978,87220,71149,86210,70325,85219,69517,84295,68763,92966,75852,93073,75939,93226,76063,93463,76255,93678,76430,93873,76588,94209,76861,94467,77070,94607,77183,94571,77153,94318,76945,93830,76547,93115,75963,92192,75210,91089,74310,89833,73284,88453,72159,86975,70952,85421,69685,83815,68375,82187,67046,80565,65723,78973,64425,77486,63211,92993,75874,93100,75961,93254,76086,93493,76279,93709,76455,93906,76615,94246,76891,94511,77105,94665,77229,94656,77221,94447,77050,94024,76704,93392,76188,92571,75518,91584,74713,90457,73793,89217,72781,87887,71696,86488,70554,85041,69374,83572,68176,82109,66982,80672,65810,79330,64715,93021,75896,93128,75983,93283,76109,93523,76304,93741,76480,93940,76641,94284,76921,94556,77141,94724,77277,94742,77291,94578,77156,94219,76863,93671,76415,92952,75828,92083,75118,91086,74305,89988,73409,88807,72445,87563,71430,86276,70380,84968,69314,83665,68250,82384,67206,81187,66229,93048,75918,93156,76006,93312,76132,93553,76328,93773,76506,93974,76669,94323,76952,94601,77177,94783,77324,94829,77361,94710,77263,94417,77023,93953,76644,93336,76141,92585,75528,91720,74822,90764,74041,89733,73200,88646,72313,87520,71394,86375,70460,85232,69528,84109,68612,83059,67755,93076,75941,93185,76029,93342,76156,93584,76353,93806,76532,94008,76696,94362,76983,94647,77214,94843,77373,94917,77432,94844,77371,94615,77184,94237,76875,93724,76456,93092,75940,92360,75342,91546,74678,90667,73961,89738,73203,88774,72416,87792,71615,86812,70815,85847,70028,84944,69292,92992,75873,93099,75960,93253,76085,93491,76278,93707,76454,93904,76613,94245,76889,94511,77105,94670,77234,94670,77233,94469,77068,94045,76722,93403,76197,92557,75507,91532,74670,90351,73706,89043,72640,87632,71488,86140,70271,84591,69007,83012,67719,81434,66431,79880,65163,78423,63975,93017,75893,93124,75980,93279,76106,93518,76300,93736,76477,93935,76638,94279,76917,94552,77138,94723,77276,94747,77295,94585,77162,94220,76864,93654,76402,92903,75789,91987,75041,90930,74178,89756,73220,88487,72185,87144,71089,85749,69951,84326,68790,82903,67629,81501,66486,80187,65413,93042,75913,93150,76000,93305,76127,93546,76322,93766,76500,93966,76662,94314,76945,94592,77170,94776,77319,94824,77357,94702,77257,94396,77007,93908,76608,93252,76072,92447,75415,91513,74653,90474,73805,89349,72887,88157,71914,86916,70902,85651,69870,84384,68836,83135,67818,81964,66862,93067,75933,93175,76021,93332,76148,93574,76345,93795,76524,93997,76687,94349,76973,94634,77204,94830,77362,94902,77420,94821,77353,94574,77151,94164,76816,93604,76358,92910,75792,92101,75131,91198,74394,90217,73595,89177,72745,88092,71861,86985,70957,85876,70052,84781,69160,83755,68322,93093,75954,93202,76042,93359,76170,93603,76367,93825,76548,94029,76713,94385,77001,94676,77237,94884,77406,94980,77483,94940,77449,94753,77296,94422,77025,93958,76646,93377,76172,92694,75614,91927,74989,91093,74308,90205,73583,89278,72827,88330,72053,87379,71277,86440,70512,85559,69793,93015,75892,93122,75979,93277,76105,93517,76299,93735,76475,93933,76636,94277,76915,94551,77137,94727,77280,94758,77304,94603,77176,94238,76878,93663,76409,92890,75778,91938,75001,90830,74096,89592,73087,88247,71989,86815,70821,85322,69602,83793,68355,82258,67103,80742,65866,79317,64704,93038,75910,93146,75997,93301,76124,93542,76319,93761,76497,93961,76659,94309,76941,94588,77167,94775,77318,94827,77359,94707,77261,94396,77006,93891,76594,93206,76035,92358,75342,91366,74533,90256,73627,89048,72641,87760,71591,86416,70494,85039,69371,83656,68243,82289,67128,81004,66079,93061,75928,93169,76016,93326,76143,93568,76339,93788,76518,93989,76681,94341,76966,94625,77197,94823,77357,94896,77415,94812,77346,94554,77135,94122,76781,93525,76294,92780,75686,91907,74973,90925,74173,89855,73299,88713,72367,87519,71393,86295,70394,85065,69391,83848,68398,82704,67465,93084,75947,93193,76035,93350,76163,93593,76360,93815,76540,94018,76704,94373,76992,94663,77227,94871,77396,94966,77472,94918,77432,94715,77265,94354,76970,93847,76556,93207,76033,92451,75417,91600,74722,90668,73962,89673,73149,88630,72299,87560,71426,86484,70548,85419,69679,84416,68861,93108,75966,93217,76055,93375,76183,93619,76381,93843,76562,94047,76728,94406,77018,94701,77258,94920,77435,95036,77529,95026,77519,94876,77396,94588,77160,94171,76819,93637,76383,93001,75864,92280,75275,91488,74629,90640,73937,89750,73211,88835,72465,87914,71713,87001,70969,86142,70267,93037,75909,93144,75996,93300,76123,93540,76318,93760,76495,93959,76657,94307,76939,94587,77166,94777,77320,94836,77367,94722,77273,94411,77019,93898,76600,93193,76024,92311,75304,91272,74457,90102,73502,88820,72456,87449,71337,86010,70163,84531,68956,83040,67740,81562,66535,80170,65399,93057,75926,93166,76013,93322,76140,93564,76336,93784,76515,93985,76678,94336,76962,94621,77193,94821,77355,94897,77417,94815,77349,94553,77134,94105,76768,93482,76259,92697,75618,91769,74861,90721,74006,89570,73067,88338,72062,87044,71006,85712,69919,84370,68824,83038,67738,81782,66714,93079,75943,93187,76031,93344,76158,93587,76355,93809,76535,94011,76699,94365,76986,94655,77221,94864,77390,94960,77467,94910,77425,94696,77250,94314,76938,93773,76496,93086,75935,92270,75269,91344,74513,90326,73683,89234,72792,88085,71855,86902,70890,85709,69916,84525,68950,83407,68038,93100,75960,93209,76048,93367,76176,93611,76374,93834,76555,94037,76720,94395,77010,94690,77248,94908,77425,95022,77517,95005,77502,94840,77367,94525,77109,94067,76735,93479,76255,92774,75680,91972,75025,91088,74303,90137,73527,89135,72710,88102,71867,87059,71016,86023,70171,85045,69373,93122,75977,93231,76066,93389,76194,93635,76393,93859,76575,94064,76741,94425,77034,94725,77276,94953,77461,95086,77569,95101,77580,94986,77485,94738,77282,94364,76976,93874,76576,93283,76094,92606,75541,91856,74929,91047,74269,90193,73572,89311,72852,88419,72124,87532,71401,86694,70718,93056,75925,93164,76012,93321,76139,93562,76335,93782,76514,93983,76677,94334,76961,94620,77192,94822,77357,94905,77423,94827,77359,94566,77144,94111,76773,93469,76249,92653,75583,91681,74790,90576,73888,89356,72893,88043,71822,86659,70692,85229,69525,83782,68345,82343,67171,80983,66061,93075,75940,93184,76028,93341,76155,93583,76352,93805,76532,94007,76696,94361,76982,94651,77217,94862,77388,94960,77467,94911,77427,94693,77248,94299,76925,93733,76463,93008,75872,92141,75164,91152,74357,90058,73465,88879,72503,87635,71487,86348,70438,85046,69375,83750,68318,82524,67318,93095,75955,93204,76044,93361,76172,93605,76369,93827,76550,94031,76715,94388,77004,94682,77242,94901,77420,95016,77512,94996,77495,94822,77353,94488,77079,93999,76680,93366,76163,92604,75541,91732,74830,90766,74041,89722,73189,88618,72289,87476,71357,86319,70413,85167,69473,84077,68584,93114,75971,93224,76060,93382,76188,93627,76387,93850,76568,94055,76734,94415,77026,94713,77268,94941,77452,95073,77558,95082,77565,94953,77458,94680,77235,94268,76898,93727,76457,93072,75922,92318,75306,91479,74622,90571,73881,89609,73096,88613,72283,87602,71459,86595,70637,85642,69859,93134,75987,93244,76076,93403,76205,93649,76404,93873,76586,94080,76754,94442,77048,94745,77293,94981,77484,95130,77604,95169,77635,95084,77565,94872,77391,94539,77118,94091,76753,93543,76305,92907,75786,92198,75207,91426,74578,90608,73910,89758,73217,88895,72513,88034,71810,87218,71144,93074,75939,93182,76027,93339,76154,93582,76351,93803,76530,94005,76694,94359,76981,94649,77216,94863,77389,94966,77472,94922,77435,94704,77257,94303,76929,93720,76453,92967,75838,92059,75097,91016,74246,89857,73301,88601,72276,87270,71190,85888,70063,84486,68918,83085,67776,81758,66693,93092,75953,93200,76041,93358,76169,93601,76366,93824,76547,94027,76712,94383,77000,94678,77239,94898,77418,95016,77512,94997,77496,94819,77350,94473,77067,93961,76649,93294,76104,92484,75444,91552,74683,90514,73836,89387,72916,88191,71941,86949,70927,85687,69897,84426,68869,83230,67894,93109,75967,93219,76056,93377,76184,93621,76382,93844,76563,94049,76729,94408,77020,94706,77262,94934,77446,95066,77553,95073,77558,94936,77445,94645,77207,94205,76847,93623,76372,92913,75793,92093,75123,91176,74375,90179,73562,89119,72697,88017,71798,86897,70884,85777,69970,84714,69103,93128,75982,93237,76071,93396,76199,93641,76398,93865,76580,94071,76747,94433,77040,94735,77285,94970,77475,95117,77594,95150,77620,95053,77540,94819,77348,94450,77047,93955,76642,93346,76145,92637,75567,91843,74918,90977,74212,90055,73459,89094,72675,88116,71878,87138,71080,86209,70322,93146,75996,93256,76085,93415,76214,93661,76414,93887,76597,94094,76765,94458,77060,94764,77308,95007,77505,95169,77635,95228,77683,95172,77636,94994,77490,94697,77248,94290,76915,93782,76500,93187,76014,92516,75466,91781,74867,90997,74227,90179,73560,89345,72879,88510,72197,87716,71550,93090,75952,93199,76040,93356,76168,93600,76365,93822,76546,94025,76710,94382,76999,94676,77238,94899,77418,95020,77516,95006,77503,94829,77358,94476,77070,93949,76639,93255,76073,92407,75381,91425,74579,90325,73681,89125,72702,87846,71659,86512,70571,85153,69462,83791,68351,82497,67295,93107,75965,93216,76053,93373,76181,93618,76380,93841,76561,94045,76726,94404,77017,94702,77259,94931,77444,95065,77552,95073,77558,94932,77442,94631,77195,94169,76819,93555,76317,92801,75702,91924,74986,90939,74182,89863,73304,88715,72367,87517,71390,86294,70392,85069,69393,83903,68442,93123,75978,93232,76067,93391,76195,93636,76394,93860,76576,94065,76742,94427,77035,94728,77280,94964,77470,95111,77589,95142,77613,95037,77527,94787,77322,94392,76999,93858,76564,93198,76025,92427,75395,91558,74686,90607,73910,89591,73081,88529,72214,87444,71330,86356,70442,85321,69597,93139,75991,93249,76080,93408,76209,93654,76409,93879,76591,94086,76759,94450,77053,94755,77301,94997,77497,95157,77626,95211,77669,95143,77613,94944,77450,94616,77181,94164,76812,93599,76351,92934,75808,92182,75195,91357,74521,90473,73800,89548,73045,88603,72274,87653,71499,86749,70761,93156,76005,93266,76094,93426,76223,93673,76424,93899,76607,94106,76775,94473,77072,94782,77322,95030,77523,95203,77663,95281,77726,95250,77699,95103,77579,94842,77365,94472,77063,94002,76679,93445,76224,92812,75707,92113,75137,91363,74525,90576,73883,89770,73225,88960,72564,88188,71934,93105,75964,93214,76052,93372,76180,93617,76379,93840,76560,94044,76725,94402,77016,94701,77257,94931,77444,95069,77555,95080,77564,94940,77448,94633,77197,94158,76809,93519,76288,92729,75643,91805,74889,90761,74037,89615,73102,88388,72101,87102,71052,85786,69978,84463,68899,83201,67869,93120,75976,93230,76065,93388,76193,93633,76392,93857,76573,94062,76740,94423,77032,94724,77276,94961,77468,95109,77588,95141,77613,95034,77524,94773,77311,94359,76973,93796,76513,93094,75940,92269,75267,91336,74505,90309,73668,89208,72769,88053,71827,86870,70862,85680,69891,84545,68965,93135,75988,93245,76077,93404,76206,93649,76405,93874,76587,94081,76754,94444,77049,94748,77296,94990,77491,95150,77621,95203,77663,95128,77601,94915,77426,94562,77138,94074,76739,93461,76239,92737,75648,91915,74977,91008,74237,90034,73442,89011,72607,87962,71752,86906,70890,85898,70068,93150,76000,93260,76089,93420,76218,93666,76418,93892,76601,94099,76769,94465,77066,94773,77315,95020,77516,95192,77654,95265,77713,95224,77678,95058,77542,94766,77304,94355,76968,93832,76540,93209,76032,92499,75452,91713,74811,90867,74120,89977,73394,89063,72649,88142,71897,87262,71179,93166,76012,93276,76102,93436,76231,93683,76432,93910,76616,94118,76785,94486,77083,94797,77335,95050,77540,95234,77688,95328,77764,95320,77756,95202,77659,94973,77472,94638,77198,94205,76844,93685,76419,93087,75931,92423,75390,91706,74804,90950,74187,90172,73552,89387,72912,88636,72299,93119,75975,93228,76064,93387,76192,93632,76391,93856,76572,94061,76739,94421,77031,94723,77275,94961,77468,95112,77590,95147,77618,95040,77530,94775,77313,94348,76964,93762,76485,93027,75885,92158,75176,91169,74369,90076,73477,88899,72517,87660,71506,86386,70467,85102,69419,83873,68417,93133,75986,93242,76075,93401,76203,93647,76403,93871,76585,94077,76752,94440,77046,94745,77293,94987,77489,95149,77620,95202,77662,95125,77598,94902,77416,94531,77113,94016,76692,93364,76160,92590,75528,91706,74807,90728,74008,89672,73147,88560,72239,87416,71306,86261,70364,85156,69462,93146,75997,93256,76086,93415,76215,93662,76415,93887,76598,94094,76766,94459,77061,94767,77310,95014,77511,95186,77649,95258,77707,95210,77667,95030,77520,94717,77263,94272,76900,93704,76436,93025,75882,92248,75248,91384,74543,90451,73782,89467,72979,88453,72151,87429,71316,86448,70516,93160,76008,93270,76097,93430,76227,93677,76427,93904,76611,94112,76779,94478,77077,94789,77328,95041,77533,95223,77679,95314,77752,95296,77737,95160,77625,94904,77415,94530,77110,94046,76715,93464,76240,92793,75692,92046,75082,91236,74422,90381,73723,89498,73003,88605,72275,87750,71577,93174,76019,93285,76109,93445,76238,93693,76440,93920,76624,94129,76793,94498,77092,94811,77346,95069,77555,95261,77710,95370,77798,95383,77807,95291,77731,95092,77568,94790,77322,94392,76996,93906,76600,93344,76140,92713,75626,92028,75066,91302,74474,90551,73861,89791,73241,89062,72646,93132,75985,93241,76074,93400,76203,93646,76402,93870,76584,94076,76751,94439,77045,94743,77292,94987,77489,95151,77621,95207,77666,95130,77602,94903,77417,94521,77104,93985,76666,93302,76109,92486,75443,91549,74679,90508,73829,89380,72908,88187,71935,86956,70931,85709,69914,84514,68939,93144,75995,93254,76084,93413,76213,93659,76413,93885,76596,94091,76763,94456,77059,94763,77307,95011,77508,95184,77648,95256,77706,95206,77664,95018,77510,94689,77241,94218,76856,93613,76363,92887,75770,92052,75088,91120,74328,90109,73503,89038,72629,87933,71727,86813,70814,85738,69937,93157,76005,93267,76094,93426,76224,93673,76424,93899,76607,94107,76776,94473,77073,94783,77323,95036,77528,95217,77675,95306,77746,95283,77726,95135,77605,94858,77378,94453,77048,93928,76618,93292,76100,92558,75500,91737,74830,90844,74101,89897,73329,88918,72530,87925,71720,86971,70942,93169,76015,93280,76105,93440,76234,93687,76435,93914,76619,94123,76788,94491,77087,94803,77340,95060,77548,95251,77702,95357,77787,95361,77789,95252,77700,95028,77517,94690,77240,94245,76876,93700,76432,93068,75916,92358,75336,91584,74705,90762,74034,89910,73339,89045,72633,88214,71955,93182,76025,93293,76115,93453,76245,93702,76447,93929,76631,94139,76801,94509,77101,94824,77356,95085,77568,95285,77729,95407,77828,95439,77853,95371,77796,95200,77656,94930,77435,94564,77136,94112,76767,93582,76334,92985,75847,92330,75313,91633,74744,90910,74154,90174,73553,89466,72975,93143,75994,93253,76083,93412,76212,93658,76412,93884,76595,94090,76762,94455,77058,94762,77306,95011,77508,95185,77649,95260,77709,95210,77668,95019,77511,94679,77233,94189,76833,93555,76315,92790,75691,91905,74968,90913,74159,89833,73277,88685,72341,87495,71371,86287,70385,85124,69436,93154,76003,93265,76093,93424,76222,93671,76422,93897,76605,94104,76773,94470,77070,94780,77321,95033,77526,95215,77673,95305,77745,95279,77723,95124,77596,94832,77357,94403,77007,93844,76550,93164,75995,92374,75351,91488,74627,90520,73838,89491,72998,88423,72127,87338,71241,86292,70389,93166,76012,93276,76102,93436,76231,93684,76432,93910,76616,94119,76785,94486,77083,94798,77336,95055,77544,95245,77697,95349,77781,95348,77780,95229,77681,94986,77482,94619,77183,94134,76786,93540,76302,92847,75736,92067,75099,91213,74402,90303,73660,89358,72889,88396,72104,87469,71348,93178,76022,93288,76111,93448,76241,93696,76443,93924,76627,94133,76796,94503,77096,94817,77350,95078,77562,95276,77722,95395,77818,95419,77836,95336,77768,95142,77609,94837,77360,94427,77025,93920,76610,93324,76124,92650,75574,91911,74971,91121,74326,90300,73656,89463,72973,88655,72315,93190,76031,93300,76121,93461,76251,93710,76453,93938,76638,94147,76808,94519,77109,94836,77365,95100,77580,95307,77747,95441,77855,95489,77893,95443,77855,95299,77736,95057,77539,94722,77265,94302,76922,93805,76515,93238,76053,92614,75544,91946,74998,91249,74430,90538,73849,89850,73289,93154,76003,93264,76092,93423,76221,93670,76421,93896,76604,94103,76773,94469,77069,94779,77320,95032,77525,95216,77674,95308,77748,95283,77726,95124,77596,94822,77349,94376,76985,93790,76506,93073,75921,92236,75238,91293,74469,90259,73625,89155,72724,88007,71788,86836,70832,85706,69911,93164,76011,93274,76100,93434,76230,93681,76430,93908,76614,94116,76783,94484,77081,94795,77333,95052,77541,95243,77696,95348,77780,95345,77777,95219,77673,94962,77463,94573,77145,94056,76723,93420,76204,92675,75596,91832,74908,90907,74153,89918,73346,88888,72505,87836,71648,86821,70819,93175,76019,93285,76109,93445,76239,93693,76440,93920,76624,94129,76793,94498,77092,94812,77347,95072,77558,95271,77718,95388,77813,95407,77827,95314,77751,95103,77577,94772,77307,94325,76942,93770,76489,93117,75956,92376,75351,91560,74685,90686,73972,89775,73228,88844,72469,87944,71735,93185,76028,93296,76118,93456,76248,93705,76449,93933,76634,94142,76804,94513,77104,94829,77360,95093,77574,95298,77740,95429,77846,95471,77878,95411,77829,95245,77693,94972,77469,94596,77162,94123,76776,93563,76318,92924,75797,92218,75221,91460,74603,90669,73957,89859,73296,89075,72657,93196,76037,93307,76126,93468,76257,93717,76459,93945,76644,94156,76814,94528,77116,94846,77374,95114,77591,95326,77762,95470,77879,95534,77930,95508,77908,95388,77809,95174,77634,94869,77384,94479,77065,94012,76684,93476,76247,92881,75761,92241,75238,91570,74691,90882,74130,90216,73586,93163,76010,93273,76100,93433,76229,93680,76430,93907,76613,94115,76782,94482,77080,94794,77332,95051,77541,95244,77697,95350,77782,95348,77779,95219,77673,94953,77456,94548,77125,94006,76682,93336,76135,92546,75491,91649,74759,90661,73952,89600,73087,88492,72183,87359,71258,86261,70363,93173,76018,93283,76107,93443,76237,93691,76438,93918,76622,94127,76791,94496,77090,94809,77344,95070,77555,95268,77716,95386,77811,95404,77824,95305,77743,95081,77559,94729,77272,94252,76883,93658,76398,92956,75824,92156,75171,91271,74450,90321,73675,89328,72864,88310,72034,87324,71229,93182,76026,93293,76115,93453,76245,93702,76447,93929,76631,94139,76801,94509,77101,94824,77357,95088,77570,95293,77736,95423,77841,95460,77870,95391,77813,95209,77664,94911,77420,94501,77085,93984,76663,93369,76161,92666,75587,91887,74951,91048,74267,90170,73550,89269,72815,88396,72103,93192,76033,93303,76123,93463,76253,93713,76455,93941,76640,94151,76811,94523,77112,94840,77369,95107,77585,95318,77756,95460,77870,95517,77916,95479,77884,95339,77769,95095,77570,94751,77288,94312,76930,93785,76500,93180,76005,92507,75456,91781,74863,91018,74241,90235,73603,89475,72983,93202,76041,93313,76131,93474,76262,93724,76464,93952,76650,94163,76820,94536,77123,94856,77382,95126,77600,95343,77776,95497,77900,95575,77963,95567,77956,95470,77875,95281,77721,95003,77493,94642,77199,94205,76841,93698,76427,93132,75965,92519,75465,91873,74938,91209,74396,90563,73869,93172,76017,93282,76107,93442,76236,93690,76437,93917,76622,94126,76791,94495,77090,94808,77343,95069,77555,95269,77717,95388,77813,95406,77826,95304,77743,95072,77553,94706,77253,94206,76845,93580,76334,92835,75726,91984,75031,91039,74260,90021,73429,88952,72558,87855,71663,86790,70793,93181,76024,93291,76114,93451,76244,93700,76445,93927,76630,94137,76799,94507,77099,94822,77355,95086,77568,95291,77734,95421,77839,95457,77867,95382,77806,95189,77647,94872,77388,94434,77030,93880,76578,93218,76038,92459,75418,91614,74729,90703,73985,89745,73204,88761,72401,87804,71620,93190,76031,93300,76121,93461,76251,93710,76453,93938,76638,94147,76808,94519,77109,94836,77366,95103,77582,95313,77752,95454,77866,95508,77909,95461,77870,95306,77742,95040,77524,94663,77217,94182,76824,93604,76352,92937,75808,92194,75201,91390,74545,90544,73855,89673,73144,88826,72453,93199,76038,93309,76128,93470,76259,93720,76461,93948,76646,94159,76817,94531,77119,94850,77377,95120,77595,95336,77770,95487,77893,95559,77950,95541,77934,95424,77839,95209,77662,94894,77405,94487,77072,93993,76668,93420,76201,92779,75677,92083,75110,91349,74511,90592,73894,89856,73292,93208,76046,93319,76136,93480,76266,93730,76469,93959,76655,94170,76826,94544,77129,94865,77389,95137,77609,95359,77789,95521,77920,95611,77992,95621,77999,95544,77935,95379,77800,95127,77594,94794,77322,94385,76988,93906,76597,93367,76157,92781,75678,92160,75172,91519,74648,90893,74138,93180,76024,93290,76113,93451,76243,93699,76445,93926,76629,94136,76799,94506,77098,94821,77354,95085,77568,95291,77735,95423,77841,95459,77869,95382,77806,95181,77641,94851,77371,94390,76995,93806,76518,93105,75946,92297,75286,91395,74550,90418,73753,89388,72913,88327,72047,87293,71204,93188,76030,93299,76120,93459,76250,93708,76452,93936,76637,94146,76806,94517,77107,94833,77364,95100,77580,95311,77751,95452,77864,95504,77906,95453,77863,95287,77727,95003,77495,94601,77166,94085,76745,93463,76237,92743,75649,91937,74992,91063,74279,90141,73526,89189,72749,88261,71992,93196,76037,93307,76126,93468,76257,93717,76459,93945,76644,94156,76814,94528,77116,94846,77374,95116,77592,95332,77767,95482,77888,95550,77943,95524,77921,95394,77814,95157,77620,94813,77338,94366,76974,93823,76530,93192,76015,92484,75437,91713,74808,90898,74144,90057,73457,89236,72787,93204,76043,93315,76133,93476,76264,93726,76466,93955,76652,94166,76823,94539,77125,94859,77385,95131,77605,95352,77783,95512,77912,95597,77981,95596,77979,95502,77902,95312,77746,95026,77512,94649,77204,94186,76826,93644,76384,93034,75885,92368,75342,91662,74766,90931,74170,90218,73587,93213,76050,93324,76140,93485,76271,93735,76473,93965,76659,94176,76831,94551,77135,94873,77395,95147,77617,95373,77800,95542,77937,95644,78019,95669,78038,95611,77990,95469,77873,95242,77687,94935,77436,94552,77124,94101,76755,93589,76337,93028,75879,92432,75393,91813,74888,91207,74394,93187,76029,93298,76119,93458,76249,93707,76451,93935,76636,94145,76806,94516,77107,94832,77363,95099,77579,95311,77751,95453,77865,95506,77907,95452,77863,95280,77722,94983,77479,94561,77134,94017,76689,93357,76151,92591,75525,91730,74823,90794,74059,89802,73250,88776,72413,87774,71595,93195,76035,93305,76125,93466,76255,93715,76458,93944,76643,94154,76813,94526,77115,94844,77372,95113,77590,95330,77765,95480,77887,95547,77941,95516,77915,95377,77800,95124,77593,94755,77292,94276,76900,93691,76423,93009,75866,92241,75239,91403,74556,90515,73831,89596,73081,88696,72347,93202,76041,93313,76131,93474,76262,93723,76464,93952,76650,94163,76820,94536,77123,94856,77382,95127,77602,95348,77780,95507,77908,95589,77974,95581,77967,95475,77879,95265,77708,94951,77451,94537,77112,94027,76697,93430,76209,92756,75659,92018,75056,91234,74417,90421,73754,89626,73105,93210,76047,93321,76137,93482,76268,93732,76471,93961,76657,94172,76828,94547,77131,94868,77391,95142,77613,95367,77795,95534,77930,95631,78008,95646,78020,95573,77959,95407,77823,95148,77611,94800,77327,94366,76973,93855,76555,93274,76081,92637,75561,91958,75007,91253,74432,90562,73868,93217,76053,93329,76144,93490,76274,93740,76477,93970,76664,94182,76835,94557,77140,94880,77401,95156,77624,95385,77810,95561,77952,95673,78042,95712,78073,95673,78040,95551,77940,95347,77773,95065,77542,94708,77251,94283,76904,93797,76507,93261,76069,92688,75602,92092,75115,91506,74637,93194,76035,93305,76125,93465,76255,93714,76457,93943,76642,94153,76812,94525,77114,94843,77372,95112,77590,95329,77765,95481,77887,95548,77942,95516,77914,95370,77795,95105,77578,94718,77262,94212,76848,93592,76342,92866,75750,92046,75081,91148,74348,90194,73569,89203,72761,88231,71968,93201,76040,93312,76130,93472,76261,93722,76463,93951,76648,94161,76819,94534,77121,94854,77380,95125,77600,95346,77779,95505,77907,95586,77972,95574,77961,95459,77867,95234,77682,94897,77407,94453,77044,93904,76596,93258,76069,92527,75472,91725,74818,90870,74121,89982,73396,89111,72685,93208,76046,93319,76136,93480,76266,93729,76469,93959,76655,94170,76826,94544,77129,94865,77389,95138,77610,95363,77792,95529,77926,95623,78002,95633,78009,95548,77939,95364,77788,95078,77555,94695,77241,94218,76852,93654,76391,93012,75867,92306,75291,91552,74676,90768,74036,89997,73407,93215,76051,93326,76141,93487,76272,93737,76475,93967,76661,94178,76832,94553,77137,94876,77397,95151,77620,95380,77806,95553,77946,95661,78033,95692,78057,95638,78012,95494,77894,95261,77703,94939,77440,94535,77110,94052,76716,93500,76265,92891,75768,92239,75236,91559,74681,90890,74135,93222,76057,93333,76147,93494,76278,93745,76481,93975,76668,94187,76839,94563,77144,94887,77406,95164,77631,95397,77819,95578,77966,95699,78064,95752,78105,95728,78085,95626,78001,95445,77852,95186,77641,94854,77369,94454,77042,93993,76666,93481,76248,92932,75800,92357,75331,91790,74869,93200,76040,93311,76130,93472,76260,93721,76462,93950,76648,94160,76818,94534,77121,94853,77379,95124,77599,95346,77778,95505,77907,95587,77973,95573,77961,95453,77862,95217,77669,94863,77380,94393,76996,93811,76521,93125,75960,92343,75323,91483,74621,90565,73872,89608,73091,88667,72323,93206,76045,93317,76135,93478,76265,93728,76468,93957,76654,94168,76824,94542,77128,94863,77387,95136,77608,95361,77790,95527,77925,95620,78000,95626,78003,95533,77927,95335,77765,95029,77514,94617,77178,94102,76758,93492,76260,92796,75692,92028,75065,91207,74395,90350,73695,89506,73007,93213,76050,93324,76140,93485,76270,93735,76473,93965,76659,94176,76831,94551,77134,94873,77395,95148,77618,95376,77803,95549,77942,95654,78027,95679,78047,95614,77993,95454,77861,95196,77650,94842,77361,94396,76996,93863,76562,93253,76063,92577,75512,91854,74922,91097,74304,90351,73696,93219,76055,93330,76145,93491,76276,93742,76479,93972,76665,94184,76837,94560,77141,94883,77403,95160,77627,95391,77815,95571,77960,95688,78055,95733,78090,95696,78059,95574,77959,95364,77787,95069,77546,94691,77237,94237,76866,93713,76438,93131,75963,92505,75452,91850,74918,91203,74390,93225,76060,93337,76150,93498,76281,93749,76485,93979,76671,94191,76843,94568,77149,94893,77411,95172,77637,95407,77827,95593,77978,95723,78083,95787,78134,95779,78127,95695,78057,95534,77925,95298,77732,94989,77480,94613,77172,94176,76816,93688,76417,93162,75987,92608,75536,92061,75089,93206,76044,93317,76134,93478,76265,93727,76467,93957,76653,94167,76824,94541,77127,94862,77387,95135,77608,95361,77790,95528,77925,95621,78000,95625,78003,95528,77922,95320,77752,94997,77489,94562,77133,94016,76687,93367,76157,92623,75550,91800,74879,90918,74159,89994,73405,89083,72662,93211,76049,93323,76139,93484,76270,93734,76472,93963,76658,94174,76829,94549,77133,94871,77394,95146,77616,95374,77801,95547,77941,95651,78025,95673,78042,95601,77982,95428,77840,95150,77613,94769,77302,94288,76909,93712,76438,93050,75898,92315,75299,91526,74655,90699,73980,89883,73314,93217,76053,93328,76144,93490,76274,93740,76477,93970,76664,94181,76835,94557,77140,94880,77401,95157,77625,95388,77812,95567,77957,95682,78050,95721,78081,95675,78042,95537,77929,95305,77739,94978,77472,94562,77132,94059,76721,93479,76248,92834,75721,92139,75154,91410,74559,90688,73970,93223,76058,93334,76148,93496,76279,93746,76482,93977,76669,94189,76841,94565,77146,94889,77408,95167,77634,95402,77823,95587,77973,95713,78075,95770,78120,95750,78103,95647,78019,95460,77865,95189,77644,94838,77356,94410,77007,93913,76601,93357,76147,92757,75658,92125,75142,91500,74632,93229,76063,93340,76153,93502,76284,93753,76488,93983,76674,94196,76847,94573,77152,94898,77416,95178,77642,95416,77835,95607,77989,95744,78100,95819,78160,95825,78164,95759,78109,95617,77993,95402,77817,95116,77583,94763,77294,94349,76956,93884,76577,93379,76165,92847,75730,92318,75298,93211,76048,93322,76138,93483,76269,93733,76472,93963,76658,94174,76829,94549,77133,94870,77393,95145,77616,95374,77801,95547,77941,95652,78026,95672,78041,95596,77978,95414,77829,95121,77589,94718,77260,94207,76843,93594,76342,92886,75765,92099,75123,91251,74431,90361,73704,89479,72985,93216,76052,93327,76143,93489,76273,93739,76476,93969,76663,94180,76834,94556,77138,94878,77400,95155,77624,95386,77811,95565,77956,95679,78048,95716,78076,95663,78032,95513,77909,95262,77704,94911,77417,94461,77049,93917,76606,93288,76092,92586,75519,91828,74901,91031,74250,90241,73606,93221,76057,93333,76147,93494,76278,93745,76481,93975,76667,94187,76839,94563,77144,94886,77406,95165,77631,95399,77821,95583,77970,95707,78070,95760,78112,95730,78087,95613,77991,95405,77820,95105,77575,94716,77257,94243,76871,93693,76422,93077,75919,92410,75374,91707,74801,91009,74232,93227,76061,93338,76151,93500,76282,93751,76486,93981,76672,94193,76844,94570,77150,94895,77413,95174,77639,95411,77831,95601,77985,95735,78093,95804,78148,95799,78142,95715,78073,95549,77937,95301,77734,94974,77467,94572,77139,94101,76754,93571,76321,92995,75852,92387,75356,91783,74862,93232,76065,93344,76156,93505,76287,93756,76490,93987,76677,94200,76850,94578,77156,94903,77420,95184,77647,95424,77841,95620,77999,95763,78115,95848,78184,95868,78198,95817,78156,95693,78055,95498,77895,95234,77679,94903,77408,94512,77089,94069,76727,93586,76333,93073,75915,92563,75498,93216,76052,93327,76142,93488,76273,93738,76476,93968,76662,94180,76834,94555,77138,94877,77399,95154,77623,95386,77811,95565,77956,95680,78048,95715,78076,95658,78028,95500,77899,95235,77682,94863,77378,94386,76988,93807,76515,93134,75967,92382,75353,91568,74688,90709,73988,89856,73292,93220,76056,93332,76146,93493,76277,93744,76480,93974,76667,94185,76838,94562,77143,94885,77405,95163,77630,95397,77820,95581,77969,95705,78068,95754,78108,95719,78078,95591,77973,95366,77788,95042,77524,94622,77181,94110,76762,93512,76275,92842,75728,92115,75134,91347,74507,90583,73884,93225,76060,93337,76150,93498,76281,93749,76484,93979,76671,94191,76843,94568,77148,94893,77411,95172,77637,95409,77829,95598,77982,95730,78088,95794,78140,95781,78128,95683,78048,95497,77895,95223,77671,94861,77375,94415,77011,93894,76585,93305,76105,92666,75583,91989,75031,91315,74481,93230,76064,93342,76154,93503,76285,93754,76489,93985,76676,94197,76848,94575,77154,94900,77417,95181,77644,95420,77838,95614,77995,95755,78109,95834,78172,95843,78179,95776,78123,95630,78003,95405,77819,95101,77571,94724,77262,94278,76898,93773,76486,93221,76036,92636,75558,92053,75082,93235,76067,93347,76158,93508,76289,93760,76493,93991,76680,94203,76853,94582,77159,94908,77423,95190,77652,95432,77847,95631,78008,95781,78129,95875,78205,95906,78229,95870,78199,95764,78112,95588,77968,95344,77768,95034,77515,94665,77213,94243,76869,93781,76492,93288,76090,92796,75688,93220,76056,93331,76146,93493,76277,93743,76480,93973,76666,94185,76838,94561,77143,94884,77404,95162,77629,95397,77819,95582,77969,95705,78068,95754,78107,95715,78074,95579,77963,95341,77768,94998,77488,94552,77124,94006,76678,93368,76157,92649,75570,91868,74933,91040,74258,90216,73585,93224,76059,93336,76149,93497,76280,93748,76484,93978,76670,94190,76842,94567,77147,94891,77410,95170,77636,95407,77828,95596,77980,95728,78087,95789,78136,95771,78120,95663,78031,95461,77866,95164,77623,94773,77303,94290,76909,93724,76447,93084,75924,92386,75355,91646,74752,90908,74149,93229,76062,93340,76153,93502,76284,93753,76487,93983,76674,94196,76846,94573,77152,94898,77415,95178,77642,95417,77836,95611,77992,95750,78105,95825,78165,95827,78165,95747,78100,95583,77965,95332,77759,94995,77485,94577,77143,94082,76739,93522,76281,92909,75781,92257,75249,91606,74718,93233,76066,93345,76157,93506,76288,93758,76491,93988,76678,94201,76851,94579,77157,94905,77421,95187,77649,95428,77844,95626,78004,95773,78123,95862,78194,95884,78211,95833,78169,95705,78064,95501,77897,95220,77667,94866,77378,94444,77033,93963,76641,93435,76210,92872,75751,92309,75291,93238,76070,93349,76160,93511,76291,93763,76495,93994,76683,94207,76855,94586,77162,94912,77427,95195,77656,95439,77853,95641,78016,95796,78142,95898,78224,95941,78258,95918,78239,95829,78165,95671,78035,95446,77851,95157,77615,94808,77330,94408,77003,93966,76642,93493,76256,93018,75868,93224,76059,93335,76149,93497,76280,93747,76483,93978,76670,94190,76842,94566,77147,94890,77409,95170,77635,95407,77827,95596,77981,95728,78087,95789,78135,95766,78116,95652,78022,95438,77847,95123,77590,94707,77250,94193,76830,93587,76336,92901,75776,92152,75164,91355,74514,90559,73864,93228,76062,93339,76152,93501,76283,93752,76487,93982,76674,94195,76846,94572,77151,94897,77414,95177,77641,95416,77835,95609,77991,95748,78103,95821,78162,95817,78158,95729,78085,95549,77938,95277,77715,94913,77418,94460,77047,93922,76609,93312,76110,92643,75565,91931,74984,91218,74402,93232,76065,93343,76155,93505,76287,93756,76490,93987,76677,94200,76850,94578,77156,94903,77419,95185,77647,95425,77842,95623,78002,95769,78120,95854,78188,95869,78199,95806,78147,95661,78029,95433,77842,95121,77587,94728,77266,94260,76884,93726,76448,93139,75968,92512,75457,91883,74944,93236,76068,93348,76159,93509,76290,93761,76494,93992,76681,94205,76853,94583,77160,94910,77425,95192,77653,95435,77850,95636,78013,95789,78136,95887,78215,95920,78241,95884,78211,95774,78121,95590,77969,95330,77757,94999,77486,94600,77161,94143,76787,93638,76375,93097,75933,92553,75490,93240,76072,93352,76162,93514,76293,93765,76497,93997,76685,94210,76858,94589,77165,94916,77430,95200,77659,95445,77858,95650,78023,95810,78153,95920,78242,95973,78284,95963,78275,95889,78213,95748,78098,95542,77929,95272,77709,94944,77440,94563,77129,94141,76785,93686,76414,93229,76040,93227,76062,93339,76152,93500,76283,93751,76486,93982,76673,94194,76845,94571,77151,94896,77414,95176,77641,95416,77834,95609,77991,95748,78103,95820,78161,95814,78155,95718,78076,95528,77920,95239,77684,94852,77368,94368,76973,93794,76504,93140,75970,92421,75384,91654,74758,90885,74130,93231,76064,93343,76155,93504,76286,93755,76490,93986,76677,94199,76849,94576,77155,94902,77419,95183,77646,95424,77841,95621,78001,95767,78118,95850,78185,95860,78192,95789,78133,95631,78004,95383,77801,95044,77524,94618,77176,94109,76761,93528,76286,92887,75763,92202,75204,91514,74643,93235,76067,93346,76158,93508,76289,93760,76493,93990,76680,94203,76852,94582,77159,94908,77423,95190,77652,95433,77848,95633,78010,95785,78133,95879,78209,95907,78230,95860,78191,95734,78088,95527,77918,95238,77682,94870,77381,94427,77020,93919,76605,93357,76146,92754,75654,92148,75159,93239,76070,93350,76161,93512,76292,93764,76496,93995,76683,94208,76856,94587,77163,94914,77428,95197,77657,95441,77855,95645,78020,95804,78148,95909,78233,95954,78268,95932,78249,95838,78173,95672,78037,95434,77841,95124,77588,94748,77281,94313,76925,93830,76531,93310,76107,92786,75679,93242,76073,93354,76164,93516,76295,93768,76499,93999,76687,94212,76860,94592,77167,94920,77433,95204,77663,95450,77862,95658,78030,95823,78163,95940,78257,96002,78307,96004,78308,95944,78258,95819,78156,95631,78001,95380,77796,95071,77544,94710,77249,94306,76919,93871,76564,93430,76204,93231,76064,93342,76155,93504,76286,93755,76489,93986,76676,94198,76848,94576,77155,94901,77418,95183,77646,95424,77841,95621,78001,95767,78118,95849,78184,95856,78189,95779,78126,95611,77988,95347,77772,94987,77478,94533,77107,93989,76662,93365,76153,92677,75592,91939,74989,91196,74384,93234,76067,93346,76157,93507,76288,93759,76492,93990,76679,94202,76852,94581,77158,94907,77422,95189,77651,95431,77847,95632,78009,95783,78132,95876,78206,95899,78224,95844,78178,95706,78065,95480,77880,95167,77624,94767,77297,94285,76904,93731,76452,93118,75951,92459,75414,91795,74872,93237,76069,93349,76160,93511,76291,93762,76495,93994,76682,94206,76855,94585,77162,94912,77427,95195,77656,95439,77853,95643,78018,95800,78145,95903,78228,95941,78258,95909,78231,95801,78142,95614,77989,95347,77771,95003,77489,94585,77148,94101,76753,93563,76314,92984,75841,92399,75364,93241,76072,93353,76163,93515,76294,93766,76498,93997,76686,94211,76858,94590,77166,94917,77431,95201,77661,95447,77860,95654,78027,95817,78159,95930,78250,95984,78293,95975,78285,95897,78221,95749,78099,95530,77919,95241,77683,94886,77393,94473,77056,94011,76679,93513,76272,93008,75860,93244,76075,93356,76166,93518,76297,93770,76501,94002,76689,94215,76862,94595,77170,94923,77435,95208,77666,95455,77866,95665,78036,95834,78173,95957,78272,96028,78328,96042,78339,95995,78300,95885,78210,95714,78069,95481,77878,95190,77641,94848,77361,94463,77047,94045,76706,93621,76360,93234,76067,93345,76157,93507,76288,93758,76492,93989,76679,94202,76851,94580,77158,94906,77422,95188,77650,95431,77847,95632,78009,95783,78132,95875,78205,95896,78221,95835,78171,95687,78050,95447,77853,95113,77580,94687,77232,94172,76811,93578,76327,92918,75789,92209,75210,91492,74625,93237,76069,93348,76160,93510,76291,93762,76495,93993,76682,94206,76854,94584,77161,94911,77426,95194,77655,95438,77852,95642,78017,95798,78144,95899,78225,95934,78252,95895,78220,95775,78121,95571,77954,95281,77717,94906,77411,94451,77039,93924,76609,93337,76130,92704,75613,92063,75090,93240,76071,93352,76162,93514,76293,93765,76497,93996,76685,94209,76857,94589,77165,94916,77430,95200,77659,95445,77858,95651,78025,95814,78156,95924,78245,95973,78284,95955,78268,95863,78193,95695,78055,95449,77854,95127,77591,94733,77268,94273,76893,93759,76473,93202,76019,92639,75559,93243,76074,93355,76165,93517,76296,93769,76500,94000,76688,94213,76860,94593,77168,94921,77434,95205,77664,95452,77864,95661,78033,95829,78168,95948,78265,96012,78316,96015,78317,95952,78265,95820,78157,95619,77992,95350,77772,95017,77500,94624,77179,94184,76819,93705,76429,93219,76032,93246,76076,93358,76167,93520,76299,93772,76503,94004,76691,94217,76864,94597,77172,94926,77437,95211,77669,95460,77870,95672,78041,95844,78181,95973,78285,96052,78348,96076,78367,96042,78338,95947,78260,95791,78132,95575,77955,95303,77732,94979,77468,94612,77168,94211,76841,93803,76508,93237,76069,93348,76159,93510,76291,93761,76494,93992,76682,94205,76854,94584,77161,94910,77425,95193,77654,95438,77852,95641,78017,95798,78144,95899,78224,95931,78250,95887,78213,95758,78107,95540,77929,95231,77676,94831,77349,94344,76952,93778,76490,93148,75975,92466,75419,91775,74855,93239,76071,93351,76162,93513,76293,93765,76497,93996,76684,94209,76857,94588,77164,94915,77429,95199,77659,95444,77857,95650,78024,95812,78155,95920,78242,95966,78279,95941,78257,95839,78173,95655,78022,95387,77803,95037,77517,94607,77166,94105,76757,93544,76298,92936,75802,92319,75298,93242,76073,93354,76164,93516,76295,93768,76499,93999,76687,94212,76860,94592,77167,94919,77432,95204,77663,95451,77862,95659,78031,95826,78166,95943,78260,96002,78307,95996,78302,95920,78239,95770,78116,95544,77931,95244,77686,94872,77382,94435,77025,93944,76624,93410,76188,92867,75745,93245,76075,93357,76166,93519,76298,93771,76502,94002,76689,94216,76862,94596,77170,94924,77436,95209,77667,95457,77868,95668,78038,95840,78177,95965,78278,96038,78336,96052,78347,96002,78306,95886,78211,95703,78060,95453,77856,95139,77599,94767,77295,94347,76952,93888,76578,93419,76195,93248,76078,93360,76169,93522,76300,93774,76504,94006,76692,94219,76865,94600,77173,94928,77440,95214,77671,95464,77873,95678,78046,95854,78188,95988,78296,96074,78366,96108,78392,96085,78373,96004,78306,95863,78191,95664,78027,95408,77819,95102,77568,94753,77283,94369,76970,93976,76649,93239,76071,93351,76161,93513,76293,93764,76497,93995,76684,94208,76857,94587,77164,94914,77429,95198,77658,95444,77857,95650,78024,95812,78155,95920,78242,95964,78276,95934,78251,95823,78160,95626,77999,95341,77765,94966,77459,94506,77084,93968,76645,93365,76152,92710,75618,92044,75074,93242,76073,93353,76163,93515,76295,93767,76499,93998,76686,94212,76859,94591,77167,94919,77432,95203,77662,95449,77862,95658,78030,95824,78165,95940,78258,95996,78302,95984,78292,95898,78221,95733,78085,95486,77884,95159,77617,94754,77285,94277,76896,93741,76458,93157,75982,92562,75497,93244,76075,93356,76166,93518,76297,93770,76501,94001,76689,94215,76862,94595,77169,94923,77435,95208,77666,95455,77866,95666,78037,95837,78175,95960,78274,96028,78329,96034,78333,95973,78282,95840,78173,95633,78003,95353,77775,95003,77488,94588,77150,94119,76767,93607,76349,93084,75922,93247,76077,93359,76168,93521,76299,93773,76503,94004,76691,94218,76864,94598,77172,94927,77438,95213,77670,95461,77871,95675,78043,95849,78185,95980,78290,96061,78355,96085,78374,96049,78343,95948,78260,95781,78124,95549,77934,95254,77693,94902,77405,94501,77078,94061,76719,93611,76351,93249,76079,93361,76170,93524,76301,93776,76506,94008,76694,94221,76867,94602,77175,94931,77442,95217,77674,95468,77876,95683,78050,95862,78195,96001,78307,96094,78382,96137,78416,96125,78405,96057,78349,95930,78245,95747,78095,95508,77899,95218,77663,94886,77391,94519,77092,94141,76783,93241,76073,93353,76163,93515,76295,93767,76499,93998,76686,94211,76859,94591,77166,94918,77431,95202,77662,95449,77861,95658,78030,95824,78164,95939,78257,95993,78300,95977,78287,95883,78209,95706,78064,95443,77849,95093,77563,94659,77208,94147,76791,93571,76320,92942,75807,92300,75283,93244,76074,93355,76165,93518,76297,93769,76501,94001,76688,94214,76861,94594,77169,94922,77434,95207,77665,95454,77866,95665,78036,95835,78173,95957,78272,96023,78324,96023,78324,95952,78265,95805,78144,95579,77959,95274,77710,94892,77398,94439,77028,93927,76610,93367,76153,92794,75685,93246,76076,93358,76167,93520,76298,93772,76503,94004,76690,94217,76863,94597,77171,94926,77437,95211,77669,95460,77870,95673,78042,95847,78183,95976,78287,96052,78348,96069,78361,96021,78321,95904,78225,95716,78071,95456,77858,95126,77589,94733,77268,94286,76902,93795,76502,93291,76091,93248,76078,93360,76169,93523,76300,93775,76505,94006,76693,94220,76866,94601,77174,94929,77440,95216,77672,95465,77874,95680,78048,95858,78192,95994,78302,96082,78372,96116,78399,96091,78378,96005,78307,95854,78183,95639,78008,95363,77781,95029,77509,94647,77197,94226,76853,93793,76500,93251,76080,93363,76171,93525,76303,93777,76507,94009,76695,94223,76868,94604,77177,94933,77443,95220,77676,95471,77879,95688,78054,95870,78201,96013,78317,96112,78397,96163,78437,96162,78435,96106,78389,95993,78296,95824,78158,95601,77976,95328,77752,95012,77494,94661,77208,94298,76911,93243,76074,93355,76165,93517,76296,93769,76501,94001,76688,94214,76861,94594,77169,94921,77434,95206,77665,95454,77865,95665,78036,95835,78173,95957,78271,96020,78322,96017,78319,95939,78254,95781,78124,95539,77926,95212,77659,94802,77325,94317,76928,93766,76479,93162,75986,92544,75482,93246,76076,93357,76167,93520,76298,93772,76502,94003,76690,94216,76863,94597,77171,94925,77437,95210,77668,95459,77869,95672,78041,95845,78182,95973,78285,96047,78344,96059,78353,96002,78306,95872,78199,95665,78030,95381,77797,95021,77503,94592,77153,94103,76754,93566,76315,93015,75865,93248,76078,93360,76168,93522,76300,93774,76504,94006,76692,94219,76865,94600,77173,94928,77440,95214,77671,95464,77873,95678,78047,95856,78190,95990,78298,96074,78366,96101,78387,96066,78358,95964,78274,95793,78134,95552,77936,95242,77683,94870,77379,94443,77031,93973,76647,93489,76252,93250,76079,93362,76170,93524,76302,93776,76506,94008,76694,94222,76867,94603,77176,94932,77442,95218,77674,95469,77877,95685,78052,95866,78198,96006,78312,96101,78388,96144,78422,96131,78410,96057,78350,95922,78238,95724,78076,95465,77865,95150,77607,94786,77310,94383,76981,93966,76641,93252,76081,93364,76172,93526,76304,93779,76508,94011,76696,94225,76869,94606,77178,94935,77445,95222,77678,95474,77881,95692,78058,95876,78207,96023,78325,96129,78410,96187,78457,96196,78463,96151,78426,96052,78344,95897,78217,95689,78047,95432,77837,95131,77591,94796,77318,94448,77033,93245,76076,93357,76166,93519,76298,93771,76502,94003,76690,94216,76863,94596,77171,94925,77437,95210,77668,95459,77869,95671,78041,95845,78181,95973,78284,96045,78342,96053,78348,95990,78296,95849,78180,95628,77999,95323,77750,94937,77435,94476,77059,93951,76629,93372,76157,92777,75671,93247,76077,93359,76168,93521,76300,93773,76504,94005,76692,94219,76865,94599,77173,94928,77439,95214,77670,95463,77872,95678,78046,95854,78189,95988,78296,96069,78362,96092,78380,96049,78343,95934,78250,95746,78095,95482,77879,95144,77603,94737,77271,94270,76890,93756,76470,93225,76036,93249,76079,93361,76170,93523,76301,93776,76506,94007,76693,94221,76867,94602,77175,94931,77441,95217,77673,95468,77876,95684,78051,95864,78196,96003,78309,96094,78382,96131,78411,96108,78391,96020,78319,95865,78192,95642,78009,95351,77772,94999,77484,94593,77152,94142,76785,93677,76405,93251,76081,93363,76171,93526,76303,93778,76507,94010,76695,94224,76869,94605,77177,94934,77444,95221,77676,95472,77880,95690,78056,95873,78204,96018,78321,96119,78402,96170,78443,96167,78440,96106,78389,95985,78290,95803,78141,95561,77943,95263,77699,94917,77417,94532,77102,94132,76775,93253,76082,93365,76173,93528,76305,93780,76509,94012,76697,94226,76871,94607,77179,94937,77447,95225,77679,95477,77883,95696,78061,95882,78211,96033,78333,96144,78422,96210,78475,96227,78489,96193,78460,96106,78388,95965,78273,95772,78115,95530,77916,95245,77684,94925,77422,94590,77149,93247,76077,93359,76168,93521,76299,93773,76504,94005,76691,94218,76865,94599,77173,94927,77439,95213,77670,95463,77872,95677,78046,95854,78189,95987,78296,96067,78360,96087,78375,96037,78334,95913,78232,95711,78067,95427,77835,95065,77538,94627,77181,94126,76772,93571,76319,92998,75851,93249,76079,93361,76169,93523,76301,93775,76505,94007,76693,94221,76866,94601,77175,94930,77441,95217,77673,95467,77875,95683,78050,95862,78195,96001,78307,96090,78378,96122,78404,96091,78378,95992,78297,95821,78156,95576,77956,95259,77697,94873,77382,94429,77019,93936,76616,93425,76200,93251,76080,93363,76171,93525,76302,93777,76507,94009,76695,94223,76868,94604,77177,94933,77443,95220,77676,95471,77879,95688,78055,95871,78202,96014,78318,96112,78397,96158,78433,96146,78422,96072,78361,95932,78247,95726,78078,95454,77856,95121,77583,94734,77268,94303,76916,93855,76550,93253,76082,93365,76172,93527,76304,93779,76509,94011,76697,94225,76870,94606,77179,94936,77446,95223,77678,95475,77882,95694,78059,95879,78209,96028,78329,96135,78415,96194,78462,96201,78467,96152,78426,96044,78338,95877,78201,95651,78016,95370,77787,95042,77518,94673,77217,94289,76904,93255,76083,93367,76174,93529,76305,93781,76510,94013,76698,94227,76872,94609,77181,94939,77448,95227,77681,95479,77885,95700,78064,95888,78216,96042,78340,96157,78433,96230,78492,96256,78512,96233,78492,96157,78430,96029,78325,95850,78178,95622,77992,95352,77771,95046,77521,94726,77260,93249,76078,93361,76169,93523,76301,93775,76505,94007,76693,94220,76866,94601,77174,94930,77441,95216,77673,95467,77875,95683,78050,95862,78195,96000,78307,96088,78377,96117,78400,96081,78370,95973,78281,95788,78130,95525,77914,95184,77636,94770,77298,94292,76907,93760,76473,93209,76023,93250,76080,93362,76171,93525,76302,93777,76507,94009,76694,94222,76868,94603,77176,94932,77443,95219,77675,95470,77878,95688,78054,95870,78201,96012,78316,96108,78393,96150,78427,96131,78410,96046,78340,95891,78213,95664,78028,95367,77785,95002,77487,94579,77141,94107,76755,93616,76355,93252,76081,93364,76172,93526,76304,93779,76508,94011,76696,94224,76869,94606,77178,94935,77445,95222,77677,95474,77881,95693,78058,95878,78208,96025,78326,96128,78410,96183,78453,96181,78451,96120,78400,95995,78298,95805,78142,95551,77934,95236,77677,94868,77377,94456,77040,94026,76689,93254,76082,93366,76173,93528,76305,93781,76510,94013,76698,94227,76871,94608,77180,94938,77447,95225,77680,95478,77884,95698,78062,95885,78214,96037,78337,96149,78427,96216,78480,96232,78492,96194,78461,96099,78383,95946,78257,95736,78085,95472,77869,95160,77614,94808,77327,94439,77026,93256,76084,93368,76175,93530,76306,93782,76511,94015,76699,94229,76873,94610,77182,94940,77449,95229,77682,95482,77887,95703,78066,95893,78220,96050,78347,96170,78443,96249,78507,96283,78534,96269,78521,96205,78468,96089,78373,95923,78237,95709,78063,95453,77853,95162,77616,94855,77365,93250,76080,93362,76170,93524,76302,93777,76506,94008,76694,94222,76868,94603,77176,94932,77443,95219,77675,95470,77878,95687,78054,95870,78201,96012,78316,96106,78392,96145,78423,96121,78402,96028,78326,95861,78188,95617,77989,95297,77728,94905,77407,94449,77035,93940,76620,93410,76187,93252,76081,93364,76172,93526,76303,93778,76508,94010,76696,94224,76869,94605,77178,94935,77445,95222,77677,95473,77881,95692,78057,95877,78207,96023,78325,96125,78407,96175,78447,96167,78440,96096,78381,95956,78266,95747,78095,95469,77867,95124,77586,94721,77257,94269,76888,93798,76503,93253,76082,93365,76173,93528,76305,93780,76509,94012,76697,94226,76870,94607,77179,94937,77446,95225,77679,95477,77883,95697,78061,95884,78212,96034,78334,96144,78422,96205,78472,96214,78478,96164,78436,96053,78345,95879,78202,95642,78008,95345,77766,94995,77480,94601,77158,94188,76821,93255,76083,93367,76174,93529,76306,93782,76510,94014,76699,94228,76872,94609,77181,94939,77448,95227,77681,95480,77886,95701,78065,95891,78218,96046,78343,96162,78437,96236,78496,96260,78515,96233,78492,96151,78424,96011,78310,95816,78150,95567,77947,95271,77705,94935,77431,94582,77142,93256,76085,93368,76175,93531,76307,93783,76512,94016,76700,94230,76874,94612,77183,94942,77450,95230,77684,95484,77889,95706,78069,95898,78224,96057,78353,96181,78453,96266,78521,96308,78554,96303,78549,96249,78504,96145,78419,95992,78294,95792,78130,95550,77932,95272,77705,94978,77465,93252,76081,93364,76172,93526,76303,93778,76508,94010,76696,94224,76869,94605,77178,94934,77444,95221,77677,95473,77880,95692,78057,95876,78207,96023,78325,96123,78406,96171,78444,96158,78433,96079,78367,95928,78243,95703,78059,95403,77814,95032,77511,94598,77156,94111,76759,93602,76343,93253,76082,93365,76173,93527,76304,93780,76509,94012,76697,94226,76870,94607,77179,94937,77446,95224,77679,95476,77883,95696,78061,95883,78212,96033,78333,96140,78420,96199,78466,96201,78467,96142,78418,96017,78316,95825,78158,95564,77945,95239,77679,94856,77366,94423,77014,93971,76644,93254,76083,93366,76174,93529,76305,93781,76510,94013,76698,94227,76872,94609,77181,94939,77448,95227,77681,95479,77885,95700,78064,95889,78217,96043,78341,96157,78433,96226,78488,96244,78502,96205,78470,96108,78389,95948,78259,95727,78078,95447,77849,95115,77578,94739,77270,94343,76947,93256,76084,93368,76175,93530,76307,93783,76511,94015,76699,94229,76873,94611,77182,94941,77450,95229,77683,95482,77888,95704,78067,95895,78222,96053,78350,96174,78447,96254,78511,96287,78537,96269,78522,96198,78463,96072,78360,95891,78212,95658,78020,95377,77791,95056,77529,94718,77253,93257,76085,93369,76176,93532,76308,93784,76513,94017,76701,94231,76874,94613,77184,94943,77451,95232,77685,95486,77890,95709,78071,95902,78227,96064,78358,96192,78461,96282,78534,96330,78572,96334,78574,96290,78538,96198,78462,96057,78346,95869,78193,95641,78006,95377,77790,95095,77560,93253,76082,93365,76173,93527,76304,93780,76509,94011,76697,94225,76870,94607,77179,94936,77446,95224,77679,95476,77883,95696,78060,95882,78212,96032,78333,96139,78418,96195,78463,96193,78461,96126,78406,95991,78294,95783,78124,95502,77895,95152,77609,94739,77271,94273,76891,93784,76492,93254,76083,93366,76174,93528,76305,93781,76510,94013,76698,94227,76871,94608,77180,94938,77447,95226,77680,95479,77885,95699,78063,95888,78216,96042,78340,96154,78431,96220,78483,96232,78492,96185,78453,96074,78362,95898,78218,95654,78019,95347,77768,94983,77470,94570,77133,94136,76778,93255,76084,93367,76175,93530,76306,93782,76511,94014,76699,94229,76873,94610,77182,94940,77449,95228,77682,95482,77887,95703,78066,95894,78221,96051,78348,96170,78443,96245,78504,96271,78524,96244,78501,96158,78431,96013,78312,95808,78143,95544,77928,95229,77671,94870,77377,94490,77067,93257,76085,93369,76176,93531,76307,93784,76512,94016,76700,94230,76874,94612,77183,94942,77451,95231,77684,95485,77889,95707,78070,95900,78225,96060,78355,96185,78456,96271,78525,96311,78557,96303,78549,96243,78500,96129,78406,95962,78269,95743,78090,95477,77873,95172,77623,94847,77358,93258,76086,93370,76177,93533,76308,93785,76513,94017,76701,94232,76875,94614,77185,94944,77452,95233,77686,95487,77892,95711,78073,95905,78230,96070,78363,96201,78469,96296,78545,96351,78589,96363,78598,96328,78569,96247,78502,96118,78396,95943,78253,95727,78076,95476,77871,95207,77651,93254,76083,93366,76173,93528,76305,93781,76510,94013,76698,94227,76871,94608,77180,94938,77447,95226,77680,95479,77885,95699,78063,95888,78216,96041,78340,96153,78430,96216,78481,96224,78486,96170,78441,96049,78342,95858,78186,95596,77971,95265,77701,94873,77380,94427,77016,93958,76633,93255,76084,93367,76174,93530,76306,93782,76511,94014,76699,94228,76872,94610,77181,94940,77449,95228,77682,95481,77887,95703,78066,95893,78220,96050,78347,96167,78441,96239,78499,96260,78516,96225,78485,96127,78405,95966,78273,95739,78088,95450,77851,95104,77569,94710,77247,94293,76906,93256,76084,93368,76175,93531,76307,93783,76512,94015,76700,94230,76873,94611,77183,94942,77450,95230,77684,95484,77889,95706,78069,95898,78224,96058,78353,96181,78453,96263,78518,96297,78545,96279,78530,96206,78469,96074,78361,95883,78205,95636,78003,95337,77759,94994,77479,94630,77181,93258,76085,93370,76176,93532,76308,93785,76513,94017,76701,94231,76875,94613,77184,94944,77452,95232,77685,95486,77891,95710,78072,95903,78228,96067,78360,96195,78464,96286,78537,96334,78575,96334,78575,96285,78533,96183,78450,96029,78323,95823,78155,95572,77950,95281,77712,94970,77458,93259,76086,93371,76177,93533,76309,93786,76514,94018,76702,94233,76876,94615,77185,94946,77453,95235,77687,95489,77893,95713,78074,95909,78233,96075,78367,96210,78476,96310,78556,96371,78605,96390,78620,96364,78598,96293,78539,96175,78442,96012,78309,95809,78143,95570,77948,95313,77738,93255,76083,93367,76174,93529,76306,93782,76511,94014,76699,94228,76872,94610,77181,94940,77449,95228,77682,95481,77887,95702,78066,95893,78220,96049,78346,96166,78440,96236,78497,96254,78510,96211,78475,96104,78387,95929,78243,95684,78043,95372,77788,94999,77483,94574,77136,94123,76768,93256,76084,93368,76175,93531,76307,93783,76512,94015,76700,94229,76873,94611,77182,94941,77450,95230,77683,95483,77889,95706,78068,95898,78224,96057,78352,96179,78451,96257,78514,96287,78537,96261,78515,96177,78446,96030,78325,95819,78153,95547,77930,95219,77662,94842,77355,94442,77028,93257,76085,93369,76176,93532,76308,93784,76513,94016,76701,94231,76874,94613,77184,94943,77451,95232,77685,95486,77890,95709,78071,95902,78227,96065,78359,96192,78461,96279,78531,96320,78564,96312,78557,96250,78505,96131,78407,95954,78263,95722,78073,95439,77842,95113,77575,94764,77290,93258,76086,93370,76177,93533,76309,93785,76513,94018,76702,94232,76875,94614,77185,94945,77453,95234,77686,95488,77892,95712,78073,95907,78231,96072,78365,96205,78472,96300,78548,96354,78592,96363,78598,96323,78565,96233,78491,96091,78374,95899,78217,95662,78023,95385,77797,95088,77554,93259,76087,93371,76178,93534,76310,93787,76514,94019,76703,94233,76876,94616,77186,94947,77454,95236,77688,95490,77894,95715,78076,95912,78235,96080,78371,96218,78482,96322,78566,96388,78619,96415,78640,96398,78625,96336,78574,96229,78486,96078,78362,95886,78206,95660,78021,95414,77820,93256,76084,93368,76175,93530,76307,93783,76511,94015,76700,94229,76873,94611,77182,94941,77450,95230,77683,95483,77888,95705,78068,95897,78224,96057,78352,96178,78450,96255,78512,96280,78532,96249,78505,96155,78428,95995,78297,95767,78111,95474,77870,95119,77581,94713,77249,94281,76896,93257,76085,93369,76176,93531,76308,93784,76512,94016,76701,94230,76874,94612,77183,94943,77451,95231,77685,95485,77890,95708,78071,95902,78227,96064,78358,96189,78459,96274,78527,96311,78557,96296,78543,96223,78483,96089,78374,95894,78214,95638,78004,95327,77750,94968,77457,94585,77144,93258,76086,93370,76177,93532,76308,93785,76513,94017,76701,94232,76875,94614,77184,94944,77452,95233,77686,95487,77892,95711,78073,95906,78230,96070,78363,96201,78469,96293,78543,96342,78582,96343,78581,96291,78539,96184,78451,96021,78317,95804,78139,95536,77921,95225,77666,94891,77394,93259,76087,93371,76177,93534,76309,93786,76514,94019,76702,94233,76876,94615,77186,94946,77453,95235,77687,95490,77894,95714,78075,95910,78234,96078,78369,96213,78478,96313,78559,96373,78607,96390,78620,96360,78594,96280,78529,96150,78422,95971,78276,95747,78092,95484,77877,95200,77645,93260,76087,93372,76178,93535,76310,93787,76515,94020,76703,94234,76877,94616,77187,94947,77455,95237,77689,95492,77895,95717,78077,95915,78237,96085,78375,96225,78488,96333,78575,96405,78633,96438,78659,96429,78651,96376,78607,96280,78528,96139,78413,95959,78265,95745,78090,95511,77898,93257,76085,93369,76176,93531,76307,93784,76512,94016,76700,94230,76874,94612,77183,94943,77451,95231,77684,95485,77890,95708,78070,95901,78227,96063,78357,96188,78459,96271,78525,96305,78552,96284,78534,96203,78467,96057,78347,95845,78174,95569,77948,95233,77674,94845,77357,94431,77019,93258,76086,93370,76176,93532,76308,93785,76513,94017,76701,94231,76875,94613,77184,94944,77452,95233,77686,95487,77891,95711,78072,95905,78230,96069,78363,96199,78467,96289,78539,96334,78575,96327,78569,96266,78518,96145,78419,95964,78271,95724,78075,95430,77834,95087,77554,94720,77255,93259,76086,93371,76177,93533,76309,93786,76514,94018,76702,94232,76876,94615,77185,94945,77453,95234,77687,95489,77893,95713,78075,95909,78233,96076,78368,96210,78476,96307,78554,96362,78598,96371,78605,96329,78570,96234,78492,96084,78369,95881,78202,95628,77995,95332,77753,95012,77492,93260,76087,93372,76178,93534,76310,93787,76515,94019,76703,94234,76877,94616,77186,94947,77454,95236,77688,95491,77895,95716,78077,95913,78236,96082,78373,96221,78485,96325,78569,96391,78622,96415,78640,96393,78622,96324,78565,96205,78467,96039,78331,95828,78158,95578,77954,95306,77732,93261,76088,93373,76179,93535,76310,93788,76515,94020,76704,94235,76877,94617,77187,94948,77455,95238,77690,95493,77896,95719,78079,95917,78239,96089,78378,96232,78493,96343,78583,96420,78645,96459,78676,96458,78674,96414,78638,96327,78566,96198,78460,96029,78322,95826,78156,95602,77973,93258,76086,93370,76176,93532,76308,93785,76513,94017,76701,94231,76875,94613,77184,94944,77452,95233,77686,95487,77891,95711,78072,95905,78230,96069,78362,96198,78466,96287,78538,96328,78571,96317,78560,96247,78503,96115,78395,95919,78234,95659,78022,95341,77761,94971,77459,94574,77135,93258,76086,93371,76177,93533,76309,93786,76514,94018,76702,94232,76876,94614,77185,94945,77453,95234,77687,95489,77893,95713,78074,95909,78233,96075,78367,96208,78474,96303,78551,96354,78592,96357,78593,96306,78551,96198,78462,96031,78325,95806,78141,95527,77913,95201,77647,94850,77360,93259,76087,93371,76178,93534,76309,93787,76514,94019,76703,94233,76876,94616,77186,94946,77454,95236,77688,95490,77894,95715,78076,95912,78235,96081,78372,96218,78482,96319,78564,96381,78613,96397,78626,96365,78599,96281,78530,96144,78417,95954,78261,95715,78066,95433,77836,95128,77586,93260,76087,93372,76178,93535,76310,93788,76515,94020,76703,94234,76877,94617,77187,94948,77455,95237,77689,95492,77896,95718,78078,95916,78238,96087,78376,96228,78490,96336,78577,96407,78635,96438,78659,96425,78647,96365,78598,96257,78510,96103,78383,95904,78220,95667,78026,95407,77814,93261,76088,93373,76179,93536,76311,93788,76516,94021,76704,94235,76878,94618,77188,94949,77456,95239,77690,95494,77897,95720,78080,95919,78241,96092,78381,96238,78498,96352,78591,96434,78657,96479,78692,96485,78696,96450,78667,96372,78603,96253,78505,96095,78376,95902,78218,95690,78044,93258,76086,93370,76177,93533,76309,93786,76514,94018,76702,94232,76875,94614,77185,94945,77453,95234,77687,95488,77893,95713,78074,95909,78232,96075,78367,96207,78474,96301,78549,96349,78588,96347,78585,96289,78537,96170,78439,95988,78290,95744,78091,95443,77845,95090,77557,94710,77246,93259,76087,93371,76178,93534,76309,93786,76514,94019,76703,94233,76876,94615,77186,94946,77454,95235,77688,95490,77894,95715,78076,95912,78235,96080,78371,96216,78481,96316,78561,96373,78607,96384,78615,96343,78581,96247,78502,96093,78376,95883,78204,95619,77988,95309,77735,94973,77460,93260,76087,93372,76178,93535,76310,93787,76515,94020,76703,94234,76877,94616,77187,94947,77455,95237,77689,95492,77895,95717,78078,95915,78238,96085,78375,96225,78488,96331,78573,96398,78627,96422,78646,96399,78626,96325,78565,96199,78462,96022,78317,95797,78133,95529,77914,95237,77676,93261,76088,93373,76179,93535,76311,93788,76516,94021,76704,94235,76878,94617,77187,94949,77456,95238,77690,95493,77897,95719,78079,95918,78240,96090,78379,96234,78495,96346,78585,96422,78647,96460,78676,96454,78671,96404,78629,96306,78549,96163,78432,95976,78279,95751,78095,95504,77893,93261,76088,93374,76179,93536,76311,93789,76516,94021,76705,94236,76878,94618,77188,94950,77457,95239,77691,95495,77898,95721,78081,95921,78243,96096,78384,96243,78503,96361,78598,96447,78667,96498,78707,96510,78717,96483,78694,96414,78637,96305,78547,96157,78426,95975,78278,95773,78112,93259,76087,93371,76178,93534,76309,93786,76514,94019,76702,94233,76876,94615,77186,94946,77454,95235,77688,95490,77894,95715,78076,95912,78235,96080,78371,96215,78480,96314,78559,96369,78604,96375,78608,96327,78568,96221,78480,96053,78343,95825,78156,95540,77923,95204,77649,94840,77352,93260,76087,93372,76178,93534,76310,93787,76515,94019,76703,94234,76877,94616,77186,94947,77454,95236,77689,95491,77895,95717,78077,95914,78237,96084,78375,96223,78487,96327,78570,96391,78622,96410,78636,96378,78609,96293,78539,96152,78424,95955,78263,95707,78059,95412,77818,95090,77556,93260,76088,93373,76179,93535,76310,93788,76515,94020,76704,94235,76877,94617,77187,94948,77455,95238,77690,95493,77896,95719,78079,95917,78240,96089,78378,96231,78493,96341,78582,96414,78640,96444,78664,96430,78651,96366,78599,96252,78505,96087,78370,95875,78197,95621,77989,95342,77761,93261,76088,93373,76179,93536,76311,93789,76516,94021,76704,94236,76878,94618,77188,94949,77456,95239,77691,95495,77897,95721,78080,95920,78242,96094,78382,96240,78500,96355,78593,96436,78658,96479,78693,96482,78694,96440,78659,96353,78587,96220,78478,96045,78335,95832,78161,95596,77968,93262,76089,93374,76180,93537,76312,93789,76517,94022,76705,94236,76879,94619,77189,94950,77457,95240,77692,95496,77899,95723,78082,95923,78244,96099,78386,96248,78507,96369,78604,96459,78677,96515,78721,96534,78736,96514,78719,96454,78670,96354,78587,96216,78474,96045,78334,95852,78177,93260,76087,93372,76178,93534,76310,93787,76515,94019,76703,94234,76877,94616,77186,94947,77454,95236,77689,95491,77895,95717,78077,95914,78237,96084,78374,96223,78486,96325,78569,96387,78618,96401,78629,96363,78597,96268,78519,96114,78393,95901,78218,95631,77998,95312,77737,94964,77453,93260,76088,93372,76179,93535,76310,93788,76515,94020,76704,94235,76877,94617,77187,94948,77455,95237,77689,95493,77896,95718,78079,95917,78239,96088,78378,96230,78492,96338,78579,96407,78635,96433,78655,96411,78636,96336,78574,96207,78469,96024,78319,95789,78127,95509,77898,95202,77647,93261,76088,93373,76179,93536,76311,93788,76516,94021,76704,94235,76878,94618,77188,94949,77456,95239,77690,95494,77897,95720,78080,95920,78241,96093,78381,96237,78498,96351,78589,96428,78652,96465,78681,96459,78675,96405,78630,96301,78545,96148,78420,95949,78257,95708,78060,95442,77842,93262,76089,93374,76180,93536,76311,93789,76516,94022,76705,94236,76879,94619,77188,94950,77457,95240,77691,95495,77898,95722,78081,95922,78243,96097,78385,96245,78504,96363,78600,96449,78669,96498,78708,96507,78714,96474,78686,96396,78622,96274,78522,96110,78388,95908,78223,95684,78039,93262,76089,93374,76180,93537,76312,93790,76517,94022,76705,94237,76879,94620,77189,94951,77458,95241,77692,95497,77899,95724,78083,95925,78246,96102,78388,96252,78510,96376,78610,96470,78686,96531,78734,96556,78754,96543,78743,96491,78700,96400,78625,96272,78520,96110,78388,95927,78238,93260,76088,93372,76178,93535,76310,93788,76515,94020,76704,94234,76877,94617,77187,94948,77455,95237,77689,95493,77896,95718,78078,95917,78239,96088,78378,96229,78492,96336,78578,96403,78632,96425,78649,96397,78624,96313,78556,96172,78440,95972,78276,95718,78069,95414,77820,95081,77548,93261,76088,93373,76179,93535,76311,93788,76516,94021,76704,94235,76878,94618,77188,94949,77456,95238,77690,95494,77897,95720,78080,95919,78241,96092,78381,96236,78497,96348,78587,96422,78647,96455,78673,96441,78661,96377,78607,96259,78511,96089,78371,95868,78191,95602,77973,95308,77733,93261,76088,93374,76179,93536,76311,93789,76516,94021,76705,94236,76878,94618,77188,94950,77457,95239,77691,95495,77898,95721,78081,95922,78243,96096,78384,96243,78503,96359,78596,96441,78662,96485,78697,96486,78697,96441,78660,96348,78583,96206,78467,96019,78314,95791,78127,95537,77920,93262,76089,93374,76180,93537,76312,93790,76517,94022,76705,94237,76879,94619,77189,94951,77457,95240,77692,95496,77899,95723,78082,95924,78245,96100,78387,96250,78508,96371,78606,96461,78678,96515,78722,96531,78734,96505,78712,96437,78655,96325,78563,96171,78438,95981,78282,95767,78107,93263,76089,93375,76180,93537,76312,93790,76517,94023,76706,94237,76880,94620,77190,94952,77458,95242,77693,95498,77900,95725,78084,95927,78247,96104,78390,96257,78514,96383,78615,96480,78694,96545,78746,96576,78771,96571,78765,96526,78728,96444,78661,96325,78563,96173,78438,95999,78296,93261,76088,93373,76179,93535,76311,93788,76516,94021,76704,94235,76878,94618,77188,94949,77456,95238,77690,95494,77897,95720,78080,95919,78241,96092,78381,96236,78497,96346,78586,96419,78644,96448,78667,96428,78650,96355,78590,96226,78484,96040,78332,95800,78135,95511,77899,95193,77640,93261,76088,93373,76179,93536,76311,93789,76516,94021,76704,94236,76878,94618,77188,94950,77456,95239,77691,95495,77898,95721,78081,95921,78243,96095,78383,96242,78502,96357,78594,96436,78658,96475,78689,96469,78684,96415,78638,96308,78551,96150,78421,95942,78251,95689,78045,95409,77816,93262,76089,93374,76180,93537,76312,93789,76517,94022,76705,94236,76879,94619,77189,94950,77457,95240,77692,95496,77899,95723,78082,95923,78244,96099,78386,96248,78507,96367,78603,96454,78673,96503,78712,96511,78718,96475,78687,96391,78618,96261,78511,96085,78368,95869,78191,95627,77993,93262,76089,93374,76180,93537,76312,93790,76517,94022,76705,94237,76879,94620,77189,94951,77458,95241,77692,95497,77900,95724,78083,95926,78246,96103,78389,96254,78512,96378,78612,96472,78687,96531,78734,96553,78752,96535,78736,96475,78687,96373,78602,96230,78485,96050,78338,95847,78172,93263,76090,93375,76181,93538,76312,93791,76517,94023,76706,94238,76880,94620,77190,94952,77458,95242,77693,95498,77901,95726,78084,95928,78248,96106,78392,96260,78517,96389,78620,96489,78701,96559,78757,96595,78786,96596,78786,96559,78755,96485,78694,96375,78604,96232,78487,96068,78352,93261,76088,93373,76179,93536,76311,93789,76516,94021,76704,94236,76878,94618,77188,94949,77456,95239,77691,95495,77898,95721,78081,95921,78243,96095,78383,96241,78501,96355,78593,96433,78656,96468,78684,96457,78674,96394,78622,96277,78525,96104,78384,95878,78199,95604,77975,95300,77727,93262,76089,93374,76180,93536,76311,93789,76516,94022,76705,94236,76879,94619,77189,94950,77457,95240,77691,95496,77898,95722,78082,95923,78244,96098,78386,96247,78506,96365,78601,96449,78669,96494,78704,96496,78705,96450,78667,96354,78588,96208,78468,96013,78308,95773,78113,95506,77894,93262,76089,93374,76180,93537,76312,93790,76517,94022,76705,94237,76879,94620,77189,94951,77458,95241,77692,95497,77899,95724,78083,95925,78246,96102,78389,96252,78510,96375,78609,96465,78682,96520,78725,96535,78737,96506,78713,96432,78652,96312,78553,96148,78419,95944,78252,95713,78063,93263,76089,93375,76180,93537,76312,93790,76517,94023,76706,94237,76880,94620,77190,94952,77458,95242,77693,95498,77900,95725,78084,95927,78247,96105,78391,96258,78515,96385,78617,96482,78695,96546,78746,96574,78769,96563,78759,96511,78716,96418,78639,96285,78530,96115,78392,95922,78234,93263,76090,93375,76181,93538,76313,93791,76518,94023,76706,94238,76880,94621,77190,94953,77459,95243,77694,95499,77901,95727,78085,95929,78249,96108,78394,96264,78519,96394,78625,96498,78708,96572,78768,96613,78801,96620,78805,96591,78781,96524,78726,96423,78642,96289,78533,96133,78405],"diffuse":[0,0,0,0,0,0,0,0,0,0,8414,6288,9180,6860,9947,7433,10714,8006,11483,8581,16110,12040,17445,13037,18782,14037,20122,15038,21464,16042,23147,17300,24892,18605,26642,19913,28396,21225,30156,22541,29578,22108,31607,23626,33644,25149,35687,26678,37738,28212,35453,26503,37668,28159,39891,29822,42123,31492,44364,33169,40820,30517,43141,32254,45472,33998,47814,35750,50166,37510,45721,34184,48087,35955,50465,37734,52854,39522,55255,41320,50196,37532,52560,39302,54937,41081,57327,42870,59729,44669,54280,40589,56608,42331,58948,44084,61302,45846,63668,47619,58007,43379,60271,45074,62549,46780,64840,48496,67144,50222,61408,45924,63590,47559,65785,49203,67994,50858,70217,52523,64510,48247,66597,49811,68697,51384,70811,52968,72938,54561,67340,50366,69323,51853,71320,53349,73329,54854,75351,56370,69920,52299,71796,53705,73683,55119,75583,56543,77496,57977,72274,54062,74039,55385,75815,56717,77604,58058,79405,59408,74420,55670,76074,56910,77740,58159,79418,59417,81107,60683,76376,57136,77922,58295,79479,59462,81047,60638,82626,61822,78159,58472,79600,59553,81051,60641,82513,61737,83985,62841,79785,59691,81124,60695,82473,61707,83832,62726,85201,63753,81267,60801,82509,61733,83761,62672,85021,63618,86291,64570,82617,61814,83767,62677,84926,63546,86093,64422,87269,65304,83848,62737,84911,63534,85982,64338,87061,65148,88148,65963,84970,63578,85951,64314,86940,65056,87936,65803,88939,66556,85992,64344,86896,65023,87808,65707,88726,66396,89651,67090,86923,65043,87756,65668,88595,66298,89441,66932,90293,67572,87772,65680,88538,66255,89310,66834,90088,67418,90871,68006,88546,66260,89249,66788,89959,67321,90673,67857,91394,68398,89250,66789,89897,67274,90548,67763,91204,68256,91865,68752,89892,67270,90485,67716,91083,68164,91685,68617,92292,69072,90478,67710,91021,68118,91569,68529,92121,68944,92678,69362,91011,68110,91509,68484,92011,68861,92517,69241,93027,69624,91496,68474,91953,68817,92413,69162,92876,69511,93343,69861,91939,68806,92357,69120,92778,69437,93202,69755,93630,70077,92342,69109,92725,69396,93110,69686,93498,69978,93890,70272,92710,69385,93059,69648,93412,69913,93767,70179,94125,70449,93045,69636,93364,69877,93687,70119,94011,70363,94339,70609,93350,69866,93642,70085,93937,70307,94233,70530,94533,70755,93628,70074,93895,70275,94164,70477,94435,70681,94709,70887,93881,70264,94125,70448,94371,70633,94619,70819,94869,71007,94112,70438,94335,70605,94559,70774,94786,70945,95014,71116,94322,70596,94525,70749,94731,70903,94938,71059,95146,71215,94513,70740,94699,70879,94887,71020,95076,71162,95266,71306,94688,70871,94858,70998,95029,71127,95201,71257,95375,71388,94847,70990,95002,71107,95158,71224,95316,71343,95475,71462,94992,71099,95134,71206,95276,71313,95420,71421,95565,71530,95124,71199,95253,71296,95383,71394,95515,71492,95647,71592,95245,71289,95362,71378,95481,71467,95601,71557,95722,71648,95354,71372,95462,71452,95570,71534,95679,71616,95790,71699,95454,71447,95552,71520,95651,71595,95751,71670,95852,71745,95545,71515,95635,71582,95725,71650,95816,71719,95908,71788,95628,71578,95710,71639,95792,71701,95875,71763,95959,71826,95704,71634,95779,71691,95854,71747,95929,71804,96006,71861,95773,71686,95841,71737,95910,71789,95979,71841,96048,71893,95836,71734,95898,71780,95960,71827,96023,71875,96087,71923,95893,71777,95950,71819,96007,71862,96064,71905,96122,71949,95946,71816,95997,71855,96049,71894,96102,71933,96154,71973,95993,71852,96040,71887,96088,71923,96135,71959,96184,71995,96037,71884,96079,71917,96123,71949,96166,71982,96210,72015,96076,71914,96115,71944,96155,71973,96194,72003,96235,72033,96112,71941,96148,71968,96184,71995,96220,72023,96257,72050,96145,71966,96178,71991,96210,72015,96244,72040,96277,72065,96175,71989,96205,72011,96235,72033,96265,72056,96295,72079,96202,72009,96229,72029,96257,72050,96284,72071,96312,72092,96227,72028,96252,72046,96277,72065,96302,72084,96327,72103,96250,72045,96272,72062,96295,72079,96318,72096,96341,72114,96271,72061,96291,72076,96312,72092,96333,72107,96354,72123,96290,72075,96308,72089,96327,72103,96346,72117,96366,72132,96307,72088,96324,72101,96341,72114,96358,72127,96376,72140,96322,72099,96338,72111,96354,72123,96370,72135,96386,72147,96337,72110,96351,72121,96365,72132,96380,72143,96394,72154,96350,72120,96363,72130,96376,72140,96389,72150,96402,72160,96362,72129,96374,72138,96385,72147,96398,72156,96410,72165,96373,72137,96383,72145,96394,72153,96405,72162,96416,72170,96382,72145,96392,72152,96402,72160,96412,72167,96422,72175,96391,72151,96400,72158,96410,72165,96419,72172,96428,72179,96400,72158,96408,72164,96416,72170,96425,72176,96433,72183,96407,72163,96415,72169,96422,72175,96430,72180,96438,72186,96414,72168,96421,72174,96428,72179,96435,72184,96442,72189,96420,72173,96427,72178,96433,72183,96439,72187,96446,72192,96426,72177,96432,72182,96437,72186,96443,72190,96449,72195]};
function rtAbsorptionMethod(method) {
  return method === 'RT_CHL' || method === 'RT_LEAF' || rtAbsIndex(method)!==null;
}
// Collapse the library to the fitted LAI histogram and red background once per date.
function rtAbsorptionResponse(histogram, background) {
  var mu = RT_PAR.mu, size = mu.length;
  var values = [];
  for (var i=0; i<2*size+2; i++) values.push(0);
  var brightness = Math.max(0, Math.min(4, 10*background));
  var lo = Math.floor(brightness), hi = Math.min(4, lo+1), fraction = brightness-lo;
  var count=histogram.reduce(function(a,b){return a+b;},0);
  histogram.forEach(function(n,lai) {
    if (!n) return;
    var weight=n/(count*100000);
    for (var j=0;j<values.length;j++) {
      var diffuse=j>=2*size;
      var stride=diffuse?2:2*size, band=diffuse?j-2*size:j;
      var table=diffuse?RT_PAR.diffuse:RT_PAR.direct;
      var low=table[(lai*5+lo)*stride+band], high=table[(lai*5+hi)*stride+band];
      values[j]+=weight*(low+(high-low)*fraction);
    }
  });
  return values;
}
// Integrate with midpoint quadrature over the sunlit half-day (symmetric about noon).
function rtSunCosines(latitude, millis) {
  var day=1+(millis-Date.UTC(new Date(millis).getUTCFullYear(),0,1))/86400000;
  var dec=.409*Math.sin(2*Math.PI*day/365-1.39), lat=latitude*Math.PI/180;
  var a=Math.sin(lat)*Math.sin(dec), b=Math.cos(lat)*Math.cos(dec);
  var sunset=Math.acos(Math.max(-1,Math.min(1,-a/Math.max(1e-12,b))));
  if (!sunset) return null;
  var points=[];
  for(var i=0;i<16;i++) points.push(Math.max(0,a+b*Math.cos(sunset*(i+.5)/16)));
  return points;
}
function rtDayAbsorption(response, latitude, millis) {
  var points=rtSunCosines(latitude,millis);
  if(!points)return null;
  var grid=RT_PAR.mu, result=[0,0,response[2*grid.length],response[2*grid.length+1]], sum=0;
  for(var i=0;i<points.length;i++) {
    var mu=points[i],h=1;
    while(h<grid.length-1 && grid[h]<mu)h++;
    var w=Math.max(0,Math.min(1,(mu-grid[h-1])/(grid[h]-grid[h-1])));
    for(var c=0;c<2;c++) result[c]+=mu*((1-w)*response[2*(h-1)+c]+w*response[2*h+c]);
    sum+=mu;
  }
  if(!sum)return null;
  result[0]/=sum; result[1]/=sum;
  return result;
}
// Spherical, randomly distributed foliage: a ray at zenith cosine mu sees LAI/mu.
function rtDayInterception(lai,points,sky) {
  if(!points)return null;
  if(lai===0)return 0;
  function intercepted(mu){return 1-Math.exp(-.5*lai/mu);}
  var beam=0, light=0, diffuse=0;
  for(var i=0;i<points.length;i++){
    beam+=points[i]*intercepted(points[i]);light+=points[i];
  }
  // Isotropic sky incident on a horizontal surface: integral of 2*mu*f(mu) dmu.
  for(var j=0;j<64;j++){
    var mu=(j+.5)/64;
    diffuse+=2*mu*intercepted(mu)/64;
  }
  return light ? (1-sky)*beam/light+sky*diffuse : null;
}
// Equal-area pixels share illumination. Weight each LAI response by its pixel count.
function rtMeanInterception(histogram,points,sky) {
  var total=0, pixels=0;
  histogram.forEach(function(count,index){
    if(count){
      total+=count*rtDayInterception(index/10,points,sky);
      pixels+=count;
    }
  });
  return total/pixels;
}
// Interpolate responses at common illumination before applying each day's solar geometry.
function rtAbsorptionPeriods(obs, periods, latitude, sky, radiation) {
  var sorted=obs.slice().sort(function(a,b){return a.millis-b.millis;});
  var curves=sorted.map(function(d){return rtAbsorptionResponse(d.lai_histogram,d.background_red);});
  var j=0, output={};
  periods.forEach(function(period){
    var totals=[0,0,0,0], energy=0, intercepted=0, observed=0, filled=0;
    for(var t=period.millis;t<period.period_end_millis;t+=86400000){
      while(j<sorted.length && sorted[j].millis<t)j++;
      var exact=j<sorted.length && sorted[j].millis===t;
      if(!exact && (!j || j===sorted.length || sorted[j].millis-sorted[j-1].millis>32*86400000))continue;
      var curve=curves[j];
      if(!exact){
        var w=(t-sorted[j-1].millis)/(sorted[j].millis-sorted[j-1].millis);
        curve=curves[j-1].map(function(v,k){return v+w*(curves[j][k]-v);});
      }
      var q=radiation[t], day=rtDayAbsorption(curve,latitude,t);
      if(typeof q!=='number'|| !isFinite(q)||q<0)continue;
      if(exact)observed++;else filled++;
      if(!day)continue;
      for(var k=0;k<4;k++)totals[k]+=q*day[k];
      var points=rtSunCosines(latitude,t);
      var value=rtMeanInterception(sorted[j].lai_histogram,points,sky);
      if(!exact)value=w*value+(1-w)*rtMeanInterception(sorted[j-1].lai_histogram,points,sky);
      intercepted+=q*value;
      energy+=q;
    }
    var row={rt_abs_observed_days:observed,rt_abs_interpolated_days:filled,
      rt_abs_missing_pct:100*(1-(observed+filled)/period.period_days),rt_diffuse_par_fraction:sky};
    if(energy>0 && (period.model_cadence==='monthly'||observed+filled===period.period_days)){
      var v=totals.map(function(n){return n/energy;});
      row.rt_fipar_daylight=intercepted/energy;
      row.rt_leaf_direct=v[0];row.rt_chl_direct=v[1];row.rt_leaf_diffuse=v[2];row.rt_chl_diffuse=v[3];
      row.FAPAR_leaf=(1-sky)*v[0]+sky*v[2];row.FAPAR_chl=(1-sky)*v[1]+sky*v[3];
    }
    output[period.date]=row;
  });
  return output;
}

function rtAbsIndex(method) {
  if(method==='RT_CHL_NDVI')return 'ndvi';
  if(method==='RT_CHL_EVI2')return 'evi2';
  return null;
}
// Pixels sharing LAI have identical modeled absorption on a given date.
// Store sufficient statistics so changing the diffuse fraction needs no new retrieval.
function rtIndexBins(pixels,covers) {
  var bins={};
  pixels.forEach(function(p,i){
    if(p.red+p.nir<=0)return;
    var lai=Math.round(-20*Math.log(1-covers[i]));
    var bin=bins[lai] || (bins[lai]={n:0,ndvi:0,ndvi_sq:0,evi2:0,evi2_sq:0});
    var ndvi=(p.nir-p.red)/(p.nir+p.red);
    var evi2=Math.max(0,Math.min(1,2.5*(p.nir-p.red)/(1+p.nir+2.4*p.red)));
    bin.n++;bin.ndvi+=ndvi;bin.ndvi_sq+=ndvi*ndvi;bin.evi2+=evi2;bin.evi2_sq+=evi2*evi2;
  });
  return bins;
}
function rtFitDaylightLine(obs,first,last,sky,index,interception) {
  var moments=[];
  obs.forEach(function(day){
    if(day.millis<first || day.millis>last)return;
    var points=rtSunCosines(day.latitude,day.millis);
    if(!points)return;
    var sx=0,sy=0,sxx=0,sxy=0,n=0;
    Object.keys(day.index_bins).forEach(function(lai){
      var bin=day.index_bins[lai], y;
      if(interception)y=rtDayInterception(Number(lai)/10,points,sky);
      else {
        var hist=[];
        hist[Number(lai)]=1;
        var response=rtAbsorptionResponse(hist,day.background_red);
        var terms=rtDayAbsorption(response,day.latitude,day.millis);
        y=(1-sky)*terms[1]+sky*terms[3];
      }
      n+=bin.n;sx+=bin[index];sxx+=bin[index+'_sq'];sy+=bin.n*y;sxy+=bin[index]*y;
    });
    if(n)moments.push({x:sx/n,y:sy/n,xx:sxx/n,xy:sxy/n,n:n});
  });
  return rtEqualDateLine(moments);
}

// Client-side deterministic fitter. No Earth Engine objects in this module.
var RT = (function () {
  var PI = Math.PI, RAD = PI / 180;
  var soils = [];
  for (var r = 1; r <= 40; r++) {
    for (var n = r; n <= 50; n++) soils.push([r, n]);
  }
  function median(values) {
    var a = values.slice().sort(function (x, y) { return x - y; });
    var m = Math.floor(a.length / 2);
    return a.length % 2 ? a[m] : (a[m - 1] + a[m]) / 2;
  }
  function geometry(sun, view, azimuth) {
    var cs = Math.cos(sun * RAD), co = Math.cos(view * RAD);
    var ss = Math.sin(sun * RAD), so = Math.sin(view * RAD);
    var psi = azimuth * RAD;
    var ks = 0, ko = 0, bf = 0, sob = 0, sof = 0;
    for (var i = 0; i < 18; i++) {
      var leaf = (i * 5 + 2.5) * RAD;
      var cl = Math.cos(leaf), sl = Math.sin(leaf);
      var a = cs * cl, b = co * cl, c = ss * sl, d = so * sl;
      var bs = Math.abs(c) > 1e-6 ? -a / c : 5;
      var bo = Math.abs(d) > 1e-6 ? -b / d : 5;
      var ts = Math.abs(bs) < 1 ? Math.acos(bs) : PI;
      var to = Math.abs(bo) < 1 ? Math.acos(bo) : PI;
      var ds = Math.abs(bs) < 1 ? c : a;
      var dv = Math.abs(bo) < 1 ? d : b;
      var chiS = 2 / PI * ((ts - PI / 2) * a + Math.sin(ts) * c);
      var chiO = 2 / PI * ((to - PI / 2) * b + Math.sin(to) * d);
      var turns = [psi, Math.abs(ts - to), PI - Math.abs(ts + to - PI)];
      turns.sort(function (x, y) { return x - y; });
      var t1 = 2 * a * b + c * d * Math.cos(psi);
      var t2 = turns[1] > 0 ?
        Math.sin(turns[1]) * (2 * ds * dv + c * d * Math.cos(turns[0]) * Math.cos(turns[2])) : 0;
      var fr = Math.max(0, ((PI - turns[1]) * t1 + t2) / (2 * PI * PI));
      var ft = Math.max(0, (-turns[1] * t1 + t2) / (2 * PI * PI));
      var weight = RT_LEAF.lidf[i];
      ks += weight * chiS / cs;
      ko += weight * chiO / co;
      bf += weight * cl * cl;
      sob += weight * fr * PI / (cs * co);
      sof += weight * ft * PI / (cs * co);
    }
    var tanS = Math.tan(sun * RAD), tanO = Math.tan(view * RAD);
    var distance = Math.sqrt(Math.max(0, tanS * tanS + tanO * tanO - 2 * tanS * tanO * Math.cos(psi)));
    return {ks: ks, ko: ko, bf: bf, sob: sob, sof: sof, distance: distance};
  }
  function j1(k, m, lai) {
    var delta = (k - m) * lai;
    if (Math.abs(delta) > .001) return (Math.exp(-m * lai) - Math.exp(-k * lai)) / (k - m);
    return .5 * lai * (Math.exp(-k * lai) + Math.exp(-m * lai)) * (1 - delta * delta / 12);
  }
  function j2(k, m, lai) { return (1 - Math.exp(-(k + m) * lai)) / (k + m); }
  function hotspot(g, lai) {
    var alpha = g.distance / .05 * 2 / (g.ks + g.ko);
    var direct = Math.exp(-g.ks * lai);
    if (alpha === 0) return [direct, (1 - direct) / (g.ks * lai)];
    var hot = lai * Math.sqrt(g.ko * g.ks);
    var interval = (1 - Math.exp(-alpha)) * .05;
    var x1 = 0, y1 = 0, f1 = 1, integral = 0;
    for (var i = 1; i <= 20; i++) {
      var x2 = i < 20 ? -Math.log(1 - i * interval) / alpha : 1;
      var y2 = -(g.ko + g.ks) * lai * x2 + hot * (1 - Math.exp(-alpha * x2)) / alpha;
      var f2 = Math.exp(y2);
      integral += y2 === y1 ? f1 * (x2 - x1) : (f2 - f1) * (x2 - x1) / (y2 - y1);
      x1 = x2; y1 = y2; f1 = f2;
    }
    return [f1, integral];
  }
  // Compute both bands for all 50 background brightness levels, once per LAI/geometry.
  function spectrum(g, lai) {
    var result = new Float64Array(100);
    if (lai === 0) {
      for (var z = 0; z < 100; z++) result[z] = (z % 50 + 1) / 100;
      return result;
    }
    var ks = g.ks, ko = g.ko, bf = g.bf;
    var tss = Math.exp(-ks * lai), too = Math.exp(-ko * lai);
    var hot = hotspot(g, lai);
    var zint = j2(ks, ko, lai);
    for (var w = 0; w < 87; w++) {
      var rho = RT_LEAF.rho[w], tau = RT_LEAF.tau[w];
      var sb = .5 * ((ks + bf) * rho + (ks - bf) * tau);
      var sf = .5 * ((ks - bf) * rho + (ks + bf) * tau);
      var vb = .5 * ((ko + bf) * rho + (ko - bf) * tau);
      var vf = .5 * ((ko - bf) * rho + (ko + bf) * tau);
      var sigb = Math.max(1e-36, .5 * ((1 + bf) * rho + (1 - bf) * tau));
      var sigf = .5 * ((1 - bf) * rho + (1 + bf) * tau);
      var att = 1 - sigf, m = Math.sqrt(att * att - sigb * sigb);
      var e1 = Math.exp(-m * lai), e2 = e1 * e1;
      var ri = (att - m) / sigb, ri2 = ri * ri, re = ri * e1;
      var den = 1 - ri2 * e2;
      var j1s = j1(ks, m, lai), j2s = j2(ks, m, lai);
      var j1o = j1(ko, m, lai), j2o = j2(ko, m, lai);
      var ps = (sf + sb * ri) * j1s, qs = (sf * ri + sb) * j2s;
      var pv = (vf + vb * ri) * j1o, qv = (vf * ri + vb) * j2o;
      var rdd = ri * (1 - e2) / den;
      var tsd = (ps - re * qs) / den;
      var tdo = (pv - re * qv) / den;
      var rdo = (qv - re * pv) / den;
      var g1 = (zint - j1s * too) / (ko + m);
      var g2 = (zint - j1o * tss) / (ks + m);
      var multiple = ((vf * ri + vb) * g1 * (sf + sb * ri) +
        (vf + vb * ri) * g2 * (sf * ri + sb) - (rdo * qs + tdo * ps) * ri) / (1 - ri2);
      var rso = (g.sob * rho + g.sof * tau) * lai * hot[1] + multiple;
      var offset = w < 51 ? 0 : 50;
      var count = w < 51 ? 51 : 36;
      for (var b = 0; b < 50; b++) {
        var background = (b + 1) / 100;
        var interaction = background * ((tss + tsd) * tdo + (tsd + tss * background * rdd) * too) /
          (1 - background * rdd);
        result[offset + b] += (rso + hot[0] * background + interaction) / count;
      }
    }
    return result;
  }
  function valid(p) {
    return ['red', 'nir', 'solar_zenith', 'view_zenith', 'relative_azimuth'].every(function (k) {
      return typeof p[k] === 'number' && isFinite(p[k]);
    }) && p.red >= 0 && p.red <= 1 && p.nir >= 0 && p.nir <= 1 &&
      p.solar_zenith >= 0 && p.solar_zenith < 90 &&
      p.view_zenith >= 0 && p.view_zenith < 90 &&
      p.relative_azimuth >= 0 && p.relative_azimuth <= 180;
  }
  // Convert an index residual to BRF units using its gradient at the observation.
  // Assumes equal independent red/NIR errors. One index term has unit weight.
  function indexTerm(red,nir,mode) {
    if(mode==='BRF')return null;
    var sum=red+nir;
    if(mode==='NDVI' && sum===0)return null;
    var den=mode==='NDVI'?sum:1+nir+2.4*red;
    var gain=mode==='NDVI'?1:2.5;
    var dRed=mode==='NDVI'?-2*nir:-2.5*(1+3.4*nir);
    var dNir=mode==='NDVI'?2*red:2.5*(1+3.4*red);
    return {observed:gain*(nir-red)/den, scale:den*den/Math.sqrt(dRed*dRed+dNir*dNir)};
  }
  function createFit(pixels,mode) {
    mode=mode||'BRF';
    if(['BRF','NDVI','EVI2'].indexOf(mode)<0)throw new Error('Unknown PROSAIL fitting objective.');
    var redCosts=new Float64Array(soils.length), nirCosts=new Float64Array(soils.length);
    var groups = {}, costs = new Float64Array(soils.length);
    var ndviCosts=new Float64Array(soils.length), evi2Costs=new Float64Array(soils.length), ndviCount=0;
    var choices = pixels.map(function () { return new Uint8Array(soils.length); });
    pixels.forEach(function (p, i) {
      if (!valid(p)) throw new Error('Invalid reflectance or geometry.');
      var key = [p.solar_zenith, p.view_zenith, p.relative_azimuth].join(',');
      if (!groups[key]) groups[key] = [];
      groups[key].push(i);
    });
    var keys = Object.keys(groups), next = 0;
    // One geometry at a time keeps memory bounded; discard its LUT after scoring its pixels.
    function step() {
      var members = groups[keys[next]], p = pixels[members[0]];
      var g = geometry(p.solar_zenith, p.view_zenith, p.relative_azimuth);
      var lut = [];
      for (var l = 0; l <= 80; l++) lut.push(spectrum(g, l / 10));
      members.forEach(function (index) {
        var pixel = pixels[index], term=indexTerm(pixel.red,pixel.nir,mode);
        var sum=pixel.red+pixel.nir;
        var obsNDVI=sum>0?(pixel.nir-pixel.red)/sum:null;
        var obsEVI2=2.5*(pixel.nir-pixel.red)/(1+pixel.nir+2.4*pixel.red);
        if(sum>0)ndviCount++;
        for (var s = 0; s < soils.length; s++) {
          var redIndex = soils[s][0] - 1, nirIndex = soils[s][1] + 49;
          var best = Infinity, selected = 0, bestRedError=0, bestNirError=0;
          for (var a = 0; a <= 80; a++) {
            var red=lut[a][redIndex], nir=lut[a][nirIndex];
            var dr=red-pixel.red, dn=nir-pixel.nir;
            var error = dr * dr + dn * dn;
            if(term){
              var predicted=mode==='NDVI'?(nir-red)/(nir+red):2.5*(nir-red)/(1+nir+2.4*red);
              var residual=(predicted-term.observed)*term.scale;
              error+=residual*residual;
            }
            if (error < best) { best=error; selected=a; bestRedError=dr; bestNirError=dn; }
          }
          var fitRed=lut[selected][redIndex], fitNir=lut[selected][nirIndex];
          if(sum>0){
            var deltaNDVI=(fitNir-fitRed)/(fitNir+fitRed)-obsNDVI;
            ndviCosts[s]+=deltaNDVI*deltaNDVI;
          }
          var deltaEVI2=2.5*(fitNir-fitRed)/(1+fitNir+2.4*fitRed)-obsEVI2;
          evi2Costs[s]+=deltaEVI2*deltaEVI2;
          costs[s] += best;
          redCosts[s]+=bestRedError*bestRedError; nirCosts[s]+=bestNirError*bestNirError;
          choices[index][s] = selected;
        }
      });
      next++;
      return next === keys.length;
    }
    function finish(tolerance) {
      var score = Array.prototype.map.call(costs, function (c) { return Math.sqrt(c / (2 * pixels.length)); });
      var best=bestBackgroundIndex(), near=[];
      var rmse=Array.prototype.map.call(redCosts,function(c,i){return Math.sqrt((c+nirCosts[i])/(2*pixels.length));});
      score.forEach(function(v,i){if(v<=score[best]+tolerance)near.push(i);});
      var covers = choices.map(function (a) { return 1 - Math.exp(-.05 * a[best]); });
      var profile = soils.map(function (soil, i) {
        return {background_red: soil[0] / 100, background_nir: soil[1] / 100,
          objective: mode, objective_score: score[i], rmse: rmse[i], fcover: median(choices.map(function (a) { return 1 - Math.exp(-.05 * a[i]); }))};
      });
      var medians = near.map(function (i) { return profile[i].fcover; });
      var records = pixels.map(function (p, i) {
        var lo = 1, hi = 0;
        near.forEach(function (s) {
          var f = 1 - Math.exp(-.05 * choices[i][s]);
          lo = Math.min(lo, f); hi = Math.max(hi, f);
        });
        return {date: p.date, x: p.x, y: p.y, red: p.red, nir: p.nir,
          solar_zenith: p.solar_zenith, view_zenith: p.view_zenith, relative_azimuth: p.relative_azimuth,
          lai: choices[i][best] / 10, fcover: covers[i], fcover_low: lo, fcover_high: hi};
      });
      return {fcover: median(covers), fcover_mean: covers.reduce(function (a, b) { return a + b; }, 0) / covers.length,
        lai: median(choices.map(function (a) { return a[best] / 10; })),
        objective:mode, objective_score:score[best], ndvi_rmse:ndviCount?Math.sqrt(ndviCosts[best]/ndviCount):null,
        evi2_rmse:Math.sqrt(evi2Costs[best]/pixels.length), red_rmse:Math.sqrt(redCosts[best]/pixels.length),
        nir_rmse:Math.sqrt(nirCosts[best]/pixels.length), rmse: rmse[best], background_red: soils[best][0] / 100, background_nir: soils[best][1] / 100,
        fcover_low: Math.min.apply(null, medians), fcover_high: Math.max.apply(null, medians),
        near_backgrounds: near.length, pixels: pixels.length, geometries: keys.length,
        profile: profile, records: records};
    }
    function bestBackgroundIndex() {
      var best = 0;
      for (var s = 1; s < costs.length; s++) if (costs[s] < costs[best]) best = s;
      return best;
    }
    function bestCovers() {
      var best=bestBackgroundIndex();
      return choices.map(function (a) { return 1 - Math.exp(-.05 * a[best]); });
    }
    return {step: step, finish: finish, bestCovers: bestCovers, brfRMSE:function(){var i=bestBackgroundIndex();return Math.sqrt((redCosts[i]+nirCosts[i])/(2*pixels.length));}, bestRed: function(){return soils[bestBackgroundIndex()][0]/100;}, geometries: keys.length};
  }
  return {createFit: createFit, valid: valid, geometry: geometry, spectrum: spectrum, median: median};
}());

function rtObjectiveSelector(initial) {
  return ui.Select({items:[{label:'Red + NIR',value:'BRF'},{label:'Red + NIR + NDVI',value:'NDVI'},
    {label:'Red + NIR + EVI2',value:'EVI2'}],value:initial||'BRF'});
}
// Earth Engine supplies screened pixels; the browser fits exact FourSAIL predictions.
// Dates are actual acquisitions retained in MYD09A1, grouped separately within composites.
function rtSample(image, region) {
  image = ee.Image(image);
  var projection = image.select('sur_refl_b01').projection();
  var state = image.select('StateQA'), quality = image.select('QA');
  var clear = state.bitwiseAnd(7).eq(0)
    .and(state.rightShift(3).bitwiseAnd(7).eq(1))
    .and(state.bitwiseAnd((3 << 8) | (1 << 10) | (1 << 12) | (1 << 13) | (1 << 15)).eq(0))
    .and(quality.bitwiseAnd(1023).eq(0));
  var reflectance = image.select(['sur_refl_b01', 'sur_refl_b02'], ['red', 'nir']).multiply(.0001);
  var angles = image.select(['SolarZenith', 'ViewZenith', 'RelativeAzimuth'],
    ['solar_zenith', 'view_zenith', 'relative_azimuth']).multiply(.01);
  angles = angles.addBands(angles.select('relative_azimuth').abs(), null, true);
  var valid = reflectance.gte(0).and(reflectance.lte(1)).reduce(ee.Reducer.min())
    .and(angles.select(['solar_zenith', 'view_zenith']).gte(0).reduce(ee.Reducer.min()))
    .and(angles.select(['solar_zenith', 'view_zenith']).lt(90).reduce(ee.Reducer.min()))
    .and(angles.select('relative_azimuth').lte(180));
  return reflectance.addBands(angles).addBands(image.select('DayOfYear').rename('doy'))
    .addBands(ee.Image.pixelCoordinates(projection)).updateMask(clear.and(valid))
    .sample({region: region, projection: projection, geometries: false, tileScale: 4});
}
function rtBox(longitude, latitude, km) {
  var zone = Math.min(60, Math.floor((longitude + 180) / 6) + 1);
  var projection = ee.Projection('EPSG:' + ((latitude >= 0 ? 32600 : 32700) + zone));
  return ee.Geometry.Point([longitude, latitude]).buffer(km * 500, ee.ErrorMargin(1, 'projected'), projection)
    .bounds(ee.ErrorMargin(1, 'projected'), projection);
}
function rtISODate(text) {
  if (!/^\d{4}-\d{2}-\d{2}$/.test(text)) return null;
  var value = Date.parse(text + 'T00:00:00Z');
  return isFinite(value) && new Date(value).toISOString().slice(0, 10) === text ? value : null;
}
function rtGroups(features, year, start, end) {
  var groups = {}, seen = {};
  features.forEach(function (feature) {
    var p = feature.properties;
    var stamp = Date.UTC(year, 0, p.doy);
    if (!RT.valid(p) || !isFinite(p.doy) || p.doy < 1 || p.doy > 366 ||
        stamp < start || stamp > end) return;
    p.date = new Date(stamp).toISOString().slice(0, 10);
    var key = p.date + '/' + p.x + '/' + p.y;
    if (seen[key]) return;
    seen[key] = true;
    if (!groups[p.date]) groups[p.date] = [];
    groups[p.date].push(p);
  });
  return groups;
}
function rtApp() {
  var previous = ui.root.widgets().getJsArray().slice();
  var map = ui.Map();
  map.setCenter(-120.966, 38.4309, 11);
  var panel = ui.Panel({style: {width: '470px', padding: '10px'}});
  var charts = ui.Panel();
  ui.root.clear();
  ui.root.add(ui.SplitPanel(panel, map));
  panel.add(ui.Button('Return to ACCLIM', function(){token++;ui.root.widgets().reset(previous);}));
  panel.add(ui.Label('ACCLIM-RT', {fontSize: '22px', fontWeight: 'bold'}));
  panel.add(ui.Label('Whole-pixel foliage cover from Aqua MODIS red and NIR. Each observation date has a shared fitted background.'));
  var source = ui.Select({items: ['Coordinates / map click', 'FeatureCollection asset center'],
    value: 'Coordinates / map click'});
  var lat = ui.Textbox({value: '38.4309'}), lon = ui.Textbox({value: '-120.966'});
  var asset = ui.Textbox({placeholder: 'projects/.../assets/...'});
  var start = ui.Textbox({value: '2012-05-01'}), end = ui.Textbox({value: '2012-05-08'});
  var box = ui.Textbox({value: '10'});
  var tolerance = ui.Textbox({value: '0.001'});
  function field(label, widget) { panel.add(ui.Label(label)); panel.add(widget); }
  field('Location source', source);
  field('Latitude', lat); field('Longitude', lon); field('Asset (accessible to app users)', asset);
  field('Start date (YYYY-MM-DD)', start); field('End date, inclusive (YYYY-MM-DD)', end);
  field('Box side (km)', box);
  var objective=rtObjectiveSelector(typeof UI!=='undefined' && UI.rtObjective ? UI.rtObjective.getValue() : 'BRF');
  field('PROSAIL fitting objective',objective);
  field('Near-best objective tolerance (BRF-equivalent)', tolerance);
  panel.add(ui.Label('MYD09A1 retains selected acquisitions within each 8-day composite. The fit uses their actual dates and angles. Leaf properties are fixed; fcover = 1 âˆ’ exp(âˆ’0.5 Ã— LAI).'));
  var status = ui.Label('Ready.');
  var run = ui.Button('Run ACCLIM-RT'), cancel = ui.Button('Cancel');
  panel.add(ui.Panel([run, cancel], ui.Panel.Layout.flow('horizontal')));
  panel.add(status); panel.add(charts);
  var token = 0, cacheKey = '', samples = {};
  var fits = {}, summaries = [], selected = null;
  map.onClick(function (point) {
    lat.setValue(point.lat.toFixed(6)); lon.setValue(point.lon.toFixed(6));
    source.setValue('Coordinates / map click');
  });
  cancel.onClick(function () {
    token++; run.setDisabled(false); render(); status.setValue('Cancelled. Completed dates remain available.');
  });
  function download(label, rows, name, holder) {
    if (!rows.length) return;
    var link = ui.Label('', {color: '#1a73e8'});
    var button = ui.Button(label);
    holder.add(button); holder.add(link);
    button.onClick(function () {
      var active = token;
      button.setDisabled(true);
      link.setValue('Preparing CSV...');
      var fields = {};
      rows.forEach(function (row) { Object.keys(row).forEach(function (key) { fields[key] = true; }); });
      ee.FeatureCollection(rows.map(function (r) { return ee.Feature(null, r); }))
        .getDownloadURL('CSV', Object.keys(fields), name, function (url, error) {
          if (active !== token) return;
          button.setDisabled(false);
          if (error || !url) { link.setValue('CSV export failed: ' + String(error)); return; }
          link.setValue('Download CSV'); link.setUrl(url);
        });
    });
  }
  function timeTable(rows, fields) {
    return {cols: [{label: 'Date', type: 'date'}].concat(fields.map(function (field) {
      return {label: field, type: 'number'};
    })), rows: rows.map(function (row) {
      var date = row.date.split('-');
      return {c: [{v: new Date(Number(date[0]), Number(date[1]) - 1, Number(date[2]))}]
        .concat(fields.map(function (field) { return {v: row[field]}; }))};
    })};
  }
  function render() {
    charts.clear();
    var complete = summaries.filter(function (r) { return r.status === 'fitted'; });
    if (!complete.length) {
      download('Download date status CSV', summaries, 'acclim_rt_dates', charts);
      return;
    }
    complete.sort(function (a, b) { return a.date.localeCompare(b.date); });
    charts.add(ui.Chart(timeTable(complete, ['fcover', 'fcover_low', 'fcover_high'])).setChartType('LineChart')
      .setOptions({title: 'Foliage cover', hAxis: {title: 'Date', format: 'yyyy-MM-dd'},
        vAxis: {title: 'Fraction', viewWindow: {min: 0, max: 1}}, pointSize: 3, lineWidth: 1}));
    charts.add(ui.Chart(timeTable(complete, ['rmse'])).setChartType('LineChart')
      .setOptions({title: 'Joint red–NIR RMSE', hAxis: {title: 'Date', format: 'yyyy-MM-dd'}, pointSize: 3}));
    charts.add(ui.Label('Cover bounds span near-best backgrounds. Coverage is the percentage of native box pixels observed on that actual date.'));
    download('Download date summary CSV', summaries, 'acclim_rt_dates', charts);
    var dates = complete.map(function (r) { return r.date; });
    var details = ui.Panel();
    var chooser = ui.Select({items: dates, value: selected && fits[selected] ? selected : dates[0],
      onChange: showDate});
    charts.add(ui.Label('Inspect date')); charts.add(chooser); charts.add(details);
    function showDate(day) {
      selected = day;
      details.clear();
      var fit = fits[day];
      var profile = {cols: [{label: 'Median fcover', type: 'number'}, {label: 'RMSE', type: 'number'}],
        rows: fit.profile.map(function (r) { return {c: [{v: r.fcover}, {v: r.rmse}]}; })};
      details.add(ui.Chart(profile).setChartType('ScatterChart')
        .setOptions({title: day + ' â€” RMSE vs median fcover',
          hAxis: {title: 'Median fcover'}, vAxis: {title: 'Joint redâ€“NIR RMSE'}, pointSize: 2, legend: 'none'}));
      details.add(ui.Label('Background red ' + fit.background_red.toFixed(2) +
        ', NIR ' + fit.background_nir.toFixed(2) + '; median LAI ' + fit.lai.toFixed(2) +
        '; ' + fit.pixels + ' pixels, ' + fit.geometries + ' distinct geometries.'));
      download('Download pixel fcover CSV', fit.records, 'acclim_rt_pixels_' + day, details);
      download('Download RMSEâ€“fcover CSV', fit.profile, 'acclim_rt_profile_' + day, details);
    }
    showDate(chooser.getValue());
  }
  run.onClick(function () {
    var first = rtISODate(start.getValue()), last = rtISODate(end.getValue());
    var km = Number(box.getValue()), tol = Number(tolerance.getValue());
    if (first === null || last === null || first > last || !box.getValue().trim() ||
        !isFinite(km) || km <= 0 || !tolerance.getValue().trim() || !isFinite(tol) || tol < 0) {
      status.setValue('Enter valid dates, a positive box size, and a nonnegative RMSE tolerance.'); return;
    }
    var active = ++token, fitMode=objective.getValue();
    run.setDisabled(true);
    fits = {}; summaries = []; selected = null; charts.clear();
    function fail(error) {
      if (active !== token) return;
      run.setDisabled(false); status.setValue('Error: ' + String(error)); render();
    }
    function begin(longitude, latitude) {
      if (active !== token) return;
      var region = rtBox(longitude, latitude, km);
      map.layers().reset([ui.Map.Layer(region, {color: 'yellow'}, 'Analysis box')]);
      map.centerObject(region);
      var key = [longitude, latitude, km, first, last].join('/');
      if (cacheKey !== key) { samples = {}; cacheKey = key; }
      var collection = ee.ImageCollection('MODIS/061/MYD09A1')
        .filterDate(ee.Date(first).advance(-7, 'day'), ee.Date(last).advance(1, 'day'))
        .filterBounds(region).sort('system:time_start');
      status.setValue('Finding Aqua MODIS composites...');
      collection.aggregate_array('system:id').evaluate(function (ids, error) {
        if (active !== token) return;
        if (error) { fail(error); return; }
        if (!ids.length) { fail('No Aqua MODIS composites in the requested range. Records begin July 2002.'); return; }
        var projection = ee.Image(ids[0]).select('sur_refl_b01').projection();
        ee.Image.pixelCoordinates(projection).sample({region: region, projection: projection,
          geometries: false, tileScale: 4}).size().evaluate(function (total, countError) {
          if (active !== token) return;
          if (countError) { fail(countError); return; }
          if (!total) { fail('No native MODIS pixel centers in this box.'); return; }
          var position = 0, failures = 0, empty = 0;
          function nextComposite() {
            if (active !== token) return;
            if (position === ids.length) {
              run.setDisabled(false); render();
              status.setValue('Finished: ' + Object.keys(fits).length + ' fitted dates; ' +
                empty + ' composites without usable dates; ' + failures + ' failed composites.');
              return;
            }
            var id = ids[position++];
            status.setValue('Loading composite ' + position + ' / ' + ids.length + '...');
            function received(data, fetchError) {
              if (active !== token) return;
              if (fetchError) {
                failures++;
                summaries.push({date: id.split('/').pop().replace(/_/g, '-'), status: 'fetch_error',
                  error: String(fetchError)});
                ui.util.setTimeout(nextComposite, 0); return;
              }
              samples[id] = data;
              var year = Number(id.split('/').pop().slice(0, 4));
              var grouped = rtGroups(data.features, year, first, last);
              var days = Object.keys(grouped).sort(), index = 0;
              if (!days.length) empty++;
              function nextDate() {
                if (active !== token) return;
                if (index === days.length) { render(); ui.util.setTimeout(nextComposite, 0); return; }
                var day = days[index++], pixels = grouped[day];
                if (pixels.length < 4) {
                  summaries.push({date: day, status: 'insufficient_pixels', pixels: pixels.length,
                    coverage_percent: 100 * pixels.length / total});
                  ui.util.setTimeout(nextDate, 0); return;
                }
                var fit = RT.createFit(pixels,fitMode), progress = 0;
                function advance() {
                  if (active !== token) return;
                  try {
                    var done = fit.step();
                    progress++;
                    status.setValue(day + ': geometry ' + progress + ' / ' + fit.geometries +
                      '; composite ' + position + ' / ' + ids.length);
                    if (!done) { ui.util.setTimeout(advance, 0); return; }
                    var result = fit.finish(tol);
                    fits[day] = result;
                    summaries.push({date: day, status: 'fitted', fcover: result.fcover,
                      fcover_mean: result.fcover_mean, fcover_low: result.fcover_low,
                      fcover_high: result.fcover_high, lai: result.lai, rmse: result.rmse, objective:fitMode, objective_score:result.objective_score, red_rmse:result.red_rmse, nir_rmse:result.nir_rmse, ndvi_rmse:result.ndvi_rmse, evi2_rmse:result.evi2_rmse,
                      background_red: result.background_red, background_nir: result.background_nir,
                      pixels: result.pixels, geometries: result.geometries,
                      coverage_percent: 100 * result.pixels / total, near_backgrounds: result.near_backgrounds,
                      rmse_tolerance: tol, latitude: latitude, longitude: longitude, box_km: km,
                      source: 'MODIS/061/MYD09A1', composite: id.split('/').pop()});
                    ui.util.setTimeout(nextDate, 0);
                  } catch (fitError) { fail(fitError); }
                }
                ui.util.setTimeout(advance, 0);
              }
              nextDate();
            }
            if (samples[id]) received(samples[id]);
            else rtSample(ee.Image(id), region).evaluate(received);
          }
          nextComposite();
        });
      });
    }
    if (source.getValue() === 'FeatureCollection asset center') {
      if (!asset.getValue().trim()) { fail('Enter an asset ID.'); return; }
      ee.FeatureCollection(asset.getValue().trim()).geometry().centroid(1).coordinates().evaluate(function (xy, error) {
        if (active !== token) return;
        if (error) { fail(error); return; }
        begin(xy[0], xy[1]);
      });
    } else {
      var latitude = Number(lat.getValue()), longitude = Number(lon.getValue());
      if (!lat.getValue().trim() || !lon.getValue().trim() || !isFinite(latitude) || !isFinite(longitude) ||
          latitude < -80 || latitude > 84 || longitude < -180 || longitude > 180) {
        fail('Enter longitude âˆ’180 to 180 and latitude âˆ’80 to 84 (UTM box).'); return;
      }
      begin(longitude, latitude);
    }
  });
}
// Retrieve the same date-specific spatial median used by the standalone RT analysis.
function rtCoverPeriods(observations, periods, gapDays) {
  var obs = observations.slice().sort(function(a,b){return a.millis-b.millis;});
  var j = 0, result = {};
  periods.forEach(function(period) {
    var values = [], measured = 0;
    for (var t = period.millis; t < period.period_end_millis; t += 86400000) {
      while (j < obs.length && obs[j].millis < t) j++;
      var right = obs[j], left = obs[j-1];
      if (right && right.millis === t) { values.push(right.cover); measured++; }
      else if (right && left && right.millis-left.millis <= gapDays*86400000) {
        values.push(left.cover + (right.cover-left.cover)*(t-left.millis)/(right.millis-left.millis));
      }
    }
    var entry = {rt_observed_days:measured, rt_interpolated_days:values.length-measured,
      rt_missing_pct:100*(1-values.length/period.period_days)};
    if (values.length && (period.model_cadence==='monthly' || values.length===period.period_days)) {
      entry.rt_fcover = period.model_cadence==='monthly' ? RT.median(values) :
        values.reduce(function(a,b){return a+b;},0)/values.length;
    }
    result[period.date]=entry;
  });
  return result;
}
function fetchRTCover(geom,start,end,km,mode,request,done) {
  var gap = SATELLITE_SOURCES.MODIS8.gap;
  var first = carbonDateMillis(start)-gap*86400000, last = carbonDateMillis(end)+gap*86400000;
  var observations = [];
  geom.centroid(1).coordinates().evaluate(function(xy,error) {
    if(request!==UI._requestId)return;
    if(error){done(null,error);return;}
    var region = rtBox(xy[0],xy[1],km);
    var collection = ee.ImageCollection('MODIS/061/MYD09A1').filterBounds(region)
      .filterDate(ee.Date(first).advance(-7,'day'),ee.Date(last).advance(1,'day')).sort('system:time_start');
    collection.aggregate_array('system:id').evaluate(function(ids,error) {
      if(request!==UI._requestId)return;
      if(error){done(null,error);return;}
      var index=0;
      function nextImage() {
        if(request!==UI._requestId)return;
        if(index===ids.length){done(observations);return;}
        var id=ids[index++];
        setStatus('Retrieving PROSAIL fcover: composite '+index+' of '+ids.length);
        rtSample(ee.Image(id),region).evaluate(function(data,error) {
          if(request!==UI._requestId)return;
          if(error){done(null,error);return;}
          var groups=rtGroups(data.features,Number(id.split('/').pop().slice(0,4)),first,last);
          var days=Object.keys(groups).filter(function(d){return groups[d].length>=4;}).sort();
          function nextDay() {
            if(request!==UI._requestId)return;
            if(!days.length){nextImage();return;}
            var date=days.shift(), fit=RT.createFit(groups[date],mode);
            setStatus('Fitting PROSAIL fcover for '+date+' ('+groups[date].length+' pixels)');
            function advance() {
              if(request!==UI._requestId)return;
              try {
                if(!fit.step()){ui.util.setTimeout(advance,0);return;}
                var covers=fit.bestCovers();
                var histogram=Array.apply(null,Array(81)).map(function(){return 0;});
                covers.forEach(function(c){histogram[Math.round(-20*Math.log(1-c))]++;});
                observations.push({millis:carbonDateMillis(date),cover:RT.median(covers),
                  index_bins:rtIndexBins(groups[date],covers),lai_histogram:histogram,background_red:fit.bestRed(),latitude:xy[1],brf_rmse:fit.brfRMSE()});
                nextDay();
              } catch(e){done(null,e.message);}
            }
            advance();
          }
          nextDay();
        });
      }
      nextImage();
    });
  });
}

// Equal total weight for each acquisition date. Keep within-date pixel covariance.
function rtEqualDateLine(moments) {
  var sx=0, sy=0, sxx=0, sxy=0, dates=moments.length, pixels=0;
  moments.forEach(function(d){
    sx+=d.x; sy+=d.y; sxx+=d.xx; sxy+=d.xy; pixels+=d.n;
  });
  var variance=dates ? sxx/dates-Math.pow(sx/dates,2) : 0;
  if(variance<=1e-12)return null;
  var slope=(sxy/dates-(sx/dates)*(sy/dates))/variance;
  return {slope:slope,intercept:(sy-slope*sx)/dates,dates:dates,pixels:pixels};
}

// END ACCLIM-RT RETRIEVAL
function runForGeometryAsync(geom,start,end,params,config,done) {
  var request=UI._requestId, cfg=mergeParams(config,{});
  var overrides = params || {};
  var autoNDVI = overrides.autoSoilNDVImin === true, autoEVI2 = overrides.autoSoilEVI2min === true;
  params=mergeParams(ACCLIM_PARAMS,overrides);
  if(params.faparMethod==='NDVI_GE_OFFSET')params.fitNDVIminGE=true;
  var gePeriods=[],ndviFit,historicalFit;
  function complete(rows){
    gePeriods=normalizeGE(gePeriods);
    var eviFit=geLinearFit(gePeriods,'EVI2_GE'),ndviDiagnosticFit=geLinearFit(gePeriods,'NDVI_GE');
    var byDate={};gePeriods.forEach(function(r){byDate[r.date]=r;});
    var spec=geSourceSpec(cfg.satelliteSource||'MODIS16');
    var enriched=rows.map(function(r){var d=mergeParams(r,byDate[r.date]||{GE:null,EVI2_GE:null,NDVI_GE:null,GE_available_pct:0});
      d.evi2_min_ge=eviFit?eviFit.intercept:null;d.ndvi_min_ge=ndviDiagnosticFit?ndviDiagnosticFit.intercept:null;
      d.NIRv=carbonFinite(d.nir)&&carbonFinite(d.NDVI)?d.nir*d.NDVI:null;
      if(params.fiparMethod==='RT_NDVI') {
        d.m_cover=params.m_cover; d.b_cover=params.b_cover;
        d.cover_fit_dates=params.cover_fit_dates; d.cover_fit_pixels=params.cover_fit_pixels;
      }
      var absIndex=rtAbsIndex(params.faparMethod);
      if(absIndex){
        d['m_abs_'+absIndex]=params['m_abs_'+absIndex];d['b_abs_'+absIndex]=params['b_abs_'+absIndex];
        d.abs_fit_dates=params.abs_fit_dates;d.abs_fit_pixels=params.abs_fit_pixels;
      }
      if(absorptionObs)d.rt_fit_objective=cfg.rtObjective||'BRF';
      if(params.faparMethod==='KNDVI_SPATIAL')d.kndvi_radius_km=sigmaRadius;
      if(params.faparMethod==='KNDVI') d.kndvi_sigma=params.kndvi_sigma;
      d.FAPARmult=params.FAPARmult;
      d.GE_scale=params.GE_scale;d.ndvi_exp=params.ndvi_exp;d.ndvi_max_used=params.ndvi_max;d.EVI2min_parameter=params.EVI2min;d.EVI2max_parameter=params.EVI2max;
      d.fge=params.fge_exp===0?1:(carbonFinite(d.GE)?Math.pow(Math.max(0,Math.min(1,d.GE/params.GE_scale)),params.fge_exp):null);
      d.ndvi_min_used=params.ndvi_min;d.ndvi_min_source=params.fitVIMinP05?'10-year observed fifth percentile':params.fitNDVIminGE?'NDVI at GE′=0':autoNDVI?'SoilGrids texture':'configured';
      if(historicalFit){d.ndvi_min_p05=historicalFit.ndvi;d.evi2_min_p05=historicalFit.evi2;d.vi_baseline_n=historicalFit.n;d.vi_baseline_start=historicalFit.start;d.vi_baseline_end=historicalFit.end;d.vi_baseline_first_observation=historicalFit.first;d.vi_baseline_last_observation=historicalFit.last;}
      if(ndviFit){d.ndvi_min_fit_n=ndviFit.n;d.ndvi_min_fit_slope=ndviFit.slope;d.ndvi_min_fit_r2=ndviFit.r2;
        d.ndvi_min_fit_extrapolated=ndviFit.extrapolated;d.ndvi_min_fit_start=start;d.ndvi_min_fit_end=end;}
      d.GE_missing_for_GPP=((params.faparMethod==='FGE' || params.faparMethod==='GE_NDVI_MIN' || params.faparMethod==='NDVI_BLEND' || params.fiparMethod==='NDVI_BLEND') && params.fge_exp>0)&&!carbonFinite(d.GE);d.fge_exp=params.fge_exp;
      d.GE_source=spec.label;d.GE_scale_m=spec.scale;d.GE_status=carbonFinite(d.GE)?'available':'unavailable';return d;
    });
    done(ee.FeatureCollection(enriched.map(function(r){return ee.Feature(geom,r);})).sort('millis'));
  }
  function materialize(fc,next) {
    fc.toList(fc.size().max(1)).map(function(f){return ee.Feature(f).toDictionary();}).evaluate(function(rows,error) {
      if(request!==UI._requestId)return;
      if(error){done(null,error);return;}
      try {next(rows,ee.FeatureCollection(rows.map(function(r){return ee.Feature(geom,r);})));
      } catch(error) {done(null,error.message);}
    });
  }
  function simulate() {
    setStatus('Aligning satellite observations with weather timesteps...');
    materialize(assembleSelectedForcing(geom,start,end,cfg),function(rows,forcing) {
      var geByDate={};gePeriods.forEach(function(r){geByDate[r.date]=r;});
      rows=rows.map(function(r){var d=mergeParams(r,{}),g=geByDate[r.date];
        d.GE=g&&carbonFinite(g.GE)?g.GE:null;
        d.GE_normalized=g&&carbonFinite(g.GE_normalized)?g.GE_normalized:null;
        d.GE_missing_for_GPP=((params.faparMethod==='FGE' || params.faparMethod==='GE_NDVI_MIN' || params.faparMethod==='NDVI_BLEND' || params.fiparMethod==='NDVI_BLEND') && params.fge_exp>0)&&!carbonFinite(d.GE);
        if(d.GE_missing_for_GPP)d.model_valid=0;
        return d;
      });
      forcing=ee.FeatureCollection(rows.map(function(r){return ee.Feature(geom,r);}));
      if(!rows.some(function(r){return r.model_valid===1||r.model_valid===true;})){complete(rows);return;}
      setStatus('Calculating native '+rows[0].model_cadence+' physiology...');
      materialize(runACCLIM(forcing,mergeParams(params,{deferBuckets:true})),function(baselineRows,baseline) {
        setStatus('Fitting soil-water buckets and bridging missing forcing...');
        // Fit once over the full trajectory; only independent final diagnostics are batched.
        var fittedRows=fitNativeBuckets(baselineRows,params), finalRows=[], offset=0;
        function finishBatch() {
          if(request!==UI._requestId)return;
          if(offset>=fittedRows.length){complete(finalRows);return;}
          var endOffset=Math.min(offset+90,fittedRows.length);
          setStatus('Calculating final flux diagnostics: '+(offset+1)+'–'+endOffset+' of '+fittedRows.length+' periods...');
          var fitted=ee.FeatureCollection(fittedRows.slice(offset,endOffset).map(function(r){return ee.Feature(geom,r);}));
          // Keep the full baseline for any fitted EVI2 bounds, consistent across batches.
          materialize(finishACCLIM(baseline,mergeParams(params,{fittedBuckets:fitted})),function(result) {
            finalRows=finalRows.concat(result);offset=endOffset;finishBatch();
          });
        }
        finishBatch();
      });
    });
  }
  cfg.timestep=cfg.timestep||(cfg.weatherPlan?cfg.weatherPlan.cadence:WEATHER_SOURCES[cfg.weatherSource||'TERRACLIMATE'].daily?'daily':'monthly');
  var periods=selectedTimestepPeriods(start,end,cfg.timestep);
  cfg.periodCalendar=carbonRowsToFC(periods);
  var first=start,last=end;
  if(cfg.weatherPlan&&cfg.weatherPlan.records.length){
    first=cfg.weatherPlan.records[0].data.date||cfg.weatherPlan.records[0].data.month+'-01';
    var lastRecord=cfg.weatherPlan.records[cfg.weatherPlan.records.length-1].data;
    if(lastRecord.period_days)last=carbonDateString(carbonDateMillis(lastRecord.date)+(lastRecord.period_days-1)*86400000);
  }
  var spec=SATELLITE_SOURCES[cfg.satelliteSource||'MODIS16'],pad=spec.gap+32;
  // Bump the version if sampling, QA, or GE aggregation semantics change.
  var sigmaRegion=null, sigmaRadius=cfg.kndviRadiusKm===undefined?5:cfg.kndviRadiusKm;
  if(params.faparMethod==='KNDVI_SPATIAL'){
    if(!carbonFinite(sigmaRadius)||sigmaRadius<=0){done(null,'Enter a positive spatial sigma radius (km).');return;}
    sigmaRegion=geom.centroid(1).buffer(sigmaRadius*1000,1);
  }
  var inputKey=JSON.stringify(['ACCLIM-inputs-v1',geom.serialize(),start,end,first,last,
    cfg.weatherPlan?'tower':(cfg.weatherSource||'TERRACLIMATE'),cfg.satelliteSource||'MODIS16',cfg.timestep]);
  if(sigmaRegion)inputKey+=':kndvi-spatial-wang-v2:'+sigmaRadius;
  var store=cfg.inputStore||sessionInputStore;
  var saved=cfg.refreshInputs?null:store.get(inputKey);
  var weatherRows=[],diagnosticRows=null;
  var absorptionObs=null;
  function alignAndSimulate(rawSatellite) {
    if(absorptionObs) {
      var radiation={};
      var tower=cfg.weatherPlan;
      var source=tower?tower.records.map(function(e){return e.data;}):weatherRows;
      var cadence=tower?tower.cadence:(WEATHER_SOURCES[cfg.weatherSource||'TERRACLIMATE'].daily?'daily':'monthly');
      source.forEach(function(r){
        var t=carbonDateMillis(r.date||r.month+'-01');
        var days=r.period_days||(cadence==='monthly'?(Date.UTC(new Date(t).getUTCFullYear(),new Date(t).getUTCMonth()+1,1)-t)/86400000:1);
        for(var i=0;i<days;i++)radiation[t+i*86400000]=r.SWRdaily;
      });
      var absorption=rtAbsorptionPeriods(absorptionObs,periods,absorptionObs[0].latitude,cfg.rtDiffuse,radiation);
      var cover=rtCoverPeriods(absorptionObs,periods,SATELLITE_SOURCES.MODIS8.gap);
      Object.keys(absorption).forEach(function(d){absorption[d]=mergeParams(cover[d],absorption[d]);});
      cfg.rtPeriods=ee.Dictionary(absorption);
    }
    cfg.alignedSatellite=carbonRowsToFC(alignSatelliteRows(rawSatellite,first,last,cfg.satelliteSource||'MODIS16'));
    simulate();
  }
  var cursor=carbonDateMillis(first)-pad*86400000, stop=carbonDateMillis(last)+pad*86400000, sampled=[];
  function satelliteBatch() {
    if(request!==UI._requestId)return;
    if(cursor>stop) {
      if (diagnosticRows!==null) store.put(inputKey,{weather:weatherRows,satellite:sampled,ge:diagnosticRows});
      alignAndSimulate(sampled);return;
    }
    var next=Math.min(cursor+90*86400000,stop+86400000),a=carbonDateString(cursor),b=carbonDateString(next-86400000);
    setStatus('Sampling satellite reflectance: '+a+' to '+b+'...');
    materialize(sampleSatellite(geom,a,b,cfg.satelliteSource||'MODIS16',true,sigmaRegion),function(rows) {
      sampled=sampled.concat(rows);cursor=next;satelliteBatch();
    });
  }
  function startInputs(){
  var weatherKey=cfg.weatherSource||'TERRACLIMATE';
  function useWeather(rows) {
    cfg.sampledWeather=carbonRowsToFC(aggregateWeatherPeriods(rows,periods,WEATHER_SOURCES[weatherKey].daily?'daily':'monthly'));
  }
  if (saved) {
    setStatus('Reusing weather and satellite inputs for this location...');
    if (!cfg.weatherPlan) {weatherRows=saved.weather;useWeather(weatherRows);}
    alignAndSimulate(saved.satellite);
    return;
  }
  if(cfg.weatherPlan)satelliteBatch();
  else {
    setStatus('Sampling selected weather data...');
    var weatherStart=WEATHER_SOURCES[weatherKey].daily?first:first.slice(0,7)+'-01';
    materialize(sampleWeather(geom,weatherStart,last,weatherKey),function(rows) {
      weatherRows=rows;useWeather(rows);satelliteBatch();
    });
  }
  }
  function useGE(diagnostic,error){
    diagnosticRows=error?null:diagnostic;
    if(request!==UI._requestId)return;
    if(error&&(params.fitNDVIminGE||((params.faparMethod==='FGE' || params.faparMethod==='GE_NDVI_MIN' || params.faparMethod==='NDVI_BLEND' || params.fiparMethod==='NDVI_BLEND') && params.fge_exp>0))){done(null,error);return;}
    try{
      gePeriods=normalizeGE(diagnostic||[]);
      cfg.gePeriods=ee.FeatureCollection(gePeriods.map(function(r){return ee.Feature(null,r);}));
      if(params.fitNDVIminGE){ndviFit=ndviMinGEFit(gePeriods,params);params.ndvi_min=ndviFit.intercept;}
      if(params.fitVIMinP05){
        fetchHistoricalVIMinima(geom,start,cfg.satelliteSource||'MODIS16',params,request,function(fit,error){
          if(error){done(null,error);return;}
          historicalFit=fit;params.ndvi_min=fit.ndvi;params.EVI2min=fit.evi2;
          cfg.modelParams=params;startInputs();
        });return;
      }
      function prepareInputs() {
        if (params.ndvi_max <= params.ndvi_min || params.EVI2max <= params.EVI2min) {
          done(null,'Each VI maximum must exceed its resolved minimum. Check vegetation parameters.');
          return;
        }
        cfg.modelParams=params;
        setStatus('NDVImin = '+params.ndvi_min.toFixed(4)+', EVI2min = '+params.EVI2min.toFixed(4)+'; preparing model inputs...');
        startInputs();
      }
      if ((autoNDVI && !params.fitNDVIminGE) || autoEVI2) {
        fetchSoilGridsTexture(geom).evaluate(function(texture, error) {
          if (request !== UI._requestId) return;
          var soil = !error && texture && soilSpectralEstimate(texture);
          if (!soil) {
            done(null,'SoilGrids texture unavailable for this location. Clear Auto: SoilGrids texture and enter numerical VI minima to run.');
            return;
          }
          if (autoNDVI && !params.fitNDVIminGE) params.ndvi_min = soil.ndvi;
          if (autoEVI2) params.EVI2min = soil.evi2;
          try { prepareInputs(); } catch (err) { done(null,err.message); }
        });
      } else {
        prepareInputs();
      }
    }catch(error){done(null,error.message);}
  }
  function startModelInputs() {
    if (saved) useGE(saved.ge);
    else fetchGEPeriods(geom,start,end,cfg.satelliteSource||'MODIS16',periods,request,useGE);
  }
  if(rtAbsorptionMethod(params.faparMethod) || params.faparMethod==='PROSAIL' || params.fiparMethod==='PROSAIL' || params.fiparMethod==='RT_NDVI') {
    if(cfg.rtDiffuse===undefined)cfg.rtDiffuse=.5;
    if(!carbonFinite(cfg.rtDiffuse)||cfg.rtDiffuse<0||cfg.rtDiffuse>1){done(null,'Enter a diffuse PAR photon fraction between 0 and 1.');return;}
    var km=cfg.rtBoxKm===undefined?10:cfg.rtBoxKm;
    if(!isFinite(km)||km<=0){done(null,'Enter a positive PROSAIL box side in km.');return;}
    var fitMode=cfg.rtObjective||'BRF';
    var rtKey=JSON.stringify(['RT-MYD09A1-v5',fitMode,geom.serialize(),km,start,end]);
    function useCover(obs,error) {
      if(request!==UI._requestId)return;
      if(error){done(null,error);return;}
      if(!obs.length){done(null,'No PROSAIL fcover fits for the selected area and date range.');return;}
      store.put(rtKey,obs);
      absorptionObs=obs;
      var absIndex=rtAbsIndex(params.faparMethod);
      if(absIndex){
        var absFit=rtFitDaylightLine(obs,carbonDateMillis(start),carbonDateMillis(end),cfg.rtDiffuse,absIndex,false);
        if(!absFit){done(null,'Chlorophyll calibration needs matched daylight pixels with varying '+absIndex.toUpperCase()+' inside the selected dates.');return;}
        params['m_abs_'+absIndex]=absFit.slope;params['b_abs_'+absIndex]=absFit.intercept;
        params.abs_fit_dates=absFit.dates;params.abs_fit_pixels=absFit.pixels;
      }
      if(params.fiparMethod==='RT_NDVI') {
        var line=rtFitDaylightLine(obs,carbonDateMillis(start),carbonDateMillis(end),cfg.rtDiffuse,'ndvi',true);
        if(!line){done(null,'PROSAIL calibration needs matched pixels with varying NDVI inside the selected dates.');return;}
        params=mergeParams(params,{m_cover:line.slope,b_cover:line.intercept,cover_fit_dates:line.dates,cover_fit_pixels:line.pixels});
      }
      startModelInputs();
    }
    var cachedCover=cfg.refreshInputs?null:store.get(rtKey);
    if(cachedCover)useCover(cachedCover);else fetchRTCover(geom,start,end,km,fitMode,request,useCover);
  } else startModelInputs();
}

function runForGeometry(geom, startDate, endDate, paramOverrides, config) {
  config = config || {};
  var forcing = assembleSelectedForcing(geom,startDate,endDate,config);
  var usable = forcing.filter(ee.Filter.eq('model_valid',1));
  return ee.FeatureCollection(ee.Algorithms.If(usable.size().gt(0),
    runACCLIM(forcing,paramOverrides),forcing));
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
// Aggregate native daily results for monthly charts without rerunning physiology.
function alignSatelliteRows(observations,start,end,key) {
  var spec=SATELLITE_SOURCES[key],groups={};
  observations.forEach(function(r){var k=String(r.millis),n=r.satellite_observation_count||1;
    if(!groups[k])groups[k]={millis:r.millis,red:0,nir:0,n:0,ge:0,geN:0,kSum:0,sigmaSum:0,kN:0};
    groups[k].red+=r.red*n;groups[k].nir+=r.nir*n;groups[k].n+=n;
    if(carbonFinite(r.kNDVI)&&carbonFinite(r.kndvi_sigma)){
      groups[k].kSum+=n*r.kNDVI;groups[k].sigmaSum+=n*r.kndvi_sigma;groups[k].kN+=n;
    }
    if(carbonFinite(r.GE)){var gn=r.GE_observation_count||n;groups[k].ge+=r.GE*gn;groups[k].geN+=gn;}
  });
  var obs=Object.keys(groups).map(function(k){var g=groups[k];return {millis:g.millis,red:g.red/g.n,nir:g.nir/g.n,kNDVI:g.kN?g.kSum/g.kN:null,kndvi_sigma:g.kN?g.sigmaSum/g.kN:null,GE:g.geN?g.ge/g.geN:null};})
    .sort(function(a,b){return a.millis-b.millis;}),rows=[],j=0;
  for(var t=carbonDateMillis(start);t<=carbonDateMillis(end);t+=86400000) {
    while(j<obs.length && obs[j].millis<t)j++;
    var r=obs[j],l=r&&r.millis===t?r:obs[j-1];
    var gap=l&&r?(r.millis-l.millis)/86400000:null,valid=l&&r&&gap<=spec.gap,same=valid&&gap===0;
    var row={date:carbonDateString(t),millis:t,GE_available_pct:0,satellite_source:key,satellite_scale_m:spec.scale,
      satellite_date_basis:spec.kind==='vi'?'pixel observation date':spec.days===1?'observation date':'composite midpoint',
      satellite_status:valid?(same?'observed':'interpolated'):'unavailable',satellite_missing_pct:same?0:100,
      satellite_unavailable_pct:valid?0:100,satellite_gap_days:valid?gap:null,satellite_max_gap_days:spec.gap,
      satellite_left_date:valid?carbonDateString(l.millis):null,satellite_right_date:valid?carbonDateString(r.millis):null};
    if(valid) {
      var w=gap?(t-l.millis)/(r.millis-l.millis):0;
      if(carbonFinite(l.kNDVI)&&carbonFinite(r.kNDVI)){
        row.kNDVI=l.kNDVI+(r.kNDVI-l.kNDVI)*w;
        row.kndvi_sigma=l.kndvi_sigma+(r.kndvi_sigma-l.kndvi_sigma)*w;
      }
      row.red=l.red+(r.red-l.red)*w;row.nir=l.nir+(r.nir-l.nir)*w;
      if(carbonFinite(l.GE)&&carbonFinite(r.GE)){row.GE=l.GE+(r.GE-l.GE)*w;row.GE_available_pct=100;}
      row.NDVI=(row.nir-row.red)/(row.nir+row.red);
      row.EVI2=Math.max(0,Math.min(1,2.5*(row.nir-row.red)/(row.nir+2.4*row.red+1)));
    }
    rows.push(row);
  }
  return rows;
}

function aggregateModelMonths(rows) {
  if (!rows.length || rows[0].model_cadence === 'monthly') return rows;
  if(rows[0].model_cadence==='weekly'){
    var expanded=[];
    rows.forEach(function(r){var days=nativePeriodSeconds(r)/86400;
      for(var i=0;i<days;i++){
        var d=mergeParams(r,{});d.date=carbonDateString(carbonDateMillis(r.date)+i*86400000);
        d.period_days=1;d.period_seconds=86400;d.daylight_seconds=nativeDaylightSeconds(r)/days;d.model_cadence='daily';
        ['precip','GPP_period_gC_m2','NPP_period_gC_m2','ET_period_mm','T_period_mm','E_period_mm','I_period_mm'].forEach(function(k){d[k]=carbonFinite(r[k])?r[k]/days:null;});
        expanded.push(d);
      }
    });
    return aggregateModelMonths(expanded).map(function(r){r.model_cadence='weekly';r.aggregation='weekly simulation allocated by calendar-day overlap';return r;});
  }
  var groups={};rows.forEach(function(r){var k=r.date.slice(0,7);if(!groups[k])groups[k]=[];groups[k].push(r);});
  return Object.keys(groups).sort().map(function(k) {
    var list=groups[k], result={}, fields={};
    list.forEach(function(r){Object.keys(r).forEach(function(f){fields[f]=true;});});
    Object.keys(fields).forEach(function(f) {
      var numeric=list.filter(function(r){return carbonFinite(r[f]);});
      if(numeric.length) {
        var total=0,weight=0;
        var daylight=/^(PredGPP|GPPmax|Amax|PredLE|PT_mean|PARday|Rnday|WUE|gs$|rs$|cica)/.test(f);
        numeric.forEach(function(r){var w=daylight?(r.daylength||0):1;total+=r[f]*w;weight+=w;});
        result[f]=weight?total/weight:null;
      } else {
        var first=list[0][f];result[f]=list.every(function(r){return r[f]===first;})?first:'mixed';
        if(first===null||first===undefined)result[f]=null;
      }
    });
    result.date=k+'-01';result.millis=carbonDateMillis(result.date);result['system:time_start']=result.millis;
    result.year=Number(k.slice(0,4));result.month=Number(k.slice(5));
    result.period_seconds=list.reduce(function(a,r){return a+nativePeriodSeconds(r);},0);
    result.daylight_seconds=list.reduce(function(a,r){return a+nativeDaylightSeconds(r);},0);
    result.period_days=result.period_seconds/86400;
    result.period_start_millis=carbonDateMillis(list[0].date);result.period_end_millis=result.period_start_millis+result.period_seconds*1000;
    result.model_cadence='daily';result.output_cadence='monthly';
    result.aggregation='daily simulation; daytime flux means weighted by daylight duration';
    result.valid_model_days=list.filter(function(r){return carbonFinite(r.PredGPPmean);}).length;
    result.model_missing_pct=100*(1-result.valid_model_days/list.length);
    result.Tdaily_valid_days=list.filter(function(r){return carbonFinite(r.Tdaily);}).length;
    ['precip','GPPdaily','NPPdaily'].forEach(function(f) {
      var valid=list.filter(function(r){return carbonFinite(r[f]);});
      var sum=valid.reduce(function(a,r){return a+r[f];},0);
      var name=f==='precip'?'precip':f.replace('daily','_period_gC_m2');
      result[name]=valid.length===list.length?sum:null;
      result[name+'_available_sum']=valid.length?sum:null;
    });
    ['ET_period_mm','T_period_mm','E_period_mm','I_period_mm'].forEach(function(f){
      var valid=list.filter(function(r){return carbonFinite(r[f]);});
      result[f]=valid.length===list.length?valid.reduce(function(a,r){return a+r[f];},0):null;
    });
    result.precip_units='mm/represented period';
    return result;
  });
}
// Chart-only monthly means; daily values remain the source for exports and totals.
// One NPP point per native weather period; totals still use the daily mass allocation.
function nativeNPPChartRows(nativeRows,daily) {
  var byDate={};daily.forEach(function(r){byDate[r.date]=r;});
  return nativeRows.map(function(r){
    var out={date:r.date},a=carbonDateMillis(r.date),b=a+nativePeriodSeconds(r)*1000,days=[];
    for(var t=a;t<b;t+=86400000){var day=byDate[carbonDateString(t)];if(day)days.push(day);}
    ['NPPdaily','NPPgrass','NPPshrub','NPPtree'].forEach(function(k){
      var valid=days.filter(function(d){return carbonFinite(d[k]);});
      out[k]=valid.length?valid.reduce(function(sum,d){return sum+d[k];},0)/valid.length:null;
    });return out;
  });
}

// Observation means use valid tower records within each model period.
function towerPeriodObservation(tower,row,start,end) {
  var a=Math.max(carbonDateMillis(row.date),carbonDateMillis(start));
  var b=Math.min(carbonDateMillis(row.date)+nativePeriodSeconds(row)*1000,carbonDateMillis(end)+86400000);
  var days=tower.daily.filter(function(d){var t=carbonDateMillis(d.date);return t>=a&&t<b;});
  var m={sum:0,n:0,secondSum:0,secondN:0,daylight:0,leSum:0,leN:0,total:0,unknownLight:0,absent:0};
  days.forEach(function(d){Object.keys(m).forEach(function(k){m[k]+=carbonFinite(d[k])?d[k]:0;});});
  m.gpp=m.n?m.sum/m.n:null;m.secondGpp=m.secondN?m.secondSum/m.secondN:null;
  m.leET=m.leN?m.leSum/m.leN*86400/2450000:null;
  m.missingPct=m.daylight?100*(1-m.n/m.daylight):null;
  m.secondMissingPct=m.daylight?100*(1-m.secondN/m.daylight):null;
  var expected=Math.max(0,(b-a)/60000/tower.intervalMinutes);m.leMissingPct=expected?100*(1-m.leN/expected):null;
  return m;
}

function nativeBiomassChartRows(nativeRows,daily) {
  var byDate={};daily.forEach(function(r){byDate[r.date]=r;});
  var conversion=0.001*CARBON_UNITS.poundsPerKilogram*CARBON_UNITS.squareMetresPerAcre;
  return nativeRows.map(function(r){
    var a=carbonDateMillis(r.date),b=a+nativePeriodSeconds(r)*1000,days=[];
    for(var t=a;t<b;t+=86400000){var d=byDate[carbonDateString(t)];if(d)days.push(d);}
    var above=days.filter(function(d){return carbonFinite(d.HerbaceousAGBdaily);});
    var below=days.map(function(d){return carbonFinite(d.NPPgrass)&&carbonFinite(d.fBNPP)?d.NPPgrass*d.fBNPP*d.carbonToBiomass*conversion:null;}).filter(carbonFinite);
    return {date:r.date,HerbaceousAGBmean:above.length?above.reduce(function(sum,d){return sum+d.HerbaceousAGBdaily;},0)/above.length:null,
      HerbaceousBGBmean:below.length?below.reduce(function(sum,n){return sum+n;},0)/below.length:null,
      HerbaceousAGBcumulative_gC_m2:days.length?days[days.length-1].HerbaceousAGBcumulative_gC_m2:null};
  });
}

function monthlyNPPChartRows(daily) {
  var groups={};
  daily.forEach(function(r){var key=r.date.slice(0,7);if(!groups[key])groups[key]=[];groups[key].push(r);});
  return Object.keys(groups).sort().map(function(key){
    var out={date:key+'-01'};
    ['NPPdaily','NPPgrass','NPPshrub','NPPtree'].forEach(function(field){
      var valid=groups[key].filter(function(r){return carbonFinite(r[field]);});
      out[field]=valid.length?valid.reduce(function(sum,r){return sum+r[field];},0)/valid.length:null;
    });
    return out;
  });
}

function carbonSummaryRows(daily,annual) {
  var groups={};daily.forEach(function(r){var key=r.date.slice(0,annual?4:7);if(!groups[key])groups[key]=[];groups[key].push(r);});
  return Object.keys(groups).sort().map(function(k) {
    var rows=groups[k],year=Number(k.slice(0,4)),month=annual?0:Number(k.slice(5))-1;
    var first=Date.UTC(year,month,1),last=annual?Date.UTC(year+1,0,1):Date.UTC(year,month+1,1);
    var expected=(last-first)/86400000, out={date:carbonDateString(first),year:year,
      expected_days:expected,requested_days:rows.length,output_cadence:annual?'annual':'monthly'};
    ['GPPdaily','NPPdaily','NPPgrass','NPPshrub','NPPtree'].forEach(function(f) {
      var valid=rows.filter(function(r){return carbonFinite(r[f]);});
      var sum=valid.reduce(function(a,r){return a+r[f];},0),name=f.replace('daily','')+'_gC_m2';
      out[name]=valid.length===expected?sum:null;out[name+'_available_sum']=valid.length?sum:null;
      out[f+'_missing_pct']=100*(1-valid.length/rows.length);out[f+'_coverage_pct']=100*valid.length/expected;
    });
    return out;
  });
}

// Numeric counterpart of soilBucketPeriod. Used by the interactive app to avoid
// expanding capacity-search trajectories into a large Earth Engine expression.
function nativeBucketPeriod(sm,capacity,precip,potential) {
  var q=potential/capacity,surplus=precip-potential,gap=Math.max(0,capacity-sm)/Math.max(surplus,1e-12),x=q*gap;
  var logRatio=x<1e-5?1-x/2+x*x/3:Math.log(1+x)/Math.max(x,1e-30);
  var t=surplus>0?Math.min(1,gap*logRatio):1,z=q*t,decay=Math.exp(-z);
  var phi=z<1e-4?1-z/2+z*z/6-z*z*z/24:(1-decay)/Math.max(z,1e-30);
  var psi=z<1e-4?.5-z/6+z*z/24-z*z*z/120:(1-phi)/Math.max(z,1e-30);
  var mean=Math.max(0,Math.min(capacity,sm*t*phi+precip*t*t*psi+capacity*(1-t)));
  var evaporation=potential*mean/capacity,overflow=Math.max(0,surplus)*(1-t);
  return {meanSM:mean,fraction:mean/capacity,evaporation:evaporation,overflow:overflow,
    endSM:Math.max(0,Math.min(capacity,sm+precip-evaporation-overflow))};
}
function nativeBucketBridge(input) {
  var fields=['_bucket_rain','_bucket_demand','_bucket_surface','PT_mean','FIPAR','fwet_mean','fSM_mean',
    'PredLEsday_mean','PredLEday_mean','PredLEcday_mean','PredLEiday_mean'];
  var rows=input.map(function(r) {
    var f=mergeParams(r,{}),valid=['PT_mean','FIPAR','fwet_mean','daylength'].every(function(k){return carbonFinite(f[k]);});
    var demand=f.PT_mean*(1-f.FIPAR);
    f.period_seconds=nativePeriodSeconds(f);f.period_days=f.period_seconds/86400;
    f.daylight_seconds=nativeDaylightSeconds(f);
    f._bucket_rain=carbonFinite(f.precip)?f.precip/f.period_seconds:null;
    f._bucket_demand=valid?Math.max(0,demand)*(1-f.fwet_mean)*f.daylight_seconds/f.period_seconds/2450000:null;
    f._bucket_surface=valid?demand*f.fwet_mean:null;
    f.bucket_forcing_missing_pct=f.bucket_forcing_missing_pct>0||fields.some(function(k){return !carbonFinite(f[k]);})?100:0;
    return f;
  }).sort(function(a,b){return a.millis-b.millis;});
  fields.forEach(function(k) {
    var valid=rows.filter(function(r){return carbonFinite(r[k]);}),j=0;
    rows.forEach(function(r) {
      while(j<valid.length&&valid[j].millis<r.millis)j++;
      if(carbonFinite(r[k]))return;
      var l=valid[j-1],right=valid[j];
      if(l&&right)r[k]=l[k]+(right[k]-l[k])*(r.millis-l.millis)/(right.millis-l.millis);
    });
  });
  rows.forEach(function(r) {
    r.precip=carbonFinite(r._bucket_rain)?r._bucket_rain*r.period_seconds:null;
    r.bucket_fit_observed=r.bucket_forcing_missing_pct?0:1;
    r.bucket_forcing_status=!r.bucket_forcing_missing_pct?'supplied':fields.every(function(k){return carbonFinite(r[k]);})?'interpolated':'unavailable';
  });
  return rows;
}
function nativeBucketTrajectory(rows,capacity,history) {
  var sm=capacity,sse=0,weight=0,out=[];
  rows.forEach(function(r) {
    var conversion=r.daylight_seconds/2450000;
    var p=nativeBucketPeriod(sm,capacity,r.precip,r._bucket_demand*r.period_seconds);
    var soilLE=p.evaporation/conversion,totalLE=r._bucket_surface+soilLE;
    var w=r.PT_mean*(1-r.FIPAR)>0 && r.bucket_fit_observed===1?r.period_seconds:0;
    sse+=Math.pow(totalLE-r.PredLEsday_mean,2)*w;weight+=w;
    if(history)out.push(mergeParams(r,{SM:sm,SM_mean:p.meanSM,SM_end:p.endSM,SMmax:capacity,
      fSM_start:sm/capacity,SM_units:'mm effective evaporative storage',fSM_mean:p.fraction,
      fSM_RHVPD:r.fSM_mean,PredLEsday_RHVPD:r.PredLEsday_mean,PredLEday_RHVPD:r.PredLEday_mean,
      PredLEs_surface_mean:r._bucket_surface,PredLEs_soil_mean:soilLE,
      Esoil_period:p.evaporation,Esurface_period:r._bucket_surface*conversion,bucket_overflow:p.overflow,
      bucket_water_limited:r._bucket_demand>0&&p.fraction<1?1:0,
      SM_integration:'analytic period mean; uniform precipitation within period',
      PredLEsday_mean:totalLE,PredLEday_mean:r.PredLEcday_mean+r.PredLEiday_mean+totalLE}));
    sm=p.endSM;
  });
  return {mse:sse/Math.max(1,weight),weight:weight,history:out};
}
function nativeFitBucket(input,options) {
  var cfg=mergeParams(SOIL_BUCKET,options||{});
  if(!(cfg.minimum>0&&cfg.maximum>cfg.minimum&&cfg.initial>=cfg.minimum&&cfg.initial<=cfg.maximum&&cfg.gridSize>=3&&cfg.gridSize%1===0))throw new Error('Invalid bucket capacity settings.');
  var rows=nativeBucketBridge(input),required=['precip','period_days','daylength','fwet_mean','PT_mean','FIPAR','fSM_mean',
    'PredLEsday_mean','PredLEday_mean','PredLEcday_mean','PredLEiday_mean','_bucket_demand','_bucket_surface'];
  function valid(r){return required.every(function(k){return carbonFinite(r[k]);})&&r.precip>=0&&r.period_seconds>0&&r.daylight_seconds>0&&r.daylight_seconds<=r.period_seconds&&r.daylength>0&&r.daylength<=24&&r.fwet_mean>=0&&r.fwet_mean<=1&&r.FIPAR>=0&&r.FIPAR<=1&&r.fSM_mean>=0&&r.fSM_mean<=1;}
  var fitRows=rows.filter(valid),duration=fitRows.reduce(function(a,r){return a+r.period_days;},0);
  var eligible=fitRows.length&&duration>=365&&Math.abs((fitRows[fitRows.length-1].millis-fitRows[0].millis)/1000+fitRows[fitRows.length-1].period_seconds-duration*86400)<1e-6;
  var best={capacity:cfg.initial,mse:0},spread=0,informative=false,apply=false,allFull=false;
  function score(cap){return {capacity:cap,mse:nativeBucketTrajectory(fitRows,cap,false).mse};}
  if(eligible) {
    var lo=Math.log(cfg.minimum),hi=Math.log(cfg.maximum),step=(hi-lo)/(cfg.gridSize-1);
    best=score(cfg.initial);var min=Infinity,max=-Infinity;
    for(var i=0;i<cfg.gridSize;i++){var candidate=score(Math.exp(lo+i*step));min=Math.min(min,candidate.mse);max=Math.max(max,candidate.mse);if(candidate.mse<best.mse)best=candidate;}
    spread=max-min;
    for(var pass=0;pass<3;pass++){
      var center=Math.log(best.capacity),challenger=best;
      for(var j=-4;j<=4;j++){var refined=score(Math.exp(Math.max(lo,Math.min(hi,center+j*step/4))));if(refined.mse<challenger.mse)challenger=refined;}
      best=challenger;step/=4;
    }
    informative=spread>1e-8&&nativeBucketTrajectory(fitRows,best.capacity,false).weight>0;
    allFull=fitRows.every(function(r){return r.precip>=Math.max(0,r.PT_mean*(1-r.FIPAR))*(1-r.fwet_mean)*r.daylight_seconds/2450000;});
    apply=informative&&!allFull;
  }
  var boundary=best.capacity<=cfg.minimum*1.001||best.capacity>=cfg.maximum/1.001;
  var status=allFull?'unconstrained_full_bucket':!eligible?'baseline_requires_complete_contiguous_year':
    informative?(boundary?'fitted_at_search_boundary':'fitted'):'baseline_unidentifiable_capacity';
  var fitted=apply?nativeBucketTrajectory(fitRows,best.capacity,true).history:fitRows.map(function(r){
    var f=mergeParams(r,{fSM_RHVPD:r.fSM_mean,PredLEsday_RHVPD:r.PredLEsday_mean,PredLEday_RHVPD:r.PredLEday_mean});
    if(allFull){var potential=r.PT_mean*(1-r.FIPAR),full=r.fwet_mean*potential+(1-r.fwet_mean)*Math.max(0,potential);
      f.fSM_mean=1;f.fSM_start=1;f.PredLEsday_mean=full;f.PredLEday_mean=full+r.PredLEcday_mean+r.PredLEiday_mean;}
    return f;
  });
  var totalSeconds=rows.reduce(function(a,r){return a+r.period_seconds;},0);
  var missing=rows.reduce(function(a,r){return a+r.period_seconds*r.bucket_forcing_missing_pct;},0)/Math.max(1e-9,totalSeconds);
  return fitted.concat(rows.filter(function(r){return !valid(r);})).sort(function(a,b){return a.millis-b.millis;}).map(function(r){return mergeParams(r,{
    SM_fit_status:status,SM_fit_RMSE:apply?Math.sqrt(best.mse):null,SM_initial_guess:cfg.initial,SM_search_min:cfg.minimum,SM_search_max:cfg.maximum,
    SM_initialization:'full at first bracketed forcing; no spin-up',precip_units:'mm/period',bucket_record_missing_pct:missing});});
}
function fitNativeBuckets(input,params) {
  var p=mergeParams(ACCLIM_PARAMS,params||{}),near=nativeFitBucket(input,p.soilBucket);
  var proxy=near.map(function(r){var f={};['date','millis','precip','period_days','period_seconds','daylight_seconds','daylength','fwet_mean','PT_mean','fSM_mean','model_valid','bucket_forcing_missing_pct','bucket_fit_observed'].forEach(function(k){f[k]=r[k];});
    return mergeParams(f,{FIPAR:0,PredLEsday_mean:r.PredLEday_RHVPD,PredLEday_mean:r.PredLEday_RHVPD,PredLEcday_mean:0,PredLEiday_mean:0});});
  var et=nativeFitBucket(proxy,p.etBucket);
  return near.map(function(r,i){var e=et[i],hasFit=carbonFinite(e.SMmax),full=e.SM_fit_status==='unconstrained_full_bucket';
    function fitted(k){return hasFit?e[k]:null;}function moisture(k){return full?e[k]:fitted(k);}
    var result=mergeParams(r,{SM_et:fitted('SM'),SM_et_mean:fitted('SM_mean'),SM_et_end:fitted('SM_end'),SMmax_et:fitted('SMmax'),
      fSM_et_mean:moisture('fSM_mean'),fSM_et_start:moisture('fSM_start'),PredLE_ETeffective:moisture('PredLEsday_mean'),PredLE_ETreference:r.PredLEday_RHVPD,
      ETbucket_period:fitted('Esoil_period'),ETsurface_period:fitted('Esurface_period'),ETbucket_overflow:fitted('bucket_overflow'),
      SM_et_fit_status:e.SM_fit_status,SM_et_fit_RMSE:e.SM_fit_RMSE,SM_et_initial_guess:e.SM_initial_guess,SM_et_search_min:e.SM_search_min,
      SM_et_search_max:e.SM_search_max,SM_et_initialization:e.SM_initialization,SM_et_integration:fitted('SM_integration'),
      SM_et_units:'mm effective ET storage',ETbucket_period_units:'mm/period'});
    Object.keys(result).forEach(function(k){if(result[k]===undefined)result[k]=null;});
    return result;
  });
}

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
  var annual={},allocations={};
  annualRows.forEach(function(r){annual[r.year]=r;});
  var byDay=periodCarbonByDay(monthlyRows,start,end);
  var daily = [], windows = {}, periods = [];
  for (var ms=start;ms<=end;ms+=dayMs) {
    var date=new Date(ms),year=date.getUTCFullYear();
    var current=byDay[carbonDateString(ms)],gpp=current?current.gpp:null;
    var mode=current?current.mode:'missing_weather_period';
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
    if (mode === 'allocated_from_weather_period') w.held++;
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
    var weather=sampleWeather(geom,yearStart,yearEnd.advance(-1,'day'),config.weatherSource||'TERRACLIMATE');
    var validWeather=weather.filter(ee.Filter.notNull(['Tdaily'])).map(function(f){return f.set('_tdays',f.getNumber('Tdaily').multiply(f.getNumber('period_days')));});
    var days=ee.Number(validWeather.aggregate_sum('period_days'));
    var mat=ee.Algorithms.If(days.eq(yearEnd.difference(yearStart,'day')),ee.Number(validWeather.aggregate_sum('_tdays')).divide(days.max(1)),null);
    rows.push(ee.Feature(null,cover.combine(ee.Dictionary({MAT:mat,year:y,temperatureMonths:12,temperature_days:days,coverImages:coverIC.size()}))));
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
function fetchSoilGridsTexture(geom) {
  function depthMean(name) {
    var layers = ee.Image('projects/soilgrids-isric/' + name + '_mean');
    var valid = layers.updateMask(layers.gte(0).and(layers.lte(1000)));
    // SoilGrids g/kg to mass %, weighted by layer thickness over 0–30 cm.
    return valid.select(name + '_0-5cm_mean').multiply(5)
      .add(valid.select(name + '_5-15cm_mean').multiply(10))
      .add(valid.select(name + '_15-30cm_mean').multiply(15))
      .divide(300).rename(name + '_0_30cm_pct');
  }
  var sand = depthMean('sand');
  return sand.addBands(depthMean('silt')).addBands(depthMean('clay')).reduceRegion({
    reducer: ee.Reducer.mean(), geometry: geom, crs: sand.projection(),
    scale: 250, maxPixels: 1e10, tileScale: 4});
}

// Approximate endmembers read from "Soil spectra in VNIR region":
// https://lcluc.umd.edu/sites/default/files/lcluc_documents/Chirag.pdf
// Sand/silt/clay at 640 nm: 0.26/0.32/0.26; at 850 nm: 0.37/0.42/0.33.
function soilSpectralEstimate(texture) {
  var sand = texture.sand_0_30cm_pct, silt = texture.silt_0_30cm_pct, clay = texture.clay_0_30cm_pct;
  var valid = [sand, silt, clay].every(function(x) { return carbonFinite(x) && x >= 0 && x <= 100; });
  var total = sand + silt + clay;
  if (!valid || total <= 0) return null;
  var red = (0.26 * sand + 0.32 * silt + 0.26 * clay) / total;
  var nir = (0.37 * sand + 0.42 * silt + 0.33 * clay) / total;
  return {red: red, nir: nir, ndvi: (nir - red) / (nir + red),
    evi2: 2.5 * (nir - red) / (nir + 2.4 * red + 1)};
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
      if(variable==='Tdaily' && row.model_cadence==='daily' && row.Tdaily_valid_days!==n)complete=false;
      sum += row[variable] * (variable === 'Tdaily' ? n : 1);
      days += n;
    }
    return year + ': ' + (complete ? (variable === 'Tdaily' ? sum / days : sum).toFixed(1) +
      (variable === 'Tdaily' ? ' °C' : ' mm/year') : 'unavailable (incomplete year/data)');
  });
  return (variable === 'Tdaily' ? 'Mean annual air temperature (24-hour): ' : 'Annual precipitation: ') + values.join('; ');
}

// Matched-pixel baseline fit; retain the selected model reflectance source.
function geSourceSpec(key) {
  var s=SATELLITE_SOURCES[key];
  if(key.indexOf('MODIS')===0){var daily=key==='MODIS1';return {modis:true,daily:daily,scale:500,gap:32,
    ids:daily?['MODIS/061/MOD09GA','MODIS/061/MYD09GA']:['MODIS/061/MOD09A1','MODIS/061/MYD09A1'],
    label:daily?'MOD09GA / MYD09GA (500 m, daily)':'MOD09A1 / MYD09A1 (500 m, 8-day; actual observation dates)'};}
  return {modis:false,scale:s.scale,gap:s.gap,label:s.label};
}
function geDiagnosticBands(red,nir,green,blue,mask) {
  [red,nir,green,blue].forEach(function(b){mask=mask.and(b.gte(0)).and(b.lte(1)).and(b.mask());});
  mask=mask.and(red.add(nir).gt(0));
  return green.multiply(2).subtract(red).subtract(blue).rename('GE')
    .addBands(nir.subtract(red).multiply(2.5).divide(nir.add(red.multiply(2.4)).add(1)).rename('EVI2_GE'))
    .addBands(nir.subtract(red).divide(nir.add(red)).rename('NDVI_GE')).updateMask(mask);
}
function sampleGEDiagnostics(geom,start,end,key) {
  var spec=geSourceSpec(key),images;
  if(spec.modis){
    images=ee.ImageCollection([]);
    spec.ids.forEach(function(id){images=images.merge(ee.ImageCollection(id).filterBounds(geom).filterDate(start,ee.Date(end).advance(1,'day')).map(function(img){
      var d=ee.Date(img.get('system:time_start')),qa=img.select(spec.daily?'QC_500m':'QA'),state=img.select(spec.daily?'state_1km':'StateQA');
      var mask=clearModisState(state).and(state.rightShift(3).bitwiseAnd(7).eq(1))
        .and(state.rightShift(6).bitwiseAnd(3).lt(3)).and(qa.bitwiseAnd(3).eq(0)).and(qa.rightShift(2).bitwiseAnd(65535).eq(0));
      var out=geDiagnosticBands(img.select('sur_refl_b01').multiply(0.0001),img.select('sur_refl_b02').multiply(0.0001),
        img.select('sur_refl_b04').multiply(0.0001),img.select('sur_refl_b03').multiply(0.0001),mask);
      var day=spec.daily?ee.Image.constant(d.millis().divide(86400000)):
        img.select('DayOfYear').subtract(1).add(ee.Date.fromYMD(d.get('year'),1,1).millis().divide(86400000));
      return out.addBands(day.rename('observation_day')).set({sample_millis:d.millis()});
    }));});
  } else {
    images=satelliteImages(geom,ee.Date(start),ee.Date(end).advance(1,'day'),key).map(function(img){
      var red=img.select('red'),nir=img.select('nir'),ge=img.select('GE');
      return ge.addBands(nir.subtract(red).multiply(2.5).divide(nir.add(red.multiply(2.4)).add(1)).rename('EVI2_GE'))
        .addBands(nir.subtract(red).divide(nir.add(red)).rename('NDVI_GE'))
        .updateMask(ge.mask().and(red.mask()).and(nir.mask())).addBands(img.select('observation_day'))
        .set({sample_millis:img.get('sample_millis')});
    });
  }
  return ee.FeatureCollection(images.toList(images.size().max(1)).map(function(item){
    var img=ee.Image(item),d=ee.Date(img.get('sample_millis'));
    var count=spec.modis&&!spec.daily?ee.Date.fromYMD(d.get('year').add(1),1,1).difference(d,'day').min(8):ee.Number(1);
    return ee.List.sequence(0,count.subtract(1)).map(function(i){
      var t=d.advance(i,'day'),values=img.select(['GE','EVI2_GE','NDVI_GE']).updateMask(img.select('observation_day').eq(t.millis().divide(86400000)))
        .reduceRegion({reducer:ee.Reducer.mean(),geometry:geom,scale:spec.scale,maxPixels:1e10,tileScale:4});
      return ee.Feature(null,values).set({date:t.format('YYYY-MM-dd'),millis:t.millis()});
    });
  }).flatten()).filter(ee.Filter.notNull(['GE','EVI2_GE','NDVI_GE']));
}
// Same weights, dates and support for all three diagnostic quantities.
function aggregateGEDiagnostics(observations,periods,gapDays) {
  var fields=['GE','EVI2_GE','NDVI_GE'],groups={};
  observations.forEach(function(r){if(!fields.every(function(k){return carbonFinite(r[k]);}))return;
    var g=groups[r.millis]||(groups[r.millis]={millis:r.millis,n:0,GE:0,EVI2_GE:0,NDVI_GE:0});
    g.n++;fields.forEach(function(k){g[k]+=r[k];});
  });
  var obs=Object.keys(groups).map(function(k){var g=groups[k];fields.forEach(function(f){g[f]/=g.n;});return g;}).sort(function(a,b){return a.millis-b.millis;}),j=0;
  return periods.map(function(p){var out={date:p.date,GE:null,EVI2_GE:null,NDVI_GE:null,GE_available_pct:0,GE_observed_pct:0},n=0,direct=0,sums={GE:0,EVI2_GE:0,NDVI_GE:0};
    var a=carbonDateMillis(p.date),b=a+nativePeriodSeconds(p)*1000;
    for(var t=a;t<b;t+=86400000){
      while(j<obs.length&&obs[j].millis<t)j++;
      var r=obs[j],l=r&&r.millis===t?r:obs[j-1];
      if(!l||!r||(r.millis-l.millis)>gapDays*86400000)continue;
      var w=r.millis===l.millis?0:(t-l.millis)/(r.millis-l.millis);n++;if(l.millis===r.millis)direct++;
      fields.forEach(function(k){sums[k]+=l[k]+(r[k]-l[k])*w;});
    }
    if(n)fields.forEach(function(k){out[k]=sums[k]/n;});
    out.GE_available_pct=100*n/((b-a)/86400000);out.GE_observed_pct=100*direct/((b-a)/86400000);return out;
  });
}
// Normalize the paired selected-timestep observations over this run only.
function normalizeGE(rows) {
  var pairs=rows.filter(function(r){return carbonFinite(r.GE)&&carbonFinite(r.EVI2_GE)&&carbonFinite(r.NDVI_GE);});
  var lo=null,hi=null;
  pairs.forEach(function(r){lo=lo===null?r.GE:Math.min(lo,r.GE);hi=hi===null?r.GE:Math.max(hi,r.GE);});
  return rows.map(function(r){var d={};Object.keys(r).forEach(function(k){d[k]=r[k];});
    d.GE_min=lo;d.GE_max=hi;
    d.GE_normalization_applied=lo!==null&&lo<0;
    var valid=lo!==null&&carbonFinite(r.GE)&&carbonFinite(r.EVI2_GE)&&carbonFinite(r.NDVI_GE);
    d.GE_normalized=valid?(lo<0?(hi-lo>1e-12?(r.GE-lo)/(hi-lo):null):r.GE):null;
    return d;
  });
}
function geLinearFit(rows,field) {
  var pairs=normalizeGE(rows).filter(function(r){return carbonFinite(r.GE_normalized)&&carbonFinite(r[field]);}),n=pairs.length;
  if(n<3)return null;
  var x=pairs.reduce(function(a,r){return a+r.GE_normalized;},0)/n,y=pairs.reduce(function(a,r){return a+r[field];},0)/n,xx=0,xy=0,yy=0;
  pairs.forEach(function(r){xx+=Math.pow(r.GE_normalized-x,2);xy+=(r.GE_normalized-x)*(r[field]-y);yy+=Math.pow(r[field]-y,2);});
  if(xx<1e-14)return null;
  return {n:n,slope:xy/xx,intercept:y-xy/xx*x,r2:yy>0?xy*xy/(xx*yy):null,
    geMin:pairs[0].GE_min,geMax:pairs[0].GE_max,extrapolated:pairs.every(function(r){return r.GE_normalized>0;})||pairs.every(function(r){return r.GE_normalized<0;})};
}
function ndviMinGEFit(diagnostic,params) {
  var fit=geLinearFit(diagnostic,'NDVI_GE');
  if(!fit)throw new Error('Cannot estimate NDVImin at GE′=0: at least three valid timestep pairs with varying GE are required. Choose a longer date range or a different timestep/source.');
  if(!carbonFinite(fit.intercept)||fit.intercept < -1||fit.intercept > 1||fit.intercept>=params.ndvi_max)
    throw new Error('NDVI at GE′=0 ('+fit.intercept+') is outside [-1,1] or is not below NDVImax ('+params.ndvi_max+'). NDVImin was not applied.');
  return fit;
}
function fetchGEPeriods(geom,start,end,key,periods,requestId,done) {
  var spec=geSourceSpec(key),cursor=carbonDateMillis(start)-(spec.gap+32)*86400000;
  var stop=carbonDateMillis(end)+(spec.gap+32)*86400000,all=[];
  function batch(){
    if(requestId!==UI._requestId)return;
    try {
      if(cursor>stop){done(aggregateGEDiagnostics(all,periods,spec.gap));return;}
      var next=Math.min(cursor+32*86400000,stop+86400000);
      setStatus('Sampling GE/vegetation-index pairs: '+carbonDateString(cursor)+' to '+carbonDateString(next-86400000)+'...');
      var sample=sampleGEDiagnostics(geom,carbonDateString(cursor),carbonDateString(next-86400000),key);
      sample.toList(sample.size().max(1)).map(function(f){return ee.Feature(f).toDictionary();}).evaluate(function(result,error){
        if(requestId!==UI._requestId)return;
        if(error){done(null,'NDVImin estimation failed: '+error);return;}
        all=all.concat(result||[]);cursor=next;batch();
      });
    }catch(error){done(null,error.message);}
  }
  batch();
}

function addGEComparison(rows,geom,start,end,requestId,regionName,key) {
  rows=normalizeGE(rows);
  function chart(data){return ui.Chart({cols:[{label:'GE′',type:'number'},{label:'EVI2',type:'number'},{label:'NDVI',type:'number'}],
    rows:data.filter(function(r){return carbonFinite(r.GE)&&carbonFinite(r.EVI2_GE)&&carbonFinite(r.NDVI_GE);}).map(function(r){return {c:[{v:r.GE_normalized},{v:r.EVI2_GE},{v:r.NDVI_GE}]};})})
    .setChartType('ScatterChart').setOptions({title:regionChartTitle('EVI2 and NDVI vs GE′',regionName),chartArea:{top:80},
      hAxis:{title:'GE′: min–max scaled if GEmin < 0; otherwise GE'},vAxis:{title:'Vegetation index'},legend:{position:'top'},
      colors:['#8E24AA','#2E7D32'],pointSize:3,lineWidth:0,
      trendlines:{0:{type:'linear',color:'#8E24AA',lineWidth:1,visibleInLegend:false},1:{type:'linear',color:'#2E7D32',lineWidth:1,visibleInLegend:false}}});}
  addResultChart(chart(rows));
  var spec=geSourceSpec(key),fits=['EVI2_GE','NDVI_GE'].map(function(k){
    var f=geLinearFit(rows,k),label=k.replace('_GE','');
    return f?label+' at GE′=0: '+f.intercept.toFixed(4)+'; R²='+(f.r2===null?'undefined':f.r2.toFixed(3))+'; n='+f.n+(f.extrapolated?' (extrapolated)':''):label+': insufficient GE variation for a fit';
  });
  var applied=rows.length&&carbonFinite(rows[0].ndvi_min_used)?'Applied NDVImin = '+rows[0].ndvi_min_used.toFixed(4)+' ('+rows[0].ndvi_min_source+') for this run.':'NDVImin unavailable.';
  addChartNote(ui.Label({value:spec.label+'; matched pixels and dates at the selected timestep. GE′ uses min–max scaling only when the paired GE minimum is negative; otherwise GE′ = GE.\n'+fits.join('\n')+'\n'+applied+(rows.length&&carbonFinite(rows[0].evi2_min_p05)?' Historical EVI2min = '+rows[0].evi2_min_p05.toFixed(4)+'; '+rows[0].vi_baseline_n+' observed dates from '+rows[0].vi_baseline_first_observation+' to '+rows[0].vi_baseline_last_observation+' (requested '+rows[0].vi_baseline_start+' to '+rows[0].vi_baseline_end+').':' Fitted GE EVI2min is diagnostic; the configured EVI2min parameter is unchanged.'),style:{whiteSpace:'pre-wrap'}}));

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
  // Paired native model values; omit missing pairs rather than manufacturing zeros.
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
      error || !texture ? {} : texture, c4Chart);
  }
  // Texture failure leaves depths missing but does not suppress other outputs.
  try { fetchSoilGridsTexture(geom).evaluate(ready); } catch (error) { ready(null, error); }
}
function loadCarbonProductsWithTexture(output, geom, start, end, requestId, regionName, texture, c4Chart) {
  var context = fetchCarbonAnnualContext(geom, start, end, {weatherPlan:UI._runWeatherPlan,weatherSource:UI._lastRun.weatherSource,satelliteSource:UI._lastRun.satelliteSource});
  var monthlyProperties = output.toList(output.size().max(1)).map(function(f) { return ee.Feature(f).toDictionary(); });
  var annualProperties = context.toList(context.size().max(1)).map(function(f) { return ee.Feature(f).toDictionary(); });
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
        var depths = effectiveDepthFields(row, texture.sand_0_30cm_pct);
        row.silt_0_30cm_pct = carbonFinite(texture.silt_0_30cm_pct) ? texture.silt_0_30cm_pct : null;
        row.clay_0_30cm_pct = carbonFinite(texture.clay_0_30cm_pct) ? texture.clay_0_30cm_pct : null;
        Object.keys(depths).forEach(function(key) { row[key] = depths[key]; });
        var tDepth = transpirationSourceDepth([row]);
        Object.keys(tDepth).forEach(function(key) { row[key] = tDepth[key]; });
      });
      var nativeRows=result.monthly;
      var products = buildCarbonProducts(nativeRows, result.annual, start, end);
      var byDate={};nativeRows.forEach(function(r){byDate[r.date]=r;});
      UI._lastNative=carbonRowsToFC(nativeRows);
      var monthlyExportRows=aggregateModelMonths(nativeRows);
      result.monthly=nativeRows;
      result.monthly.forEach(function(row){var d=transpirationSourceDepth([row]);Object.keys(d).forEach(function(k){row[k]=d[k];});});
      output = carbonRowsToFC(result.monthly);
      UI._lastOutput = carbonRowsToFC(nativeRows); // selected timestep for charts and export
      UI._monthlyExportReady = true;
      var bucketInfo = {};
      var bucketRow=result.monthly.filter(function(r){return carbonFinite(r.PredGPPmean);})[0];
      if (bucketRow) Object.keys(bucketRow).forEach(function(key) {
        bucketInfo[key] = bucketRow[key];
      });
      var selectedPeriodDepth = transpirationSourceDepth(result.monthly);
      Object.keys(selectedPeriodDepth).forEach(function(key) {
        bucketInfo[key] = selectedPeriodDepth[key];
      });
      UI._lastDaily = carbonRowsToFC(products.daily.map(function(r){return nativeRows[0].model_cadence==='daily'?mergeParams(byDate[r.date]||{},r):r;}));
      UI._lastAnnual=carbonRowsToFC(carbonSummaryRows(products.daily,true));
      UI._lastCarbonMonthly=carbonRowsToFC(carbonSummaryRows(products.daily,false));
      if (UI._annualGPPNote) UI._annualGPPNote.setValue(annualCarbonNote(products.daily, ['GPPdaily'], ['GPP']));
      UI._last16day = carbonRowsToFC(products.npp16day);
      addOptionalSeriesChart({rows: nativeNPPChartRows(nativeRows,products.daily)},
        ['NPPdaily', 'NPPgrass', 'NPPshrub', 'NPPtree'],
        ['NPPeco', 'NPPgrass', 'NPPshrub', 'NPPtree'],
        ['#000000', '#2E7D32', '#EF6C00', '#1976D2'], [1, 2, 3], {
          title: regionChartTitle('NPP ('+(nativeRows[0].model_cadence||'monthly')+'; g C per day)', regionName), chartArea: {top: 80},
          hAxis: {title: 'Date'}, vAxis: {title: 'g C/m2/day'}, lineWidth: 1.5, pointSize: 0
        }, function(rows, fields, labels) { return annualCarbonNote(products.daily, fields, labels); });
      var biomassChartRows=nativeBiomassChartRows(nativeRows,products.daily);
      var hasBiomass=biomassChartRows.some(function(r){return carbonFinite(r.HerbaceousAGBmean)||carbonFinite(r.HerbaceousBGBmean);});
      if (hasBiomass) {
        var annualChartRows=biomassChartRows;
        addOptionalSeriesChart({rows: annualChartRows},
          ['HerbaceousAGBmean', 'HerbaceousAGBcumulative_gC_m2', 'HerbaceousBGBmean'],
          ['ANPP', 'Accumulated aboveground production (period end)', 'BNPP'], ['#2E7D32', '#2E7D32', '#000000'], [2], {
            title: regionChartTitle('Herbaceous Biomass Production', regionName), chartArea: {top: 80},
            hAxis: {title: 'Model period start'},
            vAxes: {0: {title: 'Mean Dry Matter Production (lbs/acre/day)'},
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
      if (carbonFinite(bucketInfo.SMmax_et)) addChartNote(ui.Label(annualClimateNote(monthlyExportRows, start, end, 'precip')));
      if (carbonFinite(bucketInfo.SMmax)) {
        addSoilMoistureChart(output, geom, requestId, regionName, bucketInfo);
      }
      c4Chart(result.monthly,null);
      var soil = soilSpectralEstimate(texture);
      addChartNote(ui.Label('Texture-based bare-soil estimate (0–30 cm): ' + (soil ?
        'R (640 nm) = ' + soil.red.toFixed(3) + ', NIR (850 nm) = ' + soil.nir.toFixed(3) +
        ', NDVImin = ' + soil.ndvi.toFixed(3) + ', EVI2min = ' + soil.evi2.toFixed(3) : 'unavailable')));
      if (bucketInfo.FAPAR_method_requested === 'EVI2') {
        addChartNote(ui.Label(bucketInfo.FAPAR_method === 'EVI2'
          ? 'FAPAR = EVI2 scaled between EVI2min=' + bucketInfo.FAPAR_EVI2min.toFixed(3) +
            ' and EVI2max=' + bucketInfo.FAPAR_EVI2max.toFixed(3) +
            ' (' + bucketInfo.FAPAR_fit_status + ')'
          : 'FAPAR: NDVI² fallback. Fit needs at least 3 paired months, a positive slope, and minimum below maximum.'));
      }
      addWeatherCharts(output, result.monthly, regionName);
      addGEComparison(result.monthly,geom,start,end,requestId,regionName,UI._lastRun.satelliteSource);
      var valid = products.daily.filter(function(r) { return carbonFinite(r.NPPdaily); }).length;
      var biomassDays = products.daily.filter(function(r) { return carbonFinite(r.HerbaceousAGBdaily); }).length;
      setStatus('Ready: ' + valid + '/' + products.daily.length + ' days have GPP/NPP; ' +
        biomassDays + ' have biomass estimates. Bucket forcing missing before interpolation: ' +
        (carbonFinite(bucketInfo.bucket_record_missing_pct)?bucketInfo.bucket_record_missing_pct.toFixed(2)+'%':'unavailable') +
        '. Incomplete windows have blank full totals in CSV.');
    } catch (e) { setStatus('Carbon products error: ' + e.message); }
  });
}
// END CARBON ASYNC UI


// Standalone Aqua MODIS LST-NDVI diagnostic: no ACCLIM forcing or physiology.
function lstNDVISettings(yearText,boxText) {
  var year=Number(yearText),km=Number(boxText),latest=new Date().getUTCFullYear();
  if(!/^\d{4}$/.test(String(yearText))||year<2002||year>latest)throw new Error('Choose a MODIS year from 2002 through '+latest+'.');
  if(!isFinite(km)||km<0.5||km>50)throw new Error('Box side must be between 0.5 and 50 km.');
  return {year:year,boxKm:km};
}
function lstNDVIBox(geom,km) {
  var center=geom.centroid(1),ll=center.coordinates();
  // Local UTM grid (polar stereographic outside UTM coverage), in metres.
  var lon=ee.Number(ll.get(0)),lat=ee.Number(ll.get(1));
  var zone=lon.add(180).divide(6).floor().add(1).max(1).min(60);
  var code=ee.Number(ee.Algorithms.If(lat.gte(0),32600,32700)).add(zone);
  code=ee.Number(ee.Algorithms.If(lat.gte(84),32661,ee.Algorithms.If(lat.lt(-80),32761,code)));
  var proj=ee.Projection(ee.String('EPSG:').cat(code.format('%d')));
  var xy=center.transform(proj,1).coordinates(),x=ee.Number(xy.get(0)),y=ee.Number(xy.get(1)),half=km*500;
  return ee.Geometry.Rectangle(ee.List([x.subtract(half),y.subtract(half),x.add(half),y.add(half)]),proj,false);
}
function lstNDVIReflectance(image,projection) {
  var red=image.select('sur_refl_b01').multiply(0.0001),nir=image.select('sur_refl_b02').multiply(0.0001);
  var state=image.select('StateQA'),qa=image.select('QA');
  var clear=state.bitwiseAnd(3).eq(0).and(state.bitwiseAnd(1<<2).eq(0))
    .and(state.rightShift(8).bitwiseAnd(3).eq(0)).and(state.bitwiseAnd(1<<10).eq(0))
    .and(state.bitwiseAnd(1<<12).eq(0)).and(state.bitwiseAnd(1<<13).eq(0)).and(state.bitwiseAnd(1<<15).eq(0));
  var good=clear.and(qa.bitwiseAnd(3).lte(1)).and(qa.rightShift(2).bitwiseAnd(15).eq(0))
    .and(qa.rightShift(6).bitwiseAnd(15).eq(0)).and(red.gte(0)).and(red.lte(1))
    .and(nir.gte(0)).and(nir.lte(1)).and(nir.add(red).gt(0));
  // Calculate NDVI at 500 m FIRST, then area-average valid NDVI within each LST cell.
  var ndvi=nir.subtract(red).divide(nir.add(red)).updateMask(good).rename('NDVI');
  var mean=ndvi.reduceResolution({reducer:ee.Reducer.mean(),maxPixels:64}).reproject(projection);
  var fraction=ndvi.mask().unmask(0).reduceResolution({reducer:ee.Reducer.mean(),maxPixels:64})
    .reproject(projection).rename('NDVI_valid_fraction');
  var reflectance=red.rename('red').addBands(nir.rename('nir')).updateMask(good)
    .reduceResolution({reducer:ee.Reducer.mean(),maxPixels:64}).reproject(projection);
  return mean.addBands(fraction).addBands(reflectance);
}
// QA choices apply only to the standalone diagnostic, not model thermal diagnostics.
function lstNDVIQualityLabel(quality) {
  return quality==='strict'?'Strict: good retrieval, LST error ≤1 K':'Valid retrieval, LST error ≤2 K';
}
function lstNDVIQualityMask(qc,raw,quality) {
  var strict=quality==='strict';
  return qc.bitwiseAnd(3).lte(strict?0:1)
    .and(qc.rightShift(2).bitwiseAnd(3).eq(0))
    .and(qc.rightShift(6).bitwiseAnd(3).lte(strict?0:1)).and(raw.gt(0));
}
function sampleLSTNDVIPeriods(geom,year,boxKm,firstMonth,lastMonth,quality) {
  var box=lstNDVIBox(geom,boxKm),start=ee.Date.fromYMD(year,firstMonth,1);
  var end=lastMonth===12?ee.Date.fromYMD(year+1,1,1):ee.Date.fromYMD(year,lastMonth+1,1);
  var reflectance=ee.ImageCollection('MODIS/061/MYD09A1').filterDate(start,end).filterBounds(box);
  var lstSource=ee.ImageCollection('MODIS/061/MYD11A2');
  var projection=ee.Image(lstSource.first()).select('LST_Day_1km').projection();
  var lst=lstSource.filterDate(start,end).filterBounds(box);
  var joined=ee.ImageCollection(ee.Join.saveFirst('paired_lst').apply(reflectance,lst,
    ee.Filter.equals({leftField:'system:time_start',rightField:'system:time_start'})));
  var grid=ee.Image.pixelCoordinates(projection).rename(['pixel_x','pixel_y']);
  var pixels=grid.sample({region:box,projection:projection,geometries:false,tileScale:4}).size();
  var samples=ee.FeatureCollection(joined.toList(joined.size().max(1)).map(function(item){
    var image=ee.Image(item),thermal=ee.Image(image.get('paired_lst'));
    var qc=thermal.select('QC_Day'),raw=thermal.select('LST_Day_1km');
    var good=lstNDVIQualityMask(qc,raw,quality);
    var temperature=raw.multiply(0.02).subtract(273.15).updateMask(good).rename('LST_C');
    var aligned=lstNDVIReflectance(image,projection).addBands(temperature).addBands(qc.rename('LST_QC')).addBands(grid);
    var date=image.date(),stop=ee.Date(date.advance(8,'day').millis().min(ee.Date.fromYMD(year+1,1,1).millis()));
    return aligned.sample({region:box,projection:projection,geometries:false,tileScale:4,dropNulls:true}).map(function(f){
      return f.set({date:date.format('YYYY-MM-dd'),end_exclusive:stop.format('YYYY-MM-dd'),month:date.get('month'),
        pixel_id:f.getNumber('pixel_x').format('%.1f').cat('_').cat(f.getNumber('pixel_y').format('%.1f')),
        box_km:boxKm,year:year,LST_quality:quality||'valid_2k',LST_source:'MODIS/061/MYD11A2',NDVI_source:'MODIS/061/MYD09A1'});
    });
  })).flatten();
  return ee.Dictionary({rows:samples.toList(samples.size().max(1)).map(function(f){return ee.Feature(f).toDictionary();}),
    pixels:pixels,paired_periods:joined.size(),reflectance_periods:reflectance.size(),lst_periods:lst.size(),
    pixel_size_m:projection.nominalScale()});
}
function lstNDVIMonthFit(rows,month) {
  var points=rows.filter(function(r){return r.month===month&&carbonFinite(r.NDVI)&&carbonFinite(r.LST_C);});
  if(points.length<2)return null;
  var mx=0,my=0;points.forEach(function(r){mx+=r.NDVI;my+=r.LST_C;});mx/=points.length;my/=points.length;
  var xx=0,xy=0,maxX=-Infinity;
  points.forEach(function(r){xx+=(r.NDVI-mx)*(r.NDVI-mx);xy+=(r.NDVI-mx)*(r.LST_C-my);maxX=Math.max(maxX,r.NDVI);});
  if(xx<=1e-12)return null;
  var slope=xy/xx,intercept=my-slope*mx,lo=0,hi=maxX;
  if(hi<=0)return null;
  // Draw the fitted span and extrapolate toward NDVI=0, clipped at LST=60 C.
  if(Math.abs(slope)<1e-12){if(intercept>60)return null;}
  else if(slope<0)lo=Math.max(lo,(60-intercept)/slope);
  else hi=Math.min(hi,(60-intercept)/slope);
  if(hi<=lo)return null;
  return {month:month,n:points.length,slope:slope,intercept:intercept,
    endpoints:[{x:lo,y:Math.min(60,intercept+slope*lo)},{x:hi,y:Math.min(60,intercept+slope*hi)}]};
}
function lstNDVIChartTable(rows,period) {
  var months=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'],all=period==='all';
  var selected=rows.filter(function(r){return (all||r.date===period)&&carbonFinite(r.NDVI)&&carbonFinite(r.LST_C);});
  var monthIds=all?[1,2,3,4,5,6,7,8,9,10,11,12]:[Number(period.slice(5,7))];
  var cols=[{id:'NDVI',label:'NDVI (500 m values averaged within LST pixel)',type:'number'}];
  monthIds.forEach(function(m){cols.push({label:all?months[m-1]:'Daytime LST',type:'number'});});
  monthIds.forEach(function(m){cols.push({label:months[m-1]+' monthly linear fit',type:'number'});});
  function cells(x){var c=[{v:x}];for(var i=0;i<monthIds.length*2;i++)c.push({v:null});return c;}
  var result=selected.map(function(r){var c=cells(r.NDVI);c[1+monthIds.indexOf(r.month)]={v:r.LST_C};return {c:c};});
  monthIds.forEach(function(m,index){var fit=lstNDVIMonthFit(rows,m);if(!fit)return;
    fit.endpoints.forEach(function(p){var c=cells(p.x);c[1+monthIds.length+index]={v:p.y};result.push({c:c});});
  });
  result.sort(function(a,b){return a.c[0].v-b.c[0].v;});
  return {cols:cols,rows:result};
}
function lstNDVIFitSeries(period) {
  var colors=['#313695','#4575b4','#74add1','#abd9e9','#66bd63','#1a9850','#fee08b','#fdae61','#f46d43','#d73027','#a50026','#762a83'];
  var all=period==='all',n=all?12:1,series={};
  for(var i=0;i<n;i++){
    var color=colors[all?i:Number(period.slice(5,7))-1];
    series[i]={color:color,lineWidth:0,pointSize:2};
    series[n+i]={color:color,lineWidth:2,pointSize:0,visibleInLegend:false};
  }
  return series;
}

// Percentiles use linear interpolation at (n - 1) * p (R type 7).
function lstNDVIStats(values) {
  var sorted=values.filter(carbonFinite).sort(function(a,b){return a-b;});
  if(!sorted.length)return {n:0,mean:null,p05:null,p95:null};
  function quantile(p){var pos=(sorted.length-1)*p,lo=Math.floor(pos),hi=Math.ceil(pos);return sorted[lo]+(sorted[hi]-sorted[lo])*(pos-lo);}
  return {n:sorted.length,mean:sorted.reduce(function(a,b){return a+b;},0)/sorted.length,p05:quantile(.05),p95:quantile(.95)};
}
function lstNIRRedChartTable(rows,period) {
  var months=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'],all=period==='all';
  var cols=[{label:'Red reflectance',type:'number'}];
  (all?months:['NIR reflectance']).forEach(function(label){cols.push({label:label,type:'number'});});
  return {cols:cols,rows:rows.filter(function(r){return (all||r.date===period)&&carbonFinite(r.NDVI)&&carbonFinite(r.LST_C)&&carbonFinite(r.red)&&carbonFinite(r.nir);}).map(function(r){
    var cells=[{v:r.red}];for(var m=1;m<=(all?12:1);m++)cells.push({v:!all||r.month===m?r.nir:null});return {c:cells};
  })};
}
function lstNDVITimeTable(rows,year,field) {
  var groups={};rows.forEach(function(r){if(carbonFinite(r.NDVI)&&carbonFinite(r.LST_C)){if(!groups[r.date])groups[r.date]=[];groups[r.date].push(r[field]);}});
  var table={cols:[{label:'Composite start',type:'date'},{label:'Mean',type:'number'},
    {label:'5th percentile',type:'number'},{label:'95th percentile',type:'number'}],rows:[]};
  for(var t=Date.UTC(year,0,1);t<Date.UTC(year+1,0,1);t+=8*86400000){
    var key=new Date(t).toISOString().slice(0,10),stat=lstNDVIStats(groups[key]||[]),p=key.split('-');
    table.rows.push({c:[{v:new Date(Number(p[0]),Number(p[1])-1,Number(p[2]))},{v:stat.mean},{v:stat.p05},{v:stat.p95}]});
  }
  return table;
}

function runLSTNDVI(geom) {
  var settings;
  try{settings=lstNDVISettings(UI.lstNDVIYear.getValue(),UI.lstNDVIBoxKm.getValue());}
  catch(error){setStatus(error.message);return;}
  var quality=UI.lstNDVIQuality.getValue();
  var request=UI._requestId,rows=[],periods=0,available=0,pixelCount=0,month=1,pixelSize=0;
  UI.runButton.setDisabled(true);
  var box=lstNDVIBox(geom,settings.boxKm);
  if(UI.map){
    UI.map.centerObject(box);
    if(UI._lstNDVILayer)UI.map.layers().remove(UI._lstNDVILayer);
    UI._lstNDVILayer=ui.Map.Layer(box,{color:'FF8800'},'LST-NDVI sampling box');UI.map.layers().add(UI._lstNDVILayer);
  }
  function finish(){
    UI.runButton.setDisabled(false);
    if(!rows.length){
      var reason=!pixelCount?'No LST grid-cell centers inside this box; increase the box size.':
        !periods?'No paired product intervals available for this year and box. Aqua records begin in July 2002.':
        periods+' paired product intervals and '+pixelCount+' grid cells were found, but no pixels passed both reflectance and LST screening ('+lstNDVIQualityLabel(quality)+').'+
          (quality==='strict'?' Try the valid-retrieval ≤2 K option.':' Try another year or a larger box.');
      setStatus('No paired, quality-screened LST-NDVI pixels. '+reason);return;
    }
    UI.resultsPanel.style().set({shown:true});UI.resultsOverlay.setExpanded(true);
    function chart(period){return ui.Chart(lstNDVIChartTable(rows,period)).setChartType('LineChart').setOptions({
      title:'Aqua MODIS daytime LST vs NDVI — '+(period==='all'?settings.year:period)+'\n'+settings.boxKm+' km square; '+pixelCount+' LST grid cells',
      hAxis:{title:'NDVI (mean within native 1 km LST pixel)'},vAxis:{title:'Daytime LST (°C)'},
      pointSize:2,lineWidth:0,interpolateNulls:true,series:lstNDVIFitSeries(period),legend:{position:period==='all'?'bottom':'none'},
      colors:period==='all'?['#313695','#4575b4','#74add1','#abd9e9','#66bd63','#1a9850','#fee08b','#fdae61','#f46d43','#d73027','#a50026','#762a83']:['#1976D2']});}
    var dates={};rows.forEach(function(r){dates[r.date]=true;});
    var annualAllowed=rows.length<=20000,initial=annualAllowed?'all':Object.keys(dates).sort()[0];
    addResultChart(chart(initial));
    var replace=chartReplacer();
    UI._chartCard.add(ui.Select({items:(annualAllowed?[{label:'All 8-day intervals (colored by month)',value:'all'}]:[]).concat(Object.keys(dates).sort().map(function(d){return {label:d,value:d};})),
      value:initial,onChange:function(d){replace(chart(d));if(replaceReflectance)replaceReflectance(reflectanceChart(d));},style:{stretch:'horizontal'}}));
    if(!annualAllowed)addChartNote(ui.Label('This box has more than 20,000 valid pairs. Select an 8-day interval to view all its pixels, or reduce the box for an annual view. No pairs were subsampled.'));
    var partial=rows.filter(function(r){return r.NDVI_valid_fraction<0.999;}).length;
    addChartNote(ui.Label({value:rows.length+' valid pixel–interval pairs across '+Object.keys(dates).length+' intervals. '+pixelCount+' grid-cell centers in box; native spacing '+pixelSize.toFixed(1)+' m.\n'+
      'NDVI is the area-weighted mean of quality-screened 500 m NDVI, not NDVI of averaged reflectance. '+partial+' pairs have partial NDVI coverage. Reflectance cloud, shadow and snow flags are excluded; no gap filling.\n'+
      'MYD09A1 selected observations and MYD11A2 8-day means share the composite interval, not necessarily the observation day. LST QA: '+lstNDVIQualityLabel(quality)+'. Missing/cloud retrievals and nonzero data-quality flags are excluded. Dates are composite starts.',style:{whiteSpace:'pre-wrap'}}));
    function reflectanceChart(period){
      return ui.Chart(lstNIRRedChartTable(rows,period)).setChartType('ScatterChart').setOptions({
        title:'NIR vs red reflectance — '+(period==='all'?settings.year:period),
        hAxis:{title:'Red surface reflectance (fraction)'},vAxis:{title:'NIR surface reflectance (fraction)'},
        pointSize:2,lineWidth:0,legend:{position:period==='all'?'bottom':'none'},
        colors:period==='all'?['#313695','#4575b4','#74add1','#abd9e9','#66bd63','#1a9850','#fee08b','#fdae61','#f46d43','#d73027','#a50026','#762a83']:['#1976D2']});
    }
    addResultChart(reflectanceChart(initial));
    var replaceReflectance=chartReplacer();
    addChartNote(ui.Label('NIR versus red uses the same valid pixel–interval pairs and follows the LST–NDVI interval selector. Both reflectance bands are area-averaged over the same quality-screened 500 m subpixels within each LST cell.'));
    [{field:'NDVI',name:'NDVI',units:'NDVI',color:'#2E7D32'},
      {field:'LST_C',name:'Daytime LST',units:'Daytime LST (°C)',color:'#D84315'}].forEach(function(spec){
      addResultChart(ui.Chart(lstNDVITimeTable(rows,settings.year,spec.field)).setChartType('LineChart').setOptions({
        title:spec.name+' time series — '+settings.year+'\nMean and spatial 5th–95th percentiles of paired pixels',
        hAxis:{title:'8-day composite start'},vAxis:{title:spec.units},legend:{position:'bottom'},
        series:{0:{color:spec.color,lineWidth:2,pointSize:3},1:{color:'#757575',lineDashStyle:[4,4],lineWidth:1},
          2:{color:'#757575',lineDashStyle:[4,4],lineWidth:1}},interpolateNulls:false}));
    });
    addChartNote(ui.Label('Monthly ordinary least-squares fits use all valid paired pixels in each calendar month. Lines extend toward NDVI = 0 and stop at 60 °C if reached first. An individual 8-day view retains its full-month fit. Months without NDVI variation have no fit.'));
    UI._lstNDVIRows=rows;
    setStatus('LST-NDVI ready: '+rows.length+' pairs; '+periods+' paired product intervals ('+available+' reflectance intervals available). No ACCLIM model calculations were run.');
  }
  function next(){
    if(request!==UI._requestId)return;
    if(month>12){finish();return;}
    setStatus('Sampling LST-NDVI: '+settings.year+', months '+month+'–'+Math.min(12,month+2)+'…');
    sampleLSTNDVIPeriods(geom,settings.year,settings.boxKm,month,Math.min(12,month+2),quality).evaluate(function(data,error){
      if(request!==UI._requestId)return;
      if(error||!data){UI.runButton.setDisabled(false);setStatus('LST-NDVI error: '+(error||'No response'));return;}
      rows=rows.concat(data.rows||[]);periods+=data.paired_periods;available+=data.reflectance_periods;pixelCount=data.pixels;pixelSize=data.pixel_size_m;
      month+=3;next();
    });
  }
  next();
}

// Standalone soil-reference search. Empirical support does not verify soil identity.
function soilQuantile(values,q) {
  var v=values.slice().sort(function(a,b){return a-b;});
  if(!v.length)return null;
  var p=(v.length-1)*q,i=Math.floor(p);return v[i]+(v[Math.min(i+1,v.length-1)]-v[i])*(p-i);
}
function soilEVI2(red,nir){return 2.5*(nir-red)/(nir+2.4*red+1);}
function soilLineEVI2(red,fit){return soilEVI2(red,fit.a*red+fit.b);}
function soilLineNDVI(red,fit){var nir=fit.a*red+fit.b,total=nir+red;return total>0?(nir-red)/total:null;}
function soilFit(rows,q) {
  if(rows.length<30)return null;
  var counts={};rows.forEach(function(r){counts[r.month]=(counts[r.month]||0)+1;});
  var months=Object.keys(counts).length;
  function profile(a){
    var v=rows.map(function(r){return {z:r.nir-a*r.red,w:1/(months*counts[r.month])};});
    v.sort(function(x,y){return x.z-y.z;});
    var cumulative=0,b=v[v.length-1].z;
    for(var i=0;i<v.length;i++){cumulative+=v[i].w;if(cumulative>=q){b=v[i].z;break;}}
    var loss=0;v.forEach(function(r){var e=r.z-b;loss+=r.w*(e>=0?q*e:(q-1)*e);});
    return {a:a,b:b,loss:loss};
  }
  // Convex profile loss; broad numerical bracket, with boundary solutions rejected.
  var lo=-20,hi=20,g=(Math.sqrt(5)-1)/2,c=hi-g*(hi-lo),d=lo+g*(hi-lo),fc=profile(c),fd=profile(d);
  for(var j=0;j<42;j++){
    if(fc.loss<fd.loss){hi=d;d=c;fd=fc;c=hi-g*(hi-lo);fc=profile(c);}
    else{lo=c;c=d;fc=fd;d=lo+g*(hi-lo);fd=profile(d);}
  }
  var fit=profile((lo+hi)/2);return Math.abs(fit.a)>19.99?null:fit;
}
function soilPrepare(rows,cut){
  var selected=rows.filter(function(r){return carbonFinite(r.red)&&carbonFinite(r.nir)&&r.red>=0&&r.red<=1&&r.nir>=0&&r.nir<=1&&r.red+r.nir>0&&soilEVI2(r.red,r.nir)<=cut;});
  var lo=soilQuantile(selected.map(function(r){return r.red;}),.01),hi=soilQuantile(selected.map(function(r){return r.red;}),.99);
  var core=selected.filter(function(r){return r.red>=lo&&r.red<=hi;}),counts={};
  core.forEach(function(r){counts[r.month]=(counts[r.month]||0)+1;});
  var months=Object.keys(counts).filter(function(m){return counts[m]>=30;});
  core=core.filter(function(r){return months.indexOf(String(r.month))>=0;});
  // Report only the brightness span actually supported by eligible months.
  if(core.length){lo=Math.min.apply(null,core.map(function(r){return r.red;}));hi=Math.max.apply(null,core.map(function(r){return r.red;}));}
  return {rows:core,lo:lo,hi:hi,months:months,selected:selected.length};
}
function soilAssessLegacy(rows,settings){
  var p=soilPrepare(rows,settings.cut),out={total:rows.length,n:p.rows.length,months:p.months.length,lo:p.lo,hi:p.hi,reasons:[],supported:false};
  if(p.rows.length<100||p.months.length<3||p.hi-p.lo<settings.range){out.reasons.push('Insufficient samples, months or red range');return out;}
  var fit=soilFit(p.rows,.05);if(!fit){out.reasons.push('Numerical slope limit');return out;}out.fit=fit;
  if(fit.a<=0)out.reasons.push('Nonpositive slope');
  function departure(other){return Math.max(Math.abs((other.a-fit.a)*p.lo+other.b-fit.b),Math.abs((other.a-fit.a)*p.hi+other.b-fit.b));}
  var q2=soilFit(p.rows,.02),q10=soilFit(p.rows,.10);
  out.quantileShift=q2&&q10?Math.max(departure(q2),departure(q10)):Infinity;
  if(out.quantileShift>settings.tolerance)out.reasons.push('Quantile sensitivity');
  out.omissionShift=0;
  p.months.forEach(function(m){var v=p.rows.filter(function(r){return String(r.month)!==m;}),f=soilFit(v,.05);out.omissionShift=Math.max(out.omissionShift,f?departure(f):Infinity);});
  if(out.omissionShift>settings.tolerance)out.reasons.push('Month sensitivity');
  out.curvature=0;out.bins=0;
  for(var k=0;k<6;k++){
    var left=p.lo+(p.hi-p.lo)*k/6,right=p.lo+(p.hi-p.lo)*(k+1)/6;
    var residuals=p.rows.filter(function(r){return r.red>=left&&(k===5?r.red<=right:r.red<right);}).map(function(r){return r.nir-fit.a*r.red-fit.b;});
    if(residuals.length>=20){out.bins++;out.curvature=Math.max(out.curvature,Math.abs(soilQuantile(residuals,.05)));}
  }
  if(out.bins<4)out.reasons.push('Insufficient brightness-bin support');
  if(out.curvature>settings.tolerance)out.reasons.push('Curved lower boundary');
  out.thresholdShift=0;
  [settings.cut-.05,settings.cut+.05].forEach(function(cut){
    var v=soilPrepare(rows,cut),f=v.rows.length>=100&&v.months.length>=3&&v.hi-v.lo>=settings.range?soilFit(v.rows,.05):null;
    out.thresholdShift=Math.max(out.thresholdShift,f?departure(f):Infinity);
  });
  if(out.thresholdShift>settings.tolerance)out.reasons.push('EVI2 threshold sensitivity or insufficient support');
  // PVI=0 exactly: NIR=aR+b. Uniform Red grid avoids random-draw regression bias.
  var values=[],ndviValues=[];for(var t=0;t<=1000;t++){var red=p.lo+(p.hi-p.lo)*t/1000;values.push(soilLineEVI2(red,fit));ndviValues.push(soilLineNDVI(red,fit));}
  out.eviMin=soilQuantile(values,.5);out.eviLow=Math.min.apply(null,values);out.eviHigh=Math.max.apply(null,values);
  if(ndviValues.every(carbonFinite)){out.ndviMin=soilQuantile(ndviValues,.5);out.ndviLow=Math.min.apply(null,ndviValues);out.ndviHigh=Math.max.apply(null,ndviValues);}
  else{out.ndviMin=null;out.ndviLow=null;out.ndviHigh=null;out.reasons.push('Undefined NDVI along reference line');}
  if(fit.a*p.lo+fit.b<0||fit.a*p.hi+fit.b>1||out.eviMin<0||out.eviMin>settings.cut)out.reasons.push('Implausible reference reflectance or EVI2');
  out.supported=out.reasons.length===0;return out;
}
// Pooled geometric soil-edge candidates. All thresholds are research diagnostics.
function soilLandReflectance(image,projection){
  // StateQA 3–5: retain land (1), exclude water and mixed shorelines BEFORE averaging.
  return lstNDVIReflectance(image.updateMask(image.select('StateQA').rightShift(3).bitwiseAnd(7).eq(1)),projection);
}
function soilOLS(points){
  if(points.length<3)return null;
  var mx=0,my=0;points.forEach(function(p){mx+=p.x;my+=p.y;});mx/=points.length;my/=points.length;
  var xx=0,xy=0,yy=0;points.forEach(function(p){xx+=(p.x-mx)*(p.x-mx);xy+=(p.x-mx)*(p.y-my);yy+=(p.y-my)*(p.y-my);});
  if(xx<1e-14)return null;var a=xy/xx,b=my-a*mx;return {a:a,b:b,r2:yy>0?xy*xy/(xx*yy):null,n:points.length};
}
function soilRobustAnchors(points){
  var slopes=[];for(var i=0;i<points.length;i++)for(var j=i+1;j<points.length;j++)if(Math.abs(points[i].x-points[j].x)>1e-8)slopes.push((points[j].y-points[i].y)/(points[j].x-points[i].x));
  if(!slopes.length)return null;var a=soilQuantile(slopes,.5),b=soilQuantile(points.map(function(p){return p.y-a*p.x;}),.5);return {a:a,b:b};
}
function soilSimulation(fit,settings){
  var bounds=settings.simBounds||[0,1,0,1],seed=173,points=[];
  function random(){seed=(1664525*seed+1013904223)%4294967296;return seed/4294967296;}
  // Common deterministic uniform draws for all candidate geometries.
  for(var i=0;i<8192;i++){
    var red=bounds[0]+(bounds[1]-bounds[0])*random(),nir=bounds[2]+(bounds[3]-bounds[2])*random();
    var pvi=(nir-fit.a*red-fit.b)/Math.sqrt(1+fit.a*fit.a);
    if(pvi>=0&&nir+red>0)points.push({x:pvi,y:soilEVI2(red,nir),ndvi:(nir-red)/(nir+red)});
  }
  if(points.length<100)return null;
  var evi=soilOLS(points),ndvi=soilOLS(points.map(function(p){return {x:p.x,y:p.ndvi};}));
  if(!evi||!ndvi)return null;
  var fits=[evi],ndviFits=[ndvi];[.5,.75].forEach(function(q){var limit=soilQuantile(points.map(function(p){return p.x;}),q),subset=points.filter(function(p){return p.x<=limit;});fits.push(soilOLS(subset));ndviFits.push(soilOLS(subset.map(function(p){return {x:p.x,y:p.ndvi};})));});
  return {evi:evi,ndvi:ndvi,points:points,bounds:bounds,horizons:fits,ndviHorizons:ndviFits,
    horizonShift:Math.max.apply(null,fits.map(function(f){return Math.abs(f.b-evi.b);})),ndviHorizonShift:Math.max.apply(null,ndviFits.map(function(f){return Math.abs(f.b-ndvi.b);}))};
}
function soilGeometryCandidate(rows,lo,hi,bins,q,segment,tol){
  var groups=[],anchors=[],k;
  for(k=0;k<bins;k++)groups.push([]);
  rows.forEach(function(r){groups[Math.min(bins-1,Math.floor((r.red-lo)/(hi-lo)*bins))].push(r);});
  groups.forEach(function(g,index){if(g.length>=20)anchors.push({x:soilQuantile(g.map(function(r){return r.red;}),.5),y:soilQuantile(g.map(function(r){return r.nir;}),q),bin:index,n:g.length});});
  var used=anchors.filter(function(p){return segment==='all'||(segment==='bright'?p.bin>=Math.floor(bins/3):p.bin<Math.ceil(2*bins/3));});
  if(used.length<4)return null;
  var fit=soilRobustAnchors(used);if(!fit)return null;
  // Detrend before taking quantiles so bin width does not bias the NIR anchor.
  for(k=0;k<2;k++){
    used=used.map(function(p){var g=groups[p.bin];return {x:p.x,y:fit.a*p.x+fit.b+soilQuantile(g.map(function(r){return r.nir-fit.a*r.red-fit.b;}),q),bin:p.bin,n:p.n};});
    fit=soilRobustAnchors(used);
  }
  if(!fit||fit.a<=0||fit.a>20)return null;
  var support=[],supportBins=[],dates={},pixels={},belowFractions=[],edgeTol=Math.max(.002,tol/2);
  groups.forEach(function(g,index){
    if(g.length<20)return;
    var near=g.filter(function(r){return Math.abs(r.nir-fit.a*r.red-fit.b)<=edgeTol;});
    belowFractions.push(g.filter(function(r){return r.nir-fit.a*r.red-fit.b < -edgeTol;}).length/g.length);
    if(near.length>=3){supportBins.push(index);support=support.concat(near);}
  });
  support.forEach(function(r){if(r.date)dates[r.date]=true;if(r.pixel_x!==undefined&&r.pixel_y!==undefined)pixels[r.pixel_x+':'+r.pixel_y]=true;});
  var edgeLo=support.length?soilQuantile(support.map(function(r){return r.red;}),.05):lo;
  var edgeHi=support.length?soilQuantile(support.map(function(r){return r.red;}),.95):hi;
  var span=(edgeHi-edgeLo)/(hi-lo),below=belowFractions.reduce(function(a,b){return a+b;},0)/Math.max(1,belowFractions.length);
  var residuals=rows.map(function(r){return r.nir-fit.a*r.red-fit.b;}),top=soilQuantile(residuals,.95);
  var veg=rows.filter(function(r){return r.nir-fit.a*r.red-fit.b>=top;});
  var apex={x:soilQuantile(veg.map(function(r){return r.red;}),.5),y:soilQuantile(veg.map(function(r){return r.nir;}),.5)};
  var height=apex.y-fit.a*apex.x-fit.b;
  // Barycentric enclosure: a proposed soil edge plus a robust vegetation region.
  function enclosed(r){var w=height>0?(r.nir-fit.a*r.red-fit.b)/height:Infinity;
    if(w<-.05||w>1.05)return false;var baseX=(r.red-w*apex.x)/Math.max(1e-8,1-w);return w>.95||(baseX>=edgeLo-.05*(hi-lo)&&baseX<=edgeHi+.05*(hi-lo));}
  var enclosure=groups.filter(function(g){return g.length>=20;}).map(function(g){return g.filter(enclosed).length/g.length;});
  var coverage=enclosure.reduce(function(a,b){return a+b;},0)/Math.max(1,enclosure.length),reasons=[];
  if(supportBins.length<4||support.length<20||span<.5)reasons.push('Sparse or narrow edge support');
  if(below>.1)reasons.push('Line cuts through lower observations');
  if(fit.a*edgeLo+fit.b<0||fit.a*edgeHi+fit.b>1)reasons.push('Nonphysical reference reflectance');
  if(height<.02||coverage<.35)reasons.push('Incomplete triangle / weak vegetation-region support');
  var hasDates=rows.some(function(r){return !!r.date;}),hasPixels=rows.some(function(r){return r.pixel_x!==undefined&&r.pixel_y!==undefined;});
  if(hasDates&&Object.keys(dates).length<2)reasons.push('Edge supported on only one date');
  if(hasPixels&&Object.keys(pixels).length<3)reasons.push('Edge supported by too few pixels');
  return {id:bins+' bins / q='+q+' / '+segment,fit:fit,q:q,bins:bins,segment:segment,anchors:used,support:support,
    edgeLo:edgeLo,edgeHi:edgeHi,supportBins:supportBins.length,supportDates:Object.keys(dates).length,supportPixels:Object.keys(pixels).length,
    metadataComplete:hasDates&&hasPixels,span:span,below:below,triangleCoverage:coverage,apex:apex,height:height,reasons:reasons,
    score:3*below+Math.max(0,1-span)+(1-coverage)+Math.max(0,4-supportBins.length),plausible:reasons.length===0};
}
function soilAssessGeometry(input,settings){
  var rows=input.filter(function(r){return carbonFinite(r.red)&&carbonFinite(r.nir)&&r.red>=0&&r.red<=1&&r.nir>=0&&r.nir<=1&&r.red+r.nir>0;});
  var lo=soilQuantile(rows.map(function(r){return r.red;}),.01),hi=soilQuantile(rows.map(function(r){return r.red;}),.99);
  rows=rows.filter(function(r){return r.red>=lo&&r.red<=hi;});
  var months={};rows.forEach(function(r){if(r.month!==undefined)months[r.month]=true;});
  var out={method:'geometry',total:input.length,n:rows.length,months:Object.keys(months).length,lo:lo,hi:hi,reasons:[],warnings:[],supported:false,candidates:[]};
  if(rows.length<100||hi-lo<settings.range){out.reasons.push('Insufficient observations or red range');return out;}
  [12,18,24].forEach(function(bins){[.01,.02,.05].forEach(function(q){['all','bright','dark'].forEach(function(segment){
    var c=soilGeometryCandidate(rows,lo,hi,bins,q,segment,settings.tolerance);if(c)out.candidates.push(c);
  });});});
  if(!out.candidates.length){out.reasons.push('Too few populated brightness bins');return out;}
  out.candidates.sort(function(a,b){return (a.plausible?0:1)-(b.plausible?0:1)||a.score-b.score;});
  var plausible=out.candidates.filter(function(c){return c.plausible;});
  // Compare intercepts only after geometric screening; never select for a desired value.
  var compare=plausible.length?plausible:out.candidates.slice(0,1);
  compare.forEach(function(c){c.sim=soilSimulation(c.fit,settings);});
  var valid=compare.filter(function(c){return c.sim;});
  var median=soilQuantile(valid.map(function(c){return c.sim.evi.b;}),.5);
  valid.forEach(function(c){c.rank=c.score+(Math.abs(c.sim.evi.b-median)+c.sim.horizonShift)/Math.max(.001,settings.tolerance);});
  valid.sort(function(a,b){return a.rank-b.rank;});
  var best=valid.length?valid[0]:out.candidates[0];out.selected=best;out.fit=best.fit;
  out.supported=best.plausible&&!!best.sim;out.reasons=best.reasons.slice();
  if(!best.sim)out.reasons.push('Insufficient nonnegative-PVI simulation support');
  if(!best.metadataComplete)out.warnings.push('Date/pixel support cannot be verified from this input');
  if(!out.supported)out.warnings.push('Exploratory boundary only; missing soil branch remains possible');
  if(best.sim){
    out.eviMin=best.sim.evi.b;out.ndviMin=best.sim.ndvi.b;
    var evi=valid.map(function(c){return c.sim.evi.b;}),ndvi=valid.map(function(c){return c.sim.ndvi.b;});
    out.eviLow=Math.min.apply(null,evi);out.eviHigh=Math.max.apply(null,evi);out.ndviLow=Math.min.apply(null,ndvi);out.ndviHigh=Math.max.apply(null,ndvi);
    if(out.eviHigh-out.eviLow>settings.tolerance)out.warnings.push('EVI2 intercept varies across plausible geometries');
    if(best.sim.horizonShift>settings.tolerance)out.warnings.push('EVI2 intercept depends on PVI fitting range');
    if(best.sim.ndviHorizonShift>settings.tolerance)out.warnings.push('NDVI intercept depends on PVI fitting range');
    if(out.eviMin<0||out.eviMin>1)out.warnings.push('EVI2 intercept outside 0–1');
    var local=soilSimulation(best.fit,{simBounds:[lo,hi,soilQuantile(rows.map(function(r){return r.nir;}),.01),soilQuantile(rows.map(function(r){return r.nir;}),.99)]});
    out.observedDomainEvi=local?local.evi.b:null;
    if(local&&Math.abs(local.evi.b-out.eviMin)>settings.tolerance)out.warnings.push('EVI2 intercept depends on reflectance sampling bounds');
  }
  var exact=[],nd=[];for(var i=0;i<=1000;i++){var x=best.edgeLo+(best.edgeHi-best.edgeLo)*i/1000;exact.push(soilLineEVI2(x,best.fit));nd.push(soilLineNDVI(x,best.fit));}
  out.exactEviMedian=soilQuantile(exact,.5);out.exactNdviMedian=nd.every(carbonFinite)?soilQuantile(nd,.5):null;
  // Keep only the selected simulation cloud; alternatives retain regression metrics.
  out.candidates.forEach(function(c){if(c!==best&&c.sim)delete c.sim.points;});
  return out;
}
function soilAssess(rows,settings){return settings.method==='geometry'?soilAssessGeometry(rows,settings):soilAssessLegacy(rows,settings);}
function soilGeometryDisplay(holder,assessment,rows){
  var candidates=assessment.candidates;
  if(!candidates.length){holder.add(ui.Label(assessment.reasons.join('; ')));return;}
  holder.add(ui.Label((assessment.supported?'Geometrically supported candidate for review. ':'Exploratory fit: '+assessment.reasons.join('; ')+'. ')+(assessment.warnings||[]).join('; ')));
  holder.add(ui.Label('Intercept range across plausible geometries: EVI2 '+soilFormat(assessment.eviLow)+'–'+soilFormat(assessment.eviHigh)+', NDVI '+soilFormat(assessment.ndviLow)+'–'+soilFormat(assessment.ndviHigh)+'. These are sensitivity ranges, not confidence intervals.'));
  var detail=ui.Panel();
  function show(id){
    detail.clear();var c=candidates.filter(function(v){return v.id===id;})[0];if(!c.sim||!c.sim.points)c.sim=soilSimulation(c.fit,UI._soilGeometrySettings);
    candidates.forEach(function(other){if(other!==c&&other.sim)delete other.sim.points;});
    var f=c.fit;detail.add(ui.Label('NIR = '+soilFormat(f.a)+' × Red + '+soilFormat(f.b)+'. '+c.support.length+' near-edge pairs in '+c.supportBins+' brightness bins, '+c.supportDates+' dates and '+c.supportPixels+' pixels.'));
    detail.add(ui.Label('Edge span '+(100*c.span).toFixed(1)+'%; below-line fraction (bin-balanced) '+(100*c.below).toFixed(1)+'%; triangle enclosure '+(100*c.triangleCoverage).toFixed(1)+'%. '+(c.reasons.length?c.reasons.join('; '):'Geometric checks passed; soil identity unverified.')));
    var points=rows.map(function(r){return {c:[{v:r.red},{v:r.nir},{v:null},{v:null},{v:null}]};});
    c.support.forEach(function(r){points.push({c:[{v:r.red},{v:null},{v:r.nir},{v:null},{v:null}]});});
    var vertices=[{x:c.edgeLo,y:f.a*c.edgeLo+f.b},{x:c.edgeHi,y:f.a*c.edgeHi+f.b},c.apex];
    vertices.forEach(function(p){points.push({c:[{v:p.x},{v:null},{v:null},{v:p.y},{v:null}]});});
    [c.edgeLo,c.edgeHi].forEach(function(x){points.push({c:[{v:x},{v:null},{v:null},{v:null},{v:f.a*x+f.b}]});});
    points.sort(function(a,b){return a.c[0].v-b.c[0].v;});
    detail.add(ui.Chart({cols:[{label:'Red',type:'number'},{label:'Screened observations',type:'number'},{label:'Edge support',type:'number'},{label:'Dark / bright / vegetation regions',type:'number'},{label:'Candidate soil edge',type:'number'}],rows:points}).setChartType('LineChart').setOptions({title:'Pooled geometric candidate',hAxis:{title:'Red reflectance'},vAxis:{title:'NIR reflectance'},interpolateNulls:true,series:{0:{pointSize:2,lineWidth:0,color:'#B0BEC5'},1:{pointSize:3,lineWidth:0,color:'#EF6C00'},2:{pointSize:8,lineWidth:0,color:'#1565C0'},3:{pointSize:0,lineWidth:2,color:'#111111'}}}));
    var dates={};c.support.forEach(function(r){var date=r.date||'Date unavailable';dates[date]=(dates[date]||0)+1;});
    detail.add(ui.Label('Edge contributions by composite date: '+Object.keys(dates).sort().map(function(d){return d+': '+dates[d];}).join(', ')));
    if(c.sim){
      detail.add(ui.Label('Simulation bounds [Rmin, Rmax, NIRmin, NIRmax]: '+c.sim.bounds.join(', ')+'. '+c.sim.points.length+' nonnegative-PVI pairs retained from 8192 uniform draws (seed 173).'));
      detail.add(ui.Label('EVI2min intercept '+soilFormat(c.sim.evi.b)+'; slope '+soilFormat(c.sim.evi.a)+'; R² '+soilFormat(c.sim.evi.r2)+'. NDVImin intercept '+soilFormat(c.sim.ndvi.b)+'; R² '+soilFormat(c.sim.ndvi.r2)+'.'));
      var exact=[];for(var t=0;t<=1000;t++)exact.push(soilLineEVI2(c.edgeLo+(c.edgeHi-c.edgeLo)*t/1000,c.fit));
      detail.add(ui.Label('EVI2 intercepts using full / lower 50% / lower 75% PVI ranges: '+c.sim.horizons.map(function(f){return soilFormat(f.b);}).join(', ')+'. Exact soil-line EVI2 median (comparison only): '+soilFormat(soilQuantile(exact,.5))+'.'));
      ['evi','ndvi'].forEach(function(key){var reg=c.sim[key],name=key==='evi'?'EVI2':'NDVI';
        var table=c.sim.points.map(function(p){return {c:[{v:p.x},{v:key==='evi'?p.y:p.ndvi},{v:null}]};});
        [0,Math.max.apply(null,c.sim.points.map(function(p){return p.x;}))].forEach(function(x){table.push({c:[{v:x},{v:null},{v:reg.a*x+reg.b}]});});table.sort(function(a,b){return a.c[0].v-b.c[0].v;});
        detail.add(ui.Chart({cols:[{label:'PVI',type:'number'},{label:'Simulated '+name,type:'number'},{label:'OLS fit',type:'number'}],rows:table}).setChartType('LineChart').setOptions({title:name+' versus PVI — simulated reflectances',hAxis:{title:'PVI'},vAxis:{title:name},interpolateNulls:true,series:{0:{pointSize:2,lineWidth:0},1:{pointSize:0,lineWidth:2,color:'#111111'}}}));
      });
    }
  }
  holder.add(ui.Select({items:candidates.map(function(c){return {label:c.id+' — '+(c.plausible?'plausible':'exploratory')+(c.sim?' — EVI2 intercept '+soilFormat(c.sim.evi.b):''),value:c.id};}),value:assessment.selected.id,onChange:show,style:{stretch:'horizontal'}}));holder.add(detail);show(assessment.selected.id);
}
function soilFormat(v){return carbonFinite(v)?v.toFixed(4):'unavailable';}

function soilSearchSettings(){
  var s=lstNDVISettings(UI.soilYear.getValue(),UI.soilBox.getValue());
  s.method=UI.soilMethod?UI.soilMethod.getValue():'legacy';
  s.simBounds=s.method==='geometry'?[UI.soilRmin,UI.soilRmax,UI.soilNmin,UI.soilNmax].map(function(w){return Number(w.getValue());}):[0,1,0,1];
  if(s.simBounds.some(function(v){return !isFinite(v)||v<0||v>1;})||s.simBounds[0]>=s.simBounds[1]||s.simBounds[2]>=s.simBounds[3])throw new Error('Simulation bounds must be within 0–1, with each minimum smaller than its maximum.');
  s.rings=Number(UI.soilRings.getValue());s.cut=Number(UI.soilCut.getValue());s.range=Number(UI.soilRange.getValue());s.tolerance=Number(UI.soilTolerance.getValue());
  if(s.boxKm<2||s.boxKm>10)throw new Error('Reference box side must be 2–10 km.');
  if([0,1,2].indexOf(s.rings)<0)throw new Error('Choose a supported search extent.');
  if(s.method==='legacy'&&(!isFinite(s.cut)||s.cut<.05||s.cut>.8))throw new Error('EVI2 ceiling must be 0.05–0.8.');
  if(!isFinite(s.range)||s.range<.01||s.range>.3)throw new Error('Minimum red range must be 0.01–0.3.');
  if(!isFinite(s.tolerance)||s.tolerance<.001||s.tolerance>.05)throw new Error('Stability tolerance must be 0.001–0.05 reflectance.');
  var latValue=UI.soilLat.getValue(),lonValue=UI.soilLon.getValue();
  s.latText=latValue==null?'':String(latValue).trim();s.lonText=lonValue==null?'':String(lonValue).trim();
  if(s.latText||s.lonText){
    s.lat=Number(s.latText);s.lon=Number(s.lonText);
    if(!s.latText||!s.lonText||!isFinite(s.lat)||!isFinite(s.lon)||Math.abs(s.lat)>90||Math.abs(s.lon)>180)throw new Error('Enter both valid reference-center coordinates, or leave both blank.');
  }
  return s;
}
function soilSearchTiles(target,reference,settings){
  var base=lstNDVIBox(reference,settings.boxKm),proj=base.projection(),xy=reference.centroid(1).transform(proj,1).coordinates();
  var x=ee.Number(xy.get(0)),y=ee.Number(xy.get(1)),step=settings.boxKm*1000,tiles=[];
  for(var ring=0;ring<=settings.rings;ring++)for(var iy=-ring;iy<=ring;iy++)for(var ix=-ring;ix<=ring;ix++){
    if(Math.max(Math.abs(ix),Math.abs(iy))!==ring)continue;
    var cx=x.add(ix*step),cy=y.add(iy*step),box=ee.Geometry.Rectangle(ee.List([cx.subtract(step/2),cy.subtract(step/2),cx.add(step/2),cy.add(step/2)]),proj,false);
    tiles.push({id:'R'+ring+' E'+ix+' N'+iy,box:box,ring:ring,center:box.centroid(1).coordinates(),distance:box.centroid(1).distance(target.centroid(1),1).divide(1000)});
  }
  return tiles;
}
function sampleSoilReference(box,year,firstMonth,lastMonth){
  var source=ee.ImageCollection('MODIS/061/MYD09A1'),projection=ee.Image(source.first()).select('sur_refl_b01').projection().scale(2,2);
  var start=ee.Date.fromYMD(year,firstMonth,1),end=lastMonth===12?ee.Date.fromYMD(year+1,1,1):ee.Date.fromYMD(year,lastMonth+1,1);
  var images=source.filterBounds(box).filterDate(start,end),grid=ee.Image.pixelCoordinates(projection).rename(['pixel_x','pixel_y']);
  var samples=ee.FeatureCollection(images.toList(images.size().max(1)).map(function(item){
    var image=ee.Image(item);
    // Reflectance-only sampling: do not require a coincident valid LST retrieval.
    var reflectance=soilLandReflectance(image,projection).select(['red','nir','NDVI_valid_fraction']).addBands(grid);
    return reflectance.sample({region:box,projection:projection,geometries:false,tileScale:4,dropNulls:true}).map(function(f){return f.set({date:image.date().format('YYYY-MM-dd'),month:image.date().get('month')});});
  })).flatten();
  return ee.Dictionary({rows:samples.toList(samples.size().max(1)).map(function(f){return ee.Feature(f).toDictionary();}),periods:images.size()});
}
function runSoilSearch(target){
  var settings;try{settings=soilSearchSettings();}catch(e){setStatus(e.message);return;}
  var reference=settings.latText?ee.Geometry.Point([settings.lon,settings.lat]):target.centroid(1);
  UI._soilGeometrySettings=settings;
  var tiles=soilSearchTiles(target,reference,settings),request=UI._requestId,index=0,month=1,rows=[],results=[];
  UI.runButton.setDisabled(true);UI.resultsPanel.style().set({shown:true});UI.resultsOverlay.setExpanded(true);
  if(UI.map){UI.map.centerObject(lstNDVIBox(reference,settings.boxKm*(2*settings.rings+1)));UI._soilLayers=UI._soilLayers||[];
    var targetLayer=ui.Map.Layer(target.centroid(1),{color:'00FFFF'},'Original model target');UI.map.layers().add(targetLayer);UI._soilLayers.push(targetLayer);}
  var progress=ui.Label('Searching separate reference boxes…');UI.chartHolder.add(progress);
  function render(){
    var supported=results.filter(function(r){return r.assessment.supported;}).sort(function(a,b){return a.distance-b.distance;});
    var summary=supported.length?supported.length+' candidate references for review. Nearest: '+supported[0].id+' ('+supported[0].distance.toFixed(1)+' km from target).':'No candidate passed the configured support checks. Inspect the exploratory fits and diagnostics; this does not establish that soil is absent.';
    if(supported.length){var vals=supported.map(function(r){return r.assessment.eviMin;});summary+=' Candidate EVI2min median '+soilQuantile(vals,.5).toFixed(4)+', range '+Math.min.apply(null,vals).toFixed(4)+'–'+Math.max.apply(null,vals).toFixed(4)+'.';var ndviVals=supported.map(function(r){return r.assessment.ndviMin;});summary+=' Candidate NDVImin median '+soilQuantile(ndviVals,.5).toFixed(4)+', range '+Math.min.apply(null,ndviVals).toFixed(4)+'–'+Math.max.apply(null,ndviVals).toFixed(4)+'.';}
    progress.setValue(summary);
    results.forEach(function(r){var a=r.assessment;UI.chartHolder.add(ui.Label(r.id+' | '+r.distance.toFixed(1)+' km | '+(!carbonFinite(a.eviMin)?'EVI2min unavailable':'EVI2min '+a.eviMin.toFixed(4))+' | '+(carbonFinite(a.ndviMin)?'NDVImin '+a.ndviMin.toFixed(4):'NDVImin unavailable')+' | '+(a.supported?'supported candidate':a.reasons.join('; '))));});
    UI.chartHolder.add(ui.Label(settings.method==='geometry'?'Pooled geometric candidates; soil identity unverified. EVI2min and NDVImin are simulated index–PVI regression intercepts. Candidate ranges describe sensitivity, not confidence intervals. No model parameter was changed.': 'Empirical screening only: verify substrate, terrain, litter and exposed-soil identity before transferring a reference. No model parameter was changed. EVI2min and NDVImin are medians along PVI=0 on the same uniform red grid; their ranges describe brightness dependence, not confidence intervals.'));
    var holder=ui.Panel();UI.chartHolder.add(holder);
    function show(id){
      holder.clear();var r=results.filter(function(v){return v.id===id;})[0],a=r.assessment;
      holder.add(ui.Label(r.id+' | center '+r.center[1].toFixed(5)+', '+r.center[0].toFixed(5)+' | '+r.distance.toFixed(1)+' km from original target'));
      holder.add(ui.Label(a.supported?'Statistically supported candidate; soil identity unverified.':'Not supported: '+a.reasons.join('; ')));
      holder.add(ui.Label(a.n+' fitting pairs / '+a.total+' total, '+a.months+' eligible months.'));
      if(a.method==='geometry'){soilGeometryDisplay(holder,a,r.rows);return;}
      if(!a.fit)return;
      function fmt(v){return carbonFinite(v)?v.toFixed(4):'insufficient support';}
      holder.add(ui.Label('NIR = '+fmt(a.fit.a)+' × Red + '+fmt(a.fit.b)+'. Red range '+fmt(a.lo)+'–'+fmt(a.hi)+'. EVI2min '+fmt(a.eviMin)+' ('+fmt(a.eviLow)+'–'+fmt(a.eviHigh)+'). NDVImin '+fmt(a.ndviMin)+' ('+fmt(a.ndviLow)+'–'+fmt(a.ndviHigh)+').'));
      holder.add(ui.Label('Maximum NIR shifts: quantile '+fmt(a.quantileShift)+', omitted month '+fmt(a.omissionShift)+', EVI2 ceiling ±0.05 '+fmt(a.thresholdShift)+'. Binned boundary departure '+fmt(a.curvature)+'.'));
      var points=r.rows.map(function(v){return {c:[{v:v.red},{v:v.nir},{v:null}]};});
      points.push({c:[{v:a.lo},{v:null},{v:a.fit.a*a.lo+a.fit.b}]},{c:[{v:a.hi},{v:null},{v:a.fit.a*a.hi+a.fit.b}]});
      points.sort(function(x,y){return x.c[0].v-y.c[0].v;});
      holder.add(ui.Chart({cols:[{id:'red',label:'Red',type:'number'},{id:'nir',label:'All QA-valid pairs',type:'number'},{id:'fit',label:'5% boundary',type:'number'}],rows:points}).setChartType('LineChart').setOptions({title:'Reference '+r.id+' — '+settings.year,hAxis:{title:'Red reflectance'},vAxis:{title:'NIR reflectance'},interpolateNulls:true,series:{0:{lineWidth:0,pointSize:2,color:'#78909C'},1:{lineWidth:2,pointSize:0,color:'#111111'}}}));
      var line=[];for(var k=0;k<=100;k++){var red=a.lo+(a.hi-a.lo)*k/100;line.push({c:[{v:red},{v:soilLineEVI2(red,a.fit)},{v:soilLineNDVI(red,a.fit)}]});}
      holder.add(ui.Chart({cols:[{label:'Red',type:'number'},{label:'EVI2 at PVI = 0',type:'number'},{label:'NDVI at PVI = 0',type:'number'}],rows:line}).setChartType('LineChart').setOptions({title:'Bare-soil candidate indices versus brightness',hAxis:{title:'Red reflectance'},vAxis:{title:'Index at PVI = 0'},legend:{position:'bottom'}}));
    }
    UI.chartHolder.add(ui.Select({items:results.map(function(r){return {label:r.id+' — '+(r.assessment.supported?'supported candidate':'not supported')+' — '+r.distance.toFixed(1)+' km',value:r.id};}),value:supported.length?supported[0].id:results[0].id,onChange:show,style:{stretch:'horizontal'}}));
    show(supported.length?supported[0].id:results[0].id);
    UI._soilSearch={settings:settings,results:results};UI.runButton.setDisabled(false);setStatus('Soil-reference search complete. '+summary+' No ACCLIM model calculations were run.');
  }
  function next(){
    if(request!==UI._requestId)return;
    if(index>=tiles.length){render();return;}
    var tile=tiles[index];setStatus('Soil reference '+(index+1)+'/'+tiles.length+' ('+tile.id+'), months '+month+'–'+(month+2)+'…');
    sampleSoilReference(tile.box,settings.year,month,month+2).evaluate(function(data,error){
      if(request!==UI._requestId)return;
      if(error||!data){UI.runButton.setDisabled(false);setStatus('Soil-reference sampling failed: '+(error||'No response')+'. Reduce search extent or retry.');return;}
      rows=rows.concat(data.rows||[]);month+=3;
      if(month<=12){next();return;}
      var assessment=soilAssess(rows,settings),savedRows=rows;
      ee.Dictionary({center:tile.center,distance:tile.distance}).evaluate(function(meta,err){
        if(request!==UI._requestId)return;
        if(err||!meta){UI.runButton.setDisabled(false);setStatus('Reference geometry error: '+err);return;}
        results.push({id:tile.id,center:meta.center,distance:meta.distance,assessment:assessment,rows:savedRows});
        if(UI.map){var layer=ui.Map.Layer(tile.box,{color:assessment.supported?'22AA55':'CC8844'},tile.id+' '+(assessment.supported?'supported candidate':'not supported'),true);UI.map.layers().add(layer);UI._soilLayers.push(layer);}
        progress.setValue(results.length+'/'+tiles.length+' boxes evaluated. '+results.filter(function(r){return r.assessment.supported;}).length+' supported candidates.');
        rows=[];month=1;index++;ui.util.setTimeout(next,25);
      });
    });
  }
  next();
}

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

// Standalone observed LST-NDVI analysis uses the same geometry controls.
UI.analysisMode=ui.Select({items:[{label:'Run ACCLIM model',value:'model'},
  {label:'Plot MODIS LST vs NDVI (no model run)',value:'lst_ndvi'},
  {label:'Search for soil reference (no model run)',value:'soil_search'}],value:'model',
  onChange:function(mode){UI.lstNDVIControls.style().set({shown:mode==='lst_ndvi'});UI.soilControls.style().set({shown:mode==='soil_search'});UI.runButton.setLabel(mode==='soil_search'?'Search soil references':mode==='lst_ndvi'?'Plot LST vs NDVI':'Run ACCLIM');},style:{stretch:'horizontal'}});
UI.lstNDVIYear=ui.Textbox({value:String(Number(DEFAULTS.startDate.slice(0,4))),placeholder:'Year (YYYY)'});
UI.lstNDVIBoxKm=ui.Textbox({value:'10',placeholder:'Box side (km)'});
UI.lstNDVIQuality=ui.Select({items:[
  {label:lstNDVIQualityLabel('valid_2k'),value:'valid_2k'},
  {label:lstNDVIQualityLabel('strict'),value:'strict'}],value:'valid_2k',style:{stretch:'horizontal'}});
UI.lstNDVIControls=ui.Panel({widgets:[ui.Label('LST-NDVI calendar year'),UI.lstNDVIYear,
  ui.Label('Square box side (km; centered on selected geometry)'),UI.lstNDVIBoxKm,
  ui.Label('LST quality screening'),UI.lstNDVIQuality,
  ui.Label('Aqua MYD09A1 NDVI + MYD11A2 daytime LST, paired 8-day intervals. Uses the year above; model dates and weather selections are ignored.')],style:{shown:false}});
UI.inputPanel.add(ui.Label('Analysis')).add(UI.analysisMode).add(UI.lstNDVIControls);

UI.soilMethod=ui.Select({items:[{label:'Pooled edge + triangle geometry + PVI intercept',value:'geometry'},{label:'Previous month-balanced quantile method',value:'legacy'}],value:'geometry'});
UI.soilRmin=ui.Textbox({value:'0'});UI.soilRmax=ui.Textbox({value:'1'});
UI.soilNmin=ui.Textbox({value:'0'});UI.soilNmax=ui.Textbox({value:'1'});
UI.soilYear=ui.Textbox({value:String(Number(DEFAULTS.startDate.slice(0,4)))});
UI.soilBox=ui.Textbox({value:'10'});
UI.soilRings=ui.Select({items:[{label:'Original/reference box only',value:'0'},{label:'3 × 3 boxes (default)',value:'1'},{label:'5 × 5 boxes',value:'2'}],value:'1'});
UI.soilLat=ui.Textbox({value:'',placeholder:'Optional reference-center latitude'});
UI.soilLon=ui.Textbox({value:'',placeholder:'Optional reference-center longitude'});
UI.soilCut=ui.Textbox({value:'0.25'});
UI.soilRange=ui.Textbox({value:'0.04'});
UI.soilTolerance=ui.Textbox({value:'0.01'});
UI.soilControls=ui.Panel({widgets:[ui.Label('Soil-reference method'),UI.soilMethod,ui.Label('Soil-reference calendar year'),UI.soilYear,
  ui.Label('Reference-box side (2–10 km)'),UI.soilBox,ui.Label('Expansion around reference center'),UI.soilRings,
  ui.Label('Leave both coordinates blank to search around the original target. Enter both to shift only the reference search.'),UI.soilLat,UI.soilLon,
  ui.Label('EVI2 ceiling (previous method only; ignored for geometry)'),UI.soilCut,ui.Label('Minimum red-reflectance span'),UI.soilRange,
  ui.Label('Research tolerance (NIR edge checks / EVI2 intercept warnings)'),UI.soilTolerance,
  ui.Label('Simulation red minimum / maximum'),UI.soilRmin,UI.soilRmax,
  ui.Label('Simulation NIR minimum / maximum'),UI.soilNmin,UI.soilNmax,
  ui.Label('Geometry: pooled dates, no EVI2 cutoff; 12/18/24 red bins and 1/2/5% lower quantiles. Fits seek support in ≥4 brightness intervals, ≥2 dates and ≥3 pixels when metadata is available. Triangle support and intercept sensitivity are diagnostics, not proof of soil identity. Simulations: 8192 reproducible independent uniform draws; retain PVI ≥0. Cloud/shadow/snow/water screening before aggregation. No automatic model-baseline changes.')],style:{shown:false}});
UI.inputPanel.add(UI.soilControls);


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
  {label: 'FAPAR = kNDVI (spatial mean |NIR − red| sigma)', value: 'KNDVI_SPATIAL'},
  {label: 'FAPAR = kNDVI (fixed sigma)', value: 'KNDVI'},
  {label: 'FAPAR = m_abs_ndvi × NDVI + b_abs_ndvi (PROSAIL chlorophyll fit)', value: 'RT_CHL_NDVI'},
  {label: 'FAPAR = m_abs_evi2 × EVI2 + b_abs_evi2 (PROSAIL chlorophyll fit)', value: 'RT_CHL_EVI2'},
  {label: 'FAPAR = PROSAIL chlorophyll absorption', value: 'RT_CHL'},
  {label: 'FAPAR = PROSAIL green-leaf absorption', value: 'RT_LEAF'},
  {label: 'FAPAR = PROSAIL daylight interception', value: 'PROSAIL'},
  {label: 'FAPAR = FIPAR × FAPARmult', value: 'FIPAR'},
  {label: 'FAPAR = NDVI²', value: 'NDVI'},
  {label: 'FAPAR = NIRv (NIR × NDVI)', value: 'NIRV'},
  {label: 'FAPAR = NDVI^(R × 0.5 + (1 − R) × 2)', value: 'NDVI_RED_EXP'},
  {label: 'FAPAR = NDVI^ndvi_exp − fitted NDVImin (GE′ = 0)', value: 'NDVI_GE_OFFSET'},
  {label: 'FAPAR = fge', value: 'FGE'},
  {label: 'FAPAR = min(fge, NDVI^ndvi_exp)', value: 'GE_NDVI_MIN'},
  {label: 'FAPAR = w × NDVI^ndvi_exp + (1 − w) × fge', value: 'NDVI_BLEND'},
  {label: 'FAPAR = FIPAR; VI minima = 10-year fifth percentiles', value: 'FIPAR_P05'},
  {label: 'FAPAR = FIPAR; NDVImin fitted at GE′ = 0', value: 'FIPAR_GE'},
  {label: 'FAPAR = EVI2 scaled to FIPAR', value: 'EVI2_FITTED'},
  {label: 'FAPAR = EVI2 FIPAR-fitted min; fixed max 0.6', value: 'EVI2_FITTED_MIN'},
  {label: 'FAPAR = scaled EVI2 (fixed bounds, EVI2exp)', value: 'EVI2_FIXED_POWER'}], value: 'FIPAR',
  style: {stretch: 'horizontal', margin: '4px 0px'}});
UI.inputPanel.add(ui.Label('FAPAR method (applied on Run ACCLIM)'));
UI.inputPanel.add(UI.faparSelector);
UI.kndviSigma=ui.Textbox({value:String(ACCLIM_PARAMS.kndvi_sigma)});
UI.kndviPanel=ui.Panel([ui.Label('kNDVI sigma (fractional reflectance units)'),UI.kndviSigma]);
UI.inputPanel.add(UI.kndviPanel);
UI.kndviRadius=ui.Textbox({value:'5'});
UI.kndviRegionPanel=ui.Panel([ui.Label('Spatial sigma radius around site center (km)'),UI.kndviRadius]);
UI.inputPanel.add(UI.kndviRegionPanel);
UI.reuseInputs=ui.Checkbox({label:'Reuse sampled weather and satellite inputs',value:true});
UI.inputPanel.add(UI.reuseInputs);
UI.inputPanel.add(ui.Label('Vegetation parameters (applied on Run ACCLIM)'));
UI.vegetationParameters = {};
UI.soilMinimumSelectors = {};
[
  {key:'GE_scale', label:'GE_scale'},
  {key:'ndvi_exp', label:'ndvi_exp'},
  {key:'FAPARmult', label:'FAPARmult (direct FIPAR option)'},
  {key:'fDgpp_exp', label:'fDgpp_exp (GPP VPD response; 0 disables)'},
  {key:'fDtrans_exp', label:'fDtrans_exp (transpiration VPD response; 0 disables)'},
  {key:'fMtrans_exp', label:'fMtrans_exp (transpiration seasonal FAPAR response; 0 disables)'},
  {key:'fMgpp_exp', label:'fMgpp_exp (GPP seasonal FAPAR response; 0 disables)'},
  {key:'EVI2min', label:'EVI2min'},
  {key:'EVI2exp', label:'EVI2exp (scaled EVI2)'},
  {key:'EVI2max', label:'EVI2max'},
  {key:'ndvi_min', label:'NDVImin'},
  {key:'ndvi_max', label:'NDVImax'}
].forEach(function(parameter) {
  UI.inputPanel.add(ui.Label(parameter.label));
  var input = ui.Textbox({value:String(ACCLIM_PARAMS[parameter.key]), style:{stretch:'horizontal'}});
  UI.vegetationParameters[parameter.key] = input;
  UI.inputPanel.add(input);
  if (parameter.key === 'EVI2min' || parameter.key === 'ndvi_min') {
    var auto = ui.Checkbox({label:'Auto: SoilGrids texture', value:false,
      onChange:function(selected) { input.setDisabled(selected); }});
    UI.soilMinimumSelectors[parameter.key] = auto;
    UI.inputPanel.add(auto);
  }
});
UI.inputPanel.add(ui.Label('Fitted/percentile methods estimate their bounds for each run. The fixed-maximum experiment uses EVI2max = 0.6.'));
function readVegetationParameters() {
  var values = {};
  Object.keys(UI.vegetationParameters).forEach(function(key) {
    var text = UI.vegetationParameters[key].getValue().trim();
    if (UI.soilMinimumSelectors[key] && UI.soilMinimumSelectors[key].getValue()) {
      values[key === 'ndvi_min' ? 'autoSoilNDVImin' : 'autoSoilEVI2min'] = true;
      return;
    }
    var value = Number(text);
    if (!text || !isFinite(value)) throw new Error(key + ' must be a finite number.');
    values[key] = value;
  });
  if (values.GE_scale <= 0) throw new Error('GE_scale must be positive.');
  if (values.FAPARmult < 0) throw new Error('FAPARmult must be nonnegative.');
  ['fDgpp_exp','fDtrans_exp','fMgpp_exp','fMtrans_exp'].forEach(function(key) {
    if (values[key] < 0) throw new Error(key + ' must be nonnegative.');
  });
  if (values.EVI2exp < 0) throw new Error('EVI2exp must be nonnegative.');
  if (values.ndvi_exp < 0) throw new Error('ndvi_exp must be nonnegative.');
  if (values.ndvi_max <= values.ndvi_min) throw new Error('NDVImax must exceed NDVImin.');
  if (values.EVI2max <= values.EVI2min) throw new Error('EVI2max must exceed EVI2min.');
  return values;
}

UI.sharedFPAR = ui.Checkbox({label: 'Use FAPAR for both photosynthesis and water fluxes', value: false});
UI.inputPanel.add(UI.sharedFPAR);
UI.inputPanel.add(ui.Label('FIPAR method (before optional shared FAPAR override)'));
UI.fiparSelector = ui.Select({
  items: [
    {label: 'FIPAR = m_cover × NDVI + b_cover (PROSAIL daylight fit)', value: 'RT_NDVI'},
    {label: 'FIPAR = PROSAIL daylight interception', value: 'PROSAIL'},
    {label: 'FIPAR = scaled NDVI', value: 'NDVI'},
    {label: 'FIPAR = w × NDVI^ndvi_exp + (1 − w) × fge', value: 'NDVI_BLEND'}
  ],
  value: 'RT_NDVI',
  style: {stretch: 'horizontal'}
});
UI.inputPanel.add(UI.fiparSelector);
UI.rtBoxKm=ui.Textbox({value:'10',style:{stretch:'horizontal'}});
UI.rtObjective=rtObjectiveSelector('BRF');
UI.rtBoxPanel=ui.Panel([ui.Label('PROSAIL fitting objective'),UI.rtObjective,ui.Label('PROSAIL box side (km), Aqua MODIS 8-day'),UI.rtBoxKm]);
UI.inputPanel.add(UI.rtBoxPanel);
UI.rtDiffuse=ui.Textbox({value:'0.5'});
UI.rtSkyPanel=ui.Panel([ui.Label('Assumed diffuse fraction of incident PAR photons (0–1)'),UI.rtDiffuse]);
UI.inputPanel.add(UI.rtSkyPanel);
UI.inputPanel.add(ui.Button('PROSAIL retrieval diagnostics', rtApp));
function updateRTOptions() {
  UI.kndviRegionPanel.style().set('shown',UI.faparSelector.getValue()==='KNDVI_SPATIAL');
  UI.kndviPanel.style().set('shown',UI.faparSelector.getValue()==='KNDVI');
  var rtWater=['PROSAIL','RT_NDVI'].indexOf(UI.fiparSelector.getValue())>=0;
  if(rtWater)UI.sharedFPAR.setValue(false,false);
  UI.sharedFPAR.setDisabled(rtWater);
  var absorption=rtAbsorptionMethod(UI.faparSelector.getValue());
  UI.rtSkyPanel.style().set('shown',absorption||rtWater||UI.faparSelector.getValue()==='PROSAIL');
  UI.rtBoxPanel.style().set('shown',absorption||rtWater||UI.faparSelector.getValue()==='PROSAIL');
}
UI.fiparSelector.onChange(updateRTOptions);
UI.faparSelector.onChange(updateRTOptions);
updateRTOptions();
UI.inputPanel.add(ui.Label('Blend weight w = NDVI^blend_exp'));
setupTowerImport();
UI.weatherSelector=ui.Select({items:Object.keys(WEATHER_SOURCES).map(function(k){return {label:WEATHER_SOURCES[k].label,value:k};}),value:'TERRACLIMATE'});
UI.satelliteSelector=ui.Select({items:Object.keys(SATELLITE_SOURCES).map(function(k){return {label:SATELLITE_SOURCES[k].label,value:k};}),value:'MODIS16'});
UI.timestepSelector=ui.Select({items:[{label:'Daily',value:'daily'},{label:'Weekly',value:'weekly'},{label:'Monthly',value:'monthly'}],value:'monthly'});
UI.inputPanel.add(ui.Label('Timestep (all inputs, model, charts and exports)')).add(UI.timestepSelector);
UI.inputPanel.add(ui.Label('Weeks start January 1; end periods are clipped to selected dates. Coarser weather means are held constant and precipitation allocated by duration, with flags in the export.'));
UI.inputPanel.add(ui.Label('Weather data (applied on Run ACCLIM)')).add(UI.weatherSelector);
UI.inputPanel.add(ui.Label('Satellite reflectance (applied on Run ACCLIM)')).add(UI.satelliteSelector);
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
// Optional chart series use the native model values already fetched for carbon products.
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
// Normalize legacy monthly and new selectable tower weather into explicit intervals.
function towerWeatherPlan(tower,start,end) {
  if(!tower||!tower.weather||tower.weather.version!==2)throw new Error('Reprocess the source CSV/ZIP with the current importer for tower weather.');
  var cadence=tower.weather.cadence||'monthly';
  if(['daily','weekly','monthly'].indexOf(cadence)<0)throw new Error('Invalid tower weather timestep.');
  if(cadence==='monthly') {
    var plan=legacyTowerWeatherPlan(tower,start,end);plan.cadence=cadence;plan.records=[];
    var d=new Date(start.slice(0,7)+'-01T00:00:00Z');
    while(d.toISOString().slice(0,7)<=end.slice(0,7)){
      var key=d.toISOString().slice(0,7),next=Date.UTC(d.getUTCFullYear(),d.getUTCMonth()+1,1),entry=plan.byMonth[key];
      plan.records.push(entry||{valid:false,data:{month:key,date:key+'-01',period_days:(next-d.getTime())/86400000}});d=new Date(next);
    }
    return plan;
  }
  var weatherRows=tower.weather.periods;
  // Full saved/batch JSON retains daily weather; compact copy payloads carry periods.
  if(!Array.isArray(weatherRows)&&Array.isArray(tower.weather.daily)){
    if(cadence==='daily')weatherRows=tower.weather.daily;
    else{
      var groups={},seenDays={};
      tower.weather.daily.forEach(function(d){
        if(seenDays[d.date])throw new Error('Duplicate daily tower weather.');seenDays[d.date]=true;
        var t=carbonDateMillis(d.date),year=new Date(t).getUTCFullYear(),jan1=Date.UTC(year,0,1);
        var a=jan1+Math.floor((t-jan1)/604800000)*604800000,key=carbonDateString(a);
        if(!groups[key])groups[key]=[];groups[key].push(d);
      });
      weatherRows=Object.keys(groups).sort().map(function(key){
        var rows=groups[key],a=carbonDateMillis(key),days=Math.min(7,(Date.UTC(new Date(a).getUTCFullYear()+1,0,1)-a)/86400000);
        var r={date:key,period_days:days,expected:days*1440/tower.intervalMinutes,counts:{},complete:rows.length===days&&rows.every(function(d){return d.complete===true;})};
        ['ta','sw','ea','pa','rain','rn'].forEach(function(k){r.counts[k]=rows.reduce(function(n,d){return n+(d.counts&&carbonFinite(d.counts[k])?d.counts[k]:0);},0);});
        ['Tmin','Tmax','Tdaily','SWRdaily','ea','Pdaily','precip','Rndaily'].forEach(function(k){r[k]=rows.length===days&&rows.every(function(d){return carbonFinite(d[k]);})?rows.reduce(function(n,d){return n+d[k];},0)/(k==='precip'?1:days):null;});
        return r;
      });
    }
  }
  if(!Array.isArray(weatherRows))throw new Error('Use the importer to prepare '+cadence+' tower weather.');
  var supplied={},fields=['Tmin','Tmax','Tdaily','SWRdaily','ea','Pdaily','precip','Rndaily'],counts=['ta','sw','ea','pa','rain','rn'];
  if(tower.intervalMinutes!==30&&tower.intervalMinutes!==60)throw new Error('Unsupported tower recording interval.');
  weatherRows.forEach(function(r){
    if(!/^\d{4}-\d{2}-\d{2}$/.test(r.date)||!isFinite(Date.parse(r.date))||new Date(r.date).toISOString().slice(0,10)!==r.date||supplied[r.date])throw new Error('Invalid or duplicate tower-weather period.');
    var a=Date.parse(r.date+'T00:00:00Z'),yearStart=Date.UTC(new Date(a).getUTCFullYear(),0,1);
    var days=cadence==='daily'?1:Math.min(7,(Date.UTC(new Date(a).getUTCFullYear()+1,0,1)-a)/86400000);
    if(r.period_days!==days||(cadence==='weekly'&&(a-yearStart)%604800000!==0))throw new Error('Invalid tower-weather period duration or weekly anchor.');
    var expected=days*1440/tower.intervalMinutes;
    supplied[r.date]={data:r,valid:r.complete===true&&r.expected===expected&&fields.every(function(k){return carbonFinite(r[k]);})&&counts.every(function(k){return r.counts&&r.counts[k]===expected;})&&r.Tmin<=r.Tdaily&&r.Tdaily<=r.Tmax&&r.ea>=0&&r.Pdaily>0&&r.precip>=0};
  });
  var plan={cadence:cadence,records:[],segments:[],missing:[],rejected:[],byMonth:{},source:tower.weather.source||{}};
  (tower.weather.months||[]).forEach(function(r){plan.byMonth[r.month]={data:r,valid:r.complete===true};});
  var a=carbonDateMillis(start),stop=carbonDateMillis(end),segment=null;
  if(cadence==='weekly'){var y=Date.UTC(new Date(a).getUTCFullYear(),0,1);a=y+Math.floor((a-y)/604800000)*604800000;}
  while(a<=stop){
    var b=cadence==='daily'?a+86400000:Math.min(a+604800000,Date.UTC(new Date(a).getUTCFullYear()+1,0,1));
    var key=carbonDateString(a),entry=supplied[key]||{valid:false,data:{date:key,period_days:(b-a)/86400000}};
    plan.records.push(entry);
    if(entry.valid){if(!segment){segment=[];plan.segments.push(segment);}segment.push(entry.data);}
    else{segment=null;plan.missing.push(key);plan.rejected.push(key+': incomplete or missing '+cadence+' weather');}
    a=b;
  }
  return plan;
}

function legacyTowerWeatherPlan(tower, start, end) {
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
  var date=row.date||row.month+'-01',d=ee.Date(date);
  var days=row.period_days||d.advance(1,'month').difference(d,'day');
  return ee.Feature(geom,{date:date,millis:d.millis(),year:d.get('year'),month:d.get('month'),
    Tmin:(row.Tmin===undefined?null:row.Tmin),Tmax:(row.Tmax===undefined?null:row.Tmax),Tdaily:(row.Tdaily===undefined?null:row.Tdaily),SWRdaily:(row.SWRdaily===undefined?null:row.SWRdaily),ea:(row.ea===undefined?null:row.ea),
    Pdaily:(row.Pdaily===undefined?null:row.Pdaily),precip:(row.precip===undefined?null:row.precip),Rndaily:(row.Rndaily===undefined?null:row.Rndaily),period_days:days,period_seconds:ee.Number(days).multiply(86400),period_start_millis:d.millis(),period_end_millis:d.advance(days,'day').millis(),
    weather_source:'tower',weather_missing_day_pct:row.weather_missing_day_pct||0,weather_source_cadence:row.weather_source_cadence||'unknown',weather_temporally_allocated_pct:row.weather_temporally_allocated_pct||0,weather_temporal_method:row.weather_temporal_method||'source values',weather_energy_basis:'LE_CORR + H_CORR',weather_status:'complete',weather_expected_intervals:(row.expected===undefined?null:row.expected)});
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
  UI.weatherNote=ui.Label({value:'Optional: tower temperature, humidity, shortwave, pressure, precipitation and LE_CORR + H_CORR. Carbon requires complete months; bucket-only gaps are interpolated, with no gridded substitutions. CO2 and satellite vegetation remain unchanged.',style:{whiteSpace:'pre-wrap'}});
  UI.inputPanel.add(UI.weatherNote);
}
function addTowerControls(cache, tower, start, end, lat, lon, coordinateMode) {
  var note=ui.Label({value:'Import tower data in the inputs before running ACCLIM.',style:{whiteSpace:'pre-wrap'}});
  UI._chartCard.add(note); if(!tower)return;
  var months={};
  var coverage=ui.Label({value:'',style:{whiteSpace:'pre-wrap'}});
  var monthChoice=ui.Select({items:[],placeholder:'Tower coverage by timestep',onChange:function(key){
    var m=months[key]; if(!m)return;
    function pct(x){return x===null?'unavailable':x.toFixed(1)+'%';}
    coverage.setValue(key+': missing GPP '+pct(m.missingPct)+(tower.secondVariable?'; missing '+tower.secondVariable+' '+pct(m.secondMissingPct):'')+
      '; LE_CORR missing '+pct(m.leMissingPct)+'. Valid primary/daytime records: '+m.n+'/'+m.daylight+'. Unknown radiation: '+(m.total?(100*m.unknownLight/m.total).toFixed(1)+'%':'unavailable')+
      '; absent intervals: '+m.absent+'. No coverage cutoff. Missing percentages exclude records whose daytime status is unknown.');
  }});
  UI._chartCard.add(monthChoice).add(coverage);
  function update(){
    var count=0,partial=[],keys=[];
    cache.rows.concat(cache.monthlyRows||[]).forEach(function(row){
      var key=row.date,m=towerPeriodObservation(tower,row,start,end);months[key]=m;
      row.TowerLE_CORR_mm_day=m?m.leET:null;row.TowerLE_CORR_count=m?m.leN:null;row.TowerLE_CORR_missing_pct=m?m.leMissingPct:null;
      row.TowerGPP=m?m.gpp:null;row.TowerGPP2=m?m.secondGpp:null;row.TowerGPP_count=m?m.n:null;
      row.TowerGPP_missing_pct=m?m.missingPct:null;row.TowerGPP2_missing_pct=m?m.secondMissingPct:null;
      if(m&&keys.indexOf(key)<0){keys.push(key);if(months[key]&&(months[key].n||months[key].secondN))count++;}
    });
    monthChoice.items().reset(keys);if(keys.length)monthChoice.setValue(keys[0]);
    var mismatch=!coordinateMode||tower.latitude===null||tower.longitude===null||Math.abs(lat-tower.latitude)>0.00001||Math.abs(lon-tower.longitude)>0.00001;
    var provenance=tower.provenance||{};
    note.setValue(tower.site+': '+count+' periods with GPP; observations averaged within each model timestep. Primary: '+(tower.gppVariable||'GPP')+(tower.secondVariable?'; additional: '+tower.secondVariable:'')+
      '. '+(mismatch?'Verify selected location against tower coordinates. ':'')+'Whole-footprint observations; dates '+start+' through '+end+'. '+
      (partial.length?'Partial-month coverage: '+partial.join(', ')+'. Model points remain full-month estimates. ':'')+
      (provenance.citation?'Citation: '+provenance.citation+' ':'')+(provenance.productId?'Product: '+provenance.productId:''));
    UI._lastNative=carbonRowsToFC(cache.rows);
    UI._lastOutput=carbonRowsToFC(cache.rows);
    cache.refresh.forEach(function(refresh){refresh();});
  }
  if(cache.rows)update();else cache.ready.push(update);
}

// Independent observation requests: never include model graphs in LST reductions.
function loadThermalComparisons(cache, geom, start, end, requestId) {
  var dates=carbonRowsToFC(cache.rows.map(function(r){return {date:r.date,period_days:nativePeriodSeconds(r)/86400};}));
  var first = ee.Date(start);
  var last = ee.Date(end).advance(1,'day');
  function load(name, field, fetch) {
    var note = ui.Label(name + ': loading observations… Model results are ready.');
    addChartNote(note);
    function finish(result, error) {
      if (requestId !== UI._requestId) return;
      if (error || !result) { note.setValue(name + ': unavailable — ' + (error || 'empty response')); return; }
      var byDate = {};
      result.forEach(function(row) { byDate[row.date] = row; });
      cache.rows.concat(cache.monthlyRows||[]).forEach(function(row) {
        var observed = byDate[row.date];
        if (observed) Object.keys(observed).forEach(function(key) { if (key !== 'date') row[key] = observed[key]; });
      });
      var count = cache.rows.filter(function(row) { return carbonFinite(row[field]); }).length;
      var details = field === 'Landsat_LST' ? ' Valid scenes/month: ' + cache.rows.map(function(row) {
        return row.date.slice(0, 7) + ': ' + (row.Landsat_LST_scene_count || 0);
      }).join(', ') : ' Quality-screened observations.';
      note.setValue(name + ': '+UI._lastRun.timestep+' observations; '+count+' periods with quality-screened observations.');
      UI._lastNative = carbonRowsToFC(cache.rows);
      UI._lastOutput=carbonRowsToFC(cache.rows);
      cache.refresh.forEach(function(refresh) { refresh(); });
    }
    try {
      var periodRows=cache.rows.map(function(r){return {date:r.date,period_days:nativePeriodSeconds(r)/86400};});
      var offset=0,all=[];
      function batch() {
        if(requestId!==UI._requestId)return;
        if(offset>=periodRows.length){finish(all);return;}
        var subset=periodRows.slice(offset,offset+30);
        var observations=fetch(carbonRowsToFC(subset),geom,ee.Date(subset[0].date),
          ee.Date(subset[subset.length-1].date).advance(subset[subset.length-1].period_days,'day'));
        observations.toList(observations.size().max(1)).map(function(f){return ee.Feature(f).toDictionary();}).evaluate(function(result,error){
          if(error){finish(null,error);return;}all=all.concat(result||[]);offset+=subset.length;batch();
        });
      }
      batch();
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
  label: 'Prepare CSV download (selected timestep)',
  onClick: function() { onExportCsv(); },
  style: { stretch: 'horizontal', padding: '6px',
           margin: '0px 4px 0px 0px' }
});
UI.exportPanel.add(UI.exportCsv);
UI.csvDownloadLink = ui.Label({value: '', style: {shown: false, color: '#1565C0'}});
UI.exportPanel.add(UI.csvDownloadLink);
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
  UI._soilSearch=null;
  if(UI.map&&UI._soilLayers)UI._soilLayers.forEach(function(layer){UI.map.layers().remove(layer);});
  UI._soilLayers=[];
  UI._requestId++;
  UI._monthlyExportReady = false;
  UI._lastNative=null; UI._lastAnnual=null; UI._lastCarbonMonthly=null; UI._lstNDVIRows=null;
  UI._downloadRequestId = (UI._downloadRequestId || 0) + 1;
  UI.csvDownloadLink.style().set({shown: false});
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
  if(UI.analysisMode.getValue()==='lst_ndvi'){runLSTNDVI(geom);return;}
  if(UI.analysisMode.getValue()==='soil_search'){runSoilSearch(geom);return;}
  if (!startDate || !endDate) {
    setStatus('Please fill both start and end dates (YYYY-MM-DD).');
    return;
  }

  setStatus('Assembling forcing data…');

  try {
    var manualParameters = readVegetationParameters();
    if(UI.faparSelector.getValue()==='KNDVI'){
      manualParameters.kndvi_sigma=Number(UI.kndviSigma.getValue());
      if(!(isFinite(manualParameters.kndvi_sigma)&&manualParameters.kndvi_sigma>0))throw new Error('kNDVI sigma must be finite and positive.');
    }
    if (carbonDateMillis(endDate) < carbonDateMillis(startDate))
      throw new Error('End date must be on or after start date.');
    UI._runWeatherPlan=null;
    if(UI.useTowerWeather.getValue()) {
      if(!runTower || !towerCoordinateMode || !carbonFinite(runTower.latitude) || !carbonFinite(runTower.longitude) || Math.abs(towerLat-runTower.latitude)>0.00001 || Math.abs(towerLon-runTower.longitude)>0.00001)
        throw new Error('Use the imported tower coordinates for a tower-weather run.');
      UI._runWeatherPlan=selectedTowerWeatherPlan(runTower,startDate,endDate,UI.timestepSelector.getValue());
      if(!UI._runWeatherPlan.segments.length){
        UI.weatherNote.setValue('Tower weather unavailable:\n'+UI._runWeatherPlan.rejected.join('\n'));
        throw new Error('No complete tower-weather periods. '+UI._runWeatherPlan.rejected.slice(0,2).join(' | ')+'. See the weather note for all periods.');
      }
      UI.weatherNote.setValue('Tower weather selected ('+UI._runWeatherPlan.cadence+'). Unavailable periods: '+(UI._runWeatherPlan.missing.join(', ')||'none')+'. Bucket forcing is interpolated between valid periods; missing carbon remains unavailable. CO2 and satellite vegetation remain unchanged. Supplied gap-filled fields are retained.');
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
    var runDataConfig={kndviRadiusKm:Number(UI.kndviRadius.getValue()),rtObjective:UI.rtObjective.getValue(),rtDiffuse:UI.rtDiffuse.getValue().trim()===''?NaN:Number(UI.rtDiffuse.getValue()),rtBoxKm:Number(UI.rtBoxKm.getValue()),refreshInputs:!UI.reuseInputs.getValue(),timestep:UI.timestepSelector.getValue(),weatherPlan:UI._runWeatherPlan,weatherSource:UI.weatherSelector.getValue(),satelliteSource:UI.satelliteSelector.getValue()};
    runForGeometryAsync(geom, startDate, endDate, mergeParams(manualParameters, {fiparMethod: UI.fiparSelector.getValue(), useSharedFPAR: UI.sharedFPAR.getValue(), fitVIMinP05: UI.faparSelector.getValue() === 'FIPAR_P05', fitNDVIminGE: UI.faparSelector.getValue() === 'FIPAR_GE' || UI.faparSelector.getValue() === 'NDVI_GE_OFFSET', faparMethod: (UI.faparSelector.getValue()==='KNDVI'||UI.faparSelector.getValue()==='KNDVI_SPATIAL'||rtAbsorptionMethod(UI.faparSelector.getValue())) ? UI.faparSelector.getValue() : UI.faparSelector.getValue() === 'PROSAIL' ? 'PROSAIL' : UI.faparSelector.getValue() === 'NIRV' ? 'NIRV' : UI.faparSelector.getValue() === 'NDVI_RED_EXP' ? 'NDVI_RED_EXP' : UI.faparSelector.getValue() === 'NDVI_GE_OFFSET' ? 'NDVI_GE_OFFSET' : UI.faparSelector.getValue() === 'GE_NDVI_MIN' ? 'GE_NDVI_MIN' : UI.faparSelector.getValue() === 'NDVI' ? 'NDVI' : UI.faparSelector.getValue() === 'NDVI_BLEND' ? 'NDVI_BLEND' : UI.faparSelector.getValue() === 'FGE' ? 'FGE' : (UI.faparSelector.getValue().indexOf('FIPAR') === 0 ? 'FIPAR' : 'EVI2'), fitEVI2: UI.faparSelector.getValue() === 'EVI2_FITTED' || UI.faparSelector.getValue() === 'EVI2_FITTED_MIN', powerFixedEVI2: UI.faparSelector.getValue() === 'EVI2_FIXED_POWER', fitEVI2maxFixed: UI.faparSelector.getValue() === 'EVI2_FITTED_MIN'}), runDataConfig, function(output,modelError) {
    if(runRequestId!==UI._requestId)return;
    if(modelError){UI.runButton.setDisabled(false);setStatus('Model error: '+modelError);return;}
    try {
    UI._lastOutput = output;
    UI._lastRun = {timestep:runDataConfig.timestep,start: startDate, end: endDate, site: regionName, regionName: regionName, weatherSource:runDataConfig.weatherPlan?'tower':runDataConfig.weatherSource,satelliteSource:runDataConfig.satelliteSource};

    regionName += UI._runWeatherPlan ? '\nTower weather ('+UI._runWeatherPlan.cadence+')' : '\n'+WEATHER_SOURCES[runDataConfig.weatherSource].label;
    regionName+='\nTimestep: '+runDataConfig.timestep;
    // Center the map on the active geometry so the user sees where they ran.
    if (UI.map) UI.map.centerObject(geom);

    setStatus('Building chart (server-side)...');

    var conductanceChart;
    var chartCache = {rows: null, ready: [], refresh: []};
    var c4Chart = function(rows,monthlyRows) {
      chartCache.monthlyRows=monthlyRows;
      chartCache.rows = rows;
      output = carbonRowsToFC(rows);
      addRemainingDiagnosticCharts();
      if (UI._lastRun.timestep === 'monthly') {
        addResultChart(ui.Chart.feature.byFeature(output, 'date', ['satellite_missing_pct'])
          .setOptions({
            title: regionChartTitle('Monthly satellite missing data', regionName),
            hAxis: {title: 'Date'},
            vAxis: {title: 'Missing NDVI/EVI2 days (%)', viewWindow: {min: 0, max: 100}},
            legend: {position: 'none'}, lineWidth: 1, pointSize: 3, interpolateNulls: false
          }));
        addChartNote(ui.Label('Monthly NDVI/EVI2 use the median of available daily values after gap filling. Missing percentage counts days still unavailable.'));
      }

      var resistance = rows.map(function(row) { return row.rs; }).filter(carbonFinite);
      var conductanceOptions = conductanceChart.getOptions();
      conductanceOptions.vAxes[1].viewWindow.max = resistance.length ? 10 * Math.min.apply(null, resistance) : 1000;
      conductanceChart.setOptions(conductanceOptions);
      chartCache.ready.forEach(function(ready) { ready(); });
      chartCache.ready = [];
      addOptionalSeriesChart(chartCache, ['FIPAR', 'FAPAR', 'C4frac', 'NDVI', 'EVI2', 'LAI', 'GE', 'GE_normalized', 'fge', 'NIRv'],
        ['FIPAR', 'FAPAR', 'C4 fraction', 'NDVI', 'EVI2', 'Green LAI', 'GE', 'GE′', 'fge', 'NIRv'],
        ['#D55E00', '#000000', '#1976D2', '#2E7D32', '#8E24AA', '#EF6C00', '#00838F', '#795548', '#D81B60', '#43A047'], [3, 4, 5, 6, 7, 8, 9], {
          series: {0: {color: '#D55E00'}, 1: {color: '#000000', lineDashStyle: [6, 4]},
            2: {color: '#1976D2', lineDashStyle: [6, 4]},
            3: {color: '#2E7D32'}, 4: {color: '#8E24AA'},
            5: {color: '#EF6C00', targetAxisIndex: 1}},
          title: regionChartTitle('FAPAR, FIPAR & C4 Fraction', regionName), chartArea: {top: 80},
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
      'Tower LE_CORR: available full-day flux mean × 86400 / 2450000 (mm/day). Includes night and negative fluxes; missing percentages are in Tower coverage by timestep.':
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
    addChartNote(ui.Label(annualClimateNote(aggregateModelMonths(chartCache.rows), startDate, endDate, 'Tdaily')));
    loadThermalComparisons(chartCache, geom, startDate, endDate, runRequestId);

    }
    UI.resultsPanel.style().set({shown: true});
    UI.resultsOverlay.setExpanded(true);
    setStatus('Calculating daily carbon totals and loading annual RAP cover...');
    loadCarbonProducts(output, geom, startDate, endDate, runRequestId, regionName, c4Chart);
    } catch(renderError){UI.runButton.setDisabled(false);setStatus('Chart error: '+renderError.message);}
    });
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

// Direct table downloads work in published Earth Engine Apps without Tasks.
function prepareCsvDownload(collection, cadence) {
  var run = UI._lastRun;
  var requestId = UI._requestId;
  var downloadId = (UI._downloadRequestId || 0) + 1;
  UI._downloadRequestId = downloadId;
  UI.csvDownloadLink.style().set({shown: false});
  var name = ('acclim_' + run.site + '_' + cadence + '_' + run.start + '_to_' + run.end)
    .replace(/[^A-Za-z0-9_-]/g, '_');
  var tabular = collection.map(function(f) {
    return f.setGeometry(null).set({export_site: run.site,
      export_start: run.start, export_end: run.end, export_weather_source: run.weatherSource,export_satellite_source:run.satelliteSource});
  });
  setStatus('Preparing ' + cadence + ' CSV download...');
  function finish(url, error) {
    if (requestId !== UI._requestId || downloadId !== UI._downloadRequestId) return;
    if (error || !url) {
      setStatus('CSV download could not be prepared: ' + (error && error.message || error || 'No download URL returned.'));
      return;
    }
    UI.csvDownloadLink.setValue('Download ' + cadence + ' CSV');
    UI.csvDownloadLink.setUrl(url);
    UI.csvDownloadLink.style().set({shown: true});
    setStatus('CSV ready. Click Download ' + cadence + ' CSV in the results export panel.');
  }
  try {
    tabular.getDownloadURL({format: 'csv', filename: name, callback: finish});
  } catch (error) { finish(null, error); }
}

function onExportCarbon(cadence) {
  var collection = cadence === 'daily' ? UI._lastDaily : UI._last16day;
  if (!collection || !UI._lastRun) { setStatus('Run the model and wait for carbon products before exporting.'); return; }
  prepareCsvDownload(collection, cadence);
}

function onExportCsv() {
  if (!UI._lastOutput || !UI._lastRun || !UI._monthlyExportReady) {
    setStatus('Run the model and wait for results before exporting.');
    return;
  }
  prepareCsvDownload(UI._lastOutput, UI._lastRun.timestep||'monthly');
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
