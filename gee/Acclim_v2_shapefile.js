/**
 * ACCLIM research web app — selectable weather and satellite reflectance.
 * TerraClimate: monthly physiology. ERA5-Land / GRIDMET: daily physiology.
 * Native output is retained; daily fluxes are aggregated for monthly comparisons.
 * Tower-weather input currently uses the importer's monthly summaries.
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
  faparMethod: 'FIPAR',    // Direct FIPAR by default; scaled EVI2 remains selectable
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
  fDexp:          0.5,      // VPD attenuation exponent (VPDday_mean in kPa)
  GPPmean2GPPmax: 0.67,      // ratio of mean GPP to GPPmax asymptote
  Oi:             20900,    // leaf internal O2 partial pressure (Pa)
  QRC3:           12,       // C3 quantum requirement (mol photons / mol CO2)
  QRC4:           16.67,    // C4 quantum requirement
  EVI2min:        0.0,     // EVI2 of bare ground
  EVI2max:        0.83,     // EVI2 of full green canopy
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
  ndvi_min: 0,
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
function satelliteImages(geom,start,end,key) {
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
      var red,nir,mask=ee.Image.constant(1);
      if(spec.kind==='landsat') {
        red=img.select('red');nir=img.select('nir');
        dateDay=dateDay.add(spec.days/2); // Catalog composites have no selected-pixel date band.
        sampleDay=d.advance(spec.days/2,'day');
      } else if(spec.kind==='hls') {
        red=img.select('B4');nir=img.select(id.indexOf('HLSS30')>=0?'B8A':'B5');
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
      mask=mask.and(red.gte(0)).and(nir.gte(0)).and(red.lte(1)).and(nir.lte(1)).and(red.add(nir).gt(0));
      return red.rename('red').addBands(nir.rename('nir')).addBands(dateDay.rename('observation_day'))
        .updateMask(mask).set({'system:time_start':img.get('system:time_start'),sample_millis:sampleDay.millis()});
    });
    merged=merged.merge(images);
  });
  return merged;
}
function sampleSatellite(geom,start,end,key,unpadded) {
  var spec=SATELLITE_SOURCES[key], pad=spec.gap+32;
  var a=ee.Date(start).advance(unpadded?0:-pad,'day'), b=ee.Date(end).advance(unpadded?0:pad,'day');
  var images=satelliteImages(geom,a,b.advance(1,'day'),key);
  // Group pixels by their actual observation day BEFORE reducing spatially.
  // Same-day Terra/Aqua or HLS values are averaged with a common red/NIR mask.
  // Sample existing images only, rather than rebuilding every composite on every calendar day.
  var sampled=ee.FeatureCollection(images.toList(images.size().max(1)).map(function(item) {
    var img=ee.Image(item),startDay=ee.Date(img.get('sample_millis'));
    var count=spec.kind==='vi'?ee.Date.fromYMD(startDay.get('year').add(1),1,1).difference(startDay,'day').min(16):ee.Number(1);
    return ee.List.sequence(0,count.subtract(1)).map(function(i) {
      var d=startDay.advance(i,'day');
      var pair=img.select(['red','nir']).updateMask(img.select('observation_day').eq(d.millis().divide(86400000)));
      var values=pair.reduceRegion({reducer:ee.Reducer.mean(),geometry:geom,scale:spec.scale,maxPixels:1e10,tileScale:4});
      return ee.Feature(null,values).set({date:d.format('YYYY-MM-dd'),millis:d.millis()});
    });
  }).flatten()).filter(ee.Filter.notNull(['red','nir']));
  var observed=sampled.distinct(['millis']).map(function(f) {
    var same=sampled.filter(ee.Filter.eq('millis',f.get('millis')));
    return f.set({red:same.aggregate_mean('red'),nir:same.aggregate_mean('nir'),satellite_observation_count:same.size()});
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
    return f.set({red:ee.Algorithms.If(valid,red,null),nir:ee.Algorithms.If(valid,nir,null),
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
  var cadence=config.weatherPlan?config.weatherPlan.cadence:(WEATHER_SOURCES[key].daily?'daily':'monthly');
  var daily=cadence==='daily';
  var cal=config.weatherPlan?ee.FeatureCollection(config.weatherPlan.records.map(function(e){return towerWeatherFeature(e.data,geom);})):sourceCalendar(start,end,daily);
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
  var satellite=config.alignedSatellite||alignedSatellite(geom,first,lastDate,sat,config.sampledSatellite), p=ACCLIM_PARAMS;
  var lat=ee.Number(geom.centroid({maxError:1}).coordinates().get(1));
  var raw=cal.map(function(f) {
    var t=f.getNumber('millis'), days=f.getNumber('period_days');
    var matches=weather.filter(ee.Filter.eq('millis',t));
    var wf=ee.Feature(ee.Algorithms.If(matches.size(),matches.first(),ee.Feature(null,{})));
    var vi=satellite.filter(ee.Filter.gte('millis',t)).filter(ee.Filter.lt('millis',t.add(days.multiply(86400000))));
    var valid=vi.filter(ee.Filter.notNull(['red','nir']));
    var props=ee.Dictionary({});
    // Unavailable satellite periods have null means; Dictionary.set rejects null.
    // Omit those properties so the existing completeness check flags the period.
    ['red','nir','NDVI','EVI2','satellite_missing_pct','satellite_unavailable_pct'].forEach(function(k){
      var mean=vi.aggregate_mean(k);
      props=ee.Dictionary(ee.Algorithms.If(ee.Algorithms.IsEqual(mean,null),props,props.set(k,mean)));
    });
    var dl=ee.Number(ee.List.sequence(0,days.subtract(1)).map(function(i){
      return astronomicalDaylength(lat,doyFromDate(ee.Date(t).advance(i,'day')));
    }).reduce(ee.Reducer.mean()));
    var complete=wf.toDictionary(['Tmin','Tmax','Tdaily','SWRdaily','ea','Pdaily']).size().eq(6).and(valid.size().eq(days)).and(dl.gt(0));
    if(config.weatherPlan)complete=complete.and(ee.Number(ee.Algorithms.If(ee.Algorithms.IsEqual(wf.get('tower_weather_valid'),null),0,wf.get('tower_weather_valid'))).eq(1));
    return f.setGeometry(geom).set(wf.toDictionary()).set(props).set({
      daylength:dl,daylight_seconds:dl.multiply(3600).multiply(days),
      doy:doyFromDate(ee.Date(t)),dt:days,ca_daily:co2ForDate(ee.Date(t)),
      model_cadence:cadence,weather_source:config.weatherPlan?'tower':key,
      model_valid:complete,weather_missing_pct:ee.Number(7).subtract(wf.toDictionary(['Tmin','Tmax','Tdaily','SWRdaily','ea','Pdaily','precip']).size()).multiply(100/7),
      precipitation_missing_pct:ee.Algorithms.If(ee.Algorithms.IsEqual(wf.get('precip'),null),100,0),
      satellite_source:sat,satellite_scale_m:SATELLITE_SOURCES[sat].scale,
      satellite_date_basis:ee.Feature(vi.first()).get('satellite_date_basis'),
      satellite_status:ee.Algorithms.If(valid.size().eq(days),ee.Algorithms.If(vi.aggregate_sum('satellite_missing_pct'),'interpolated','observed'),'unavailable'),
      satellite_gap_days:vi.aggregate_max('satellite_gap_days'),satellite_left_date:ee.Feature(vi.sort('millis').first()).get('satellite_left_date'),
      satellite_right_date:ee.Feature(vi.sort('millis',false).first()).get('satellite_right_date'),satellite_max_gap_days:SATELLITE_SOURCES[sat].gap});
  });
  var ready=raw.filter(ee.Filter.eq('model_valid',1)).map(function(f) {
    var tday=f.getNumber('Tmax').multiply(0.75).add(f.getNumber('Tmin').multiply(0.25));
    var ea=f.getNumber('ea'), albedo=f.getNumber('red').add(f.getNumber('nir')).divide(2).max(0).min(1);
    var fipar=computeFIPAR(f.getNumber('NDVI'),p);
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
    var FIPAR = computeFIPAR(NDVI, p);
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
function computeFIPAR(ndvi, p) {
  return ee.Number(ndvi).subtract(p.ndvi_min).max(0)
    .divide(p.ndvi_max - p.ndvi_min).pow(p.ndvi_exp);
}

/**
 * Compute mean-state-dependent diagnostics from a state dictionary +
 * today's instantaneous values. Pure function of inputs; no side effects.
 * Returns a dictionary of all derived quantities to attach to the output feature.
 */
function photosyntheticFAPAR(instant, p) {
  if (p.faparMethod === 'FIPAR') return ee.Number(instant.get('FIPAR'));
  var FAPAR_in = ee.Number(instant.get('EVI2'));
  var FAPAR = ee.Number(instant.get('NDVI')).pow(2);
  if (p.faparMethod === 'EVI2' && !p.fitEVI2) {
    FAPAR = FAPAR_in.subtract(p.EVI2min).divide(p.EVI2max - p.EVI2min).max(0).min(1);
  } else if (p.faparMethod === 'EVI2') {
    var fitRange = ee.Number(ee.Algorithms.If(p.faparFit.valid, p.faparFit.max.subtract(p.faparFit.min), 1));
    var scaledEVI2 = FAPAR_in.subtract(p.faparFit.min)
      .divide(fitRange).max(0).min(1);
    FAPAR = ee.Number(ee.Algorithms.If(p.faparFit.valid, scaledEVI2, FAPAR));
  }

  return FAPAR;
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
  var faparFiparValid = FIPAR.gt(1e-9);
  var faparFiparRatio = (p.useSharedFPAR || p.faparMethod === 'FIPAR') ? ee.Number(1) : FAPAR.divide(FIPAR.max(1e-9));
  var zCoef = ee.Number(1.5).multiply(p.GPPmean2GPPmax).multiply(2500).multiply(18)
    .divide(ee.Number(p.QRC3).multiply(p.PT_coef).multiply(ee.Number(1).add(p.fresp)));
  var zC3 = zCoef.multiply(
    s_mean.add(p.gamma_psy).multiply(VPDday_mean).multiply(PARday_mean).multiply(1000)
  ).divide(s_mean.multiply(Rnday_mean).multiply(Pmean)).multiply(faparFiparRatio);

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

  // GPPmax (umolCO2/m2s) = QY * FAPAR * PARday_mean
  // fD = min(1, 1 / VPDday_mean^fDexp); max(1) avoids division by zero.
  var fD = ee.Number(1).divide(VPDday_mean.pow(p.fDexp).max(1));
  var GPPmax_C3 = QYC3_mean.multiply(FAPAR).multiply(PARday_mean).multiply(fD);
  var GPPmax_C4 = QYC4_mean.multiply(FAPAR).multiply(PARday_mean).multiply(fD);

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
    FAPAR:            FAPAR,
    FIPAR:            FIPAR,
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
  var instant = ee.Dictionary({EVI2: EVI2, NDVI: NDVI, FIPAR: f.getNumber('FIPAR')});

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
function fitEVI2FAPAR(forcingFC, params) {
  var p = mergeParams(ACCLIM_PARAMS, params || {});
  var pairs = forcingFC.filter(ee.Filter.notNull(['NDVI', 'EVI2']))
    .filter(ee.Filter.gte('NDVI', -1)).filter(ee.Filter.lte('NDVI', 1))
    .map(function(f) { return f.set('FIPAR_fit_target', computeFIPAR(f.getNumber('NDVI'), p)); });
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
  var shared = photosyntheticFAPAR(f, p);
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
  forcingFC=forcingFC.map(withPeriodTiming);
  var missingForcing = forcingFC.filter(ee.Filter.eq('model_valid',0));
  forcingFC = forcingFC.filter(ee.Filter.notNull(['Tmin','Tmax','SWRdaily','ea','NDVI','EVI2','Rnday']));
  if (p.faparMethod !== 'NDVI' && p.faparMethod !== 'EVI2' && p.faparMethod !== 'FIPAR') throw new Error('Unknown FAPAR method.');
  if (p.faparMethod === 'EVI2' && p.fitEVI2 && !p.faparFit) p.faparFit = fitEVI2FAPAR(forcingFC, p);
  if (p.faparMethod === 'EVI2' && p.fitEVI2 && p.fitEVI2maxFixed) p.faparFit = fixedUpperEVI2Fit(p.faparFit);
  if (p.faparMethod === 'EVI2' && !p.fitEVI2 &&
      (!isFinite(p.EVI2min) || !isFinite(p.EVI2max) || p.EVI2max <= p.EVI2min)) {
    throw new Error('Fixed EVI2 bounds must be finite with EVI2max > EVI2min.');
  }
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
        'precipitation_missing_pct','satellite_source','satellite_scale_m','satellite_missing_pct','satellite_unavailable_pct','satellite_status',
        'satellite_date_basis','satellite_gap_days','satellite_left_date','satellite_right_date','satellite_max_gap_days',
        'Tmin','Tmax','Tdaily','SWRdaily','ea','Pdaily','NDVI','EVI2','red','nir',
        'precip','precip_units','fSM_mean','fSM_et_mean']);
    return ee.Feature(f.geometry(),f.toDictionary(keep));
  });
  return bucketOutput.filter(ee.Filter.notNull(['PredGPPmean'])).map(function(f) {
    f = f.set('FAPAR_method_requested', p.faparMethod);
    if (p.faparMethod === 'EVI2' && !p.fitEVI2) {
      f = f.set({FAPAR_method: 'EVI2', FAPAR_fit_status: 'fixed bounds',
        FAPAR_EVI2min: p.EVI2min, FAPAR_EVI2max: p.EVI2max});
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

function runForGeometryAsync(geom,start,end,params,config,done) {
  var request=UI._requestId, cfg=mergeParams(config,{});
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
      if(!rows.some(function(r){return r.model_valid===1||r.model_valid===true;})){done(forcing);return;}
      setStatus('Calculating native '+rows[0].model_cadence+' physiology...');
      materialize(runACCLIM(forcing,mergeParams(params,{deferBuckets:true})),function(baselineRows,baseline) {
        setStatus('Fitting soil-water buckets and bridging missing forcing...');
        // Fit once over the full trajectory; only independent final diagnostics are batched.
        var fittedRows=fitNativeBuckets(baselineRows,params), finalRows=[], offset=0;
        function finishBatch() {
          if(request!==UI._requestId)return;
          if(offset>=fittedRows.length){done(ee.FeatureCollection(finalRows.map(function(r){return ee.Feature(geom,r);})).sort('millis'));return;}
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
  var monthly=cfg.weatherPlan?cfg.weatherPlan.cadence==='monthly':!WEATHER_SOURCES[cfg.weatherSource||'TERRACLIMATE'].daily;
  var first=monthly?start.slice(0,7)+'-01':start;
  var last=monthly?carbonDateString(Date.UTC(Number(end.slice(0,4)),Number(end.slice(5,7)),1)-86400000):end;
  if(cfg.weatherPlan&&cfg.weatherPlan.records.length){
    first=cfg.weatherPlan.records[0].data.date||cfg.weatherPlan.records[0].data.month+'-01';
    var lastRecord=cfg.weatherPlan.records[cfg.weatherPlan.records.length-1].data;
    if(lastRecord.period_days)last=carbonDateString(carbonDateMillis(lastRecord.date)+(lastRecord.period_days-1)*86400000);
  }
  var spec=SATELLITE_SOURCES[cfg.satelliteSource||'MODIS16'],pad=spec.gap+32;
  var cursor=carbonDateMillis(first)-pad*86400000, stop=carbonDateMillis(last)+pad*86400000, sampled=[];
  function satelliteBatch() {
    if(request!==UI._requestId)return;
    if(cursor>stop) {
      cfg.alignedSatellite=carbonRowsToFC(alignSatelliteRows(sampled,first,last,cfg.satelliteSource||'MODIS16'));
      simulate();return;
    }
    var next=Math.min(cursor+90*86400000,stop+86400000),a=carbonDateString(cursor),b=carbonDateString(next-86400000);
    setStatus('Sampling satellite reflectance: '+a+' to '+b+'...');
    materialize(sampleSatellite(geom,a,b,cfg.satelliteSource||'MODIS16',true),function(rows) {
      sampled=sampled.concat(rows);cursor=next;satelliteBatch();
    });
  }
  if(cfg.weatherPlan)satelliteBatch();
  else {
    setStatus('Sampling selected weather data...');
    materialize(sampleWeather(geom,first,last,cfg.weatherSource||'TERRACLIMATE'),function(rows,weather) {
      cfg.sampledWeather=weather;satelliteBatch();
    });
  }
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
    if(!groups[k])groups[k]={millis:r.millis,red:0,nir:0,n:0};
    groups[k].red+=r.red*n;groups[k].nir+=r.nir*n;groups[k].n+=n;
  });
  var obs=Object.keys(groups).map(function(k){var g=groups[k];return {millis:g.millis,red:g.red/g.n,nir:g.nir/g.n};})
    .sort(function(a,b){return a.millis-b.millis;}),rows=[],j=0;
  for(var t=carbonDateMillis(start);t<=carbonDateMillis(end);t+=86400000) {
    while(j<obs.length && obs[j].millis<t)j++;
    var r=obs[j],l=r&&r.millis===t?r:obs[j-1];
    var gap=l&&r?(r.millis-l.millis)/86400000:null,valid=l&&r&&gap<=spec.gap,same=valid&&gap===0;
    var row={date:carbonDateString(t),millis:t,satellite_source:key,satellite_scale_m:spec.scale,
      satellite_date_basis:spec.kind==='vi'?'pixel observation date':spec.days===1?'observation date':'composite midpoint',
      satellite_status:valid?(same?'observed':'interpolated'):'unavailable',satellite_missing_pct:same?0:100,
      satellite_unavailable_pct:valid?0:100,satellite_gap_days:valid?gap:null,satellite_max_gap_days:spec.gap,
      satellite_left_date:valid?carbonDateString(l.millis):null,satellite_right_date:valid?carbonDateString(r.millis):null};
    if(valid) {
      var w=gap?(t-l.millis)/(r.millis-l.millis):0;
      row.red=l.red+(r.red-l.red)*w;row.nir=l.nir+(r.nir-l.nir)*w;
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
  var m={sum:0,n:0,secondSum:0,secondN:0,daylight:0,leSum:0,leN:0};
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
      if(variable==='Tdaily' && row.model_cadence==='daily' && row.Tdaily_valid_days!==n)complete=false;
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
  // EVI2 paired with native-timestep Priestley-Taylor daytime energy flux.
  var eviPTPairs = rows.filter(function(r) {
    return carbonFinite(r.PT_mean) && carbonFinite(r.EVI2);
  });
  var eviPTTable = {cols: [
    {id: 'PT_mean', label: 'PT_mean (W/m²)', type: 'number'},
    {id: 'EVI2', label: 'EVI2', type: 'number'}],
    rows: eviPTPairs.map(function(r) { return {c: [{v: r.PT_mean}, {v: r.EVI2}]}; })};
  addResultChart(ui.Chart(eviPTTable).setChartType('ScatterChart').setOptions({
    title: regionChartTitle('EVI2 vs PT_mean', regionName), chartArea: {top: 80},
    hAxis: {title: 'PT_mean (W/m²; daytime)'}, vAxis: {title: 'EVI2'},
    legend: {position: 'none'}, colors: ['#8E24AA'],
    pointSize: 4, pointShape: 'circle', lineWidth: 0
  }));
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
        var depths = effectiveDepthFields(row, sandPercent);
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
      UI._lastOutput = carbonRowsToFC(monthlyExportRows); // monthly export; charts use native rows
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
      c4Chart(result.monthly,monthlyExportRows);
      if (bucketInfo.FAPAR_method_requested === 'EVI2') {
        addChartNote(ui.Label(bucketInfo.FAPAR_method === 'EVI2'
          ? 'FAPAR = EVI2 scaled between EVI2min=' + bucketInfo.FAPAR_EVI2min.toFixed(3) +
            ' and EVI2max=' + bucketInfo.FAPAR_EVI2max.toFixed(3) +
            ' (' + bucketInfo.FAPAR_fit_status + ')'
          : 'FAPAR: NDVI² fallback. Fit needs at least 3 paired months, a positive slope, and minimum below maximum.'));
      }
      addWeatherCharts(output, result.monthly, regionName);
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
  {label: 'FAPAR = FIPAR', value: 'FIPAR'},
  {label: 'FAPAR = EVI2 scaled to FIPAR', value: 'EVI2_FITTED'},
  {label: 'FAPAR = EVI2 FIPAR-fitted min; fixed max 0.6', value: 'EVI2_FITTED_MIN'},
  {label: 'FAPAR = EVI2 scaled to fixed bounds', value: 'EVI2_FIXED'}], value: 'FIPAR',
  style: {stretch: 'horizontal', margin: '4px 0px'}});
UI.inputPanel.add(ui.Label('FAPAR method (applied on Run ACCLIM)'));
UI.inputPanel.add(UI.faparSelector);
UI.sharedFPAR = ui.Checkbox({label: 'Use FAPAR for both photosynthesis and water fluxes', value: false});
UI.inputPanel.add(UI.sharedFPAR);
setupTowerImport();
UI.weatherSelector=ui.Select({items:Object.keys(WEATHER_SOURCES).map(function(k){return {label:WEATHER_SOURCES[k].label,value:k};}),value:'TERRACLIMATE'});
UI.satelliteSelector=ui.Select({items:Object.keys(SATELLITE_SOURCES).map(function(k){return {label:SATELLITE_SOURCES[k].label,value:k};}),value:'MODIS16'});
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
    weather_source:'tower',weather_energy_basis:'LE_CORR + H_CORR',weather_status:'complete',weather_expected_intervals:(row.expected===undefined?null:row.expected)});
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
    cache.rows.concat(cache.monthlyRows||[]).forEach(function(row){
      var key=row.date.slice(0,7),m=towerPeriodObservation(tower,row,start,end);
      row.TowerLE_CORR_mm_day=m?m.leET:null;row.TowerLE_CORR_count=m?m.leN:null;row.TowerLE_CORR_missing_pct=m?m.leMissingPct:null;
      row.TowerGPP=m?m.gpp:null;row.TowerGPP2=m?m.secondGpp:null;row.TowerGPP_count=m?m.n:null;
      row.TowerGPP_missing_pct=m?m.missingPct:null;row.TowerGPP2_missing_pct=m?m.secondMissingPct:null;
      if(m&&keys.indexOf(key)<0){keys.push(key);if(months[key]&&(months[key].n||months[key].secondN))count++;}
    });
    monthChoice.items().reset(keys);if(keys.length)monthChoice.setValue(keys[0]);
    var mismatch=!coordinateMode||tower.latitude===null||tower.longitude===null||Math.abs(lat-tower.latitude)>0.00001||Math.abs(lon-tower.longitude)>0.00001;
    var provenance=tower.provenance||{};
    note.setValue(tower.site+': '+count+' months with GPP; observations averaged within each model timestep. Primary: '+(tower.gppVariable||'GPP')+(tower.secondVariable?'; additional: '+tower.secondVariable:'')+
      '. '+(mismatch?'Verify selected location against tower coordinates. ':'')+'Whole-footprint observations; dates '+start+' through '+end+'. '+
      (partial.length?'Partial-month coverage: '+partial.join(', ')+'. Model points remain full-month estimates. ':'')+
      (provenance.citation?'Citation: '+provenance.citation+' ':'')+(provenance.productId?'Product: '+provenance.productId:''));
    UI._lastNative=carbonRowsToFC(cache.rows);
    if(cache.monthlyRows)UI._lastOutput=carbonRowsToFC(cache.monthlyRows);
    cache.refresh.forEach(function(refresh){refresh();});
  }
  if(cache.rows)update();else cache.ready.push(update);
}

// Independent observation requests: never include model graphs in LST reductions.
function loadThermalComparisons(cache, geom, start, end, requestId) {
  var dates=sourceCalendar(start,end,false).select(['date']);
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
      cache.rows.concat(cache.monthlyRows||[]).forEach(function(row) {
        var observed = byDate[row.date];
        if (observed) Object.keys(observed).forEach(function(key) { if (key !== 'date') row[key] = observed[key]; });
      });
      var count = cache.rows.filter(function(row) { return carbonFinite(row[field]); }).length;
      var details = field === 'Landsat_LST' ? ' Valid scenes/month: ' + cache.rows.map(function(row) {
        return row.date.slice(0, 7) + ': ' + (row.Landsat_LST_scene_count || 0);
      }).join(', ') : ' Quality-screened observations.';
      note.setValue(name + ': monthly observations; '+count+' timestamps match native model points. Quality-screened observations.');
      UI._lastNative = carbonRowsToFC(cache.rows);
      if(cache.monthlyRows)UI._lastOutput=carbonRowsToFC(cache.monthlyRows);
      cache.refresh.forEach(function(refresh) { refresh(); });
    }
    try {
      var observations = fetch(dates, geom, first, last);
      observations.toList(observations.size().max(1)).map(function(f) { return ee.Feature(f).toDictionary(); }).evaluate(finish);
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
  label: 'Prepare monthly CSV download',
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
UI.exportPanel.add(ui.Button({label:'Export native model CSV',onClick:function(){if(UI._lastNative)prepareCsvDownload(UI._lastNative,'native');}}));
UI.exportPanel.add(ui.Button({label:'Export annual carbon CSV',onClick:function(){if(UI._lastAnnual)prepareCsvDownload(UI._lastAnnual,'annual_carbon');}}));
UI.exportPanel.add(ui.Button({label:'Export monthly carbon CSV',onClick:function(){if(UI._lastCarbonMonthly)prepareCsvDownload(UI._lastCarbonMonthly,'monthly_carbon');}}));
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
    if (carbonDateMillis(endDate) < carbonDateMillis(startDate))
      throw new Error('End date must be on or after start date.');
    UI._runWeatherPlan=null;
    if(UI.useTowerWeather.getValue()) {
      if(!runTower || !towerCoordinateMode || !carbonFinite(runTower.latitude) || !carbonFinite(runTower.longitude) || Math.abs(towerLat-runTower.latitude)>0.00001 || Math.abs(towerLon-runTower.longitude)>0.00001)
        throw new Error('Use the imported tower coordinates for a tower-weather run.');
      UI._runWeatherPlan=towerWeatherPlan(runTower,startDate,endDate);
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
    var runDataConfig={weatherPlan:UI._runWeatherPlan,weatherSource:UI.weatherSelector.getValue(),satelliteSource:UI.satelliteSelector.getValue()};
    runForGeometryAsync(geom, startDate, endDate, {useSharedFPAR: UI.sharedFPAR.getValue(), faparMethod: UI.faparSelector.getValue() === 'FIPAR' ? 'FIPAR' : 'EVI2', fitEVI2: UI.faparSelector.getValue() === 'EVI2_FITTED' || UI.faparSelector.getValue() === 'EVI2_FITTED_MIN', fitEVI2maxFixed: UI.faparSelector.getValue() === 'EVI2_FITTED_MIN'}, runDataConfig, function(output,modelError) {
    if(runRequestId!==UI._requestId)return;
    if(modelError){UI.runButton.setDisabled(false);setStatus('Model error: '+modelError);return;}
    try {
    UI._lastOutput = output;
    UI._lastRun = {start: startDate, end: endDate, site: regionName, regionName: regionName, weatherSource:runDataConfig.weatherPlan?'tower':runDataConfig.weatherSource,satelliteSource:runDataConfig.satelliteSource};

    regionName += UI._runWeatherPlan ? '\nTower weather ('+UI._runWeatherPlan.cadence+')' : '\n'+WEATHER_SOURCES[runDataConfig.weatherSource].label;
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
    setStatus('Run the model and wait for monthly results before exporting.');
    return;
  }
  prepareCsvDownload(UI._lastOutput, 'monthly');
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
