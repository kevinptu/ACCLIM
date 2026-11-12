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
if (typeof ui !== 'undefined') rtApp();
