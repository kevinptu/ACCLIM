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
  function createFit(pixels) {
    var groups = {}, costs = new Float64Array(soils.length);
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
        var pixel = pixels[index];
        for (var s = 0; s < soils.length; s++) {
          var redIndex = soils[s][0] - 1, nirIndex = soils[s][1] + 49;
          var best = Infinity, selected = 0;
          for (var a = 0; a <= 80; a++) {
            var dr = lut[a][redIndex] - pixel.red, dn = lut[a][nirIndex] - pixel.nir;
            var error = dr * dr + dn * dn;
            if (error < best) { best = error; selected = a; }
          }
          costs[s] += best;
          choices[index][s] = selected;
        }
      });
      next++;
      return next === keys.length;
    }
    function finish(tolerance) {
      var rmse = Array.prototype.map.call(costs, function (c) { return Math.sqrt(c / (2 * pixels.length)); });
      var best = 0;
      for (var s = 1; s < soils.length; s++) if (rmse[s] < rmse[best]) best = s;
      var near = [];
      rmse.forEach(function (v, i) { if (v <= rmse[best] + tolerance) near.push(i); });
      var covers = choices.map(function (a) { return 1 - Math.exp(-.05 * a[best]); });
      var profile = soils.map(function (soil, i) {
        return {background_red: soil[0] / 100, background_nir: soil[1] / 100,
          rmse: rmse[i], fcover: median(choices.map(function (a) { return 1 - Math.exp(-.05 * a[i]); }))};
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
        rmse: rmse[best], background_red: soils[best][0] / 100, background_nir: soils[best][1] / 100,
        fcover_low: Math.min.apply(null, medians), fcover_high: Math.max.apply(null, medians),
        near_backgrounds: near.length, pixels: pixels.length, geometries: keys.length,
        profile: profile, records: records};
    }
    return {step: step, finish: finish, geometries: keys.length};
  }
  return {createFit: createFit, valid: valid, geometry: geometry, spectrum: spectrum, median: median};
}());
