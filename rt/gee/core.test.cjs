/* Run with node rt/gee/core.test.cjs from the ACCLIM-RT checkout.
 * reference_brf.json comes from Python prosail 2.0.5 FourSAIL, at 1 nm.
 */
const fs = require('fs');
const path = require('path');
const vm = require('vm');
const assert = require('assert');
const context = vm.createContext({console});
vm.runInContext(fs.readFileSync(path.join(__dirname, 'ACCLIM_RT_GEE.js'), 'utf8').split('// BEGIN ACCLIM-RT RETRIEVAL')[1].split('// END ACCLIM-RT RETRIEVAL')[0], context);
const {RT, rtGroups, rtISODate} = context;
let maximumError = 0;
for (const reference of JSON.parse(fs.readFileSync(path.join(__dirname, 'reference_brf.json')))) {
  const prediction = RT.spectrum(RT.geometry(...reference.angles), reference.lai);
  const actual = [
    prediction[Math.round(reference.soil[0] * 100) - 1],
    prediction[Math.round(reference.soil[1] * 100) + 49]
  ];
  actual.forEach((value, band) => {
    assert(Number.isFinite(value));
    const error = Math.abs(value - reference.brf[band]);
    maximumError = Math.max(maximumError, error);
    assert(error < 1e-10, JSON.stringify({reference, actual}));
  });
}
assert.equal(rtISODate('2023-02-29'), null);
assert.equal(rtISODate('2024-02-29'), Date.UTC(2024, 1, 29));
assert.equal(rtISODate('2012-5-01'), null);
const properties = {red: .19, nir: .19, solar_zenith: 30, view_zenith: 20,
  relative_azimuth: 60, x: 1, y: 2, doy: 122};
const groups = rtGroups([{properties}, {properties: {...properties}},
  {properties: {...properties, x: 2, doy: 121}},
  {properties: {...properties, x: 3, red: -1}}],
  2012, Date.UTC(2012, 4, 1), Date.UTC(2012, 4, 1));
assert.equal(Object.keys(groups).length, 1);
assert.equal(groups['2012-05-01'].length, 1);
const fit = RT.createFit([0, 1, 2, 3].map(x => ({...properties, x})));
while (!fit.step()) {}
const result = fit.finish(0);
assert.equal(result.geometries, 1);
assert.equal(result.rmse, 0);
assert.equal(result.fcover, 0);
assert.equal(result.background_red, .19);
assert.equal(result.background_nir, .19);
assert.equal(result.records.length, 4);
console.log('PASS BRF reference, hotspot, date boundaries, duplicate screening, zero cover, and geometry reuse.');
console.log('Maximum reference BRF error:', maximumError);
const day = 86400000;
const obs = [{millis:0,cover:0},{millis:2*day,cover:.6},{millis:40*day,cover:.8}];
function period(cadence, days) {
  return context.rtCoverPeriods(obs,[{date:'test',millis:0,period_end_millis:days*day,
    period_days:days,model_cadence:cadence}],32).test;
}
assert.equal(period('daily',1).rt_fcover,0);
assert.equal(period('weekly',3).rt_interpolated_days,1);
assert(Math.abs(period('weekly',3).rt_fcover-.3)<1e-12);
assert.equal(period('weekly',7).rt_fcover,undefined);
assert.equal(period('monthly',31).rt_fcover,.3);
assert.equal(RT.median(fit.bestCovers()), result.fcover);

const lineDays=[{millis:0,calibration:context.rtDayMoments([{red:1,nir:1},{red:0,nir:1}],[0,1])},
  {millis:1,calibration:context.rtDayMoments([{red:1,nir:1}],[1])},
  {millis:2,calibration:null}];
const line=context.rtFitCoverLine(lineDays,0,1);
assert(Math.abs(line.m_cover-1/3)<1e-12);
assert(Math.abs(line.b_cover-2/3)<1e-12);
assert.equal(line.cover_fit_pixels,3);
assert.equal(context.rtFitCoverLine(lineDays,1,2),null);
assert.equal(context.rtFitCoverLine(lineDays,0,0).m_cover,1);
