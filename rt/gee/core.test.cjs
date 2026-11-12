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

const moments=[{x:.5,y:.5,xx:.5,xy:.5,n:2},{x:0,y:1,xx:0,xy:0,n:1}];
const line=context.rtEqualDateLine(moments);
assert(Math.abs(line.slope-1/3)<1e-12);
assert(Math.abs(line.intercept-2/3)<1e-12);
assert.equal(context.rtEqualDateLine(moments.slice(1)),null);
assert.equal(context.rtDayInterception(1,[1],0),1-Math.exp(-.5));
assert.equal(context.rtDayInterception(1,[.5],0),1-Math.exp(-1));
assert.equal(context.rtDayInterception(0,[.5],.5),0);
assert.equal(context.rtDayInterception(1,null,.5),null);
let skyRef=0;
for(let i=0;i<10000;i++){const mu=(i+.5)/10000;skyRef+=2*mu*(1-Math.exp(-.5/mu))/10000;}
assert(Math.abs(context.rtDayInterception(1,[.5],1)-skyRef)<1e-4);

// Match independently evaluated spectral FourSAIL fluxes for LAI 1.3, mu .55, flat background .2.
const hist = Array(81).fill(0); hist[13]=4;
const response=context.rtAbsorptionResponse(hist,.2), angles=context.RT_PAR.mu;
const j=angles.indexOf(.55)*2;
[response[j],response[j+1],...response.slice(-2)].forEach((v,i)=>{
  assert(Math.abs(v-[.684791895,.558815381,.713196543,.533485125][i])<.000006);
});
const bare=Array(81).fill(0);bare[0]=4;
assert(context.rtAbsorptionResponse(bare,.4).every(v=>v===0));
const first=Date.UTC(2012,4,1);
const observations=[{millis:first,background_red:.2,lai_histogram:bare},
  {millis:first+2*day,background_red:.2,lai_histogram:hist}];
const periods=[0,1,2].map(n=>({date:String(n),millis:first+n*day,period_end_millis:first+(n+1)*day,period_days:1,model_cadence:'daily'}));
const radiation={[first]:50,[first+day]:150,[first+2*day]:250};
const daily=context.rtAbsorptionPeriods(observations,periods,38.43,.5,radiation);
assert.equal(daily['0'].FAPAR_chl,0);
assert.equal(daily['1'].rt_abs_interpolated_days,1);
assert(daily['1'].FAPAR_leaf>daily['1'].FAPAR_chl);
const aggregate={date:'week',millis:first,period_end_millis:first+3*day,period_days:3,model_cadence:'weekly'};
const week=context.rtAbsorptionPeriods(observations,[aggregate],38.43,.5,radiation).week;
assert(Math.abs(week.FAPAR_chl-(150*daily['1'].FAPAR_chl+250*daily['2'].FAPAR_chl)/450)<1e-12);
assert.equal(context.rtAbsorptionPeriods(observations,[aggregate],38.43,1,{}).week.FAPAR_leaf,undefined);
console.log('PASS FourSAIL spectral absorption, bare ground, bracketed interpolation and photon-weighted aggregation.');

// Compare the compressed per-LAI fit with explicit pixel-by-pixel weighted least squares.
for(const index of ['ndvi','evi2']){
  const paired=[];
  const scenes=[[[.1,.2,1],[.1,.4,12],[.2,.2,1]],[[.03,.26,18],[.08,.31,12]]].map((scene,t)=>{
    const pixels=scene.map(([red,nir])=>({red,nir}));
    const covers=scene.map(([,,lai])=>1-Math.exp(-.05*lai));
    for(const [r,n,lai] of scene){
      const h=Array(81).fill(0);h[lai]=1;
      const v=context.rtDayAbsorption(context.rtAbsorptionResponse(h,.19),38,first+t*day);
      paired.push({x:index==='ndvi'?(n-r)/(n+r):2.5*(n-r)/(1+n+2.4*r),
        y:.7*v[1]+.3*v[3],w:1/scene.length});
    }
    return {millis:first+t*day,latitude:38,background_red:.19,index_bins:context.rtIndexBins(pixels,covers)};
  });
  const xbar=paired.reduce((s,p)=>s+p.w*p.x,0)/2,ybar=paired.reduce((s,p)=>s+p.w*p.y,0)/2;
  const m=paired.reduce((s,p)=>s+p.w*(p.x-xbar)*(p.y-ybar),0)/paired.reduce((s,p)=>s+p.w*(p.x-xbar)**2,0);
  const fit=context.rtFitDaylightLine(scenes,first,first+day,.3,index);
  assert(Math.abs(fit.slope-m)<1e-12 && Math.abs(fit.intercept-(ybar-m*xbar))<1e-12);
  assert.equal(fit.pixels,5);assert.equal(fit.dates,2);
}
console.log('PASS date-weighted pixel regressions for chlorophyll absorption against NDVI and EVI2.');

// Reconstruct the winning objective and band errors from the selected LAIs and background.
const spectra=[[.06,.26],[.13,.24],[0,0],[.17,.32]].map(([red,nir],i)=>({
  red,nir,x:i,y:0,solar_zenith:32,view_zenith:14,relative_azimuth:110
}));
for(const mode of ['BRF','NDVI','EVI2']){
  const f=RT.createFit(spectra,mode);while(!f.step()){}
  const res=f.finish(0);
  const I=(r,n)=>mode==='NDVI'?(n-r)/(n+r):2.5*(n-r)/(1+n+2.4*r);
  let bands=0,extra=0;
  for(const p of res.records){
    const s=RT.spectrum(RT.geometry(32,14,110),p.lai);
    const r=s[Math.round(100*res.background_red)-1],n=s[Math.round(100*res.background_nir)+49];
    bands+=(r-p.red)**2+(n-p.nir)**2;
    if(mode==='BRF'||(mode==='NDVI'&&p.red+p.nir===0))continue;
    const h=1e-6;
    const grad=Math.hypot(I(p.red+h,p.nir)-I(p.red-h,p.nir),I(p.red,p.nir+h)-I(p.red,p.nir-h))/(2*h);
    extra+=((I(r,n)-I(p.red,p.nir))/grad)**2;
  }
  assert(Math.abs(Math.sqrt(bands/8)-res.rmse)<1e-10);
  assert(Math.abs(Math.sqrt((bands+extra)/8)-res.objective_score)<1e-9);
  assert.equal(res.objective_score,Math.min(...res.profile.map(r=>r.objective_score)));
}
console.log('PASS scaled index objectives and independent BRF diagnostics.');
