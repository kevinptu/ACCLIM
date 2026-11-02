const fs=require('fs'),vm=require('vm'),assert=require('assert/strict');
for(const file of ['gee/Acclim_v2_shapefile.js','../outputs/ACCLIM_GEE_web_app_script_Private.txt']){
 const s=fs.readFileSync(file,'utf8'),c={};vm.createContext(c);new vm.Script(s);
 vm.runInContext(s.match(/var SATELLITE_SOURCES = \{[\s\S]*?\n\};/)[0]+'\n'+s.slice(s.indexOf('function nativePeriodSeconds('),s.indexOf('function sourceCalendar('))+'\n'+s.slice(s.indexOf('var CARBON_UNITS ='),s.indexOf('/** Annual RAP'))+'\n'+s.slice(s.indexOf('function geSourceSpec('),s.indexOf('function addWeatherCharts(')),c);
 assert.equal(c.geSourceSpec('MODIS1').ids[0],'MODIS/061/MOD09GA');
 assert.equal(c.geSourceSpec('MODIS16').ids[0],'MODIS/061/MOD09A1');assert.equal(c.geSourceSpec('MODIS8').scale,500);
 const t=Date.UTC(2024,0,1),obs=[{millis:t,GE:-.1,EVI2_GE:.1,NDVI_GE:.2},{millis:t+2*86400000,GE:.1,EVI2_GE:.3,NDVI_GE:.6}];
 const rows=c.aggregateGEDiagnostics(obs,c.selectedTimestepPeriods('2024-01-01','2024-01-03','daily'),32);
 assert.equal(rows[1].GE,0);assert(Math.abs(rows[1].EVI2_GE-.2)<1e-10);assert(Math.abs(rows[1].NDVI_GE-.4)<1e-10);
 assert.equal(rows[1].GE_observed_pct,0);assert.equal(rows[1].GE_available_pct,100);
 const f=c.geLinearFit(rows,'NDVI_GE');assert(Math.abs(f.intercept-.4)<1e-10);assert(Math.abs(f.slope-2)<1e-10);assert(Math.abs(f.r2-1)<1e-10);assert.equal(f.extrapolated,false);
 const week=c.aggregateGEDiagnostics(obs,c.selectedTimestepPeriods('2024-01-01','2024-01-03','weekly'),32)[0];assert.equal(week.GE,0);assert(Math.abs(week.NDVI_GE-.4)<1e-10);
 // A missing index rejects the entire triplet; no independent support or fabricated zeros.
 const missing=c.aggregateGEDiagnostics([obs[0],{...obs[1],NDVI_GE:null}],c.selectedTimestepPeriods('2024-01-01','2024-01-03','daily'),32);
 assert.equal(missing[1].GE,null);assert.equal(missing[1].EVI2_GE,null);assert.equal(missing[1].NDVI_GE,null);
 assert.equal(c.geLinearFit([{GE:1,NDVI_GE:.2},{GE:1,NDVI_GE:.3},{GE:1,NDVI_GE:.4}],'NDVI_GE'),null);
 const negative=rows.map(r=>({...r,GE:r.GE-1}));assert.equal(c.geLinearFit(negative,'NDVI_GE').extrapolated,true);
 const fit=c.ndviMinGEFit(rows,{ndvi_max:.95});assert(Math.abs(fit.intercept-.4)<1e-10);
 assert.throws(()=>c.ndviMinGEFit(rows.slice(0,2),{ndvi_max:.95}),/at least three/);
 assert.throws(()=>c.ndviMinGEFit(rows,{ndvi_max:.3}),/not below NDVImax/);
 assert.throws(()=>c.ndviMinGEFit(rows.map(r=>({...r,NDVI_GE:r.NDVI_GE-2})),{ndvi_max:.95}),/outside/);

}
console.log('PASS MODIS source mapping, common-support interpolation/aggregation, zero/negative GE, regression intercepts, missing data and degenerate fits');

