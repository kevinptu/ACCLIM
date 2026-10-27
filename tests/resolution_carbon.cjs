const fs=require('fs'),vm=require('vm'),assert=require('assert/strict');
for(const file of ['gee/Acclim_v2_shapefile.js','../outputs/ACCLIM_GEE_web_app_script_Private.txt']) {
 const s=fs.readFileSync(file,'utf8'),ctx={};vm.createContext(ctx);
 vm.runInContext(s.match(/var ACCLIM_PARAMS = \{[\s\S]*?\n\};/)[0]+'\n'+
 s.slice(s.indexOf('function mergeParams('),s.indexOf('// SECTION 6'))+'\n'+
 s.slice(s.indexOf('function nativePeriodSeconds('),s.indexOf('function sourceCalendar('))+'\n'+
 s.slice(s.indexOf('var CARBON_UNITS ='),s.indexOf('/** Annual RAP')),ctx);
 const climate=[{year:2024,MAT:10,temperatureMonths:12,AFG:20,PFG:40,SHR:20,TRE:20}];
 const rows=Array.from({length:29},(_,i)=>({date:'2024-02-'+String(i+1).padStart(2,'0'),
   model_cadence:'daily',period_days:1,PredGPPmean:i%2?20:10,daylength:i%2?12:8,precip:1}));
 let p=ctx.buildCarbonProducts(rows,climate,'2024-02-01','2024-02-29');
 assert.equal(p.daily.length,29);assert(p.daily.every(r=>r.interpolation==='native_daily'));
 const total=rows.reduce((a,r)=>a+r.PredGPPmean*r.daylength*3600*12e-6,0);
 assert(Math.abs(p.daily.reduce((a,r)=>a+r.GPPdaily,0)-total)<1e-10);
 let sum=ctx.carbonSummaryRows(p.daily,false)[0];assert(Math.abs(sum.GPP_gC_m2-total)<1e-10);assert.equal(sum.expected_days,29);
 let month=ctx.aggregateModelMonths(rows)[0];
 const expected=rows.reduce((a,r)=>a+r.PredGPPmean*r.daylength,0)/rows.reduce((a,r)=>a+r.daylength,0);
 assert(Math.abs(month.PredGPPmean-expected)<1e-10);assert.equal(month.precip,29);
 delete rows[10].PredGPPmean;
 p=ctx.buildCarbonProducts(rows,climate,'2024-02-01','2024-02-29');
 assert.equal(p.daily[10].GPPdaily,null);assert.equal(p.daily[11].interpolation,'native_daily');
 sum=ctx.carbonSummaryRows(p.daily,false)[0];assert.equal(sum.GPP_gC_m2,null);
 assert(sum.GPP_gC_m2_available_sum>0);assert(Math.abs(sum.GPPdaily_missing_pct-100/29)<1e-10);
 const partial=ctx.carbonSummaryRows(p.daily.slice(0,10),false)[0];assert.equal(partial.GPP_gC_m2,null);assert.equal(partial.requested_days,10);
 const monthly=ctx.buildCarbonProducts([{date:'2024-02-01',PredGPPmean:10,daylength:12}],climate,'2024-02-01','2024-02-29');
 assert(monthly.daily.every(r=>r.interpolation==='allocated_from_weather_period'));
 assert(Math.abs(monthly.daily[0].GPPdaily-5.184)<1e-10);
}
console.log('PASS native daily carbon, leap-year totals, daylight-weighted monthly means, missing/partial coverage, and retained monthly mode');
