const assert=require('assert/strict'),fs=require('fs'),vm=require('vm');
const core=require('../docs/tower-core.js'),weather=require('../docs/tower-weather.js');
const header='TIMESTAMP_START,TIMESTAMP_END,GPP_DT_VUT_REF,SW_IN_F,TA_F,VPD_F,PA_F,P_F,LE_CORR,H_CORR';
const stamp=t=>new Date(t).toISOString().replace(/[-:T]/g,'').slice(0,12);
const start=Date.UTC(2024,0,1),end=Date.UTC(2025,0,1),lines=[header];
for(let t=start;t<end;t+=1800000){
 const h=new Date(t).getUTCHours(),day=h>=6&&h<18,d=Math.floor((t-start)/86400000);
 lines.push([stamp(t),stamp(t+1800000),day?10+d%5:0,day?200:0,day?20+d%3:10,10,90,.1,day?100:-30,day?50:-20].join(','));
}
const data=core.parse(lines.join('\n'),'FLX_US-Test_FULLSET_HH_2024.csv',{metadata:{latitude:38.4309,longitude:-120.966,site:'Test'}});
assert.equal(data.weather.daily.length,366);assert.equal(data.weather.daily[0].Tmin,10);assert.equal(data.weather.daily[0].Tmax,20);
assert(Math.abs(data.weather.daily[0].precip-4.8)<1e-10);
assert.equal(data.weather.daily[0].Rndaily,50);
const restored=core.restoreSummary(JSON.stringify(data));assert.equal(restored.weather.daily.length,366);
const transfers={};
for(const cadence of ['daily','weekly','monthly']){
 const transfer=core.transferSummary(restored,'2024-01-01','2024-12-31',cadence);transfers[cadence]=transfer;
 assert.equal(transfer.weather.periods.length,{daily:366,weekly:53,monthly:12}[cadence]);
 assert.equal(transfer.weather.periods.reduce((n,r)=>n+r.period_days,0),366);
 assert(Math.abs(transfer.weather.periods.reduce((n,r)=>n+r.precip,0)-366*4.8)<1e-8);
 assert(transfer.weather.periods.every(r=>r.complete));assert(!transfer.weather.daily);
 console.log(cadence+': '+transfer.weather.periods.length+' weather rows, '+JSON.stringify(transfer).length+' transfer characters');
}
assert.equal(transfers.weekly.weather.periods.at(-1).period_days,2);
const partial=weather.periods(data.weather,'weekly','2024-01-30','2024-02-02',30);assert.equal(partial[0].date,'2024-01-29');assert.equal(partial[0].period_days,7);
const missing=JSON.parse(JSON.stringify(data));missing.weather.daily[3].counts.ta=47;missing.weather.daily[3].complete=false;missing.weather.daily[3].Tdaily=null;
assert.equal(weather.periods(missing.weather,'weekly','2024-01-01','2024-01-14',30)[0].complete,false);
assert.equal(weather.periods(missing.weather,'weekly','2024-01-01','2024-01-14',30)[1].complete,true);
missing.weather.daily.splice(4,1);assert.equal(weather.periods(missing.weather,'daily','2024-01-05','2024-01-05',30)[0].complete,false);
assert.throws(()=>weather.periods({...data.weather,daily:undefined},'daily','2024-01-01','2024-01-02',30),/Reprocess/);
for(const file of ['gee/Acclim_v2_shapefile.js','../outputs/ACCLIM_GEE_web_app_script_Private.txt']){
 const s=fs.readFileSync(file,'utf8'),c={};vm.createContext(c);new vm.Script(s);
 vm.runInContext(s.match(/var ACCLIM_PARAMS = \{[\s\S]*?\n\};/)[0]+'\n'+s.slice(s.indexOf('function mergeParams('),s.indexOf('// SECTION 6'))+'\n'+s.slice(s.indexOf('function nativePeriodSeconds('),s.indexOf('function sourceCalendar('))+'\n'+s.slice(s.indexOf('var CARBON_UNITS ='),s.indexOf('/** Annual RAP'))+'\n'+s.slice(s.indexOf('function towerWeatherPlan('),s.indexOf('function setupTowerImport(')),c);
 for(const cadence of ['daily','weekly','monthly']){
  const plan=c.towerWeatherPlan(transfers[cadence],'2024-01-01','2024-12-31');assert.equal(plan.records.length,transfers[cadence].weather.periods.length);assert.equal(plan.missing.length,0);
  assert.equal(c.towerAnnualTemperature(plan,2024).temperatureMonths,12);
 }
 assert.equal(c.towerWeatherPlan(data,'2024-01-01','2024-12-31').records.length,366);
 const fullWeekly=JSON.parse(JSON.stringify(data));fullWeekly.weather.cadence='weekly';assert.equal(c.towerWeatherPlan(fullWeekly,'2024-01-01','2024-12-31').records.length,53);
 const absent=c.towerWeatherPlan(transfers.daily,'2023-12-31','2024-01-02');assert.equal(absent.records.length,3);assert.equal(absent.missing[0],'2023-12-31');
 const forged=JSON.parse(JSON.stringify(transfers.daily));forged.weather.periods[0].counts.ta=47;assert.equal(c.towerWeatherPlan(forged,'2024-01-01','2024-01-02').missing.length,1);
 const rows=[{date:'2024-01-29',model_cadence:'weekly',period_seconds:604800,daylight_seconds:7*43200,period_days:7,daylength:12,PredGPPmean:10,GPPdaily:5.184,NPPdaily:2.592,precip:7}];
 const products=c.buildCarbonProducts(rows,[],'2024-01-29','2024-02-04');
 const npp=c.nativeNPPChartRows(rows,products.daily);assert.equal(npp.length,1);assert(Math.abs(npp[0].NPPdaily-2.592)<1e-10);
 const biomass=c.nativeBiomassChartRows(rows,products.daily);assert.equal(biomass.length,1);
 const monthly=c.aggregateModelMonths(rows);assert.equal(monthly.length,2);assert.equal(monthly[0].precip,3);assert.equal(monthly[1].precip,4);
 const observation=c.towerPeriodObservation(data,rows[0],'2024-01-29','2024-02-04');assert.equal(observation.n,7*24);assert.equal(observation.leN,7*48);
 assert(s.includes('result.monthly=nativeRows;'));assert(s.includes('nativeNPPChartRows(nativeRows,products.daily)'));
 assert(!s.includes('UI._lastOutput=carbonRowsToFC(cache.rows)'));
}
if(process.argv[2])fs.writeFileSync(process.argv[2],JSON.stringify(data));
console.log('PASS daily weather persistence, three cadences, leap/year/week boundaries, missing periods, legacy monthly compatibility, native chart cadence and observation aggregation');
