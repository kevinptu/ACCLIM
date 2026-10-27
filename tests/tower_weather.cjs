const assert=require('assert/strict'),fs=require('fs'),vm=require('vm');
const core=require('../docs/tower-core.js');
const header='TIMESTAMP_START,TIMESTAMP_END,GPP_NT_VUT_REF,SW_IN_F,TA_F,VPD_F,PA_F,P_F,NETRAD,LE_CORR,H_CORR,PA,P';
const stamp=ms=>new Date(ms).toISOString().replace(/[-:T]/g,'').slice(0,12);
let rows=[];const first=Date.UTC(2024,0,1),last=Date.UTC(2024,1,1);
for(let t=first;t<last;t+=1800000){let hour=new Date(t).getUTCHours(),day=hour>=6&&hour<18;rows.push([stamp(t),stamp(t+1800000),day?10:0,day?200:0,day?20:10,10,90,0.1,9999,day?100:-30,day?50:-20,80,0].join(','));}
const text=header+'\n'+rows.join('\n'),name='FLX_US-Test_FULLSET_HH_2024.csv';
const d=core.parse(text,name),m=d.weather.months[0];d.weather.cadence='monthly';
assert.equal(m.complete,true);assert.equal(d.weather.source.pa,'PA_F');assert.equal(d.weather.source.rain,'P_F');assert.equal(d.weather.version,2);assert.equal(m.Tmin,10);assert.equal(m.Tmax,20);assert.equal(m.Tdaily,15);assert.equal(m.SWRdaily,100);assert.equal(m.Rndaily,50);assert.equal(m.Pdaily,90);assert.ok(Math.abs(m.precip-148.8)<1e-8);
const es=t=>0.61121*Math.exp(17.502*t/(t+240.97));assert.ok(Math.abs(m.ea-((es(20)+es(10))/2-1))<1e-12);assert.equal(d.months[0].gpp,10);
let fluxMissing=rows.slice();let cells=fluxMissing[0].split(',');cells[9]='-9999';fluxMissing[0]=cells.join(',');const fluxBad=core.parse(header+'\n'+fluxMissing.join('\n'),name);fluxBad.weather.cadence='monthly';assert.equal(fluxBad.weather.months[0].counts.rn,1487);assert.equal(fluxBad.weather.months[0].Rndaily,null);
let missing=rows.slice();missing[0]=missing[0].replace(',90,',',-9999,');let bad=core.parse(header+'\n'+missing.join('\n'),name);bad.weather.cadence='monthly';assert.equal(bad.weather.months[0].complete,false);assert.equal(bad.weather.months[0].Pdaily,null);assert.equal(bad.months[0].gpp,10);
assert.equal(core.parse(header+'\n'+rows.slice(1).join('\n'),name).weather.months[0].complete,false);
const copied=core.transferSummary(d,'2024-01-10','2024-01-15');assert.equal(copied.weather.months[0].expected,1488);assert.equal(copied.daily.length,6);assert.equal(core.restoreSummary(JSON.stringify(d)).weather.months[0].complete,true);
const wt=require('../docs/tower-weather.js'),qhead=['AIR_TEMPERATURE','RELATIVE_HUMIDITY','AIR_PRESSURE','RAINFALL','NET_RADIATION','INCOMING_SHORTWAVE_RADIATION'];
let q=wt.create(qhead,null,qhead[5],core.num,true);q.checkUnits(['deg C','%','kPa','mm','W/m2','W/m2']);assert.throws(()=>q.checkUnits(['K','%','kPa','mm','W/m2','W/m2']),/units/);
for(let t=first;t<last;t+=1800000)q.add(new Date(t).toISOString().slice(0,10),['20','50','90','0','100','200']);assert.ok(Math.abs(q.finish(1800000).months[0].ea-es(20)*0.5)<1e-12);
for(const file of ['gee/Acclim_v2_shapefile.js',...(process.env.ACCLIM_PRIVATE_SCRIPT?[process.env.ACCLIM_PRIVATE_SCRIPT]:[])]){
 const s=fs.readFileSync(file.startsWith('gee/')?__dirname+'/../'+file:file,'utf8');new vm.Script(s);const ctx={};vm.createContext(ctx);
 vm.runInContext('function carbonFinite(x){return typeof x==="number"&&isFinite(x);}\n'+s.slice(s.indexOf('function towerWeatherPlan('),s.indexOf('function setupTowerImport(')),ctx);
 vm.runInContext(s.slice(s.indexOf('function parseTowerMonthly('),s.indexOf('// Validate imported monthly weather')),ctx);
 const leMonth=ctx.towerMonthsInRange(ctx.parseTowerMonthly(JSON.stringify(d)),'2024-01-01','2024-01-31')['2024-01'];
 assert.ok(Math.abs(leMonth.leET-35*86400/2450000)<1e-12);assert.equal(leMonth.leMissingPct,0);
 const lePartial=ctx.towerMonthsInRange(ctx.parseTowerMonthly(JSON.stringify(fluxBad)),'2024-01-01','2024-01-31')['2024-01'];
 assert.equal(lePartial.leN,1487);assert.ok(lePartial.leET!==null);assert.ok(Math.abs(lePartial.leMissingPct-100/1488)<1e-12);
 const leCopy=ctx.parseTowerMonthly(JSON.stringify(core.transferSummary(d,'2024-01-01','2024-01-31')));assert.equal(leCopy.leVariable,'LE_CORR');assert.equal(leCopy.daily[0].leN,48);
 let plan=ctx.towerWeatherPlan(d,'2024-01-01','2024-03-31');assert.equal(plan.segments.length,1);assert.equal(plan.missing.join(','),'2024-02,2024-03');assert.equal(ctx.towerAnnualTemperature(plan,2024).MAT,null);
 let separated={...d,weather:{...d.weather,months:[m,{...m,month:'2024-03'}]}};plan=ctx.towerWeatherPlan(separated,'2024-01-01','2024-03-31');assert.equal(plan.segments.length,2);assert.equal(plan.missing[0],'2024-02');
 plan=ctx.towerWeatherPlan(bad,'2024-01-01','2024-01-31');assert.equal(plan.segments.length,0);assert.ok(plan.rejected[0].includes('pressure: 1487/1488'));
 assert.throws(()=>ctx.towerWeatherPlan({...d,weather:{...d.weather,version:1}},'2024-01-01','2024-01-31'),/Reprocess/);
 assert.throws(()=>ctx.towerWeatherPlan({version:3},'2024-01-01','2024-02-01'),/Reprocess/);

 // Incomplete monthly weather must remain missing in allocated daily carbon.
 const params=s.match(/var ACCLIM_PARAMS = \{[\s\S]*?\n\};/)[0];
 const merge=s.slice(s.indexOf('function mergeParams('),s.indexOf('// SECTION 6'));
 const carbon=s.slice(s.indexOf('var CARBON_UNITS ='),s.indexOf('// BEGIN CARBON ASYNC UI'));
 const cc={};vm.createContext(cc);vm.runInContext(params+'\n'+merge+'\n'+s.slice(s.indexOf('function nativePeriodSeconds('),s.indexOf('function sourceCalendar('))+'\n'+carbon,cc);
 const products=cc.buildCarbonProducts([{date:'2024-01-01',PredGPPmean:10,daylength:12,weather_source:'tower'}, {date:'2024-02-01',weather_source:'tower',weather_status:'missing_weather'}, {date:'2024-03-01',PredGPPmean:20,daylength:12,weather_source:'tower'}],[],'2024-01-01','2024-03-31');
 assert.ok(products.daily.filter(r=>r.date.startsWith('2024-02')).every(r=>r.GPPdaily===null));
 assert.equal(products.daily[0].weather_source,'tower');
 const widgets=[];const ui={Label:o=>o,Textbox:o=>o,Button:o=>o,Checkbox:o=>o};
 const uc={ui,UI:{inputPanel:{add:w=>widgets.push(w)}}};vm.createContext(uc);
 vm.runInContext(s.slice(s.indexOf('function setupTowerImport('),s.indexOf('function addTowerControls(')),uc);uc.setupTowerImport();
 assert.equal(widgets.find(w=>w.label==='Use Tower Weather Data').value,false);

 // A false complete flag or forged count never allows an incomplete month.
 assert.equal(ctx.towerWeatherPlan({...bad,weather:{...bad.weather,months:[{...bad.weather.months[0],complete:true}]}},'2024-01-01','2024-01-31').segments.length,0);
}
console.log('Weather units, full-day means, daily extrema, rainfall accumulation, missing intervals, JSON roundtrip and GEE gap segmentation passed.');
