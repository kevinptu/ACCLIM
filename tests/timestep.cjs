const fs=require('fs'),vm=require('vm'),assert=require('assert/strict');
for(const file of ['gee/Acclim_v2_shapefile.js','../outputs/ACCLIM_GEE_web_app_script_Private.txt']){
 const s=fs.readFileSync(file,'utf8'),c={};vm.createContext(c);
 vm.runInContext(s.match(/var ACCLIM_PARAMS = \{[\s\S]*?\n\};/)[0]+'\n'+s.match(/var SOIL_BUCKET = .*?;/)[0]+'\n'+
 s.slice(s.indexOf('function mergeParams('),s.indexOf('// SECTION 6'))+'\n'+
 s.slice(s.indexOf('function nativePeriodSeconds('),s.indexOf('function sourceCalendar('))+'\n'+
 s.slice(s.indexOf('var CARBON_UNITS ='),s.indexOf('/** Annual RAP')),c);
 const start=Date.UTC(2024,0,1),end=Date.UTC(2024,0,31),seconds=31*86400;
 const monthly={date:'2024-01-01',millis:start,period_seconds:seconds,daylight_seconds:seconds/2,daylength:12,
   PredGPPmean:10,precip:31,PT_mean:120,FIPAR:.5,fwet_mean:.1,fSM_mean:.5,
   PredLEsday_mean:33,PredLEday_mean:73,PredLEcday_mean:30,PredLEiday_mean:10};
 const daily=Array.from({length:31},(_,i)=>({...monthly,date:new Date(start+i*86400000).toISOString().slice(0,10),millis:start+i*86400000,
   period_seconds:86400,daylight_seconds:43200,precip:1}));
 const m=c.periodCarbonByDay([monthly],start,end),d=c.periodCarbonByDay(daily,start,end);
 const total=x=>Object.values(x).reduce((a,r)=>a+r.gpp,0);
 assert(Math.abs(total(m)-total(d))<1e-10);assert(Math.abs(total(m)-10*seconds/2*12e-6)<1e-10);
 const one=c.nativeBucketTrajectory(c.nativeBucketBridge([monthly]),100,true).history[0];
 const many=c.nativeBucketTrajectory(c.nativeBucketBridge(daily),100,true).history;
 assert(Math.abs(one.SM_end-many[30].SM_end)<1e-9);
 assert(Math.abs(one.Esoil_period-many.reduce((a,r)=>a+r.Esoil_period,0))<1e-9);
 assert(Math.abs(Math.exp(-seconds/(7*86400))-Math.pow(Math.exp(-86400/(7*86400)),31))<1e-12);
 const leap={...monthly,date:'2024-02-01',period_seconds:29*86400,daylight_seconds:29*10*3600};
 const bins=c.periodCarbonByDay([leap],Date.UTC(2024,1,1),Date.UTC(2024,1,29));
 assert.equal(Object.keys(bins).length,29);assert(Math.abs(total(bins)-10*29*10*3600*12e-6)<1e-9);
 // Explicit interval metadata supports arbitrary durations, independent of cadence labels.
 const short=Array.from({length:4},(_,i)=>({...monthly,period_start_millis:start+i*21600000,
   period_seconds:21600,daylight_seconds:10800})).reverse();
 assert(Math.abs(c.periodCarbonByDay(short,start,start)['2024-01-01'].gpp-5.184)<1e-10);
 short.pop();assert.equal(c.periodCarbonByDay(short,start,start)['2024-01-01'].gpp,null);
 assert.throws(()=>c.periodCarbonByDay([monthly,monthly],start,end),/Overlapping/);
 assert(s.includes("var w=elapsedSeconds.multiply(-1).divide(ee.Number(p.tc).multiply(86400)).exp()"));
 assert(s.includes("var initialSeconds=(params.initialization_days===undefined?5:params.initialization_days)*86400"));
 assert(!s.includes("var nativeDaily = monthlyRows"));
}
console.log('PASS one period integration for daily/monthly/subdaily metadata, carbon and bucket subdivision invariance, leap days, gaps, and second-based acclimation');
