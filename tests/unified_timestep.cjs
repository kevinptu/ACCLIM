const assert=require('assert/strict'),fs=require('fs'),vm=require('vm');
const core=require('../docs/tower-core.js');
for(const file of ['gee/Acclim_v2_shapefile.js','../outputs/ACCLIM_GEE_web_app_script_Private.txt']){
 const s=fs.readFileSync(file,'utf8'),c={};vm.createContext(c);new vm.Script(s);
 vm.runInContext(s.slice(s.indexOf('function mergeParams('),s.indexOf('// SECTION 6'))+'\n'+s.slice(s.indexOf('function nativePeriodSeconds('),s.indexOf('function sourceCalendar('))+'\n'+s.slice(s.indexOf('var CARBON_UNITS ='),s.indexOf('/** Annual RAP')),c);
 const periods=c.selectedTimestepPeriods('2024-01-30','2024-03-02','weekly');
 assert.equal(periods[0].date,'2024-01-30');assert.equal(periods[0].period_days,6);assert.equal(periods.at(-1).period_end_exclusive,'2024-03-03');
 for(let i=1;i<periods.length;i++)assert.equal(periods[i].millis,periods[i-1].period_end_millis);
 assert.equal(periods.reduce((n,p)=>n+p.period_days,0),33);
 const dates=c.selectedTimestepPeriods('2024-01-01','2024-12-31','daily');
 const fields={Tmin:10,Tmax:30,Tdaily:20,SWRdaily:200,ea:1,Pdaily:95,precip:2,Rndaily:100};
 const daily=dates.map(p=>({...p,...fields,complete:true,expected:48,counts:{ta:48,sw:48,ea:48,pa:48,rain:48,rn:48}}));
 for(const cadence of ['daily','weekly','monthly']){
  const ps=c.selectedTimestepPeriods('2024-01-01','2024-12-31',cadence),r=c.aggregateWeatherPeriods(daily,ps,'daily',30);
  assert.equal(r.length,{daily:366,weekly:53,monthly:12}[cadence]);assert.equal(r.reduce((n,x)=>n+x.precip,0),732);
  assert(r.every(x=>Math.abs(x.Tdaily-20)<1e-10&&x.complete&&x.model_cadence===cadence));
  const tower={intervalMinutes:30,weather:{version:2,daily,months:[]}};
  const plan=c.selectedTowerWeatherPlan(tower,'2024-01-01','2024-12-31',cadence);assert.equal(plan.records.length,r.length);assert.equal(plan.missing.length,0);
 }
 const gap=daily.map(r=>({...r}));gap[1].Tdaily=null;gap[1].complete=false;
 const g=c.aggregateWeatherPeriods(gap,c.selectedTimestepPeriods('2024-01-01','2024-01-14','weekly'),'daily',30);
 assert.equal(g[0].Tdaily,null);assert(!g[0].complete);assert(g[1].complete);assert(Math.abs(g[0].weather_missing_pct-100/7)<1e-8);
 // Coarse source allocation conserves totals, also across week/month boundaries.
 const coarse=[{date:'2024-01-01',period_days:31,...fields,precip:62},{date:'2024-02-01',period_days:29,...fields,Tdaily:22,precip:87}];
 const fine=c.aggregateWeatherPeriods(coarse,c.selectedTimestepPeriods('2024-01-01','2024-02-29','daily'),'monthly');
 assert.equal(fine.reduce((n,r)=>n+r.precip,0),149);assert(fine.every(r=>r.weather_temporally_allocated_pct===100));
 const week=c.aggregateWeatherPeriods(coarse,c.selectedTimestepPeriods('2024-01-29','2024-02-04','weekly'),'monthly');
 assert.equal(week[0].precip,18);assert(Math.abs(week[0].Tdaily-(20*3+22*4)/7)<1e-8);
 assert.throws(()=>c.aggregateWeatherPeriods([coarse[0],coarse[0]],dates,'monthly'),/overlapping/);
 const tower={intervalMinutes:30,daily:[{date:'2024-01-01',sum:24,n:24,secondSum:48,secondN:24,daylight:24,total:48,unknownLight:0,absent:0,leSum:4800,leN:48},{date:'2024-01-02',sum:30,n:10,secondSum:60,secondN:10,daylight:24,total:48,unknownLight:0,absent:0,leSum:2400,leN:24}]};
 const obs=c.towerPeriodObservation(tower,{date:'2024-01-01',period_days:2},'2024-01-01','2024-01-02');assert.equal(obs.gpp,54/34);assert.equal(obs.secondGpp,108/34);assert.equal(obs.leMissingPct,25);
}
assert(!fs.readFileSync('docs/tower-csv-importer.html','utf8').includes('id="weatherCadence"'));
assert(!fs.readFileSync('docs/tower-ui.js','utf8').includes("el('weatherCadence')"));
console.log('PASS unified intervals, leap years, clipped weeks, all weather cadences, coarse-source conservation/flags, missingness and count-weighted tower fluxes');

