const fs=require('fs'),vm=require('vm'),assert=require('assert/strict');
for(const file of ['gee/Acclim_v2_shapefile.js','../outputs/ACCLIM_GEE_web_app_script_Private.txt']) {
 const s=fs.readFileSync(file,'utf8'),ctx={};vm.createContext(ctx);
 vm.runInContext(s.match(/var ACCLIM_PARAMS = \{[\s\S]*?\n\};/)[0]+'\n'+s.match(/var SOIL_BUCKET = .*?;/)[0]+'\n'+
 s.slice(s.indexOf('function mergeParams('),s.indexOf('// SECTION 6'))+'\n'+
 'function carbonFinite(x){return typeof x==="number"&&isFinite(x);}\n'+
 s.slice(s.indexOf('function nativePeriodSeconds('),s.indexOf('function sourceCalendar('))+'\n'+
 s.slice(s.indexOf('function nativeBucketPeriod('),s.indexOf('function buildCarbonProducts(')),ctx);
 const rows=Array.from({length:365},(_,i)=>({date:new Date(Date.UTC(2023,0,1)+i*86400000).toISOString().slice(0,10),millis:Date.UTC(2023,0,1)+i*86400000,
   period_days:1,daylength:12,precip:i%8===0?4:0,PT_mean:120,FIPAR:.5,fwet_mean:.1,fSM_mean:.5,
   PredLEsday_mean:33,PredLEday_mean:73,PredLEcday_mean:30,PredLEiday_mean:10,model_valid:1}));
 for(const i of [10,11]){delete rows[i].precip;delete rows[i].PT_mean;rows[i].model_valid=0;}
 const result=ctx.nativeFitBucket(rows);
 // Independent live Earth Engine fit of this same fixture: 17.038796819296145 mm.
 assert(Math.abs(result[0].SMmax-17.038796819296145)<1e-9);
 assert.equal(result[10].bucket_forcing_status,'interpolated');assert.equal(result[10].bucket_fit_observed,0);
 assert(Math.abs(result[0].bucket_record_missing_pct-200/365)<1e-12);
 for(let i=0;i<result.length;i++) {
   const r=result[i];assert(Math.abs(r.SM+r.precip-r.Esoil_period-r.bucket_overflow-r.SM_end)<1e-8);
   if(i)assert.equal(r.SM,result[i-1].SM_end);
 }
 assert.equal(rows[10].precip,undefined); // No mutation of supplied inputs.
 const et=ctx.fitNativeBuckets(rows,{});assert.equal(et.length,365);assert(et.every(r=>Object.values(r).every(v=>v!==undefined)));
 assert.equal(et[10].PredGPPmean,undefined); // No invented carbon.
 const edges=[{date:'2022-12-31',millis:Date.UTC(2022,11,31),period_days:1,daylength:12},...rows];
 const bridged=ctx.nativeBucketBridge(edges);assert.equal(bridged[0].bucket_forcing_status,'unavailable');
 for(const input of [[100,100,0,300],[0,100,80,40],[20,100,300,40],[50,100,0,0]]) {
   const p=ctx.nativeBucketPeriod(...input);let state=input[0],evap=0;
   for(let i=0;i<100;i++){const q=ctx.nativeBucketPeriod(state,input[1],input[2]/100,input[3]/100);state=q.endSM;evap+=q.evaporation;}
   assert(Math.abs(p.endSM-state)<1e-8);assert(Math.abs(p.evaporation-evap)<1e-8);
 }
}
console.log('PASS numerical fit matches live Earth Engine, gap continuity, mass balance, excluded gap targets, no extrapolation, and no invented carbon');
