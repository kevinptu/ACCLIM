const fs=require('fs'),vm=require('vm'),assert=require('assert/strict');
for(const file of ['gee/Acclim_v2_shapefile.js','../outputs/ACCLIM_GEE_web_app_script_Private.txt']){
 const s=fs.readFileSync(file,'utf8'),c={carbonFinite:v=>typeof v==='number'&&Number.isFinite(v),carbonDateMillis:s=>Date.parse(s+'T00:00:00Z'),carbonDateString:t=>new Date(t).toISOString().slice(0,10)};
 new vm.Script(s);vm.createContext(c);vm.runInContext(s.slice(s.indexOf('function historicalVIWindow('),s.indexOf('function runForGeometryAsync(')),c);
 const w=c.historicalVIWindow('2024-02-29');assert.equal(w.start,'2014-02-28');assert.equal(w.end,'2024-02-28');
 const row=(date,r,n)=>({millis:c.carbonDateMillis(date),red:r,nir:n});
 const rows=[row('2014-01-01',.2,.2),row('2019-01-01',.2,.4),row('2023-12-31',.2,.6)];
 const fit=c.historicalVIMinima(rows,'2024-01-01',{ndvi_max:1});
 assert.equal(fit.n,3);assert(Math.abs(fit.ndvi-1/30)<1e-12);assert(Math.abs(fit.evi2-(.5/1.88)*.1)<1e-12);
 const filtered=c.historicalVIMinima([...rows,row('2013-12-31',.5,0),row('2024-01-01',.5,0),row('2015-01-01',null,.2),rows[0]],'2024-01-01',{ndvi_max:1});
 assert.equal(filtered.n,3);assert.equal(filtered.ndvi,fit.ndvi);assert.equal(filtered.first,'2014-01-01');assert.equal(filtered.last,'2023-12-31');
 assert.throws(()=>c.historicalVIMinima([],'2024-01-01',{ndvi_max:1}),/at least three/);
 assert.throws(()=>c.historicalVIMinima(rows,'2024-01-01',{ndvi_max:0}),/below NDVImax/);
 assert(s.includes("value: 'FIPAR_P05'"));assert(s.includes("fitVIMinP05: UI.faparSelector.getValue() === 'FIPAR_P05'"));assert(s.includes('params.ndvi_min=fit.ndvi;params.EVI2min=fit.evi2;'));
}
console.log('PASS ten-year boundaries, leap day, independent fifth percentiles, excluded out-of-window/missing data, same-day deduplication, validation and menu wiring');
