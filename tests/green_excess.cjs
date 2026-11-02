const fs=require('fs'),vm=require('vm'),assert=require('assert/strict');
for(const file of ['gee/Acclim_v2_shapefile.js','../outputs/ACCLIM_GEE_web_app_script_Private.txt']){
 const s=fs.readFileSync(file,'utf8'),c={};vm.createContext(c);new vm.Script(s);
 vm.runInContext(s.match(/var SATELLITE_SOURCES = \{[\s\S]*?\n\};/)[0]+'\n'+s.slice(s.indexOf('var CARBON_UNITS ='),s.indexOf('/** Annual RAP')),c);
 const t=Date.UTC(2024,0,1),obs=[{millis:t,red:.1,nir:.3,GE:-.1},{millis:t+2*86400000,red:.2,nir:.4,GE:.1}];
 const rows=c.alignSatelliteRows(obs,'2024-01-01','2024-01-03','HLS');
 assert.equal(rows[0].GE,-.1);assert.equal(rows[1].GE,0);assert.equal(rows[2].GE,.1);assert(rows.every(r=>r.GE_available_pct===100));
 const noGE=obs.map(({GE,...r})=>r),old=c.alignSatelliteRows(noGE,'2024-01-01','2024-01-03','HLS');
 rows.forEach((r,i)=>{assert.equal(r.NDVI,old[i].NDVI);assert.equal(r.EVI2,old[i].EVI2);assert.equal(old[i].GE_available_pct,0);assert.equal(old[i].GE,undefined);});
 const partial=c.alignSatelliteRows([obs[0],noGE[1]],'2024-01-01','2024-01-03','HLS');assert.equal(partial[1].GE,undefined);assert(Number.isFinite(partial[1].EVI2));
 const gap=c.alignSatelliteRows([obs[0],{...obs[1],millis:t+40*86400000}],'2024-01-02','2024-01-02','HLS');assert.equal(gap[0].GE,undefined);
}
console.log('PASS negative/zero GE, bounded interpolation, missing RGB without changing vegetation indices, and both script syntaxes');
