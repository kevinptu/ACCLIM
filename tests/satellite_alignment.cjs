const fs=require('fs'),vm=require('vm'),assert=require('assert/strict');
for(const file of ['gee/Acclim_v2_shapefile.js','../outputs/ACCLIM_GEE_web_app_script_Private.txt']) {
 const s=fs.readFileSync(file,'utf8'),ctx={};vm.createContext(ctx);
 vm.runInContext(s.slice(s.indexOf('var SATELLITE_SOURCES ='),s.indexOf('function sourceCalendar('))+
 s.slice(s.indexOf('function carbonDateMillis('),s.indexOf('function carbonAllocation('))+
 s.slice(s.indexOf('function alignSatelliteRows('),s.indexOf('function aggregateModelMonths(')),ctx);
 const start=Date.UTC(2024,0,1),obs=[{millis:start,red:.1,nir:.3},{millis:start+8*86400000,red:.2,nir:.4}];
 let rows=ctx.alignSatelliteRows(obs,'2023-12-31','2024-01-10','LANDSAT8');
 assert.equal(rows[0].satellite_status,'unavailable');assert.equal(rows[1].satellite_status,'observed');
 assert.equal(rows[5].satellite_status,'interpolated');assert(Math.abs(rows[5].red-.15)<1e-12);
 assert(Math.abs(rows[5].NDVI-.4)<1e-12);assert.equal(rows[10].satellite_status,'unavailable');
 assert.equal(rows[5].satellite_date_basis,'composite midpoint');
 assert.equal(rows[5].satellite_gap_days,8);assert.equal(rows[5].satellite_missing_pct,100);assert.equal(rows[5].satellite_unavailable_pct,0);
 const gap=[obs[0],{millis:start+33*86400000,red:.2,nir:.4}];
 assert.equal(ctx.alignSatelliteRows(gap,'2024-01-02','2024-01-02','HLS')[0].satellite_status,'unavailable');
 const duplicates=[{...obs[0],satellite_observation_count:2},{millis:start,red:.4,nir:.6,satellite_observation_count:1}];
 const merged=ctx.alignSatelliteRows(duplicates,'2024-01-01','2024-01-01','MODIS16')[0];
 assert(Math.abs(merged.red-.2)<1e-12);assert.equal(merged.satellite_date_basis,'pixel observation date');
 assert.equal(ctx.alignSatelliteRows([],'2024-01-01','2024-01-01','MODIS16')[0].satellite_status,'unavailable');
}
console.log('PASS paired-reflectance interpolation, derived indices, limits, no extrapolation, duplicate-date weighting, provenance and coverage');
