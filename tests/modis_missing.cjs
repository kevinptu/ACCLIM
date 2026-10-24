const fs=require('fs'),vm=require('vm'),assert=require('assert/strict');
class N {constructor(x){if(x==null)throw Error('null arithmetic');this.x=+x;} divide(v){return new N(this.x/+v);} multiply(v){return new N(this.x*+v);} add(v){return new N(this.x+(+v));} subtract(v){return new N(this.x-(+v));} max(v){return new N(Math.max(this.x,+v));} min(v){return new N(Math.min(this.x,+v));} valueOf(){return this.x;}}
const feature=p=>({p,toDictionary(){return {get:k=>p[k]};},set(q){return feature({...p,...q});}});
class C {constructor(rows){this.rows=rows;} filterDate(){return this;} filterBounds(){return this;} map(fn){return new C(this.rows.map(fn));} filter(fn){return new C(this.rows.filter(fn));}}
const raw=[{NDVI:null,sur_refl_b01:1000,sur_refl_b02:3000},{NDVI:5000,sur_refl_b01:null,sur_refl_b02:3000},{NDVI:5000,sur_refl_b01:1000,sur_refl_b02:null},{NDVI:0,sur_refl_b01:1000,sur_refl_b02:3000}];
for(const path of ['ACCLIM/gee/Acclim_v2_shapefile.js','outputs/ACCLIM_GEE_web_app_script_Private.txt']){
 const source=fs.readFileSync(path,'utf8');new vm.Script(source);
 let rows=raw;
 const ee={Number:x=>new N(x),Date:x=>({advance(){return this;},format(){return x;},millis(){return 1;}}),Reducer:{mean:()=>0},Filter:{notNull:keys=>f=>keys.every(k=>f.p[k]!=null)},Feature:(g,p)=>p?feature(p):g,ImageCollection:()=>new C(rows.map(p=>({get:()=> '2020-03-21',select(){return this;},reduceRegion:()=>p})))};
 const ctx={ee,DATASETS:{MODIS_VI:'test'}};vm.createContext(ctx);vm.runInContext(source.slice(source.indexOf('function fetchModisVI('),source.indexOf('// SECTION 4')),ctx);
 let result=ctx.fetchModisVI({},'2020-01-01','2020-12-31',250);assert.equal(result.rows.length,1);assert.equal(+result.rows[0].p.NDVI,0);assert.equal(+result.rows[0].p.red,0.1);assert.ok(Math.abs(+result.rows[0].p.EVI2-0.5/1.54)<1e-12);
 rows=raw.slice(0,3);assert.equal(ctx.fetchModisVI({},'2020-01-01','2020-12-31',250).rows.length,0);
}
console.log('MODIS missing-band, valid-zero and all-masked adapter regressions passed for both scripts.');
