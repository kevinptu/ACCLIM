const fs=require('fs'),vm=require('vm'),assert=require('assert');
const raw=x=>x instanceof N?x.v:x;
class N {constructor(v){this.v=raw(v)}add(x){return new N(this.v+raw(x))}subtract(x){return new N(this.v-raw(x))}multiply(x){return new N(this.v*raw(x))}divide(x){return new N(this.v/raw(x))}pow(x){return new N(this.v**raw(x))}max(x){return new N(Math.max(this.v,raw(x)))}gt(x){return new N(this.v>raw(x))}gte(x){return new N(this.v>=raw(x))}and(x){return new N(this.v&&raw(x))}}
class F{constructor(d){this.d=d}getNumber(k){return new N(this.d[k])}set(k,v){return new F({...this.d,[k]:raw(v)})}}
class FC{constructor(rows){this.rows=rows}filter(){return this}map(fn){return new FC(this.rows.map(fn))}size(){return new N(this.rows.length)}reduceColumns(r,keys){assert.deepEqual(keys,['EVI2','FIPAR_fit_target']);const x=this.rows.map(f=>f.d[keys[0]]),y=this.rows.map(f=>f.d[keys[1]]),mx=x.reduce((a,b)=>a+b)/x.length,my=y.reduce((a,b)=>a+b)/y.length,b=x.reduce((v,z,i)=>v+(z-mx)*(y[i]-my),0)/x.reduce((v,z)=>v+(z-mx)**2,0);return {offset:my-b*mx,scale:b}}}
for(const file of ['gee/Acclim_v2_shapefile.js','../outputs/ACCLIM_GEE_web_app_script_Private.txt']){
const s=fs.readFileSync(file,'utf8');const ctx={ACCLIM_PARAMS:{ndvi_min:.15,ndvi_max:.97,ndvi_exp:1.25},mergeParams:(a,b)=>({...a,...b}),ee:{Number:x=>new N(x),Dictionary:x=>({get:k=>x[k]}),Filter:{notNull:()=>{},gte:()=>{},lte:()=>{}},Reducer:{linearFit:()=>{}},Algorithms:{If:(c,a,b)=>raw(c)?a:b,IsEqual:(a,b)=>raw(a)===raw(b)}}};vm.createContext(ctx);
vm.runInContext(s.slice(s.indexOf('function computeFIPAR('),s.indexOf('/**',s.indexOf('function computeFIPAR('))),ctx);
vm.runInContext(s.slice(s.indexOf('function fitEVI2FAPAR('),s.indexOf('// Preserve the regression-derived')),ctx);
for(const exponent of [1.25,1.5]){const rows=[.1,.3,.6,.8].map(y=>new F({NDVI:.15+.82*y**(1/exponent),EVI2:.12+.55*y}));const fit=ctx.fitEVI2FAPAR(new FC(rows),{ndvi_exp:exponent});assert(Math.abs(fit.min.v-.12)<1e-10);assert(Math.abs(fit.max.v-.67)<1e-10);assert(fit.valid.v);}
}
console.log('PASS FIPAR regression recovers known endpoints with default and overridden NDVI parameters in both scripts');
