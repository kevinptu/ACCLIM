const fs=require('fs'),vm=require('vm'),assert=require('assert');
class N{constructor(x){this.v=x instanceof N?x.v:x}add(x){return new N(this.v+new N(x).v)}subtract(x){return new N(this.v-new N(x).v)}multiply(x){return new N(this.v*new N(x).v)}divide(x){return new N(this.v/new N(x).v)}pow(x){return new N(this.v**x)}}
for(const file of ['gee/Acclim_v2_shapefile.js','../outputs/ACCLIM_GEE_web_app_script_Private.txt']){
const s=fs.readFileSync(file,'utf8');new vm.Script(s);const ctx={ee:{Number:x=>new N(x),Algorithms:{IsEqual:(a,b)=>a===b,If:(c,a,b)=>c?a:b}},photosyntheticFAPAR:()=>new N(.7)};vm.createContext(ctx);
vm.runInContext(s.slice(s.indexOf('function griddedDaytimeRn('),s.indexOf('function assembleForcing(')),ctx);
vm.runInContext(s.slice(s.indexOf('function applySharedFPAR('),s.indexOf('function runACCLIM(')),ctx);
for(const source of ['gridded','tower']){
 const data={FIPAR:.4,Tdaily:20,Tday:25,ea:1.5,albedo:.2,SWRdaily:200,daylength:12,Rndaily:100,Rnday:200,weather_source:source};
 const f={get:k=>data[k],getNumber:k=>new N(data[k]),set:o=>o};
 const result=ctx.applySharedFPAR(f,{em_v:.97,em_nv:.92,sb:5.67e-8});assert.equal(result.FIPAR.v,.7);assert.equal(result.FIPAR_NDVI.v,.4);
 if(source==='tower'){assert.equal(result.Rnday,200);assert.equal(result.Rndaily,100);}else{assert(Number.isFinite(result.Rnday.v));assert.notEqual(result.Rnday.v,200);}
}
assert(s.includes('p.useSharedFPAR ? ee.Number(1) : FAPAR.divide'));
assert(s.indexOf('p.faparFit = fitEVI2FAPAR')<s.indexOf('if (p.useSharedFPAR) forcingFC'));
}
console.log('PASS shared fraction, preserved raw FIPAR, gridded radiation update, unchanged tower energy and fitting order');
