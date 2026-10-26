const fs=require('fs'),vm=require('vm'),assert=require('assert');
class N{constructor(v){this.v=v instanceof N?v.v:v}subtract(x){return new N(this.v-new N(x).v)}gt(x){return new N(this.v>new N(x).v)}and(x){return new N(Boolean(this.v)&&Boolean(new N(x).v))}}
for(const path of ['gee/Acclim_v2_shapefile.js','../outputs/ACCLIM_GEE_web_app_script_Private.txt']){
 const s=fs.readFileSync(path,'utf8');new vm.Script(s);
 const ctx={ee:{Number:x=>new N(x)}};vm.createContext(ctx);
 vm.runInContext(s.slice(s.indexOf('function fixedUpperEVI2Fit('),s.indexOf('function runACCLIM(')),ctx);
 for(const [min,valid,expected] of [[.15,true,true],[.6,true,false],[.7,true,false],[.15,false,false]]){
 const original={min:new N(min),max:new N(.9),valid,a:1,b:2,count:12};
 const result=ctx.fixedUpperEVI2Fit(original);assert.equal(result.min.v,min);assert.equal(result.max.v,.6);assert.equal(result.valid.v,expected);assert.equal(original.max.v,.9);
 }
 assert(s.includes("value: 'EVI2_FITTED_MIN'"));assert(s.includes("fitEVI2: UI.faparSelector.getValue() !== 'EVI2_FIXED'"));
}
console.log('PASS fitted lower endpoint, fixed upper endpoint, invalid ranges, preserved fit, UI wiring, both script syntaxes');
