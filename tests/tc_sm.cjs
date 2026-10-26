const fs=require('fs'),vm=require('vm'),assert=require('assert');
class N{constructor(x){this.v=x instanceof N?x.v:x}multiply(x){return new N(this.v*new N(x).v)}subtract(x){return new N(this.v-new N(x).v)}add(x){return new N(this.v+new N(x).v)}}
for(const file of ['gee/Acclim_v2_shapefile.js','../outputs/ACCLIM_GEE_web_app_script_Private.txt']){
const s=fs.readFileSync(file,'utf8');new vm.Script(s);const c={ee:{Number:x=>new N(x),Algorithms:{If:(x,a,b)=>x?a:b}},state:{get:()=>.8},w:new N(Math.exp(-30/7)),isFirst:false,lowLight:false};vm.createContext(c);
vm.runInContext(s.slice(s.indexOf('  function updateMean('),s.indexOf('  // R script multiplies Pdaily')),c);
assert(Math.abs(c.updateMean('x',.2).v-(.8*Math.exp(-30/7)+.2*(1-Math.exp(-30/7))))<1e-12);
assert(Math.abs(c.updateMean('x',.2,new N(Math.exp(-1))).v-(.8/Math.E+.2*(1-1/Math.E)))<1e-12);
c.lowLight=true;assert.equal(c.updateMean('x',.2,new N(.5)).v,.8);c.lowLight=false;c.isFirst=true;assert.equal(c.updateMean('x',.2,new N(.5)).v,.8);
assert(s.includes('var fwet_mean = RHmin_mean.pow(10)'));assert(s.includes("running.get('RHmin_sm_mean')"));assert(s.includes(".pow(ee.Number(running.get('VPDmax_sm_mean')))"));
assert.equal((s.match(/RHmin_sm_mean: RHmin_sm_mean/g)||[]).length,2);
}
console.log('PASS separate 7/30-day weights, unchanged first-step/low-light handling, fwet isolation and persistence');
