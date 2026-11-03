const fs=require('fs'),vm=require('vm'),assert=require('assert/strict');
class N{constructor(v){this.v=v instanceof N?v.v:v;if(typeof this.v!=='number')throw Error('missing GE');}divide(v){return new N(this.v/new N(v).v);}multiply(v){return new N(this.v*new N(v).v);}max(v){return new N(Math.max(this.v,v));}min(v){return new N(Math.min(this.v,v));}pow(v){return new N(this.v**v);}}
for(const file of ['gee/Acclim_v2_shapefile.js','../outputs/ACCLIM_GEE_web_app_script_Private.txt']){
 const s=fs.readFileSync(file,'utf8');new vm.Script(s);const c={ee:{Number:v=>new N(v)}};vm.createContext(c);
 vm.runInContext(s.slice(s.indexOf('function greenExcessMultiplier('),s.indexOf('function computeDiagnostics(')),c);
 for(const [ge,exp,want] of [[-.1,1,0],[0,1,0],[.05,1,.5],[.1,1,1],[.2,1,1],[.025,.5,.5],[-.1,.5,0],[null,0,1]])assert.equal(c.greenExcessMultiplier(ge,{fge_exp:exp}).v,want);
 // Execute the actual two capacity expressions, checking identical attenuation.
 for(const ge of [0,.05,.1]){
  Object.assign(c,{QYC3_mean:new N(.04),QYC4_mean:new N(.06),FAPAR:new N(.5),PARday_mean:new N(1000),fD:new N(.8),fM:new N(.5),fge:c.greenExcessMultiplier(ge,{fge_exp:1})});
  vm.runInContext(s.match(/var GPPmax_C3 =[^;]+;/)[0]+s.match(/var GPPmax_C4 =[^;]+;/)[0],c);
  assert(Math.abs(c.GPPmax_C3.v-8*ge/.1)<1e-12);assert(Math.abs(c.GPPmax_C4.v-12*ge/.1)<1e-12);
 }
 assert(s.includes("if(p.fge_exp>0)forcingFC=forcingFC.filter(ee.Filter.notNull(['GE']))"));
 assert(s.includes('d.GE=g&&carbonFinite(g.GE)?g.GE:null;'));
}
console.log('PASS raw GE attenuation, clipping, fractional exponent, disabled/missing GE, both GPP capacities and forcing wiring');
