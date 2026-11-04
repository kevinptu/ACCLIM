const fs=require('fs'),vm=require('vm'),assert=require('assert/strict');
class N{constructor(v){this.v=v instanceof N?v.v:v;if(typeof this.v!=='number')throw Error('missing GE');}add(v){return new N(this.v+new N(v).v);}subtract(v){return new N(this.v-new N(v).v);}divide(v){return new N(this.v/new N(v).v);}multiply(v){return new N(this.v*new N(v).v);}max(v){return new N(Math.max(this.v,v));}min(v){return new N(Math.min(this.v,new N(v).v));}pow(v){return new N(this.v**v);}}
for(const file of ['gee/Acclim_v2_shapefile.js','../outputs/ACCLIM_GEE_web_app_script_Private.txt']){
 const s=fs.readFileSync(file,'utf8');const scale=Number(s.match(/return ee.Number\(ge\).divide\(([\d.]+)\)/)[1]);new vm.Script(s);const c={ee:{Number:v=>new N(v)}};vm.createContext(c);
 vm.runInContext(s.slice(s.indexOf('function ndviBlendWeight('),s.indexOf('function computeDiagnostics(')),c);
 for(const [ge,exp,want] of [[-.1,1,0],[0,1,0],[.05,1,.5],[.1,1,1],[.2,1,1],[.025,.5,.5],[-.1,.5,0],[null,0,1]])assert(Math.abs(c.greenExcessMultiplier(ge===null?null:ge*scale/.1,{fge_exp:exp}).v-want)<1e-12);
 vm.runInContext(s.slice(s.indexOf('function photosyntheticFAPAR('),s.indexOf('function greenExcessMultiplier(')),c);
 for(const [ge,fexp,expected] of [[.05,1,.5],[.05,2,.25],[.05,0,1],[.2,1,1]]){
  assert(Math.abs(c.photosyntheticFAPAR({get:k=>ge*scale/.1},{faparMethod:'FGE',fge_exp:1,fexp}).v-expected)<1e-12);
 }
 for(const [gePrime,ndvi,fexp,want] of [[.2,.6,1,.2],[.7,.3,1,.3],[.5,.5,2,.25],[.2,-.1,.5,0],[0,.6,1,0]]){
  assert(Math.abs(c.photosyntheticFAPAR({get:k=>k==='GE_normalized'?gePrime:ndvi},{faparMethod:'GE_NDVI_MIN',fexp}).v-want)<1e-12);
 }
 for(const shared of [false,true]){
  for(const weight of [0,1,1.6,2]) for(const exponent of [1,1.6,2]) for(const ndvi of [-.1,0,.5,1]){
   const data={NDVI:ndvi,FIPAR:shared?.99:.6,FIPAR_NDVI:.6,GE:scale*.2};
   const result=c.photosyntheticFAPAR({get:k=>data[k]},{faparMethod:'NDVI_BLEND',useSharedFPAR:shared,fge_exp:1,fexp:2,ndvi_exp:exponent,blend_exp:weight}).v;
   assert(Math.abs(result-(Math.max(0,ndvi)**weight*Math.max(0,ndvi)**exponent+(1-Math.max(0,ndvi)**weight)*.2)**2)<1e-12);
  }
 }
 vm.runInContext(s.slice(s.indexOf('function computeFIPAR('),s.indexOf('/**',s.indexOf('function computeFIPAR('))),c);
 for(const weight of [0,1,1.6,2]) for(const exponent of [1,1.6,2]) for(const n of [-.1,0,.5,1]){
  const expected=Math.max(0,n)**weight*Math.max(0,n)**exponent+(1-Math.max(0,n)**weight)*.2;
  assert(Math.abs(c.computeFIPAR(n,{fiparMethod:'NDVI_BLEND',fge_exp:1,ndvi_exp:exponent,blend_exp:weight},scale*.2).v-expected)<1e-12);
 }
 assert.throws(()=>c.ndviBlendWeight(.5,{blend_exp:-1}),/nonnegative/);
 for(const key of ['GPPmax_C3','GPPmax_C4'])assert(!s.match(new RegExp('var '+key+' =[^;]+;'))[0].includes('fge'));
 assert(s.includes("if((p.faparMethod==='FGE' || p.faparMethod==='NDVI_BLEND' || p.fiparMethod==='NDVI_BLEND') && p.fge_exp>0)forcingFC=forcingFC.filter(ee.Filter.notNull(['GE']))"));
 assert(s.includes('d.GE=g&&carbonFinite(g.GE)?g.GE:null;'));
}
console.log('PASS raw GE attenuation, clipping, fractional exponent, disabled/missing GE, both GPP capacities and forcing wiring');
