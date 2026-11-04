const fs=require('fs'),vm=require('vm'),assert=require('assert/strict');
class N{constructor(v){this.v=v instanceof N?v.v:v;if(typeof this.v!=='number')throw Error('missing GE');}add(v){return new N(this.v+new N(v).v);}subtract(v){return new N(this.v-new N(v).v);}divide(v){return new N(this.v/new N(v).v);}multiply(v){return new N(this.v*new N(v).v);}max(v){return new N(Math.max(this.v,v));}min(v){return new N(Math.min(this.v,v));}pow(v){return new N(this.v**v);}}
for(const file of ['gee/Acclim_v2_shapefile.js','../outputs/ACCLIM_GEE_web_app_script_Private.txt']){
 const s=fs.readFileSync(file,'utf8');const scale=Number(s.match(/return ee.Number\(ge\).divide\(([\d.]+)\)/)[1]);new vm.Script(s);const c={ee:{Number:v=>new N(v)}};vm.createContext(c);
 vm.runInContext(s.slice(s.indexOf('function greenExcessMultiplier('),s.indexOf('function computeDiagnostics(')),c);
 for(const [ge,exp,want] of [[-.1,1,0],[0,1,0],[.05,1,.5],[.1,1,1],[.2,1,1],[.025,.5,.5],[-.1,.5,0],[null,0,1]])assert(Math.abs(c.greenExcessMultiplier(ge===null?null:ge*scale/.1,{fge_exp:exp}).v-want)<1e-12);
 vm.runInContext(s.slice(s.indexOf('function photosyntheticFAPAR('),s.indexOf('function greenExcessMultiplier(')),c);
 for(const [ge,fexp,expected] of [[.05,1,.5],[.05,2,.25],[.05,0,1],[.2,1,1]]){
  assert(Math.abs(c.photosyntheticFAPAR({get:k=>ge*scale/.1},{faparMethod:'FGE',fge_exp:1,fexp}).v-expected)<1e-12);
 }
 for(const shared of [false,true]){
  for(const ndvi of [0,.5,1]){
   const data={NDVI:ndvi,FIPAR:shared?.99:.6,FIPAR_NDVI:.6,GE:scale*.2};
   const result=c.photosyntheticFAPAR({get:k=>data[k]},{faparMethod:'NDVI_BLEND',useSharedFPAR:shared,fge_exp:1,fexp:2}).v;
   assert(Math.abs(result-(ndvi*.6+(1-ndvi)*.2)**2)<1e-12);
  }
 }
 for(const key of ['GPPmax_C3','GPPmax_C4'])assert(!s.match(new RegExp('var '+key+' =[^;]+;'))[0].includes('fge'));
 assert(s.includes("if((p.faparMethod==='FGE' || p.faparMethod==='NDVI_BLEND') && p.fge_exp>0)forcingFC=forcingFC.filter(ee.Filter.notNull(['GE']))"));
 assert(s.includes('d.GE=g&&carbonFinite(g.GE)?g.GE:null;'));
}
console.log('PASS raw GE attenuation, clipping, fractional exponent, disabled/missing GE, both GPP capacities and forcing wiring');
