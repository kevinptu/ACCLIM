const fs=require('fs'),vm=require('vm'),assert=require('assert/strict');
class N{constructor(v){this.v=v instanceof N?v.v:v;}multiply(v){return new N(this.v*new N(v).v);}add(v){return new N(this.v+new N(v).v);}subtract(v){return new N(this.v-new N(v).v);}divide(v){return new N(this.v/new N(v).v);}max(v){return new N(Math.max(this.v,v));}min(v){return new N(Math.min(this.v,v));}pow(v){return new N(this.v**new N(v).v);}}
for(const file of ['gee/Acclim_v2_shapefile.js','../outputs/ACCLIM_GEE_web_app_script_Private.txt']){
 const s=fs.readFileSync(file,'utf8'),c={ee:{Number:v=>new N(v),Algorithms:{If:(b,y,n)=>b?y:n}}};new vm.Script(s);vm.createContext(c);
 vm.runInContext(s.slice(s.indexOf('function photosyntheticFAPAR('),s.indexOf('function greenExcessMultiplier(')),c);
 const values={red:.2,nir:.4,NDVI:.5,EVI2:.3,FIPAR:.4,FIPAR_NDVI:.4},instant={get:k=>values[k]};
 const modes=[{faparMethod:'NIRV'},{faparMethod:'NDVI_RED_EXP'},{faparMethod:'FIPAR'}, {faparMethod:'NDVI'}, {faparMethod:'EVI2',EVI2min:.1,EVI2max:.6}, {faparMethod:'EVI2',fitEVI2:true,faparFit:{valid:true,min:new N(.1),max:new N(.6)}}, {faparMethod:'EVI2',fitEVI2:true,faparFit:{valid:false,min:new N(.1),max:new N(.6)}}];
 for(const p of modes){const base=c.photosyntheticFAPAR(instant,{...p,fexp:1}).v;for(const fexp of [.5,1,2])assert(Math.abs(c.photosyntheticFAPAR(instant,{...p,fexp}).v-base**fexp)<1e-12);}
 assert(Math.abs(c.photosyntheticFAPAR(instant,{faparMethod:'NIRV',fexp:1}).v-.2)<1e-12);
 assert(Math.abs(c.photosyntheticFAPAR(instant,{faparMethod:'NDVI_RED_EXP',fexp:1}).v-.5**1.7)<1e-12);
 values.NDVI=-.1;assert.equal(c.photosyntheticFAPAR(instant,{faparMethod:'NDVI_RED_EXP',fexp:.5}).v,0);values.NDVI=.5;
 const mult={faparMethod:'FIPAR',FAPARmult:.5,fexp:1};
 assert.equal(c.photosyntheticFAPAR(instant,mult).v,.2);
 assert.equal(c.photosyntheticFAPAR(instant,{...mult,fitVIMinP05:true}).v,.4);
 values.FIPAR=.2;assert.equal(c.photosyntheticFAPAR(instant,{...mult,useSharedFPAR:true}).v,.2);
 values.FIPAR=.16;assert(Math.abs(c.photosyntheticFAPAR(instant,{faparMethod:'FIPAR',useSharedFPAR:true,fexp:2}).v-.16)<1e-12);
 assert(s.includes('photosyntheticFAPAR(f, mergeParams(p,{useSharedFPAR:false}))'));
 assert(s.includes("p.faparMethod === 'FIPAR' && p.fexp === 1"));
}
console.log('PASS exponent across FIPAR, NDVI fallback, fixed/fitted EVI2, default identity and shared-FPAR single application');
