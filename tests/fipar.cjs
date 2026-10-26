const fs = require('fs'), vm = require('vm'), assert = require('assert');
const paths = ['gee/Acclim_v2_shapefile.js', '../outputs/ACCLIM_GEE_web_app_script_Private.txt'];
class N {
 constructor(x) {this.v = x instanceof N ? x.v : x;}
 add(x){return new N(this.v+new N(x).v);} subtract(x){return new N(this.v-new N(x).v);}
 multiply(x){return new N(this.v*new N(x).v);} divide(x){return new N(this.v/new N(x).v);}
 min(x){return new N(Math.min(this.v,new N(x).v));}
 max(x){return new N(Math.max(this.v,new N(x).v));} pow(x){return new N(this.v**new N(x).v);}
}
for (const path of paths) {
 const s=fs.readFileSync(path,'utf8').replace(/\r/g,''); new vm.Script(s);
 const ctx={ee:{Number:x=>new N(x)},p:{ndvi_min:.15,ndvi_max:.97,ndvi_exp:1.25}};
 vm.createContext(ctx);
 vm.runInContext(s.slice(s.indexOf('function computeFIPAR('),s.indexOf('/**\n * Compute mean-state-dependent diagnostics')).replace(/\r/g,''),ctx);
 for(const ndvi of [-.2,.15,.5,.97,1]) assert(Math.abs(ctx.computeFIPAR(ndvi,ctx.p).v-(Math.max(0,ndvi-.15)/.82)**1.25)<1e-12);
 assert(!s.includes('greencover'));
 assert(!s.includes('fipar_ndvi_'));
 assert.strictEqual((s.match(/computeFIPAR\(/g)||[]).length,3); // Shared definition, forcing calculation, and calibration target.
 assert(s.includes("FIPAR: f.getNumber('FIPAR')"));
 assert(s.includes("var FIPAR = ee.Number(instant.get('FIPAR'));"));
 assert(s.includes('var em_sur = FIPAR.multiply(p.em_v)'));
 vm.runInContext(s.slice(s.indexOf('function photosyntheticFAPAR('),s.indexOf('function computeDiagnostics(')),ctx);
 ctx.ee.Algorithms={If:(c,a,b)=>c?a:b};
 function fapar(method,fit,evi,valid=true){
   Object.assign(ctx.p,{faparMethod:method,fitEVI2:fit,EVI2min:.1,EVI2max:.8});
   ctx.p.faparFit=fit?{min:new N(.2),max:new N(.6),valid}:undefined;
   ctx.instant={get:k=>k==='EVI2'?evi:.5}; ctx.FAPAR_in=new N(evi);
   return ctx.photosyntheticFAPAR(ctx.instant,ctx.p).v;
 }
 assert(Math.abs(fapar('EVI2',false,.45)-.5)<1e-12);
 assert.equal(fapar('EVI2',false,0),0); assert.equal(fapar('EVI2',false,1),1);
 assert(Math.abs(fapar('EVI2',true,.4)-.5)<1e-12);
 assert.equal(fapar('EVI2',true,.4,false),.25);
 assert.equal(fapar('NDVI',false,.4),.25);
 assert(s.includes("p.fitEVI2 && !p.faparFit"));
 assert(s.includes("overrides.faparMethod==='EVI2' && overrides.fitEVI2"));
 const block=s.slice(s.indexOf('  var PredLEcday_mean ='),s.indexOf('  return ee.Dictionary',s.indexOf('  var PredLEcday_mean =')));
 function water(fapar,fipar){Object.assign(ctx,{FAPAR:new N(fapar),FIPAR:new N(fipar),PT_mean:new N(100),fcold_mean:new N(.8),fwet_mean:new N(.1),fSM_mean:new N(.5)});vm.runInContext(block,ctx);return ['PredLEcday_mean','PredLEiday_mean','PredLEsday_mean','PredLEday_mean'].map(k=>ctx[k].v);}
 assert.deepStrictEqual(water(.2,.6),water(.9,.6));
 const expected=[43.2,6,22,71.2]; water(.2,.6).forEach((v,i)=>assert(Math.abs(v-expected[i])<1e-10));
 assert.notDeepStrictEqual(water(.2,.3),water(.2,.6));
 assert(s.includes('QYC3_mean.multiply(FAPAR).multiply(PARday_mean)'));
 assert(s.includes('QYC4_mean.multiply(FAPAR).multiply(PARday_mean)'));
 const waterTail=s.slice(s.indexOf("    var potential = r.getNumber('PT_mean')"),s.indexOf('function fitEVI2FAPAR'));
 assert(!/get(?:Number)?\('FAPAR'\)/.test(waterTail));
 assert(waterTail.includes('ee.Dictionary({FIPAR: 0,'));
}
console.log('PASS FIPAR endpoints, exact formula, water partition independence from FAPAR, photosynthesis preservation, and both script syntaxes');


