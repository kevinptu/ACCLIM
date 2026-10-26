const fs=require('fs'),vm=require('vm'),assert=require('assert');
class N{constructor(x){this.v=x instanceof N?x.v:x}add(x){return new N(this.v+new N(x).v)}subtract(x){return new N(this.v-new N(x).v)}multiply(x){return new N(this.v*new N(x).v)}divide(x){return new N(this.v/new N(x).v)}pow(x){return new N(this.v**x)}}
for(const file of ['gee/Acclim_v2_shapefile.js','../outputs/ACCLIM_GEE_web_app_script_Private.txt']){
const s=fs.readFileSync(file,'utf8');new vm.Script(s);const ctx={ee:{Number:x=>new N(x)}};vm.createContext(ctx);vm.runInContext(s.slice(s.indexOf('function griddedDaytimeRn('),s.indexOf('function assembleForcing(')),ctx);
for(const hours of [8,12,16]){
const kelvin=25+273.15,longwave=5.67e-8*kelvin**4*(1.31*(10*1.5/kelvin)**(1/7)-.97),shortwave=.8*200*24/hours;
const result=ctx.griddedDaytimeRn(...[200,hours,25,1.5,.2,.97].map(x=>new N(x)),{sb:5.67e-8}).v;
assert(Math.abs(result-shortwave-longwave)<1e-10);assert(result>shortwave+longwave*24/hours);
}
assert.equal((s.match(/var Rnday = f.getNumber\('Rnday'\)/g)||[]).length,2);
assert(s.includes('var Rnday = towerRows ? Rndaily.multiply(ee.Number(24).divide(dl))'));
}
console.log('PASS daytime radiation at three daylengths, no longwave daylength multiplier, initialization/step wiring and unchanged tower branch');
