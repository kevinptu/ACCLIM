const fs=require('fs'),vm=require('vm'),assert=require('assert');
class N{constructor(v){this.v=v instanceof N?v.v:v}add(x){return new N(this.v+new N(x).v)}multiply(x){return new N(this.v*new N(x).v)}divide(x){return new N(this.v/new N(x).v)}max(x){return new N(Math.max(this.v,x))}gt(x){return this.v>x}}
for(const path of ['gee/Acclim_v2_shapefile.js','../outputs/ACCLIM_GEE_web_app_script_Private.txt']){
const s=fs.readFileSync(path,'utf8');const block=s.slice(s.indexOf('  var faparFiparValid'),s.indexOf('  // Quadratic solves'));
for(const ratio of [0,.5,1,2]){let c={ee:{Number:x=>new N(x)},p:{GPPmean2GPPmax:.67,QRC3:12,PT_coef:1.26,fresp:.1,gamma_psy:.066}};for(const [k,v]of Object.entries({FAPAR:ratio*.4,FIPAR:.4,s_mean:.15,VPDday_mean:1.5,PARday_mean:1000,Rnday_mean:250,Pmean:101325}))c[k]=new N(v);vm.runInNewContext(block,c);
const zMol=1.5*.67*(2500000*.018)/(12*1.26*1.1)*(.15+.066)/.15*(1500/101325)*(.001/250)*ratio;
assert(Math.abs(c.zC3.v-zMol*1e6)<1e-9);assert(c.faparFiparValid);}}
console.log('PASS z ppm versus independent mol/mol calculation and FAPAR/FIPAR scaling');
