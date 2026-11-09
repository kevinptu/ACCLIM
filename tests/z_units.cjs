const fs = require('fs'), vm = require('vm'), assert = require('assert/strict');
class N {
  constructor(x) { this.v = Number(x); }
  valueOf() { return this.v; }
  add(x) { return new N(this.v + Number(x)); }
  subtract(x) { return new N(this.v - x); }
  multiply(x) { return new N(this.v * x); }
  divide(x) { return new N(this.v / x); }
  pow(x) { return new N(this.v ** x); }
  sqrt() { return new N(Math.sqrt(this.v)); }
  max(x) { return new N(Math.max(this.v, x)); }
  gt(x) { return this.v > x; }
}
for (const path of ['gee/Acclim_v2_shapefile.js', '../outputs/ACCLIM_GEE_web_app_script_Private.txt']) {
  const script = fs.readFileSync(path, 'utf8');
  new vm.Script(script);
  const block = script.slice(script.indexOf('  var fD ='), script.indexOf('  var cicaC4_mean'));
  for (const ratio of [0, .5, 1, 2]) {
    for (const fDexp of [0, .5]) {
      for (const [fDexp_T, fMexp_T] of [[0,0],[1,0],[0,1],[1,1],[.5,2]]) for (const seasonal of [0, .6, 1]) {
        const c = {
          ee: {Number:x=>new N(x), Algorithms:{If:(condition,a,b)=>condition?a:b}},
          instant: {contains:()=>true, get:()=>seasonal},
          p: {GPPmean2GPPmax:.67, QRC3:12, PT_coef:1.26, fresp:.1, gamma_psy:.066, fDexp, fDexp_T, fMexp_T}
        };
        for (const [key,value] of Object.entries({FAPAR:ratio*.4, FIPAR:.4, s_mean:.15,
          VPDday_mean:1.5, PARday_mean:1000, Rnday_mean:250, Pmean:101325,
          gamma_star_mean:4, ca_mean:410})) c[key] = new N(value);
        vm.runInNewContext(block,c);
        const zMol = 1.5*.67*(2500000*.018)/(12*1.26*1.1)*(.15+.066)/.15*(1500/101325)*(.001/250)*ratio;
        const factor = (1.5**(-fDexp))**fDexp_T * seasonal**fMexp_T;
        assert(Math.abs(c.transpirationFactor-factor)<1e-12);
        assert.equal(c.validTransRatio,factor>0);
        Object.assign(c,{PT_mean:new N(100),fcold_mean:new N(.8),fwet_mean:new N(.25)});
        vm.runInNewContext(script.slice(script.indexOf('  var PredLEcday_mean ='),script.indexOf('  var PredLEiday_mean =')),c);
        assert(Math.abs(c.PredLEcday_mean-24*factor)<1e-10);
        const z = zMol*1e6*(factor>0 ? 1.5**(-fDexp)*seasonal/factor : 1);
        assert(Math.abs(c.zC3-z)<1e-9);
        assert(c.faparFiparValid);
        const ci = c.cicaC3_mean*410, gamma = 4e6/101325;
        assert(ci>=gamma && ci<=410);
        assert(Math.abs((410-ci)-z*(ci-gamma)/(ci+2*gamma))<1e-9);
      }
    }
  }
}
console.log('PASS C3 z units and ci/ca against the WUE equation for both scripts');
