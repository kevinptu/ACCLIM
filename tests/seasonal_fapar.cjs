const fs=require('fs'),vm=require('vm'),assert=require('assert/strict');
const val=x=>x instanceof Num?x.v:x;
class Num{constructor(x){this.v=Number(val(x));}add(x){return new Num(this.v+val(x));}subtract(x){return new Num(this.v-val(x));}divide(x){return new Num(this.v/val(x));}pow(x){return new Num(this.v**val(x));}max(x){return new Num(Math.max(this.v,val(x)));}eq(x){return new Num(this.v===val(x));}neq(x){return new Num(this.v!==val(x));}lt(x){return new Num(this.v<val(x));}gt(x){return new Num(this.v>val(x));}and(x){return new Num(this.v&&val(x));}or(x){return new Num(this.v||val(x));}}
const dict=o=>({get:k=>o[k],contains:k=>Object.hasOwn(o,k)});
const date=t=>({advance:(n,u)=>{assert.equal(u,'second');return date(val(t)+n*1000);},get:k=>{assert.equal(k,'year');return new Date(val(t)).getUTCFullYear();},getRelative:(u,v)=>{assert.equal(u,'day');assert.equal(v,'year');const d=new Date(val(t));return Math.floor((val(t)-Date.UTC(d.getUTCFullYear(),0,1))/86400000);}});
const ee={Number:x=>new Num(x),Date:date,Dictionary:dict,Algorithms:{If:(c,a,b)=>val(c)?a:b}};
for(const file of ['gee/Acclim_v2_shapefile.js','../outputs/ACCLIM_GEE_web_app_script_Private.txt']){
 const s=fs.readFileSync(file,'utf8'),c={ee};vm.createContext(c);new vm.Script(s);
 vm.runInContext(s.slice(s.indexOf('function seasonalFAPARState('),s.indexOf('function acclimStep(')),c);
 function step(f,lat,end,state={step:0}){return c.seasonalFAPARState(new Num(f),new Num(lat),new Num(Date.parse(end+'T00:00:00Z')),dict(state));}
 let r=step(.8,35,'2023-12-31');r=step(.4,35,'2024-01-01',{step:1,FAPARmax:r.get('FAPARmax'),FAPAR_season:r.get('FAPAR_season')});assert.equal(val(r.get('FAPARmax')),.8);assert(Math.abs(val(r.get('seasonal_ratio'))-.5)<1e-12);
 r=step(.2,35,'2024-01-02',{step:2,FAPARmax:.8,FAPAR_season:2023});assert.equal(val(r.get('FAPARmax')),.2);assert.equal(val(r.get('seasonal_ratio')),1);
 // Southern Jan 1 does not reset. Leap-year June monthly interval contains DOY 182.
 r=step(.4,-35,'2024-01-02',{step:1,FAPARmax:.8,FAPAR_season:2023});assert.equal(val(r.get('FAPARmax')),.8);
 r=step(.2,-35,'2024-07-01',{step:2,FAPARmax:.8,FAPAR_season:2023});assert.equal(val(r.get('FAPARmax')),.2);
 // Non-leap June ends before DOY 182; July 1 resets.
 r=step(.2,-35,'2023-07-01',{step:2,FAPARmax:.8,FAPAR_season:2022});assert.equal(val(r.get('FAPARmax')),.8);
 r=step(.2,-35,'2023-07-02',{step:3,FAPARmax:.8,FAPAR_season:2022});assert.equal(val(r.get('FAPARmax')),.2);
 assert.equal(val(step(0,35,'2024-01-02').get('seasonal_ratio')),0);
 assert.equal(val(step(1e-15,35,'2024-01-02').get('seasonal_ratio')),1);
 assert(s.includes("value: 'FIPAR_GE'"));if(file.startsWith('gee/')){
   assert(/\], value: 'FIPAR',\s*style:/.test(s));
   assert(/\],\s*value: 'RT_NDVI',\s*style:/.test(s));
 }else assert(s.includes("value: 'EVI2_FIXED_POWER',\n  style:"));assert(s.includes('fitNDVIminGE: false'));
}
console.log('PASS hemisphere resets, leap/non-leap boundaries, containing-month reset, zero/tiny FAPAR and default/optional menu wiring');
