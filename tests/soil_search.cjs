const fs=require('fs'),vm=require('vm'),assert=require('assert/strict');
for(const path of ['gee/Acclim_v2_shapefile.js','../outputs/ACCLIM_GEE_web_app_script_Private.txt']){
 const source=fs.readFileSync(path,'utf8');new vm.Script(source);
 const c={carbonFinite:v=>typeof v==='number'&&isFinite(v)};vm.createContext(c);
 vm.runInContext(source.slice(source.indexOf('function lstNDVISettings('),source.indexOf('var UI =')),c);
 const rows=[];
 for(let month=1;month<=6;month++)for(let i=0;i<100;i++){const red=.1+.15*i/99;rows.push({red,nir:1.2*red+.025,month});}
 const settings={cut:.25,range:.04,tolerance:.01};
 const a=c.soilAssess(rows,settings);assert(a.supported);assert(Math.abs(a.fit.a-1.2)<1e-5);assert(Math.abs(a.fit.b-.025)<1e-5);
 assert.equal(c.soilAssess([],settings).supported,false);
 assert.equal(c.soilAssess(rows.filter(r=>r.month===1),settings).supported,false);
 assert.equal(c.soilAssess(rows.map(r=>({...r,red:.12})),settings).supported,false);
 const contaminated=rows.map((r,i)=>({...r,nir:r.nir+(i%3===0?.06:0)}));assert(Math.abs(c.soilAssess(contaminated,settings).fit.a-1.2)<.01);
 const curved=rows.map(r=>({...r,nir:.02+1.2*r.red+10*(r.red-.17)**2}));assert.equal(c.soilAssess(curved,{cut:.8,range:.04,tolerance:.003}).supported,false);
 const field=value=>({getValue:()=>value});
 c.UI={soilYear:field('2014'),soilBox:field('10'),soilRings:field('1'),soilCut:field('.25'),soilRange:field('.04'),soilTolerance:field('.01'),soilLat:field(''),soilLon:field('')};
 assert.equal(c.soilSearchSettings().rings,1);c.UI.soilLat=field('32');assert.throws(()=>c.soilSearchSettings());c.UI.soilLon=field('-106');assert.equal(c.soilSearchSettings().lat,32);
 c.UI.soilBox=field('20');assert.throws(()=>c.soilSearchSettings());
 // End-to-end standalone asynchronous callback path without model execution.
 let sampling=0,status='',disabled=false;const panels=[];
 const panel=()=>({add(v){panels.push(v);return this;},clear(){},style(){return {set(){}};}});
 c.UI={_requestId:1,runButton:{setDisabled:v=>disabled=v},resultsPanel:panel(),resultsOverlay:{setExpanded(){}},chartHolder:panel()};
 c.soilSearchSettings=()=>({...settings,year:2014,boxKm:10,rings:0,latText:''});
 const target={centroid:()=>target};c.soilSearchTiles=()=>[{id:'center',center:[-106,32],distance:0,box:{}}];
 c.sampleSoilReference=()=>({evaluate:done=>{sampling++;done({rows:rows.filter(r=>r.month<=3)});}});
 c.ee={Dictionary:v=>({evaluate:done=>done(v)})};
 c.ui={Label:v=>({value:v,setValue(x){this.value=x;}}),Panel:panel,Select:v=>v,Chart:v=>({setChartType(){return this;},setOptions(){return this;}}),util:{setTimeout:fn=>fn()}};
 c.setStatus=v=>status=v;c.runSoilSearch(target);
 assert.equal(sampling,4);assert.equal(disabled,false);assert(c.UI._soilSearch.results[0].assessment.supported);assert(status.includes('No ACCLIM model calculations'));
 // Sampling errors release the button; cancellation ignores late callbacks.
 c.sampleSoilReference=()=>({evaluate:done=>done(null,'server failure')});c.runSoilSearch(target);assert.equal(disabled,false);assert(status.includes('server failure'));
 let pending;c.sampleSoilReference=()=>({evaluate:done=>pending=done});c.runSoilSearch(target);c.UI._requestId++;pending({rows});assert.equal(sampling,4);
 let run=0,model=0;c.clearResults=()=>{};c.resolveGeometry=()=>target;c.runSoilSearch=()=>run++;c.runForGeometryAsync=()=>model++;
 c.UI={analysisMode:field('soil_search'),startI:{input:field('')},endI:{input:field('')},latitude:field('32'),longitude:field('-106'),useCoordinates:field(true)};
 vm.runInContext(source.slice(source.indexOf('function onRun()'),source.indexOf('function onReset()')),c);c.onRun();assert.equal(run,1);assert.equal(model,0);
}
console.log('PASS soil fit, unsupported cases, controls, standalone routing, callbacks, errors and cancellation (public/private)');
