const fs=require('fs'),vm=require('vm'),assert=require('assert/strict');
for(const path of ['gee/Acclim_v2_shapefile.js','../outputs/ACCLIM_GEE_web_app_script_Private.txt']){
 const source=fs.readFileSync(path,'utf8');new vm.Script(source);
 const c={carbonFinite:v=>typeof v==='number'&&Number.isFinite(v)};vm.createContext(c);
 vm.runInContext(source.slice(source.indexOf('function lstNDVISettings('),source.indexOf('var UI =')),c);
 let seed=9;const rand=()=>{seed=(1664525*seed+1013904223)>>>0;return seed/4294967296;};const rows=[];
 for(let bin=0;bin<18;bin++)for(let i=0;i<70;i++){
   const red=.05+.2*(bin+rand())/18;
   // Sparse soil edge: different dates supply different brightness portions.
   const bare=i<5,nir=1.2*red+.02+(bare?rand()*.0005:.04+.15*rand());
   rows.push({red,nir,date:'2020-'+String(1+bin%12).padStart(2,'0')+'-01',month:1+bin%12,pixel_x:bin,pixel_y:i%5});
 }
 const settings={method:'geometry',range:.04,tolerance:.01,simBounds:[0,1,0,1]};
 const a=c.soilAssess(rows,settings);assert(a.candidates.length>1);assert(a.supported,JSON.stringify(a.reasons));
 assert(Math.abs(a.fit.a-1.2)<.12);assert(a.selected.supportDates>=2);assert(a.selected.supportPixels>=3);assert(a.selected.supportBins>=4);
 assert(a.selected.sim.points.length>1000);assert(Number.isFinite(a.eviMin));assert.equal(a.eviMin,a.selected.sim.evi.b);
 const again=c.soilSimulation(a.fit,settings);assert.equal(again.evi.b,a.eviMin);assert(again.points.every(p=>p.x>=0));
 assert.equal(c.soilAssess([],settings).supported,false);
 assert.equal(c.soilAssess(rows.map(r=>({...r,red:.1})),settings).supported,false);
 const oneDate=c.soilAssess(rows.map(r=>({...r,date:'2020-01-01'})),settings);assert(!oneDate.supported);assert(oneDate.candidates.every(x=>!x.plausible));
 const onePixel=c.soilAssess(rows.map(r=>({...r,pixel_x:0,pixel_y:0})),settings);assert(!onePixel.supported);
 const unknown=c.soilAssess(rows.map(({red,nir})=>({red,nir})),settings);assert(unknown.warnings.some(x=>x.includes('Date/pixel')));
 const regression=c.soilOLS([{x:0,y:.08},{x:.1,y:.28},{x:.2,y:.48}]);assert(Math.abs(regression.b-.08)<1e-12);assert(Math.abs(regression.a-2)<1e-12);
 assert.equal(c.soilOLS([{x:1,y:1},{x:1,y:2},{x:1,y:3}]),null);
 const field=value=>({getValue:()=>value});
 c.UI={soilYear:field('2020'),soilBox:field('10'),soilRings:field('0'),soilCut:field(''),soilRange:field('.04'),soilTolerance:field('.01'),soilLat:field(null),soilLon:field(null),soilMethod:field('geometry'),soilRmin:field('0'),soilRmax:field('1'),soilNmin:field('0'),soilNmax:field('1')};
 assert.equal(c.soilSearchSettings().method,'geometry');c.UI.soilRmax=field('0');assert.throws(()=>c.soilSearchSettings());
 const displayed=[],charts=[];const panel=()=>({add(v){displayed.push(v);return this;},clear(){}});
 c.UI._soilGeometrySettings=settings;c.ui={Label:value=>({value}),Panel:panel,Select:x=>x,Chart:table=>{charts.push(table);return {setChartType(){return this;},setOptions(){return this;}};}};
 c.soilGeometryDisplay(panel(),a,rows);assert.equal(charts.length,3);assert(displayed.some(x=>x.value&&x.value.includes('8192 uniform draws')));
 let status='',disabled=false,samples=0;
 const target={centroid:()=>target};
 c.UI={_requestId:1,_soilGeometrySettings:settings,runButton:{setDisabled:v=>disabled=v},resultsPanel:{style:()=>({set(){}})},resultsOverlay:{setExpanded(){}},chartHolder:panel()};
 c.soilSearchSettings=()=>({...settings,year:2020,boxKm:10,rings:0,latText:''});
 c.soilSearchTiles=()=>[{id:'center',center:[-106,32],distance:0,box:{}}];
 c.sampleSoilReference=()=>({evaluate:done=>{samples++;done({rows});}});
 c.ee={Dictionary:v=>({evaluate:done=>done(v)})};c.ui.Label=value=>({value,setValue(v){this.value=v;}});c.ui.util={setTimeout:fn=>fn()};c.setStatus=v=>status=v;
 c.runSoilSearch(target);assert.equal(samples,4);assert.equal(disabled,false);assert.equal(c.UI._soilSearch.results[0].assessment.method,'geometry');assert(status.includes('No ACCLIM model calculations'));
 assert(source.includes("image.select('StateQA').rightShift(3).bitwiseAnd(7).eq(1)"));
 assert(source.includes("value:'geometry'"));
}
console.log('PASS pooled multi-date geometry, sparse edges, metadata support, reproducible PVI regressions, controls and charts (public/private)');
