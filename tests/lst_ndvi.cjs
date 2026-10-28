const assert=require('assert/strict'),fs=require('fs'),vm=require('vm');
for(const path of ['gee/Acclim_v2_shapefile.js','../outputs/ACCLIM_GEE_web_app_script_Private.txt']) {
 const s=fs.readFileSync(path,'utf8');new vm.Script(s);
 const c={carbonFinite:x=>typeof x==='number'&&isFinite(x)};vm.createContext(c);
 vm.runInContext(s.slice(s.indexOf('function lstNDVISettings('),s.indexOf('var UI =')),c);
 assert.equal(c.lstNDVISettings('2014','10').boxKm,10);assert.throws(()=>c.lstNDVISettings('0020','10'));assert.throws(()=>c.lstNDVISettings('2014','0'));assert.throws(()=>c.lstNDVISettings('2014','Infinity'));
 const rows=[{date:'2014-01-01',month:1,NDVI:0,LST_C:20},{date:'2014-02-02',month:2,NDVI:.5,LST_C:30},{date:'2014-01-01',month:1,NDVI:null,LST_C:20}];
 assert.equal(c.lstNDVIChartTable(rows,'all').rows.length,2);assert.equal(c.lstNDVIChartTable(rows,'all').cols.length,13);assert.equal(c.lstNDVIChartTable(rows,'2014-01-01').rows[0].c[0].v,0);
 let samples=0,charts=0,note='',disabled=false;
 c.lstNDVIBox=()=>({box:true});c.sampleLSTNDVIPeriods=(g,y,k,m)=>({evaluate:done=>{samples++;done({rows:[{date:y+'-'+String(m).padStart(2,'0')+'-01',month:m,NDVI:.5,LST_C:25,NDVI_valid_fraction:1}],paired_periods:1,reflectance_periods:1,pixels:117,pixel_size_m:926.6});}});
 c.UI={_requestId:1,lstNDVIYear:{getValue:()=> '2014'},lstNDVIBoxKm:{getValue:()=> '10'},runButton:{setDisabled:v=>disabled=v},resultsPanel:{style:()=>({set(){}})},resultsOverlay:{setExpanded(){}},_chartCard:{add(){}}};
 c.ui={Chart:table=>({setChartType(){return this;},setOptions(){return this;}}),Select:x=>x,Label:x=>x};
 c.addResultChart=()=>charts++;c.chartReplacer=()=>()=>{};c.addChartNote=()=>{};c.setStatus=x=>note=x;
 c.runLSTNDVI({});assert.equal(samples,4);assert.equal(charts,1);assert.equal(disabled,false);assert(note.includes('No ACCLIM model calculations'));
 c.sampleLSTNDVIPeriods=()=>({evaluate:done=>done({rows:[],paired_periods:0,reflectance_periods:0,pixels:117,pixel_size_m:926.6})});
 c.runLSTNDVI({});assert(note.includes('No paired'));assert.equal(disabled,false);
 let plotted=0,cleared=0,modeled=0;
 c.clearResults=()=>cleared++;c.setStatus=()=>{};c.resolveGeometry=()=>({point:true});c.runLSTNDVI=()=>plotted++;c.runForGeometryAsync=()=>modeled++;
 c.UI={analysisMode:{getValue:()=> 'lst_ndvi'},startI:{input:{getValue:()=>''}},endI:{input:{getValue:()=>''}},latitude:{getValue:()=> '38'},longitude:{getValue:()=> '-120'},useCoordinates:{getValue:()=>true}};
 vm.runInContext(s.slice(s.indexOf('function onRun()'),s.indexOf('function onReset()')),c);c.onRun();
 assert.equal(plotted,1);assert.equal(cleared,1);assert.equal(modeled,0);
 assert(s.includes("ee.Filter.equals({leftField:'system:time_start',rightField:'system:time_start'})"));
 assert(s.includes("ndvi.reduceResolution({reducer:ee.Reducer.mean()"));
}
console.log('PASS settings, annual/interval scatter tables, valid zero NDVI, and standalone mode bypassing model/date requirements');
