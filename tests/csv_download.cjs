const fs=require('fs'),vm=require('vm'),assert=require('assert');
for(const file of ['gee/Acclim_v2_shapefile.js','../outputs/ACCLIM_GEE_web_app_script_Private.txt']){
 const source=fs.readFileSync(file,'utf8');new vm.Script(source);
 let pending,shown=false,url='',status='',value='';
 const link={style:()=>({set:x=>{shown=x.shown}}),setValue:x=>value=x,setUrl:x=>url=x};
 const row={setGeometry:g=>{assert.equal(g,null);return row},set:x=>{assert.equal(x.export_site,'Test site');return row}};
 const collection={map:fn=>{fn(row);return {getDownloadURL:args=>{assert.equal(args.format,'csv');assert.equal(args.filename,'acclim_Test_site_monthly_2021-01-01_to_2023-12-31');pending=args.callback}}}};
 const UI={_lastRun:{site:'Test site',start:'2021-01-01',end:'2023-12-31',weatherSource:'gridded'},_requestId:1,csvDownloadLink:link,_lastOutput:collection};
 const ctx={UI,setStatus:x=>status=x};vm.createContext(ctx);
 vm.runInContext(source.slice(source.indexOf('function prepareCsvDownload('),source.indexOf('// SECTION 9')),ctx);
 ctx.onExportCsv();assert(status.includes('wait'));assert(!pending);
 UI._monthlyExportReady=true;ctx.onExportCsv();pending('https://example.com/download');assert(shown);assert.equal(value,'Download monthly CSV');
 ctx.onExportCsv();pending(null,'Server error');assert(!shown);assert(status.includes('Server error'));
 ctx.onExportCsv();UI._requestId++;pending('stale');assert(!shown);assert.notEqual(url,'stale');
}
console.log('PASS CSV readiness, download callback, errors, stale-run protection and both script syntaxes');
