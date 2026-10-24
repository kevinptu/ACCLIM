const assert=require('node:assert/strict'),fs=require('fs'),vm=require('vm'),core=require('../docs/tower-core.js'),zip=require('../docs/vendor/fflate-0.8.2.js');
const header='TIMESTAMP_START,TIMESTAMP_END,GPP_NT_VUT_REF,GPP_DT_VUT_REF,SW_IN_F,SW_IN_F_QC,NEE_VUT_REF_QC\n';
const text=header+[
 '202401010000,202401010030,99,99,0,0,0',
 '202401010030,202401010100,2,4,10,0,0',
 '202401010100,202401010130,-9999,-9999.9,20,1,1',
 '202401010130,202401010200,6,-2,30,1,2',
 '202401010200,202401010230,100,100,-9999,2,3',
 '202402010000,202402010030,-9999,-9999,10,0,0',
 '202402010030,202402010100,8,12,20,1,1'].join('\n');
const name='AMF_US-Test_FLUXNET_FLUXMET_HH_2024-2024_v1.3_r1.csv',options={secondVariable:'GPP_DT_VUT_REF'};
const data=core.parse(text,name,options),jan=data.months[0];
assert.equal(jan.gpp,4);assert.equal(jan.secondGpp,1);assert.ok(Math.abs(jan.missingPct-100/3)<1e-10);assert.equal(jan.unknownLight,1);assert.ok(jan.absent>0);
assert.equal(core.monthly(data,'2024-02-01','2024-02-01')[0].missingPct,50);
assert.equal(core.num('-9999.999'),NaN);assert.equal(core.num('-2'),-2);
assert.throws(()=>core.parse(text+'\n202402010030,202402010100,8,12,20,1,1',name),/Duplicate/);
assert.throws(()=>core.timestamp('202402300000'),/Invalid date/);
assert.throws(()=>core.parse(header+'202401010000,202401020000,1,2,3,0,0',name),/30- or 60/);
const badm='SITE_ID,VARIABLE,DATAVALUE\nUS-Test,LOCATION_LAT,34\nUS-Test,LOCATION_LONG,-104';
const encoded=zip.zipSync({[name]:zip.strToU8(text),'US-Test_BADM.csv':zip.strToU8(badm),'ERA5_HH.csv':zip.strToU8('not flux')});
assert.equal(core.fluxFiles(core.archiveList(encoded,zip)).length,1);
const unzipped=core.archiveRead(encoded,name,zip);assert.equal(new TextDecoder().decode(unzipped.text),text);
assert.equal(core.metadata(unzipped.files,'US-Test').latitude,34);
const snapshot='site_id,location_lat,location_long,product_id,product_citation\nUS-Other,0,0,x,wrong\nUS-Test,34,-104,doi,test citation';
assert.equal(core.metadata({'snapshot.csv':snapshot},'US-Test').citation,'test citation');
// Quanterra compatibility: metadata, explicit units, local midpoint, negative GPP.
const q='"Location Name: US_NM_001"\n"Location Latitude: 34.43319 deg"\n"Location Longitude: -104.1172 deg"\nMIDPOINT_DATETIME,GPP,INCOMING_SHORTWAVE_RADIATION\nlocal,micromol/m2/s,W/m2\n2024-01-01 12:15:00,-2,50\n2024-01-01 12:45:00,-9999,60';
const qd=core.parse(q,'quanterra.csv');assert.equal(qd.latitude,34.43319);assert.equal(qd.months[0].gpp,-2);assert.equal(qd.months[0].missingPct,50);
// UI script syntax and both Earth Engine consumers, including old v2 payloads.
for(const file of ['tower-core.js','tower-worker.js','tower-ui.js'])new vm.Script(fs.readFileSync(__dirname+'/../docs/'+file,'utf8'));
for(const scriptPath of [__dirname+'/../gee/Acclim_v2_shapefile.js',...(process.env.ACCLIM_PRIVATE_SCRIPT?[process.env.ACCLIM_PRIVATE_SCRIPT]:[])]){
 const script=fs.readFileSync(scriptPath,'utf8');new vm.Script(script);
 const ctx={};vm.createContext(ctx);vm.runInContext('function carbonFinite(x){return typeof x==="number"&&isFinite(x);}\n'+script.slice(script.indexOf('function parseTowerMonthly('),script.indexOf('function setupTowerImport(')),ctx);
 const parsed=ctx.parseTowerMonthly(JSON.stringify(data)),m=ctx.towerMonthsInRange(parsed,'2024-01-01','2024-02-01');assert.equal(m['2024-01'].gpp,4);assert.equal(m['2024-01'].secondGpp,1);assert.equal(m['2024-02'].missingPct,50);
 const old={version:2,site:'test',latitude:null,longitude:null,units:'micromol/m2/s',filter:'INCOMING_SHORTWAVE_RADIATION > 0',daily:[{date:'2024-01-01',sum:4,n:2}]};assert.equal(ctx.parseTowerMonthly(JSON.stringify(old)).start,'2024-01-01');
 assert.throws(()=>ctx.parseTowerMonthly(JSON.stringify({...data,daily:[{...data.daily[0],n:9999}]})),/coverage/);
}
if(process.argv[2]){fs.mkdirSync(process.argv[2],{recursive:true});fs.writeFileSync(process.argv[2]+'/fixture.zip',encoded);fs.writeFileSync(process.argv[2]+'/snapshot.csv',snapshot);}
const messages=[],workerContext={self:{},TowerCore:core,fflate:zip,Uint8Array,importScripts:()=>{},postMessage:m=>messages.push(m)};
vm.createContext(workerContext);vm.runInContext(fs.readFileSync(__dirname+'/../docs/tower-worker.js','utf8'),workerContext);
(async()=>{
 await workerContext.self.onmessage({data:{action:'open',request:7,file:{name:'fixture.zip',size:encoded.length,arrayBuffer:async()=>encoded}}});
 assert.equal(messages.at(-1).type,'files');assert.equal(messages.at(-1).request,7);
 await workerContext.self.onmessage({data:{action:'choose',name,request:7}});assert.equal(messages.at(-1).type,'fields');
 await workerContext.self.onmessage({data:{action:'parse',options,request:8}});assert.equal(messages.at(-1).type,'result');assert.equal(messages.at(-1).data.months[0].gpp,4);assert.equal(messages.at(-1).request,8);
 console.log('Parser, ZIP, metadata, coverage, worker, Quanterra compatibility and GEE checks passed.');
})().catch(e=>{console.error(e);process.exitCode=1;});
