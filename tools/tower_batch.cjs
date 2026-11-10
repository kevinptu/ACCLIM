const fs = require('fs'), path = require('path'), vm = require('vm');
const core = require('../docs/tower-core.js'), zip = require('../docs/vendor/fflate-0.8.2.js');
const [archive, output, project, parameters] = process.argv.slice(2);
const input = fs.readFileSync(archive);
async function main() {
  const entries = core.fluxFiles(core.archiveList(input, zip));
  if (!entries.length) throw Error('No subdaily FLUXNET CSVs found.');
  let ee, ctx;
  if (output) {
    ee = require('@google/earthengine');
    ee.data.setAuthToken('', 'Bearer', process.env.ACCLIM_EE_TOKEN, 3600, [], undefined, false);
    ee.data.setProject(project);
    await new Promise((ok, fail) => ee.initialize(null, null, ok, fail));
    const source = fs.readFileSync(path.join(__dirname, '../gee/Acclim_v2_shapefile.js'), 'utf8');
    ctx = vm.createContext({ee, print: ()=>{}, console});
    vm.runInContext(source.slice(0, source.indexOf('var UI =')), ctx);
    vm.runInContext(source.slice(source.indexOf('function towerWeatherPlan('), source.indexOf('function setupTowerImport(')), ctx);
    ctx.UI = {_requestId:1}; ctx.setStatus = ()=>{};
  }
  const records = [];
  for (const {name} of entries) {
    const {text, files} = core.archiveRead(input, name, zip);
    const info = core.inspect(text, name), metadata = core.metadata(files, info.site);
    if (![metadata.latitude, metadata.longitude].every(Number.isFinite)) throw Error('Missing coordinates: '+name);
    const tower = core.parse(text, name, {variable:'GPP_DT_VUT_REF', secondVariable:'GPP_NT_CUT_REF', metadata});
    const start=tower.daily[0].date, end=tower.daily.at(-1).date;
    console.log(info.site, start, end, tower.daily.length+' days', metadata.latitude, metadata.longitude);
    if (!output) continue;
    const options = {weatherSource:'TERRACLIMATE', satelliteSource:'MODIS16', timestep:'monthly'};
    const geom = ee.Geometry.Point([metadata.longitude, metadata.latitude]);
    const fc = await new Promise((ok, fail)=>ctx.runForGeometryAsync(geom,start,end,JSON.parse(parameters),options,(v,e)=>e?fail(Error(e)):ok(v)));
    const data = await new Promise((ok,fail)=>fc.evaluate((v,e)=>e?fail(Error(e)):ok(v)));
    for (const f of data.features) {
      const r=f.properties, m=tower.months.find(m=>m.month===r.date.slice(0,7));
      const days = tower.daily.filter(d => d.date.slice(0,7) === r.date.slice(0,7));
      const leN = days.reduce((sum,d) => sum+d.leN,0);
      const towerET = leN ? days.reduce((sum,d) => sum+d.leSum,0)/leN*86400/2450000 : null;
      records.push([info.site,r.date,metadata.latitude,metadata.longitude,r.PredGPPmean,m?.gpp,m?.secondGpp,m?.missingPct,r.ET_mm_day,towerET]);
    }
  }
  if(output) fs.writeFileSync(output,'site,date,latitude,longitude,modeled_daytime_gpp,tower_daytime_gpp_dt,tower_daytime_gpp_nt_cut,tower_gpp_missing_pct,modeled_et_mm_day,tower_le_corr_mm_day\n'+records.map(r=>r.join(',')).join('\n')+'\n');
}
main().catch(e=>{console.error(e.message);process.exitCode=1;});

