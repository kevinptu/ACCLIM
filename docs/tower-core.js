/* Shared browser/Node parser. Input files stay local; no network or storage calls. */
(function(root) {
'use strict';
const weatherTools=typeof module!=='undefined'&&module.exports?require('./tower-weather.js'):root.TowerWeather;
function* chunks(source){
  if(typeof source==='string'){yield source;return;}
  const decoder=new TextDecoder();for(let i=0;i<source.length;i+=65536)yield decoder.decode(source.subarray(i,i+65536),{stream:true});yield decoder.decode();
}
function* csv(source) {
  let row=[], value='', quoted=false,afterQuote=false,skipLF=false,first=true;
  for(const text of chunks(source))for(let i=0;i<text.length;i++) {
    const c=text[i];if(first){first=false;if(c==='\uFEFF')continue;}
    if(skipLF){skipLF=false;if(c==='\n')continue;}
    if(c==='"') {if(quoted){quoted=false;afterQuote=true;}else if(afterQuote){value+='"';quoted=true;afterQuote=false;}else quoted=true;continue;}
    afterQuote=false;
    if(c===','&&!quoted){row.push(value);value='';}
    else if((c==='\n'||c==='\r')&&!quoted){skipLF=c==='\r';row.push(value);yield row;row=[];value='';}
    else value+=c;
  }
  if(quoted)throw Error('Unclosed CSV quote.');
  if(value||row.length){row.push(value);yield row;}
}
function num(value) {
  if(value==null||!String(value).trim())return NaN;
  const n=Number(value);
  return !Number.isFinite(n)||n===-999|| /^-9999(?:\.[09]+)?$/.test(String(value).trim()) ? NaN:n;
}
function timestamp(value) {
  let s=String(value||'').trim();
  if(/^\d{12}$/.test(s))s=s.slice(0,4)+'-'+s.slice(4,6)+'-'+s.slice(6,8)+'T'+s.slice(8,10)+':'+s.slice(10,12)+':00';
  s=s.replace(' ','T');if(/^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}$/.test(s))s+=':00';
  if(!/^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}$/.test(s))throw Error('Expected a subdaily local timestamp: '+value);
  const t=Date.parse(s+'Z');
  if(!Number.isFinite(t)||new Date(t).toISOString().slice(0,19)!==s)throw Error('Invalid date: '+value);
  // UTC arithmetic preserves the supplied local clock; no timezone/DST conversion.
  return t;
}
function inspect(text,name) {
  let pre=[],header=null;
  for(const row of csv(text)) {
    const h=row.map(v=>v.trim());
    if((h.includes('TIMESTAMP_START')&&h.includes('TIMESTAMP_END'))||h.includes('MIDPOINT_DATETIME')) {header=h;break;}
    pre.push(row.join(' '));if(pre.length>500)break;
  }
  if(!header)throw Error('No supported subdaily timestamp columns. Select FLUXMET/FULLSET HH or HR data, or a Quanterra CSV.');
  const variables=header.filter(h=>['GPP','GPP_NT_VUT_REF','GPP_DT_VUT_REF','GPP_NT_CUT_REF','GPP_DT_CUT_REF'].includes(h));
  const radiation=['INCOMING_SHORTWAVE_RADIATION','SW_IN_F','SW_IN_F_MDS','SW_IN'].filter(h=>header.includes(h));
  if(!variables.length||!radiation.length)throw Error('This file has no supported GPP/incoming-shortwave columns. Avoid ERA5-only, NEE-only, daily and monthly files.');
  const meta=pre.join('\n');
  function coord(label,limit){const m=meta.match(new RegExp('Location '+label+':\\s*([-+]?\\d+(?:\\.\\d+)?)','i'));return m&&Math.abs(Number(m[1]))<=limit?Number(m[1]):null;}
  const site=meta.match(/Location Name:\s*([^\r\n]+)/i), id=String(name).match(/(?:^|_)([A-Z]{2}-[A-Za-z0-9]+)(?:_|\.)/);
  return {weatherOptions:weatherTools.options(header),header,variables,radiation,quanterra:header.includes('MIDPOINT_DATETIME'),site:site?site[1].trim():id?id[1]:name,latitude:coord('Latitude',90),longitude:coord('Longitude',180)};
}
function metadata(files,site) {
  const result={site,latitude:null,longitude:null,citation:'',productId:'',sourceFiles:[]};
  // Accept a Shuttle snapshot/manifest row only when its site ID matches.
  for(const [name,text] of Object.entries(files||{})) {
    if(text.length>5000000)continue;
    if(/\.csv$/i.test(name)) {
      const rows=csv(text.replace(/^\uFEFF/,''));let h=rows.next().value||[];h=h.map(s=>s.trim());
      const si=h.indexOf('site_id');
      if(si>=0)for(const row of rows)if(row[si]===site){
        function val(k){return row[h.indexOf(k)];}
        const lat=num(val('location_lat')),lon=num(val('location_long'));
        if(Number.isFinite(lat))result.latitude=lat;if(Number.isFinite(lon))result.longitude=lon;
        result.citation=val('product_citation')||result.citation;result.productId=val('product_id')||result.productId;result.sourceFiles.push(name);
      }
      // BADM long format: VARIABLE and DATAVALUE, scoped to this site's package.
      const vi=h.indexOf('VARIABLE'),di=h.indexOf('DATAVALUE');
      if(vi>=0&&di>=0)for(const row of rows){
        const sid=h.indexOf('SITE_ID');if(sid>=0&&row[sid]!==site)continue;
        if(sid<0&&!name.includes(site))continue;
        if(row[vi]==='LOCATION_LAT')result.latitude=num(row[di]);
        if(row[vi]==='LOCATION_LONG')result.longitude=num(row[di]);
        if(row[vi]==='PRODUCT_CITATION')result.citation=row[di];
        if(row[vi]==='PRODUCT_ID'||row[vi]==='DOI')result.productId=row[di];
      }
    }
  }
  if(!Number.isFinite(result.latitude)||Math.abs(result.latitude)>90)result.latitude=null;
  if(!Number.isFinite(result.longitude)||Math.abs(result.longitude)>180)result.longitude=null;
  return result;
}
function monthly(data,start,end) {
  const out={};
  for(const d of data.daily){if(d.date<start||d.date>end)continue;const key=d.date.slice(0,7);
    const m=out[key]||(out[key]={month:key,sum:0,n:0,secondSum:0,secondN:0,daylight:0,total:0,unknownLight:0,absent:0,days:0,radiationQC:{},neeQC:{}});
    ['sum','n','secondSum','secondN','daylight','total','unknownLight','absent'].forEach(k=>m[k]+=d[k]||0);if(d.total>0)m.days++;
    ['radiationQC','neeQC'].forEach(k=>Object.keys(d[k]||{}).forEach(q=>m[k][q]=(m[k][q]||0)+d[k][q]));
  }
  return Object.keys(out).sort().map(key=>{const m=out[key];m.gpp=m.n?m.sum/m.n:null;m.secondGpp=m.secondN?m.secondSum/m.secondN:null;m.missingPct=m.daylight?100*(m.daylight-m.n)/m.daylight:null;m.secondMissingPct=m.daylight?100*(m.daylight-m.secondN)/m.daylight:null;return m;});
}
function parse(text,name,options={}) {
  const info=inspect(text,name),h=info.header;
  const variable=options.variable||info.variables[0],rad=options.radiation||info.radiation[0],other=options.secondVariable||'';
  if(!info.variables.includes(variable)||!info.radiation.includes(rad)||(other&&!info.variables.includes(other))||other===variable)throw Error('Selected fields are not available or duplicated.');
  const gi=h.indexOf(variable),ri=h.indexOf(rad),oi=h.indexOf(other),mi=h.indexOf('MIDPOINT_DATETIME'),si=h.indexOf('TIMESTAMP_START'),ei=h.indexOf('TIMESTAMP_END');
  const rqi=h.indexOf(rad+'_QC'),nqi=h.indexOf('NEE_VUT_REF_QC');
  const weather=weatherTools.create(h,options.weather,rad,num,info.quanterra);
  const days={},seen=new Set();let started=false,units=false,previous=null,interval=null,validStart=null,validEnd=null;
  function bucket(t){const date=new Date(t).toISOString().slice(0,10);return days[date]||(days[date]={date,sum:0,n:0,secondSum:0,secondN:0,daylight:0,total:0,unknownLight:0,absent:0,radiationQC:{},neeQC:{}});}
  for(const row of csv(text)) {
    if(!started){if(row.map(v=>v.trim()).join('|')===h.join('|'))started=true;continue;}
    if(info.quanterra&&!units){if(row[gi]!=='micromol/m2/s')throw Error('Expected micromol/m2/s GPP units.');weather.checkUnits(row);units=true;continue;}
    if(row.every(v=>!v.trim()))continue;
    if(row.length!==h.length)throw Error('CSV row has a different number of columns than the header.');
    let t,dt;
    if(info.quanterra){t=timestamp(row[mi]);const a=h.indexOf('STARTING_DATETIME'),b=h.indexOf('ENDING_DATETIME');if(a>=0&&b>=0)dt=timestamp(row[b])-timestamp(row[a]);}
    else{const a=timestamp(row[si]),b=timestamp(row[ei]);dt=b-a;t=(a+b)/2;}
    if(dt!=null){if(dt!==1800000&&dt!==3600000)throw Error('Only 30- or 60-minute intervals are supported.');if(interval&&interval!==dt)throw Error('Mixed interval lengths: split the file before import.');interval=dt;}
    if(seen.has(t))throw Error('Duplicate timestamp: '+new Date(t).toISOString());seen.add(t);
    if(previous!==null&&t<previous)throw Error('Timestamps must be sorted.');
    if(previous!==null&&!interval){const gap=t-previous;if(gap===1800000||gap===3600000)interval=gap;}
    if(previous!==null&&interval){const gap=t-previous;if(gap%interval!==0)throw Error('Timestamp cadence does not match interval length.');if(gap>interval){if(gap/interval>2000000)throw Error('Unreasonable timestamp gap.');for(let missing=previous+interval;missing<t;missing+=interval)bucket(missing).absent++;}}
    previous=t;
    const d=bucket(t),g=num(row[gi]),sw=num(row[ri]),g2=other?num(row[oi]):NaN;d.total++;weather.add(d.date,row);
    if(!Number.isFinite(sw)){d.unknownLight++;continue;}
    if(sw<=0)continue;
    d.daylight++;
    if(Number.isFinite(g)){d.sum+=g;d.n++;}
    if(Number.isFinite(g2)){d.secondSum+=g2;d.secondN++;}
    if(Number.isFinite(g)||Number.isFinite(g2)){if(!validStart)validStart=d.date;validEnd=d.date;}
    for(const [key,index] of [['radiationQC',rqi],['neeQC',nqi]])if(index>=0){const q=num(row[index]);if(Number.isFinite(q))d[key][String(q)]=(d[key][String(q)]||0)+1;}
  }
  if(!interval)throw Error('Cannot establish a subdaily 30/60-minute cadence.');
  if(!validStart)throw Error('No valid daytime GPP in the selected field(s).');
  const meta=options.metadata||{};
  const result={version:3,site:meta.site||info.site,latitude:meta.latitude??info.latitude,longitude:meta.longitude??info.longitude,start:validStart,end:validEnd,units:'micromol/m2/s',filter:'INCOMING_SHORTWAVE_RADIATION > 0',gppVariable:variable,secondVariable:other,radiationVariable:rad,intervalMinutes:interval/60000,provenance:{filename:name,citation:meta.citation||'',productId:meta.productId||'',qualityNote:'All finite GPP retained, including gap-filled and negative values. NEE QC is not a direct GPP quality flag.'},daily:Object.values(days).sort((a,b)=>a.date.localeCompare(b.date))};
  result.weather=weather.finish(interval);
  result.months=monthly(result,result.daily[0].date,result.daily.at(-1).date);return result;
}
function archiveList(bytes,zip) {
  const list=[];
  zip.unzipSync(bytes,{filter:f=>{list.push({name:f.name,size:f.originalSize});return false;}});
  return list.filter(f=>!f.name.includes('__MACOSX')&&!f.name.endsWith('/'));
}
function fluxFiles(list){return list.filter(f=>/\.csv$/i.test(f.name)&&/(?:FLUXMET|FULLSET|SUBSET)_(?:HH|HR)_/i.test(f.name));}
function archiveRead(bytes,name,zip) {
  const list=archiveList(bytes,zip),target=list.find(f=>f.name===name);
  if(!target||target.size>1200000000)throw Error('Missing or oversized CSV (limit 1.2 GB expanded). Extract/select a smaller CSV.');
  let total=0;
  const extracted=zip.unzipSync(bytes,{filter:f=>{const keep=f.name===name||(/(?:BADM|_BIF_|SITE|META|MANIFEST|README|LICENSE|CITATION)/i.test(f.name)&&f.originalSize<5000000&&/\.(csv|txt|json)$/i.test(f.name));if(keep)total+=f.originalSize;if(total>1250000000)throw Error('Selected archive contents exceed size limit.');return keep;}});
  const decoder=new TextDecoder(),files={};let text='';
  for(const [key,value] of Object.entries(extracted)){if(key===name)text=value;else files[key]=decoder.decode(value);}
  return {text,files};
}
function restoreSummary(text) {
  const d=JSON.parse(text);
  if(d.version!==3||d.units!=='micromol/m2/s'||d.filter!=='INCOMING_SHORTWAVE_RADIATION > 0'||!Array.isArray(d.daily)||!d.daily.length||d.daily.length>40000)throw Error('Select a saved version 3 ACCLIM summary JSON.');
  if(d.weather&&((d.weather.version!==1&&d.weather.version!==2)||!Array.isArray(d.weather.months)))throw Error('Unsupported weather summary. Reprocess the source CSV/ZIP.');
  const seen=new Set();
  for(const r of d.daily){
    if(!/^\d{4}-\d{2}-\d{2}$/.test(r.date)||!Number.isFinite(Date.parse(r.date))||new Date(r.date).toISOString().slice(0,10)!==r.date||seen.has(r.date))throw Error('Invalid or duplicate summary date.');
    seen.add(r.date);
    for(const k of ['n','secondN','daylight','total','unknownLight','absent'])if(!Number.isInteger(r[k])||r[k]<0)throw Error('Invalid summary coverage.');
    if(!Number.isFinite(r.sum)||!Number.isFinite(r.secondSum)||r.n>r.daylight||r.secondN>r.daylight||r.daylight+r.unknownLight>r.total||(!r.n&&r.sum!==0)||(!r.secondN&&r.secondSum!==0))throw Error('Inconsistent summary coverage.');
  }
  for(const [k,limit] of [['latitude',90],['longitude',180]])if(d[k]!==null&&(!Number.isFinite(d[k])||Math.abs(d[k])>limit))throw Error('Invalid '+k);
  d.daily.sort((a,b)=>a.date.localeCompare(b.date));
  const valid=d.daily.filter(r=>r.n||r.secondN);if(!valid.length)throw Error('No valid daytime GPP.');
  d.site=String(d.site||'Tower');d.start=valid[0].date;d.end=valid.at(-1).date;
  d.months=monthly(d,d.daily[0].date,d.daily.at(-1).date);return d;
}
function transferSummary(data,start,end){
  if(!start||!end||start>end)throw Error('Choose a valid comparison date range.');
  const daily=data.daily.filter(r=>r.date>=start&&r.date<=end).map(r=>Object.fromEntries(['date','sum','n','secondSum','secondN','daylight','total','unknownLight','absent'].map(k=>[k,r[k]])));
  if(!daily.some(r=>r.n||r.secondN))throw Error('No valid daytime GPP in the selected date range.');
  const out={};for(const k of ['version','site','latitude','longitude','units','filter','gppVariable','secondVariable','radiationVariable','intervalMinutes','provenance'])out[k]=data[k];
  if(data.weather)out.weather={...data.weather,months:data.weather.months.filter(m=>m.month>=start.slice(0,7)&&m.month<=end.slice(0,7))};
  out.daily=daily;out.start=start;out.end=end;return out;
}
const api={restoreSummary,transferSummary,csv,num,timestamp,inspect,metadata,monthly,parse,archiveList,fluxFiles,archiveRead};
if(typeof module!=='undefined'&&module.exports)module.exports=api;else root.TowerCore=api;
})(typeof self!=='undefined'?self:globalThis);
