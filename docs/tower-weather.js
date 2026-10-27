/* Tower meteorology: full-day inputs; no filling or gridded substitutions. */
(function(root){
'use strict';
const choices={ta:['TA_F_MDS','TA_F','TA','AIR_TEMPERATURE'],humidity:['VPD_F_MDS','VPD_F','VPD','VAPOUR_PRESSURE_DEFICIT','RH','RELATIVE_HUMIDITY'],pa:['PA_F','PA','AIR_PRESSURE'],rain:['P_F','P','RAINFALL'],le:['LE_CORR'],h:['H_CORR']};
const fields=['ta','sw','ea','pa','rain','rn'];
function options(header){const out={};for(const k in choices)out[k]=choices[k].filter(v=>header.includes(v));return out;}
function create(header,selected,sw,num,quanterra){
 const available=options(header),source={sw};for(const k in choices)source[k]=selected&&Object.hasOwn(selected,k)?selected[k]:(available[k][0]||'');
 for(const k in source)if(source[k]&&!header.includes(source[k]))throw Error('Unknown weather field: '+source[k]);
 const index=Object.fromEntries(Object.entries(source).map(([k,v])=>[k,header.indexOf(v)]));const days={};const qcMonths={};const qcColumns=Object.entries(source).map(([k,v])=>[k,header.indexOf(v+'_QC')]).filter(([,i])=>i>=0);
 function add(date,row){const d=days[date]||(days[date]={total:0,min:null,max:null});d.total++;const qc=qcMonths[date.slice(0,7)]||(qcMonths[date.slice(0,7)]={});
  for(const [k,qi] of qcColumns){const q=num(row[qi]);if(Number.isFinite(q)){const hist=qc[k]||(qc[k]={});hist[q]=(hist[q]||0)+1;}}
  const val=k=>num(row[index[k]]),t=val('ta'),h=val('humidity');
  const es=Number.isFinite(t)?0.61121*Math.exp(17.502*t/(t+240.97)):NaN;
  const isRH=/^(RH|RELATIVE_HUMIDITY)$/.test(source.humidity);
  let ea=isRH?(h>=0&&h<=100?es*h/100:NaN):es-h*(source.humidity==='VAPOUR_PRESSURE_DEFICIT'?1:0.1);
  if(!Number.isFinite(h)||h<0||ea<0)ea=NaN;
  const v={ta:t,sw:val('sw'),ea,pa:val('pa'),rain:val('rain'),rn:Number.isFinite(val('le'))&&Number.isFinite(val('h'))?val('le')+val('h'):NaN};
  if(v.pa<=0)v.pa=NaN;if(v.rain<0)v.rain=NaN;
  for(const k of fields){const a=d[k]||(d[k]={sum:0,n:0});if(Number.isFinite(v[k])){a.sum+=v[k];a.n++;}}
  if(Number.isFinite(t)){d.min=d.min===null?t:Math.min(d.min,t);d.max=d.max===null?t:Math.max(d.max,t);}
 }
 function checkUnits(row){if(!quanterra)return;const expected={ta:'deg C',sw:'W/m2',pa:'kPa',rain:'mm',le:'W/m2',h:'W/m2',humidity:source.humidity==='RELATIVE_HUMIDITY'?'%':'kPa'};for(const k in source)if(source[k]&&row[index[k]]!==expected[k])throw Error('Unsupported weather units for '+source[k]+': '+row[index[k]]);}
 function finish(interval){const months={};const perDay=86400000/interval;
  for(const [date,d] of Object.entries(days)){const key=date.slice(0,7),m=months[key]||(months[key]={month:key,counts:{},sums:{},tmin:0,tmax:0,tempDays:0});
   for(const k of fields){m.counts[k]=(m.counts[k]||0)+d[k].n;m.sums[k]=(m.sums[k]||0)+d[k].sum;}
   if(d.ta.n===perDay){m.tmin+=d.min;m.tmax+=d.max;m.tempDays++;}
  }
  const daily=Object.entries(days).sort(([a],[b])=>a.localeCompare(b)).map(([date,d])=>{
   const counts=Object.fromEntries(fields.map(k=>[k,d[k].n]));
   const mean=k=>d[k].n===perDay?d[k].sum/perDay:null;
   return {date,period_days:1,expected:perDay,counts,complete:fields.every(k=>counts[k]===perDay),
    Tmin:counts.ta===perDay?d.min:null,Tmax:counts.ta===perDay?d.max:null,Tdaily:mean('ta'),
    SWRdaily:mean('sw'),ea:mean('ea'),Pdaily:mean('pa'),precip:counts.rain===perDay?d.rain.sum:null,Rndaily:mean('rn')};
  });
  return {daily,cadence:'daily',version:2,energyBasis:'LE_CORR + H_CORR',source,policy:'complete_calendar_period',units:{Tmin:'degC',Tmax:'degC',Tdaily:'degC',SWRdaily:'W/m2',ea:'kPa',Pdaily:'kPa',precip:'mm/period',Rndaily:'W/m2'},months:Object.values(months).sort((a,b)=>a.month.localeCompare(b.month)).map(m=>{
   const [y,mo]=m.month.split('-').map(Number),nDays=new Date(Date.UTC(y,mo,0)).getUTCDate(),expected=nDays*perDay;
   const complete=fields.every(k=>m.counts[k]===expected)&&m.tempDays===nDays;
   const mean=k=>m.counts[k]===expected?m.sums[k]/expected:null;
   return {month:m.month,expected,counts:m.counts,qc:qcMonths[m.month]||{},complete,Tmin:m.tempDays===nDays?m.tmin/nDays:null,Tmax:m.tempDays===nDays?m.tmax/nDays:null,Tdaily:mean('ta'),SWRdaily:mean('sw'),ea:mean('ea'),Pdaily:mean('pa'),precip:m.counts.rain===expected?m.sums.rain:null,Rndaily:mean('rn')};
  })};
 }
 return {add,finish,checkUnits,source};
}
// Validate before allocating any date-indexed arrays, including during year editing.
function validateRange(start,end){
 const parse=value=>{if(!/^\d{4}-\d{2}-\d{2}$/.test(value))throw Error('Enter complete dates (YYYY-MM-DD).');const t=Date.parse(value+'T00:00:00Z');if(!Number.isFinite(t)||new Date(t).toISOString().slice(0,10)!==value)throw Error('Enter valid calendar dates.');return t;};
 const a=parse(start),b=parse(end);
 if(b<a)throw Error('Start date must be on or before end date.');
 if((b-a)/86400000+1>40000)throw Error('The selected range exceeds 40,000 days. Check the year fields.');
 return {start:a,end:b};
}
// Weekly blocks start January 1; the last block ends at the year boundary.
function periods(weather,cadence,start,end,intervalMinutes){
 validateRange(start,end);
 if(!['daily','weekly','monthly'].includes(cadence))throw Error('Choose daily, weekly or monthly tower weather.');
 if(cadence==='monthly')return weather.months.filter(m=>m.month>=start.slice(0,7)&&m.month<=end.slice(0,7)).map(m=>{
  const date=m.month+'-01',a=Date.parse(date+'T00:00:00Z'),d=new Date(a),b=Date.UTC(d.getUTCFullYear(),d.getUTCMonth()+1,1);
  return {...m,date,period_days:(b-a)/86400000};
 });
 if(!Array.isArray(weather.daily))throw Error('This saved JSON contains only monthly weather. Reprocess the original CSV/ZIP for daily or weekly inputs.');
 const byDate=Object.fromEntries(weather.daily.map(d=>[d.date,d])),out=[];
 let a=Date.parse(start+'T00:00:00Z'),stop=Date.parse(end+'T00:00:00Z');
 if(cadence==='weekly'){const y=Date.UTC(new Date(a).getUTCFullYear(),0,1);a=y+Math.floor((a-y)/604800000)*604800000;}
 for(;a<=stop;){
  const b=cadence==='daily'?a+86400000:Math.min(a+604800000,Date.UTC(new Date(a).getUTCFullYear()+1,0,1));
  const dates=[];for(let t=a;t<b;t+=86400000)dates.push(new Date(t).toISOString().slice(0,10));
  const rows=dates.map(d=>byDate[d]),n=dates.length,expected=n*1440/intervalMinutes,counts={};
  for(const k of fields)counts[k]=rows.reduce((sum,r)=>sum+(r&&r.counts?r.counts[k]||0:0),0);
  const r={date:dates[0],period_days:n,expected,counts,complete:rows.every(d=>d&&d.complete)&&fields.every(k=>counts[k]===expected)};
  for(const k of ['Tmin','Tmax','Tdaily','SWRdaily','ea','Pdaily','precip','Rndaily']){
   r[k]=rows.every(d=>d&&Number.isFinite(d[k]))?rows.reduce((sum,d)=>sum+d[k],0)/(k==='precip'?1:n):null;
  }
  out.push(r);a=b;
 }
 return out;
}
const api={options,create,periods,validateRange};if(typeof module!=='undefined'&&module.exports)module.exports=api;else root.TowerWeather=api;
})(typeof self!=='undefined'?self:globalThis);
