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
  return {version:2,energyBasis:'LE_CORR + H_CORR',source,policy:'complete_calendar_month',units:{Tmin:'degC',Tmax:'degC',Tdaily:'degC',SWRdaily:'W/m2',ea:'kPa',Pdaily:'kPa',precip:'mm/month',Rndaily:'W/m2'},months:Object.values(months).sort((a,b)=>a.month.localeCompare(b.month)).map(m=>{
   const [y,mo]=m.month.split('-').map(Number),nDays=new Date(Date.UTC(y,mo,0)).getUTCDate(),expected=nDays*perDay;
   const complete=fields.every(k=>m.counts[k]===expected)&&m.tempDays===nDays;
   const mean=k=>m.counts[k]===expected?m.sums[k]/expected:null;
   return {month:m.month,expected,counts:m.counts,qc:qcMonths[m.month]||{},complete,Tmin:m.tempDays===nDays?m.tmin/nDays:null,Tmax:m.tempDays===nDays?m.tmax/nDays:null,Tdaily:mean('ta'),SWRdaily:mean('sw'),ea:mean('ea'),Pdaily:mean('pa'),precip:m.counts.rain===expected?m.sums.rain:null,Rndaily:mean('rn')};
  })};
 }
 return {add,finish,checkUnits,source};
}
const api={options,create};if(typeof module!=='undefined'&&module.exports)module.exports=api;else root.TowerWeather=api;
})(typeof self!=='undefined'?self:globalThis);
