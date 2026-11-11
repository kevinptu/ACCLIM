function rtAbsorptionMethod(method) {
  return method === 'RT_CHL' || method === 'RT_LEAF' || rtAbsIndex(method)!==null;
}
// Collapse the library to the fitted LAI histogram and red background once per date.
function rtAbsorptionResponse(histogram, background) {
  var mu = RT_PAR.mu, size = mu.length;
  var values = [];
  for (var i=0; i<2*size+2; i++) values.push(0);
  var brightness = Math.max(0, Math.min(4, 10*background));
  var lo = Math.floor(brightness), hi = Math.min(4, lo+1), fraction = brightness-lo;
  var count=histogram.reduce(function(a,b){return a+b;},0);
  histogram.forEach(function(n,lai) {
    if (!n) return;
    var weight=n/(count*100000);
    for (var j=0;j<values.length;j++) {
      var diffuse=j>=2*size;
      var stride=diffuse?2:2*size, band=diffuse?j-2*size:j;
      var table=diffuse?RT_PAR.diffuse:RT_PAR.direct;
      var low=table[(lai*5+lo)*stride+band], high=table[(lai*5+hi)*stride+band];
      values[j]+=weight*(low+(high-low)*fraction);
    }
  });
  return values;
}
// Integrate with midpoint quadrature over the sunlit half-day (symmetric about noon).
function rtDayAbsorption(response, latitude, millis) {
  var day=1+(millis-Date.UTC(new Date(millis).getUTCFullYear(),0,1))/86400000;
  var dec=.409*Math.sin(2*Math.PI*day/365-1.39), lat=latitude*Math.PI/180;
  var a=Math.sin(lat)*Math.sin(dec), b=Math.cos(lat)*Math.cos(dec);
  var sunset=Math.acos(Math.max(-1,Math.min(1,-a/Math.max(1e-12,b))));
  if (!sunset) return null;
  var grid=RT_PAR.mu, result=[0,0,response[2*grid.length],response[2*grid.length+1]], sum=0;
  for(var i=0;i<16;i++) {
    var mu=Math.max(0,a+b*Math.cos(sunset*(i+.5)/16)), h=1;
    while(h<grid.length-1 && grid[h]<mu)h++;
    var w=Math.max(0,Math.min(1,(mu-grid[h-1])/(grid[h]-grid[h-1])));
    for(var c=0;c<2;c++) result[c]+=mu*((1-w)*response[2*(h-1)+c]+w*response[2*h+c]);
    sum+=mu;
  }
  if(!sum)return null;
  result[0]/=sum; result[1]/=sum;
  return result;
}
// Interpolate responses at common illumination before applying each day's solar geometry.
function rtAbsorptionPeriods(obs, periods, latitude, sky, radiation) {
  var sorted=obs.slice().sort(function(a,b){return a.millis-b.millis;});
  var curves=sorted.map(function(d){return rtAbsorptionResponse(d.lai_histogram,d.background_red);});
  var j=0, output={};
  periods.forEach(function(period){
    var totals=[0,0,0,0], energy=0, observed=0, filled=0;
    for(var t=period.millis;t<period.period_end_millis;t+=86400000){
      while(j<sorted.length && sorted[j].millis<t)j++;
      var exact=j<sorted.length && sorted[j].millis===t;
      if(!exact && (!j || j===sorted.length || sorted[j].millis-sorted[j-1].millis>32*86400000))continue;
      var curve=curves[j];
      if(!exact){
        var w=(t-sorted[j-1].millis)/(sorted[j].millis-sorted[j-1].millis);
        curve=curves[j-1].map(function(v,k){return v+w*(curves[j][k]-v);});
      }
      var q=radiation[t], day=rtDayAbsorption(curve,latitude,t);
      if(typeof q!=='number'|| !isFinite(q)||q<0)continue;
      if(exact)observed++;else filled++;
      if(!day)continue;
      for(var k=0;k<4;k++)totals[k]+=q*day[k];
      energy+=q;
    }
    var row={rt_abs_observed_days:observed,rt_abs_interpolated_days:filled,
      rt_abs_missing_pct:100*(1-(observed+filled)/period.period_days),rt_diffuse_par_fraction:sky};
    if(energy>0 && (period.model_cadence==='monthly'||observed+filled===period.period_days)){
      var v=totals.map(function(n){return n/energy;});
      row.rt_leaf_direct=v[0];row.rt_chl_direct=v[1];row.rt_leaf_diffuse=v[2];row.rt_chl_diffuse=v[3];
      row.FAPAR_leaf=(1-sky)*v[0]+sky*v[2];row.FAPAR_chl=(1-sky)*v[1]+sky*v[3];
    }
    output[period.date]=row;
  });
  return output;
}

function rtAbsIndex(method) {
  if(method==='RT_CHL_NDVI')return 'ndvi';
  if(method==='RT_CHL_EVI2')return 'evi2';
  return null;
}
// Pixels sharing LAI have identical modeled absorption on a given date.
// Store sufficient statistics so changing the diffuse fraction needs no new retrieval.
function rtIndexBins(pixels,covers) {
  var bins={};
  pixels.forEach(function(p,i){
    if(p.red+p.nir<=0)return;
    var lai=Math.round(-20*Math.log(1-covers[i]));
    var bin=bins[lai] || (bins[lai]={n:0,ndvi:0,ndvi_sq:0,evi2:0,evi2_sq:0});
    var ndvi=(p.nir-p.red)/(p.nir+p.red);
    var evi2=Math.max(0,Math.min(1,2.5*(p.nir-p.red)/(1+p.nir+2.4*p.red)));
    bin.n++;bin.ndvi+=ndvi;bin.ndvi_sq+=ndvi*ndvi;bin.evi2+=evi2;bin.evi2_sq+=evi2*evi2;
  });
  return bins;
}
function rtFitAbsorptionLine(obs,first,last,sky,index) {
  var moments=[];
  obs.forEach(function(day){
    if(day.millis<first || day.millis>last)return;
    var sx=0,sy=0,sxx=0,sxy=0,n=0;
    Object.keys(day.index_bins).forEach(function(lai){
      var hist=[];
      hist[Number(lai)]=1;
      var response=rtAbsorptionResponse(hist,day.background_red);
      var terms=rtDayAbsorption(response,day.latitude,day.millis);
      if(!terms)return;
      var y=(1-sky)*terms[1]+sky*terms[3], bin=day.index_bins[lai];
      n+=bin.n;sx+=bin[index];sxx+=bin[index+'_sq'];sy+=bin.n*y;sxy+=bin[index]*y;
    });
    if(n)moments.push({x:sx/n,y:sy/n,xx:sxx/n,xy:sxy/n,n:n});
  });
  return rtEqualDateLine(moments);
}
