// Retrieve the same date-specific spatial median used by the standalone RT analysis.
function rtCoverPeriods(observations, periods, gapDays) {
  var obs = observations.slice().sort(function(a,b){return a.millis-b.millis;});
  var j = 0, result = {};
  periods.forEach(function(period) {
    var values = [], measured = 0;
    for (var t = period.millis; t < period.period_end_millis; t += 86400000) {
      while (j < obs.length && obs[j].millis < t) j++;
      var right = obs[j], left = obs[j-1];
      if (right && right.millis === t) { values.push(right.cover); measured++; }
      else if (right && left && right.millis-left.millis <= gapDays*86400000) {
        values.push(left.cover + (right.cover-left.cover)*(t-left.millis)/(right.millis-left.millis));
      }
    }
    var entry = {rt_observed_days:measured, rt_interpolated_days:values.length-measured,
      rt_missing_pct:100*(1-values.length/period.period_days)};
    if (values.length && (period.model_cadence==='monthly' || values.length===period.period_days)) {
      entry.rt_fcover = period.model_cadence==='monthly' ? RT.median(values) :
        values.reduce(function(a,b){return a+b;},0)/values.length;
    }
    result[period.date]=entry;
  });
  return result;
}
function fetchRTCover(geom,start,end,km,request,done) {
  var gap = SATELLITE_SOURCES.MODIS8.gap;
  var first = carbonDateMillis(start)-gap*86400000, last = carbonDateMillis(end)+gap*86400000;
  var observations = [];
  geom.centroid(1).coordinates().evaluate(function(xy,error) {
    if(request!==UI._requestId)return;
    if(error){done(null,error);return;}
    var region = rtBox(xy[0],xy[1],km);
    var collection = ee.ImageCollection('MODIS/061/MYD09A1').filterBounds(region)
      .filterDate(ee.Date(first).advance(-7,'day'),ee.Date(last).advance(1,'day')).sort('system:time_start');
    collection.aggregate_array('system:id').evaluate(function(ids,error) {
      if(request!==UI._requestId)return;
      if(error){done(null,error);return;}
      var index=0;
      function nextImage() {
        if(request!==UI._requestId)return;
        if(index===ids.length){done(observations);return;}
        var id=ids[index++];
        setStatus('Retrieving PROSAIL fcover: composite '+index+' of '+ids.length);
        rtSample(ee.Image(id),region).evaluate(function(data,error) {
          if(request!==UI._requestId)return;
          if(error){done(null,error);return;}
          var groups=rtGroups(data.features,Number(id.split('/').pop().slice(0,4)),first,last);
          var days=Object.keys(groups).filter(function(d){return groups[d].length>=4;}).sort();
          function nextDay() {
            if(request!==UI._requestId)return;
            if(!days.length){nextImage();return;}
            var date=days.shift(), fit=RT.createFit(groups[date]);
            setStatus('Fitting PROSAIL fcover for '+date+' ('+groups[date].length+' pixels)');
            function advance() {
              if(request!==UI._requestId)return;
              try {
                if(!fit.step()){ui.util.setTimeout(advance,0);return;}
                var covers=fit.bestCovers();
                var histogram=Array.apply(null,Array(81)).map(function(){return 0;});
                covers.forEach(function(c){histogram[Math.round(-20*Math.log(1-c))]++;});
                observations.push({millis:carbonDateMillis(date),cover:RT.median(covers),
                  index_bins:rtIndexBins(groups[date],covers),lai_histogram:histogram,background_red:fit.bestRed(),latitude:xy[1],
                  calibration:rtDayMoments(groups[date],covers)});
                nextDay();
              } catch(e){done(null,e.message);}
            }
            advance();
          }
          nextDay();
        });
      }
      nextImage();
    });
  });
}

// Equal total weight for each acquisition date. Keep within-date pixel covariance.
function rtDayMoments(pixels,covers) {
  var sx=0, sy=0, sxx=0, sxy=0, n=0;
  pixels.forEach(function(p,i) {
    var sum=p.red+p.nir;
    if(sum<=0)return;
    var x=(p.nir-p.red)/sum, y=covers[i];
    sx+=x; sy+=y; sxx+=x*x; sxy+=x*y; n++;
  });
  return n ? {x:sx/n,y:sy/n,xx:sxx/n,xy:sxy/n,n:n} : null;
}
function rtFitCoverLine(obs,first,last) {
  var days=obs.filter(function(d){return d.calibration && d.millis>=first && d.millis<=last;});
  var fit=rtEqualDateLine(days.map(function(d){return d.calibration;}));
  return fit ? {m_cover:fit.slope,b_cover:fit.intercept,cover_fit_dates:fit.dates,cover_fit_pixels:fit.pixels} : null;
}
function rtEqualDateLine(moments) {
  var sx=0, sy=0, sxx=0, sxy=0, dates=moments.length, pixels=0;
  moments.forEach(function(d){
    sx+=d.x; sy+=d.y; sxx+=d.xx; sxy+=d.xy; pixels+=d.n;
  });
  var variance=dates ? sxx/dates-Math.pow(sx/dates,2) : 0;
  if(variance<=1e-12)return null;
  var slope=(sxy/dates-(sx/dates)*(sy/dates))/variance;
  return {slope:slope,intercept:(sy-slope*sx)/dates,dates:dates,pixels:pixels};
}
