#!/usr/bin/env node
// Shared parser entry point for local batch conversion; never sends input data.
const fs=require('fs'),path=require('path'),core=require('../docs/tower-core.js'),zip=require('../docs/vendor/fflate-0.8.2.js');
function convert(config){
 const input=path.resolve(config.input),name=path.basename(input);let text,files={},selected=name;
 if(/\.zip$/i.test(name)){
  const bytes=fs.readFileSync(input),candidates=core.fluxFiles(core.archiveList(bytes,zip));
  selected=config.archiveFile|| (candidates.length===1?candidates[0].name:null);
  if(!selected)throw Error('Choose archiveFile; candidates: '+candidates.map(f=>f.name).join(', '));
  const result=core.archiveRead(bytes,selected,zip);text=result.text;files=result.files;
 }else text=fs.readFileSync(input);
 const info=core.inspect(text,selected);
 if(config.snapshot)files[path.basename(config.snapshot)]=fs.readFileSync(config.snapshot,'utf8');
 const meta={...core.metadata(files,info.site),...(config.metadata||{})};
 return core.parse(text,selected,{...config.options,metadata:meta});
}
module.exports={convert};
if(require.main===module){try{
 const config=JSON.parse(fs.readFileSync(process.argv[2]==='-'?0:process.argv[2],'utf8')),data=convert(config);
 fs.writeFileSync(config.output,JSON.stringify(data));
 console.log(JSON.stringify({site:data.site,months:data.months.length,start:data.start,end:data.end,output:config.output}));
}catch(e){console.error(e.message);process.exitCode=1;}}
