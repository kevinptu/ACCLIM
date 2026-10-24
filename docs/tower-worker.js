importScripts('vendor/fflate-0.8.2.js','tower-core.js');
let bytes=null,text=null,name='';
self.onmessage=async function(event){const m=event.data;const send=value=>postMessage({...value,request:m.request});try{
 if(m.action==='open'){const f=m.file;bytes=null;text=null;name=f.name;if(f.size>600000000)throw Error('Limit: 600 MB input. Extract a smaller CSV.');if(/\.zip$/i.test(name)){bytes=new Uint8Array(await f.arrayBuffer());const files=TowerCore.fluxFiles(TowerCore.archiveList(bytes,fflate));if(!files.length)throw Error('No FLUXMET/FULLSET/SUBSET HH/HR CSV found.');send({type:'files',files});}else{text=new Uint8Array(await f.arrayBuffer());send({type:'fields',info:TowerCore.inspect(text,name)});}}
 if(m.action==='choose'){name=m.name;const r=TowerCore.archiveRead(bytes,name,fflate);text=r.text;const info=TowerCore.inspect(text,name);send({type:'fields',info,meta:TowerCore.metadata(r.files,info.site)});}
 if(m.action==='parse'){if(text===null)throw Error('Choose a file first.');send({type:'result',data:TowerCore.parse(text,name,m.options)});}
 }catch(e){send({type:'error',message:e.message});}};
