"""Download selected Shuttle sites and prepare ACCLIM summaries using the browser parser.

Requires Python 3.11+, Node.js 18+, and fluxnet_shuttle for online mode.
Offline --input mode needs only Python and Node. Nothing is uploaded.
"""
import argparse
import asyncio
import csv
import json
from pathlib import Path
import re
import subprocess

def convert(source, output, args, metadata=None):
    config = dict(input=str(source.resolve()), output=str(output.resolve()),
                  options=dict(variable=args.gpp, secondVariable=args.other or ''),
                  metadata=metadata or {})
    if args.radiation:
        config['options']['radiation'] = args.radiation
    if args.snapshot:
        config['snapshot'] = str(Path(args.snapshot).resolve())
    if args.archive_file:
        config['archiveFile'] = args.archive_file
    subprocess.run([args.node, str(Path(__file__).with_name('convert_tower.cjs')), '-'],
                   input=json.dumps(config), text=True, encoding='utf-8', check=True)

async def online(args, out):
    from fluxnet_shuttle import listall, download
    snapshot = args.snapshot or await listall(data_hubs=args.hubs, output_dir=str(out))
    with open(snapshot, encoding='utf-8-sig', newline='') as stream:
        catalog = list(csv.DictReader(stream))
    for site in args.sites:
        rows = [r for r in catalog if r.get('site_id') == site and (not args.hubs or r.get('data_hub','').lower() in args.hubs)]
        if len(rows) != 1:
            raise ValueError(f'{site}: expected one catalog record, found {len(rows)}. Narrow --hubs or supply a filtered --snapshot.')
        row = rows[0]
        # Pass a single verified catalog record to avoid duplicate-site downloads from other hubs.
        selected_snapshot = out / f'{site}_snapshot.csv'
        with selected_snapshot.open('w', encoding='utf-8', newline='') as f:
            writer=csv.DictWriter(f, fieldnames=list(row));writer.writeheader();writer.writerow(row)
        files = await download(site_ids=[site], snapshot_file=str(selected_snapshot), output_dir=str(out), user_info={})
        if len(files) != 1:
            raise ValueError(f'{site}: expected one downloaded archive; received {len(files)}.')
        metadata = dict(site=site, latitude=float(row['location_lat']), longitude=float(row['location_long']),
                        citation=row.get('product_citation',''), productId=row.get('product_id',''))
        convert(Path(files[0]), out/f'{site}_acclim.json', args, metadata)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--sites',nargs='+',help='Explicit site IDs; never downloads all sites.')
    mode.add_argument('--input',type=Path,help='Convert an already downloaded CSV or ZIP offline.')
    parser.add_argument('--snapshot',help='Existing Shuttle catalog CSV (also supplies metadata in offline mode).')
    parser.add_argument('--hubs',nargs='+',choices=['ameriflux','icos','tern'])
    parser.add_argument('--output-dir',type=Path,default=Path('outputs/shuttle'))
    parser.add_argument('--gpp',default='GPP_NT_VUT_REF')
    parser.add_argument('--other',default='GPP_DT_VUT_REF',help='Additional GPP field; use --other none to omit.')
    parser.add_argument('--radiation')
    parser.add_argument('--archive-file',help='Exact CSV path inside ZIP when multiple subdaily files exist.')
    parser.add_argument('--node',default='node')
    args=parser.parse_args();args.other='' if args.other=='none' else args.other
    args.output_dir.mkdir(parents=True,exist_ok=True)
    if args.input:
        convert(args.input,args.output_dir/(args.input.stem+'_acclim.json'),args)
    else:
        if any(not re.fullmatch(r'[A-Z]{2}-[A-Za-z0-9]+',s) for s in args.sites):parser.error('Invalid site ID.')
        asyncio.run(online(args,args.output_dir))

if __name__=='__main__':
    main()
