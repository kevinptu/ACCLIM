"""Run one ACCLIM comparison per FLUXMET site-year in a ZIP (same settings across sites)."""
import argparse
import os
from pathlib import Path
import subprocess
import ee

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("archive", type=Path)
p.add_argument("--output", type=Path, default=Path("five_site_comparison.csv"))
p.add_argument("--project", default="ee-landflux")
p.add_argument("--params", default="{}", help="JSON object with ACCLIM parameter overrides")
p.add_argument("--cache-dir", type=Path, help="Persistent input directory (default: beside the ZIP)")
p.add_argument("--refresh-inputs", action="store_true", help="Fetch and replace the saved samples")
args = p.parse_args()
cache_dir = args.cache_dir or args.archive.with_name(args.archive.stem + "_inputs")
ee.Initialize(project=args.project)
env = dict(os.environ, ACCLIM_EE_TOKEN=ee.data.get_persistent_credentials().token)
local_modules = Path(__file__).resolve().parents[2] / "work/ee-check/node_modules"
if local_modules.is_dir():
    env["NODE_PATH"] = os.pathsep.join(filter(None, (str(local_modules), env.get("NODE_PATH"))))
subprocess.run(["node", str(Path(__file__).with_suffix(".cjs")), str(args.archive.resolve()),
                str(args.output.resolve()), args.project, args.params, str(cache_dir.resolve()),
                "refresh" if args.refresh_inputs else "reuse"], env=env, check=True)
