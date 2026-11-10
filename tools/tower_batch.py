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
args = p.parse_args()
ee.Initialize(project=args.project)
env = dict(os.environ, ACCLIM_EE_TOKEN=ee.data.get_persistent_credentials().token)
local_modules = Path(__file__).resolve().parents[2] / "work/ee-check/node_modules"
if local_modules.is_dir():
    env["NODE_PATH"] = os.pathsep.join(filter(None, (str(local_modules), env.get("NODE_PATH"))))
subprocess.run(["node", str(Path(__file__).with_suffix(".cjs")), str(args.archive.resolve()),
                str(args.output.resolve()), args.project, args.params], env=env, check=True)
