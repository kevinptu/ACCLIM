"""Export one MODIS Aqua composite box with per-pixel observation day and grid IDs."""
import argparse
import csv
from datetime import date, timedelta
from pathlib import Path
import ee


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("project")
    p.add_argument("latitude", type=float)
    p.add_argument("longitude", type=float)
    p.add_argument("composite_start", type=date.fromisoformat)
    p.add_argument("output", type=Path)
    p.add_argument("--km", type=float, default=10)
    a = p.parse_args()
    if not 0 < a.km <= 20:
        p.error("Use a box side in (0,20] km for this interactive sampler.")
    ee.Initialize(project=a.project)
    collection = ee.ImageCollection("MODIS/061/MYD09A1").filterDate(
        a.composite_start.isoformat(), (a.composite_start + timedelta(days=1)).isoformat())
    if not collection.size().getInfo():
        p.error("No composite starts on that date.")
    image = ee.Image(collection.first())
    proj = image.select("sur_refl_b01").projection()
    point = ee.Geometry.Point([a.longitude, a.latitude])
    utm = ee.Projection("EPSG:" + str((32600 if a.latitude >= 0 else 32700) + min(60, int((a.longitude + 180) // 6) + 1)))
    region = point.buffer(a.km * 500, ee.ErrorMargin(1, "projected"), utm).bounds(ee.ErrorMargin(1, "projected"), utm)
    state, qc = image.select("StateQA"), image.select("QA")
    clear = (state.bitwiseAnd(7).eq(0)
             .And(state.rightShift(3).bitwiseAnd(7).eq(1))
             .And(state.bitwiseAnd((3 << 8) | (1 << 10) | (1 << 12) | (1 << 13) | (1 << 15)).eq(0))
             .And(qc.bitwiseAnd(1023).eq(0)))
    bands = image.select(["sur_refl_b01", "sur_refl_b02"], ["red", "nir"]).multiply(.0001)
    samples = (bands.addBands(image.select(["DayOfYear", "SolarZenith", "ViewZenith", "RelativeAzimuth"]))
               .addBands(ee.Image.pixelCoordinates(proj)).updateMask(clear)
               .sample(region=region, projection=proj, geometries=False))
    rows = [f["properties"] for f in samples.getInfo()["features"]]
    if not rows:
        p.error("No clear land pixels in this box.")
    with a.output.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["date", "x", "y", "red", "nir", "solar_zenith", "view_zenith", "relative_azimuth"])
        for r in rows:
            day = date(a.composite_start.year, 1, 1) + timedelta(days=r["DayOfYear"] - 1)
            writer.writerow([day, r["x"], r["y"], r["red"], r["nir"],
                             r["SolarZenith"]*.01, r["ViewZenith"]*.01, abs(r["RelativeAzimuth"]*.01)])


if __name__ == "__main__":
    main()

