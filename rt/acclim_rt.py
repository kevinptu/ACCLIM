"""ACCLIM-RT: test a common soil spectrum and variable canopy patch cover."""
import argparse
import csv
from collections import Counter
from itertools import product
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import prosail


def bands(spectrum):
    # Equal-weight MODIS red (620-670 nm) and NIR (841-876 nm), starting at 400 nm.
    return np.array([spectrum[220:271].mean(), spectrum[441:477].mean()])


def spectra(angles):
    """One Lambertian spectrum for exposed soil and the ground beneath each patch."""
    for brightness, dry in product([0.5, 0.75, 1, 1.25, 1.5], [0, 0.5, 1]):
        soil = brightness * (dry * prosail.spectral_lib.soil.rsoil1
                             + (1 - dry) * prosail.spectral_lib.soil.rsoil2)
        for lai, cab, leaf_angle in product([1, 3, 6], [15, 40, 65], [30, 60]):
            canopy = prosail.run_prosail(1.5, cab, 8, 0, 0.015, 0.009,
                                         lai, leaf_angle, 0.05, *angles,
                                         rsoil0=soil, factor="SDR", prospect_version="D")
            yield (brightness, dry, lai, cab, leaf_angle, *bands(soil), *bands(canopy))


def read_pixels(path, day):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        fields = ("red", "nir", "x", "y", "solar_zenith", "view_zenith", "relative_azimuth")
        if not set(("date", *fields)).issubset(reader.fieldnames or []):
            raise ValueError("Expected one area with dated grid pixels and observation angles.")
        rows = [(row["date"], [float(row[k]) for k in fields]) for row in reader]
    valid = [(d, r) for d, r in rows if np.isfinite(r).all() and
             0 <= r[0] <= 1 and 0 <= r[1] <= 1]
    if not valid:
        raise ValueError("No valid pixels.")
    counts = Counter(d for d, _ in valid)
    day = day or min(counts, key=lambda d: (-counts[d], d))
    pixels = np.array([r for d, r in valid if d == day])
    if len(pixels) < 4 or len(np.unique(pixels[:, 2:4], axis=0)) != len(pixels):
        raise ValueError("Need at least four distinct pixels and one record per pixel on the requested date.")
    return day, pixels


def predict(red, model):
    sr, sn, cr, cn = model[-4:]
    cover = np.clip((red - sr) / (cr - sr), 0, 1)
    return cover, sn + cover * (cn - sn)


def fit(train, library):
    def loss(model):
        cover, nir = predict(train[:, 0], model)
        red = model[-4] + cover * (model[-2] - model[-4])
        return np.mean((red - train[:, 0]) ** 2 + (nir - train[:, 1]) ** 2)
    return min((m for m in library if abs(m[-2] - m[-4]) > 1e-6), key=loss)


def save(path, header, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(header.split())
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--angles", nargs=3, type=float, default=None,
                        metavar=("SOLAR_ZENITH", "VIEW_ZENITH", "REL_AZIMUTH"))
    parser.add_argument("--date", help="Actual observation date (default: day with most valid pixels)")
    args = parser.parse_args()
    day, pairs = read_pixels(args.csv, args.date)
    if args.angles is None:
        args.angles = np.median(pairs[:, 4:7], axis=0)
    sun, view, azimuth = args.angles
    if not (0 <= sun < 90 and 0 <= view < 90 and 0 <= azimuth <= 360):
        parser.error("Zenith angles must be in [0,90) and relative azimuth in [0,360].")
    # Leave out contiguous quadrants, keeping all training pixels outside the test quadrant.
    mid_x, mid_y = np.median(pairs[:, 2:4], axis=0)
    groups = (pairs[:, 2] > mid_x).astype(int) + 2 * (pairs[:, 3] > mid_y)
    if len(set(groups)) < 2:
        parser.error("Need spatial variation along the pixel grid.")
    library = list(spectra(args.angles))
    args.output.mkdir(exist_ok=True, parents=True)
    fields = "brightness dry_fraction patch_LAI Cab leaf_angle soil_red soil_nir patch_red patch_nir"
    save(args.output / "library.csv", fields + " solar_zenith view_zenith relative_azimuth",
         [(*m, *args.angles) for m in library])
    scores, predictions = [], []
    squared_error = baseline_error = 0
    for group in sorted(set(groups)):
        train, test = pairs[groups != group], pairs[groups == group]
        model = fit(train, library)
        cover, nir = predict(test[:, 0], model)
        err = nir - test[:, 1]
        baseline = train[:, 1].mean() - test[:, 1]
        squared_error += np.sum(err**2)
        baseline_error += np.sum(baseline**2)
        scores.append((group, len(test), np.mean(err**2)**0.5,
                       np.mean(baseline**2)**0.5, *model))
        predictions.extend(((day, int(x), int(y), int(group), r, n, f, pred)
                           for (r, n, x, y, *_), f, pred in zip(test, cover, nir)))
    save(args.output / "folds.csv", "fold n rmse_nir rmse_constant " + fields, scores)
    save(args.output / "held_out.csv", "date x y fold red nir cover_from_red predicted_nir", predictions)
    model = fit(pairs, library)
    save(args.output / "full_fit.csv", fields, [model])
    fig, ax = plt.subplots()
    ax.scatter(*pairs[:, :2].T, s=3, alpha=.18)
    ax.plot(model[-4::2], model[-3::2], "ko-", label="Full-data fit: soil to canopy patch")
    ax.set(xlabel="Red BRF", ylabel="NIR BRF", title=f"{args.csv.stem}, {day}")
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.output / "fit.png", dpi=160)
    plt.close(fig)
    print(f"{args.csv.stem} {day}: {len(pairs)} pairs, held-out NIR RMSE "
          f"{np.sqrt(squared_error/len(pairs)):.4f}, constant "
          f"{np.sqrt(baseline_error/len(pairs)):.4f}")


if __name__ == "__main__":
    main()

