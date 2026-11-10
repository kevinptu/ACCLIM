"""Joint red/NIR cover sensitivity to one common soil pair per area and date."""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import prosail
from prosail.FourSAIL import foursail
from acclim_rt import read_pixels, save


def fit_soils(pixels, lai, soils):
    _, rho, tau = prosail.run_prospect(1.5, 40, 8, 0, .015, .009, prospect_version="D")
    idx = np.r_[220:271, 441:477]
    # FourSAIL has independent spectral channels. Batch all band-constant soil
    # levels into one call per geometry; keep PROSPECT's 1-nm spectral detail.
    levels, lookup = np.unique(soils, return_inverse=True)
    lookup = lookup.reshape(soils.shape)
    leaf_r, leaf_t = (np.tile(arr[idx], len(levels)) for arr in (rho, tau))
    ground = np.repeat(levels, len(idx))
    cover, errors, dense_errors, patch_means = [], [], [], []
    for red, nir, _, _, sun, view, az in pixels:
        sdr = foursail(leaf_r, leaf_t, -.35, -.15, 1, lai, .05, sun, view, az, ground)[17]
        spec = sdr.reshape(len(levels), -1)
        patch = np.column_stack((spec[:, :51].mean(1)[lookup[:, 0]],
                                 spec[:, 51:].mean(1)[lookup[:, 1]]))
        delta = patch - soils
        fraction = np.clip(np.sum(([red, nir] - soils)*delta, axis=1) /
                           np.sum(delta**2, axis=1), 0, 1)
        cover.append(fraction)
        errors.append((soils + fraction[:, None]*delta - [red, nir])**2)
        dense_errors.append((soils + np.maximum(fraction, .8)[:, None]*delta - [red, nir])**2)
        patch_means.append(patch)
    return np.array(cover), np.sqrt(np.mean(errors, axis=(0, 2))), np.mean(patch_means, axis=0), np.sqrt(np.mean(dense_errors, axis=(0, 2)))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("csv", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument("--date")
    args = p.parse_args()
    day, pixels = read_pixels(args.csv, args.date)
    args.out.mkdir(exist_ok=True, parents=True)
    soils = np.array([(r/100, n/100) for r in range(1, 41) for n in range(r, 51)])
    lib = prosail.spectral_lib.soil
    templates = []
    for brightness in np.arange(.5, 1.51, .05):
        for dry in np.linspace(0, 1, 21):
            spectrum = brightness*(dry*lib.rsoil1 + (1-dry)*lib.rsoil2)
            templates.append([spectrum[220:271].mean(), spectrum[441:477].mean()])
    # Grid cells within one grid spacing of the PROSAIL wet/dry soil family.
    family = np.min(np.linalg.norm(soils[:, None] - templates, axis=2), axis=1) <= .01
    fig, axes = plt.subplots(1, 4, figsize=(16, 4))
    summary = []
    for ax, lai in zip(axes, [1, 3, 5, 7]):
        cover, rmse, patches, dense_rmse = fit_soils(pixels, lai, soils)
        save(args.out/f"soils_LAI{lai}.csv", "soil_red soil_nir patch_red patch_nir joint_rmse in_template_family",
             np.column_stack((soils, patches, rmse, family)))
        for name, domain in [('broad', np.ones(len(soils), bool)), ('template', family)]:
            best = np.argmin(np.where(domain, rmse, np.inf))
            for tol in [.001, .002, .005]:
                near = domain & (rmse <= rmse[best]+tol)
                lo, hi = cover[:, near].min(1), cover[:, near].max(1)
                summary.append((lai, name, tol, *soils[best], rmse[best], near.sum(),
                                *soils[near].min(0), *soils[near].max(0),
                                np.median(cover[:, best]), np.median(hi-lo),
                                np.percentile(hi-lo, 90), np.mean(lo >= .8), np.min(dense_rmse[domain]),
                                np.ptp(patches[near, 0]), np.ptp(patches[near, 1])))
                if tol == .002:
                    save(args.out/f"pixels_{name}_LAI{lai}.csv",
                         "x y red nir cover_best cover_min cover_max",
                         np.column_stack((pixels[:, [2,3,0,1]], cover[:,best], lo, hi)))
                    if name == 'broad':
                        order = np.argsort(cover[:, best])
                        ax.fill_between(range(len(lo)), lo[order], hi[order], color='#99c9d3',
                                        label='Near-best soil pairs')
                        ax.plot(cover[order, best], color='#145563', label='Best fit')
                        ax.set(title=f'Patch LAI = {lai}', ylim=(0,1), xlabel='Pixels ranked by fitted cover')
    save(args.out/'summary.csv', 'LAI soil_domain rmse_tolerance best_soil_red best_soil_nir best_joint_rmse n_soils near_red_min near_nir_min near_red_max near_nir_max median_best_cover median_cover_width p90_cover_width fraction_min_cover_ge_0.8 best_rmse_cover_ge_0.8 patch_red_span patch_nir_span', summary)
    axes[0].set_ylabel('Fractional cover')
    axes[-1].legend()
    fig.suptitle(f'{args.csv.stem}, {day}: cover range within 0.002 BRF of minimum joint RMSE')
    fig.tight_layout()
    fig.savefig(args.out/'cover_sensitivity.png', dpi=180)
    print(day, len(pixels), 'pixels;',len(soils),'soil pairs;',sum(family),'near soil templates')


if __name__ == '__main__':
    main()
