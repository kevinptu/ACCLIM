"""Retrieve whole-pixel LAI and nadir foliage cover from directional MODIS BRF."""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import prosail
from numba import njit
from prosail.FourSAIL import foursail
from acclim_rt import read_pixels, save


@njit
def band_mean(levels, start, end, rso, gap, rdd, tsd, tdo, tss, too):
    # FourSAIL rsot soil-interaction terms, averaged at 1 nm within each band.
    out = np.empty(len(levels))
    for i, soil in enumerate(levels):
        total = 0.
        for j in range(start, end):
            total += rso[j] + gap*soil + soil*((tss+tsd[j])*tdo[j] +
                (tsd[j]+tss*soil*rdd[j])*too)/(1-soil*rdd[j])
        out[i] = total/(end-start)
    return out


def profiles(pixels, soils, lais):
    _, rho, tau = prosail.run_prospect(1.5, 40, 8, 0, .015, .009, prospect_version="D")
    idx = np.r_[220:271, 441:477]
    rho, tau = rho[idx], tau[idx]
    levels, inverse = np.unique(soils, return_inverse=True)
    ir, inn = inverse.reshape(soils.shape).T
    covers, errors, predictions = [], [], []
    for red, nir, _, _, sun, view, az in pixels:
        # LAI zero is bare soil. Every candidate shares its soil across pixels.
        best = soils.copy()
        cost = ((best-[red,nir])**2).sum(axis=1)
        selected = np.zeros(len(soils))
        for lai in lais[1:]:
            terms = foursail(rho,tau,-.35,-.15,1,lai,.05,sun,view,az,np.zeros_like(rho))
            optical = (terms[9],terms[2],terms[3],terms[6],terms[8],terms[0],terms[1])
            r = band_mean(levels,0,51,*optical)[ir]
            n = band_mean(levels,51,len(idx),*optical)[inn]
            loss = (r-red)**2+(n-nir)**2
            better = loss < cost
            selected[better], cost[better] = lai, loss[better]
            best[better,0], best[better,1] = r[better], n[better]
        covers.append(selected)
        errors.append(cost)
        predictions.append(best)
    return np.array(covers), np.sqrt(np.mean(errors,axis=0)/2), np.array(predictions)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("csv",type=Path)
    p.add_argument("out",type=Path)
    p.add_argument("--date")
    args = p.parse_args()
    day, pixels = read_pixels(args.csv,args.date)
    args.out.mkdir(parents=True,exist_ok=True)
    soils = np.array([(r/100,n/100) for r in range(1,41) for n in range(r,51)])
    lib = prosail.spectral_lib.soil
    templates = []
    for brightness in np.arange(.5,1.51,.05):
        for dry in np.linspace(0,1,21):
            s = brightness*(dry*lib.rsoil1+(1-dry)*lib.rsoil2)
            templates.append([s[220:271].mean(),s[441:477].mean()])
    family = np.linalg.norm(soils[:,None]-templates,axis=2).min(axis=1)<=.01
    lais = np.arange(81)/10
    lai, rmse, pred = profiles(pixels,soils,lais)
    cover = -np.expm1(-.5*lai)
    save(args.out/'soils.csv','red nir template joint_rmse',np.column_stack((soils,family,rmse)))
    fig, axes = plt.subplots(1,2,figsize=(10,4),sharey=True)
    results = []
    for ax,(name,domain) in zip(axes,[('broad',np.ones(len(soils),bool)),('template',family)]):
        best = np.argmin(np.where(domain,rmse,np.inf))
        for tol in [.001,.002,.005]:
            near = domain & (rmse <= rmse[best]+tol)
            lo, hi = cover[:,near].min(1), cover[:,near].max(1)
            results.append((day,len(pixels),name,tol,*soils[best],rmse[best],near.sum(),
                            np.median(lai[:,best]),np.median(cover[:,best]),np.median(hi-lo),
                            np.mean(lai[:,best]==0),np.mean(lai[:,best]==lais[-1])))
            if tol==.002:
                save(args.out/f'pixels_{name}.csv','x y red nir fitted_red fitted_nir lai nadir_cover cover_min cover_max',
                     np.column_stack((pixels[:,[2,3,0,1]],pred[:,best],lai[:,best],cover[:,best],lo,hi)))
                order=np.argsort(cover[:,best])
                ax.fill_between(np.arange(len(order)),lo[order],hi[order],color='#a7d3dc',label='Near-best soils')
                ax.plot(cover[order,best],color='#175463',label='Best soil')
                ax.set(title=name.capitalize()+' soil grid',xlabel='Pixels ranked by fitted cover',ylim=(0,1))
    axes[0].set_ylabel('Nadir foliage cover')
    axes[1].legend()
    fig.suptitle(f'{args.csv.stem}, {day} | Whole-pixel PROSAIL')
    fig.tight_layout()
    fig.savefig(args.out/'cover.png',dpi=170)
    save(args.out/'summary.csv','date pixels domain tolerance soil_red soil_nir joint_rmse near_soils median_lai median_cover median_width fraction_lai_zero fraction_lai_upper',results)


if __name__=='__main__':
    main()
