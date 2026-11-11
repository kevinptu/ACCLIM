"""Photon-weighted FourSAIL foliage and chlorophyll absorption for the GEE lookup table."""
import numpy as np
import prosail
from prosail.FourSAIL import foursail

# Node 0 approaches the horizon. No weight is assigned to sunlight below it.
MU = np.array([.0001, .01, .025, .05, .075, *np.arange(.1, 1.01, .05)])
BG = np.arange(5) / 10
WL = np.arange(400, 701)
PAR = slice(0, len(WL))


def optics(cab=40):
    pigments = prosail.spectral_lib.prospectd
    _, rho, tau = prosail.run_prospect(1.5, cab, 8, 0, .015, .009, prospect_version='D')
    terms = np.array([cab * pigments.kab, 8 * pigments.kcar, .015 * pigments.kw,
                      .009 * pigments.km])[:, PAR]
    total = terms.sum(axis=0)
    shares = np.divide(terms, total, out=np.zeros_like(terms), where=total > 0)
    return rho[PAR], tau[PAR], shares


def photon_weights(spectrum):
    # Spectral irradiance is energy based; lambda/(hc) converts to photons.
    weights = spectrum[PAR] * WL
    weights[[0, -1]] *= .5
    return weights / weights.sum()


DIRECT = photon_weights(prosail.spectral_lib.light.es)
DIFFUSE = photon_weights(prosail.spectral_lib.light.ed)


def fluxes(lai, mu, background, leaf):
    rho, tau, shares = leaf
    solar = np.degrees(np.arccos(np.clip(mu, .0001, 1)))
    background = np.broadcast_to(background, rho.shape).copy()
    t = foursail(rho, tau, -.35, -.15, 1, lai, .05, solar, 0, 0, background)
    # Total downward flux at the substrate includes canopy-background reflections.
    denominator = 1 - background * t[3]
    ground_direct = (1 - background) * (t[0] + t[6]) / denominator
    ground_diffuse = (1 - background) * t[4] / denominator
    direct = 1 - t[13] - ground_direct
    diffuse = 1 - t[12] - ground_diffuse
    return direct, diffuse, shares


def integrated(lai, mu, background, leaf):
    beam, sky, shares = fluxes(lai, mu, background, leaf)
    return [DIRECT @ beam, DIRECT @ (beam * shares[0]),
            DIFFUSE @ sky, DIFFUSE @ (sky * shares[0])]


def table():
    leaf = optics()
    direct = []
    diffuse = []
    for lai in np.arange(81) / 10:
        for bg in BG:
            for mu in MU:
                result = integrated(lai, mu, bg, leaf)
                direct.extend(result[:2])
            diffuse.extend(result[2:])
    def encode(values):
        return np.rint(np.array(values) * 100000).astype(int).tolist()
    return dict(mu=MU.round(6).tolist(), direct=encode(direct), diffuse=encode(diffuse))
