#!/usr/bin/env python3
"""
make_socrates_handoff.py — export the exact outer-HZ maximum-greenhouse column
used in the left panel of the manuscript's outer-edge radiative-transfer
benchmark figure (reference/max_greenhouse/lbl_olr_benchmark_ohz.pdf) as plain text,
so an external radiation model (e.g. SOCRATES) can be run on the identical
atmosphere.

Source of truth : reference/max_greenhouse/exocol_maxgh_8.87bar.nc  (PVER=200)
Band comparison : lbl_olr_co2_maxgh.npz, clima_band_olr_maxgh.txt

Run:  python3 reference/max_greenhouse/make_socrates_handoff.py
"""
import os
import numpy as np
import netCDF4 as nc

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'socrates_handoff')
os.makedirs(OUT, exist_ok=True)

MW = dict(H2O=18.01528, CO2=44.0095, N2=28.0134)     # g/mol
G = 9.80616                                          # m/s2 (ExoRT exo_g, Earth)

# ------------------------------------------------- LBL / Clima reference totals
def _band_integrate(wn, oc, edges):
    out = np.full(len(edges) - 1, np.nan)
    for i in range(len(out)):
        m = (wn >= edges[i]) & (wn < edges[i + 1])
        if m.sum() > 2:
            out[i] = np.trapezoid(oc[m], wn[m])
    return out


def _refs():
    """Totals for the comparison table, read from the cached benchmark products
    rather than hardcoded, so they cannot go stale when the LBL is regenerated."""
    d = np.load(os.path.join(HERE, 'lbl_olr_co2_maxgh.npz'))
    wn = d['wn'].astype(float)
    r = {'lbl_chi': float(np.trapezoid(d['olr_nu_full'].astype(float), wn)),
         'lbl_lor': float(np.trapezoid(d['olr_nu_full_nochi'].astype(float), wn))}
    c = np.loadtxt(os.path.join(HERE, 'clima_band_olr_maxgh.txt'))
    r['clima13'], r['clima16'] = float(c[:, 2].sum()), float(c[:, 3].sum())
    return r


REF = _refs()

# ---------------------------------------------------------------- read column
with nc.Dataset(os.path.join(HERE, 'exocol_maxgh_8.87bar.nc')) as ds:
    tmid = np.array(ds['tmid'][:]);  tint = np.array(ds['tint'][:])
    pmid = np.array(ds['pmid'][:]);  pint = np.array(ds['pint'][:])
    zint = np.array(ds['zint'][:])
    q    = np.array(ds['h2ommr'][:])         # specific humidity, kg/kg MOIST air
    wco2 = np.array(ds['co2mmr'][:])         # kg/kg DRY air
    wn2  = np.array(ds['n2mmr'][:])          # kg/kg DRY air
    ts   = float(ds['ts'][:]);  ps = float(ds['ps'][:])
    mwdry = float(ds['mw'][:]); cpdry = float(ds['cp'][:])
    lwup = np.array(ds['LWUP'][:]); lwdn = np.array(ds['LWDN'][:])
    swup = np.array(ds['SWUP'][:]); swdn = np.array(ds['SWDN'][:])
    band_edges = np.array(ds['wavenum_edge'][:])
    band_olr   = np.array(ds['band_lwup_toa'][:])

nlay = len(tmid)

# mass fractions of TOTAL (moist) air -> mole fractions
f = {'H2O': q, 'CO2': (1.0 - q) * wco2, 'N2': (1.0 - q) * wn2}
moles = {k: v / MW[k] for k, v in f.items()}
tot = sum(moles.values())
x = {k: v / tot for k, v in moles.items()}
mw_moist = 1000.0 / tot                                      # g/mol, per layer

# layer column amounts (molec/cm2) as a cross-check for the recipient
pdel = np.diff(pint)                                          # Pa, TOA->sfc
NA = 6.02214076e23
Ncol = {k: (f[k] * pdel / G) / (MW[k] * 1e-3) * NA * 1e-4 for k in f}

# ------------------------------------------------------------------ profile
hdr = f"""# ============================================================================
# ExoColumn outer-HZ MAXIMUM-GREENHOUSE column
# The exact atmosphere used in the left-hand panel of the outer-edge
# radiative-transfer benchmark figure of the ExoColumn paper (ExoRT n68equiv vs
# Clima vs line-by-line).  Clear sky, no aerosol, no clouds.
#
# Source file : reference/max_greenhouse/exocol_maxgh_8.87bar.nc
# Generated   : make_socrates_handoff.py
#
# This is the S_eff-minimum point of the outer-HZ sweep (Kopparapu et al. 2013
# Fig. 5 analogue): an inverse (prescribed-profile) calculation, NOT an RCE
# solution -- the profile below is an input, and only the fluxes are diagnosed.
#
# SETUP
#   Surface temperature   Ts   = {ts:.4f} K            (fixed)
#   Surface pressure      ps   = {ps:.4f} Pa  = {ps/1e5:.6f} bar
#                              = {x['CO2'][-1]*ps/1e5:.4f} bar CO2 + {x['N2'][-1]*ps/1e5:.4f} bar N2 (+ trace H2O)
#   Gravity               g    = {G} m/s^2          (Earth)
#   Model top             p    = {pint[0]:.4g} Pa
#   Layers                     = {nlay} (index 1 = TOA, {nlay} = surface)
#   Dry-air mean mol. wt.      = {mwdry:.5f} g/mol
#   Dry-air cp                 = {cpdry:.3f} J/kg/K
#   Stratosphere               = isothermal at 154 K (CO2-condensation cold trap)
#   Troposphere                = H2O-saturated (RH=1) moist adiabat from Ts,
#                                pinned to the CO2 saturation curve wherever the
#                                ascent supersaturates in CO2 (Kasting 1991),
#                                with cp_CO2 evaluated at the local temperature.
#                                CO2 is kept WELL MIXED (CO2 clouds neglected,
#                                as in Clima).
#   Composition                = CO2 + N2 + trace H2O only.  No O2, O3, CH4.
#
# SHORTWAVE (only if you also want the S_eff comparison; the figure is LW-only)
#   Stellar spectrum      = present-day Sun (G2V), S0 = 1360.0 W/m^2 at 1 AU
#   Surface albedo        = 0.32 (Lambertian, spectrally grey; direct & diffuse)
#   Zenith-angle handling = 6-point Gauss-Legendre hemispheric average, with the
#                           insolation normalised to S0/4.
#
# EXOCOLUMN / ExoRT n68equiv RESULT FOR THIS COLUMN (what SOCRATES should match)
#   OLR   (TOA LW up)     = {lwup[0]:.3f}  W/m^2
#   TOA SW down           = {swdn[0]:.3f}  W/m^2
#   TOA SW up             = {swup[0]:.3f}  W/m^2   -> planetary albedo = {swup[0]/swdn[0]:.4f}
#   Surface LW down       = {lwdn[-1]:.3f}  W/m^2
#   Surface SW down       = {swdn[-1]:.3f}  W/m^2
#   S_eff = OLR / (SWdn-SWup) = {lwup[0]/(swdn[0]-swup[0]):.4f}
#   For reference, other models on this same column (0-2000 cm^-1 OLR):
#     Clima, Kopparapu-2013 k-coeffs   {REF['clima13']:.1f} W/m^2
#     Clima, Wolf HITRAN-2016 k-coeffs {REF['clima16']:.1f} W/m^2
#     ExoRT n68equiv (this work)       {lwup[0]:.1f} W/m^2
#     LBL, PH89 sub-Lorentzian chi     {REF['lbl_chi']:.1f} W/m^2   (transparent far-wing bound)
#     LBL, pure-Lorentz wings          {REF['lbl_lor']:.1f} W/m^2   (opaque far-wing bound)
#
# NOTE ON MIXING RATIOS: vmr = mole fraction of TOTAL (moist) air, so the three
# columns sum to 1 at every layer.  H2O is a trace species here (<= 2.2e-4 by
# mass), so CO2/N2 are effectively constant with height.
# ============================================================================
"""

# --- levels (interfaces) ---
lines = [hdr, "\n# TABLE 1 -- LEVELS (interfaces), index 1 = model top, %d = surface\n" % (nlay + 1)]
lines.append("# %5s %18s %14s %14s\n" % ("i", "p_level[Pa]", "T_level[K]", "z_level[m]"))
for i in range(nlay + 1):
    lines.append("  %5d %18.8e %14.6f %14.4f\n" % (i + 1, pint[i], tint[i], zint[i]))

# --- layers (midpoints) ---
lines.append("\n\n# TABLE 2 -- LAYERS (midpoints), index 1 = top layer, %d = surface layer\n" % nlay)
lines.append("# %5s %18s %14s %16s %16s %16s %14s\n" %
             ("k", "p_layer[Pa]", "T_layer[K]", "vmr_CO2", "vmr_N2", "vmr_H2O", "mw[g/mol]"))
for k in range(nlay):
    lines.append("  %5d %18.8e %14.6f %16.9e %16.9e %16.9e %14.6f\n" %
                 (k + 1, pmid[k], tmid[k], x['CO2'][k], x['N2'][k], x['H2O'][k], mw_moist[k]))

# --- column amounts cross-check ---
lines.append("\n\n# TABLE 3 -- LAYER COLUMN AMOUNTS [molec/cm^2] (cross-check only;\n"
             "#            derived from Table 2 with g = %.5f m/s^2)\n" % G)
lines.append("# %5s %18s %18s %18s\n" % ("k", "N_CO2", "N_N2", "N_H2O"))
for k in range(nlay):
    lines.append("  %5d %18.8e %18.8e %18.8e\n" %
                 (k + 1, Ncol['CO2'][k], Ncol['N2'][k], Ncol['H2O'][k]))
lines.append("# totals: CO2 %.6e   N2 %.6e   H2O %.6e  molec/cm^2\n" %
             (Ncol['CO2'].sum(), Ncol['N2'].sum(), Ncol['H2O'].sum()))

with open(os.path.join(OUT, 'ohz_maxgh_profile.txt'), 'w') as fh:
    fh.writelines(lines)

# --------------------------------------------------------------- CSV (layers)
with open(os.path.join(OUT, 'ohz_maxgh_profile.csv'), 'w') as fh:
    fh.write("# ExoColumn outer-HZ max-greenhouse column; Ts=%.2f K, ps=%.4f Pa, g=%.5f m/s2\n"
             % (ts, ps, G))
    fh.write("# k=1 is the top layer, k=%d the surface layer; vmr = mole fraction of moist air\n" % nlay)
    fh.write("k,p_layer_Pa,T_layer_K,p_level_top_Pa,p_level_bot_Pa,"
             "T_level_top_K,T_level_bot_K,z_level_top_m,z_level_bot_m,vmr_CO2,vmr_N2,vmr_H2O\n")
    for k in range(nlay):
        fh.write("%d,%.8e,%.6f,%.8e,%.8e,%.6f,%.6f,%.4f,%.4f,%.9e,%.9e,%.9e\n" %
                 (k + 1, pmid[k], tmid[k], pint[k], pint[k+1], tint[k], tint[k+1],
                  zint[k], zint[k+1], x['CO2'][k], x['N2'][k], x['H2O'][k]))

# ------------------------------------------------- band-resolved OLR (bonus)
def band_integrate(wn, oc, edges):
    out = np.full(len(edges) - 1, np.nan)
    for i in range(len(out)):
        m = (wn >= edges[i]) & (wn < edges[i + 1])
        if m.sum() > 2:
            out[i] = np.trapezoid(oc[m], wn[m])
    return out

d = np.load(os.path.join(HERE, 'lbl_olr_co2_maxgh.npz'))
wn = d['wn'].astype(float)
lbl_chi = band_integrate(wn, d['olr_nu_full'].astype(float), band_edges)
lbl_lor = band_integrate(wn, d['olr_nu_full_nochi'].astype(float), band_edges)

with open(os.path.join(OUT, 'ohz_maxgh_band_olr.txt'), 'w') as fh:
    fh.write("# Band-resolved TOA OLR [W/m^2 per band] for the column in\n"
             "# ohz_maxgh_profile.txt, on the ExoRT n68equiv band grid.\n"
             "# This is exactly what is plotted in the right panel of the benchmark\n"
             "# figure (there divided by the band width to give W/m^2/cm^-1).\n"
             "#   exort_n68   : ExoRT n68equiv, this work        (total %.2f W/m^2)\n"
             "#   lbl_ph89chi : line-by-line, Perrin-Hartmann 1989 sub-Lorentzian\n"
             "#                 chi-factor on the CO2 wings      (total %.2f W/m^2)\n"
             "#   lbl_lorentz : line-by-line, pure-Lorentz wings (total %.2f W/m^2)\n"
             "# 'nan' in the two LBL columns = band lies outside the 10-2000 cm^-1\n"
             "# range of the line-by-line calculation (the ExoRT column is complete).\n"
             "# The Clima band OLR is on its own 55-interval grid: see\n"
             "# clima_band_olr_maxgh.txt (copied alongside this file).\n"
             "# columns: wn_lo[cm-1]  wn_hi[cm-1]  exort_n68  lbl_ph89chi  lbl_lorentz\n"
             % (np.nansum(band_olr), np.nansum(lbl_chi), np.nansum(lbl_lor)))
    for i in range(len(band_olr)):
        fh.write("%10.2f %10.2f %14.6e %14.6e %14.6e\n" %
                 (band_edges[i], band_edges[i+1], band_olr[i], lbl_chi[i], lbl_lor[i]))

# ship the Clima band OLR alongside so the package is self-contained
import shutil
shutil.copy(os.path.join(HERE, 'clima_band_olr_maxgh.txt'), OUT)

print("wrote:", OUT)
for fn in sorted(os.listdir(OUT)):
    print("  ", fn, os.path.getsize(os.path.join(OUT, fn)), "bytes")
print("\nTs=%.2f K  ps=%.4f bar  = %.4f bar CO2 + %.4f bar N2"
      % (ts, ps/1e5, x['CO2'][-1]*ps/1e5, x['N2'][-1]*ps/1e5))
print("OLR(ExoRT n68) = %.3f W/m2 ; band sum = %.3f" % (lwup[0], np.nansum(band_olr)))
print("vmr_CO2 sfc=%.6f top=%.6f ; vmr_N2 sfc=%.6f ; vmr_H2O sfc=%.4e"
      % (x['CO2'][-1], x['CO2'][0], x['N2'][-1], x['H2O'][-1]))
