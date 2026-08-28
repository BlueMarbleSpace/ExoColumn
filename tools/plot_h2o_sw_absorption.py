#!/usr/bin/env python3
"""
plot_h2o_sw_absorption.py — ExoRT n68 H2O shortwave absorption, for explaining the
IHZ planetary-albedo offset vs Kopparapu+2013 (see project_hz_albedo_offset).

The hemispheric-zenith fix (sw_zenith_quad) closed ~1/3 of the ExoColumn-vs-
Kopparapu albedo gap; the residual ~0.015-0.02 is the near-IR H2O SHORTWAVE
absorption being weaker in ExoRT's spectroscopy than in Kopparapu's.  This script
extracts ExoRT's ACTUAL H2O absorption from its correlated-k tables and plots the
band-mean cross-section vs wavelength, against the G2V solar spectrum, annotated
with the provenance difference:

  ExoColumn (ExoRT n68) : HITRAN-2016 line list + MT_CKD 3.3 continuum, 39 bands
                          across 0.25-5 um.
  Kopparapu+2013        : HITEMP-2010 line list + BPS (Paynter & Ramaswamy 2011)
                          continuum, 38 solar intervals across 0.2-4.5 um.

FINDING (with the RADIS HITEMP overlay, tools/fetch_hitemp_h2o.py):
The two use COMPARABLE spectral resolution (~36-39 vs 38 SW intervals), so the
offset is NOT a band-count effect.  And — counter to the initial hypothesis — it is
NOT the line list either: HITEMP-2010 and HITRAN-2016 H2O agree to ~1% across the
near-IR SW bands at BOTH 300 K and 500 K (median HITEMP/HITRAN ~0.99-1.00; the n68
correlated-k band-mean lands on top of both).  The strong 300->500 K growth of the
window absorption (factors 5-40x at 1.6, 2.1, 3.8-4.7 um) is a TEMPERATURE effect
captured equally by both line lists.  So the residual albedo offset vs Kopparapu is
NOT H2O line spectroscopy.  The remaining spectroscopic candidate is the water-vapour
CONTINUUM (MT_CKD 3.3 here vs BPS / Paynter & Ramaswamy 2011 in Kopparapu, which they
adopted specifically because it adds dimer absorption in the windows); non-spectro-
scopic SW-RT / moist-adiabat-profile differences may also contribute.
"""
import os
import numpy as np
import netCDF4 as nc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ExoRT tree.  Override with the EXORT_ROOT environment variable (the same
# name the build uses in config.mk) if ExoRT lives elsewhere.
EXORT = os.environ.get('EXORT_ROOT', '/models/ExoRT')
KFILE = os.path.join(EXORT, 'data/kdist/n68h2o/hitran2016/'
                     'n68_8gpt_h2o_hitran16_Nnu1e4_c25_voigt_noplinth_q0_grrtm.nc')
CFILE = os.path.join(EXORT, 'data/continuum/KH2O_MTCKD3.3_SELF.FRGN_n68_ngauss.nc')
SFILE = os.path.join(EXORT, 'data/solar/G2V_SUN_n68.nc')
OUT   = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'h2o_sw_absorption.png')

# k-coeff grids (radgrid.F90: tgrid, log10pgrid, g_weight_8gpt)
TGRID = np.array([100,125,150,175,200,225,250,275,300,325,350,375,400,425,450,475,500.])
LOG10P = np.round(np.arange(-2.0, 4.0001, 0.1), 4)         # 61 pts, pgrid in hPa
GW8 = np.array([0.30192,0.27379,0.22012,0.14595,0.04712,0.00686,0.00363,0.00061])

# Named near-IR H2O bands (um) for annotation
H2O_BANDS = [(0.94,'0.94'), (1.13,'1.13'), (1.38,'1.38'),
             (1.87,'1.87'), (2.7,'2.7'), (3.2,'3.2')]


def band_mean_sigma(k4, iT, iP):
    """k4: (NT,NP,Ng,NB) -> band-mean sigma at (iT,iP) using the 8-pt g-weights."""
    return np.tensordot(GW8, k4[iT, iP, :, :], axes=([0], [0]))


def main():
    # --- line k-table ---
    kds = nc.Dataset(KFILE)
    k4 = np.array(kds['data'][:])                      # (T,P,g,band) cm2/molec
    # --- continuum (self) ---
    cds = nc.Dataset(CFILE)
    nu_lo = np.array(cds['NU_LOW'][:]); nu_hi = np.array(cds['NU_HIGH'][:])
    Tc = np.array(cds['TEMPERATURES'][:])
    kself = np.array(cds['KSELF'][:])                  # (Tc,band,g) cm2/molec
    # --- solar ---
    sds = nc.Dataset(SFILE)
    swl = np.array(sds['wav_low'][:]); swh = np.array(sds['wav_high'][:])
    sflux = np.array(sds['solarflux'][:])              # W/m2 per band

    # --- optional RADIS HITEMP/HITRAN band-averaged sigma (from fetch_hitemp_h2o.py) ---
    radis_npz = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             'h2o_radis_bands.npz')
    radis = np.load(radis_npz) if os.path.exists(radis_npz) else None

    nu_mid = 0.5 * (nu_lo + nu_hi)
    lam = 1e4 / nu_mid                                 # um band centre
    with np.errstate(divide='ignore'):
        dlam = np.abs(1e4 / np.maximum(nu_lo, 1e-6) - 1e4 / nu_hi)   # um band width

    iP = int(np.argmin(np.abs(LOG10P - 2.0)))          # 100 hPa
    iT3 = int(np.argmin(np.abs(TGRID - 300)))
    iT5 = int(np.argmin(np.abs(TGRID - 500)))
    sig300 = band_mean_sigma(k4, iT3, iP)
    sig500 = band_mean_sigma(k4, iT5, iP)
    # self-continuum band-mean at 300 K (per molecule), same units
    iTc3 = int(np.argmin(np.abs(Tc - 300)))
    cself300 = np.tensordot(GW8, kself[iTc3, :, :], axes=([0], [1]))

    # restrict to the shortwave 0.3-5 um
    sw = (lam > 0.3) & (lam < 5.0)
    order = np.argsort(lam[sw])
    L = lam[sw][order]
    s3 = sig300[sw][order]; s5 = sig500[sw][order]; cs = cself300[sw][order]
    # solar flux per um (for the context panel), on the same SW bands
    sden = (sflux / np.maximum(dlam, 1e-9))[sw][order]

    n_sw = int(sw.sum())

    # ---------------- figure ----------------
    fig, (axS, axK) = plt.subplots(2, 1, figsize=(7.5, 6.2), dpi=300,
                                   sharex=True, height_ratios=[1, 2])
    fig.patch.set_facecolor('white')

    # (a) solar spectrum
    axS.fill_between(L, sden, step='mid', color='gold', alpha=0.65, lw=0)
    axS.plot(L, sden, drawstyle='steps-mid', color='goldenrod', lw=0.8)
    axS.set_ylabel('G2V solar flux\n(W m$^{-2}$ µm$^{-1}$)')
    axS.set_facecolor('white')
    _ttl = ('Near-IR H$_2$O shortwave absorption: ExoColumn n68 vs Kopparapu HITEMP-2010'
            if radis is not None else
            'ExoRT n68 H$_2$O shortwave absorption (HITRAN-2016 + MT_CKD 3.3) '
            'vs the solar spectrum')
    axS.set_title(_ttl, fontsize=9.5)
    axS.text(0.015, 0.9, '(a) incident sunlight', transform=axS.transAxes,
             va='top', fontsize=8, fontweight='bold')

    # (b) H2O absorption cross-section: ExoColumn n68 vs (if available) HITEMP/HITRAN
    axK.plot(L, s3, drawstyle='steps-mid', color='C0', lw=1.7,
             label='ExoColumn n68 (HITRAN-2016), 300 K')
    ratio_txt = ''
    if radis is not None:
        he3 = radis['sig_hitemp_300'][sw][order] if 'sig_hitemp_300' in radis.files else None
        he5 = radis['sig_hitemp_500'][sw][order] if 'sig_hitemp_500' in radis.files else None
        hr3 = radis['sig_hitran_300'][sw][order] if 'sig_hitran_300' in radis.files else None
        hr5 = radis['sig_hitran_500'][sw][order] if 'sig_hitran_500' in radis.files else None
        # HITEMP/HITRAN agreement statistic (the key result: line list is NOT it)
        rr = []
        for a, b in [(he3, hr3), (he5, hr5)]:
            if a is not None and b is not None:
                m = np.isfinite(a) & np.isfinite(b) & (b > 0)
                rr += list(a[m] / b[m])
        if rr:
            ratio_txt = (f'HITEMP/HITRAN median = {np.median(rr):.2f} '
                         f'(300 & 500 K)')
        if he3 is not None:   # overlaps the n68 (HITRAN-2016) curve -> agreement
            axK.plot(L, he3, drawstyle='steps-mid', color='k', lw=1.2, ls='--',
                     label='HITEMP-2010 (Kopparapu line list), 300 K')
        if he5 is not None:
            axK.plot(L, he5, drawstyle='steps-mid', color='C3', lw=1.3,
                     label='HITEMP-2010, 500 K (moist greenhouse)')
    else:
        axK.plot(L, s5, drawstyle='steps-mid', color='C3', lw=1.2, alpha=0.85,
                 label='ExoColumn n68, 500 K (table ceiling)')
    axK.plot(L, cs, drawstyle='steps-mid', color='C2', lw=1.1, ls=':',
             label='MT_CKD self-continuum (ExoColumn), 300 K')
    axK.set_yscale('log')
    axK.set_ylim(1e-28, 3e-20)
    axK.set_xlim(0.3, 5.0)
    axK.set_xlabel('Wavelength (µm)')
    axK.set_ylabel('Band-mean H$_2$O cross-section\n$\\sigma$ (cm$^2$ molecule$^{-1}$), '
                   '$p$=100 hPa')
    axK.set_facecolor('white')
    axK.text(0.015, 0.96, '(b) ExoRT n68 H$_2$O absorption', transform=axK.transAxes,
             va='top', fontsize=8, fontweight='bold')
    # band-centre labels
    for lc, name in H2O_BANDS:
        axK.axvline(lc, color='0.8', lw=0.6, zorder=0)
        axK.text(lc, 4e-20, name, rotation=90, fontsize=6.5, color='0.45',
                 ha='right', va='top')
    axK.legend(fontsize=7.5, loc='lower right', framealpha=0.92)

    if radis is not None:
        note = ('ExoColumn n68 (HITRAN-2016) and Kopparapu\'s HITEMP-2010 H$_2$O lines\n'
                f'coincide in the near-IR SW ({ratio_txt}): the LINE LIST is NOT the\n'
                'albedo offset.  Both grow strongly 300→500 K (windows fill 5–40×).\n'
                'Remaining spectroscopic difference = continuum: MT_CKD 3.3 (ours)\n'
                'vs BPS / Paynter-Ramaswamy 2011 (Kopparapu).')
    else:
        note = ('Kopparapu+2013: HITEMP-2010 + BPS continuum, 38 SW intervals\n'
                f'ExoRT n68: HITRAN-2016 + MT_CKD 3.3, {n_sw} bands in 0.3–5 µm')
    axK.text(0.015, 0.045, note, transform=axK.transAxes, va='bottom', ha='left',
             fontsize=6.6, family='monospace',
             bbox=dict(boxstyle='round', fc='white', ec='0.7', alpha=0.92))

    fig.tight_layout()
    fig.savefig(OUT, dpi=150, facecolor='white')
    print(f'  n68 SW bands (0.3-5um): {n_sw}')
    print(f'Saved: {OUT}')


if __name__ == '__main__':
    main()
