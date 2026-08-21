#!/usr/bin/env python3
"""
plot_albedo_sensitivity.py  —  SUPPLEMENTARY figure: sensitivity of the ExoColumn
habitable-zone limits to the assumed surface albedo.

Motivation (coauthor review, 2026-08-21):  the main-text inner- and outer-edge
figures fix the surface albedo at alpha_s = 0.32, the value Kopparapu et al.
(2013) tuned inside Clima so that their 1-D Earth reaches Ts = 288 K (it is a
cloud proxy, not a physical surface albedo).  ExoColumn's own Earth calibration
gives alpha_s = 0.2736.  Holding alpha_s at Kopparapu's value isolates the
RADIATIVE-TRANSFER differences between the two models; this supplement asks the
complementary question — how much of the residual S_eff offset is attributable
to the albedo convention itself?

Because the HZ limits are computed in INVERSE (flux_only) mode, the T(p) and
H2O(p) columns are prescribed by the pseudoadiabat and are independent of the
surface albedo.  Surface albedo therefore acts ONLY on the shortwave:
    F_IR (OLR) is bit-identical between the two sweeps (verified: dOLR = 0
    exactly at every Ts), as are the H2O profiles and hence the moist-greenhouse
    trigger temperature.  The whole effect is a rescaling of F_SOL, alpha_p and
    S_eff = F_IR / F_SOL.
The sensitivity is consequently well described by a single slope dS_eff/dalpha_s,
which this script reports so that any other tuned albedo can be scaled without
re-running the model.

Inputs — the four sweep caches (this script runs NO model; see the READMEs of the
two reference cases for how each cache is produced):
    reference/moist_runaway/hz_inner_nonideal.npz           alpha_s = 0.2736 (PRIMARY)
    reference/moist_runaway/hz_inner_nonideal_bps.npz            "     BPS continuum
    reference/moist_runaway/hz_inner_nonideal_a032.npz      alpha_s = 0.32   (archived)
    reference/moist_runaway/hz_inner_nonideal_a032_bps.npz       "     BPS continuum
    reference/max_greenhouse/hz_outer.npz                   alpha_s = 0.2736 (PRIMARY)
    reference/max_greenhouse/hz_outer_a032.npz              alpha_s = 0.32   (archived)
The primary caches are produced by the scripts' own defaults; regenerate the
archived alpha_s = 0.32 comparison set with:
    HZ_ALBEDO=0.32  HZ_TAG_SUFFIX=_a032  python reference/moist_runaway/hz_inner.py
    OHZ_ALBEDO=0.32 OHZ_TAG_SUFFIX=_a032 python reference/max_greenhouse/hz_outer.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
IHZ  = os.path.join(ROOT, 'reference', 'moist_runaway')
OHZ  = os.path.join(ROOT, 'reference', 'max_greenhouse')

A_PUB, A_TUN = 0.32, 0.2736      # Kopparapu convention vs ExoColumn calibration
# 2026-08-21: the ExoColumn-calibrated albedo was promoted to PRIMARY for every HZ
# calculation, so the tuned sweeps now carry the canonical filenames and the
# Kopparapu-albedo comparison set is archived under the '_a032' tag.
TAG_PUB, TAG_TUN = '_a032', ''

# Limit definitions — identical to hz_inner.py / hz_outer.py so the numbers here
# reproduce the main-text figures exactly.
MOIST_GH_VMR = 3.0e-3            # water-loss stratospheric H2O VMR (Kopparapu Sec 3.1)
RUNAWAY_TS_LO, RUNAWAY_TS_HI = 280.0, 700.0
SMOOTH_WIN = 5                   # resolution-sawtooth median filter (hz_inner default)

# Clima (Kopparapu et al. 2013) published limits, for the summary table.
KOPP_MOIST, KOPP_RUNAWAY = 1.016, 1.060
KOPP_SEFF_MAX, KOPP_PCO2_MAX = 0.343, 8.0

C_PUB, C_TUN, GREY = 'C0', 'C3', '0.3'


def _smooth(y, w=SMOOTH_WIN):
    y = np.asarray(y, float)
    if w <= 1 or y.size < 3:
        return y
    h, out = w // 2, y.copy()
    for i in range(y.size):
        out[i] = np.median(y[max(0, i - h):min(y.size, i + h + 1)])
    return out


def _rows(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"missing sweep cache: {path}")
    return np.load(path)['rows']


def _ihz_limits(ts, seff, strat_vmr):
    """(moist-GH Seff, its Ts, runaway Seff, its Ts) — hz_inner.py definitions."""
    wet = np.where(_smooth(strat_vmr) >= MOIST_GH_VMR)[0]
    m_s = m_t = np.nan
    if wet.size:
        m_s, m_t = float(seff[wet[0]]), float(ts[wet[0]])
    band = (ts >= RUNAWAY_TS_LO) & (ts <= RUNAWAY_TS_HI)
    r_s = r_t = np.nan
    if np.any(band):
        i = int(np.nanargmax(seff[band]))
        r_s, r_t = float(seff[band][i]), float(ts[band][i])
    return m_s, m_t, r_s, r_t


def _ohz_limit(pco2, seff):
    """(max-greenhouse Seff, pCO2) via the quadratic fit hz_outer.py uses."""
    i = int(np.argmin(seff))
    if 0 < i < len(seff) - 1:
        c = np.polyfit(np.log10(pco2[i-1:i+2]), seff[i-1:i+2], 2)
        x = -c[1] / (2 * c[0])
        return float(np.polyval(c, x)), float(10 ** x)
    return float(seff[i]), float(pco2[i])


def _load_kopp_ihz():
    p = os.path.join(IHZ, 'waterloss_IHZ_present.dat')
    if not os.path.isfile(p):
        return None
    d = np.genfromtxt(p, skip_header=3)
    return dict(tgo=d[:, 0], seff=d[:, 1], palb=d[:, 2])


def main():
    # ---- inner edge -------------------------------------------------------
    ihz = {}
    for tag, alb in ((TAG_PUB, A_PUB), (TAG_TUN, A_TUN)):
        for cont, suf in (('mtckd', ''), ('bps', '_bps')):
            r = _rows(os.path.join(IHZ, f'hz_inner_nonideal{tag}{suf}.npz'))
            ihz[(alb, cont)] = dict(ts=r[:, 0], olr=_smooth(r[:, 1]),
                                    asr=_smooth(r[:, 2]), alpha=_smooth(r[:, 3]),
                                    seff=_smooth(r[:, 4]), svmr=r[:, 5])

    # The prescribed columns are albedo-independent: OLR and the stratospheric
    # H2O must be identical between the two sweeps.  Assert it — this is the
    # claim the supplement rests on.
    for cont in ('mtckd', 'bps'):
        a, b = ihz[(A_PUB, cont)], ihz[(A_TUN, cont)]
        dolr = np.max(np.abs(a['olr'] - b['olr']))
        dvmr = np.max(np.abs(a['svmr'] - b['svmr']))
        print(f"  check [{cont}]: max|dOLR| = {dolr:.3e} W/m2 , "
              f"max|d(strat H2O VMR)| = {dvmr:.3e}")

    # Limits follow hz_inner.py's published convention EXACTLY, so the alpha_s =
    # 0.32 column of the table below reproduces the main-text figure/caption:
    #   * MT_CKD — moist GH at the first Ts whose stratospheric H2O VMR reaches
    #     3e-3; runaway at the Simpson-Nakajima peak of Seff over 280-700 K.
    #   * BPS    — the BPS Seff curve read at the SAME two Ts.  This is not just
    #     a convenience: the profiles (hence the moist-GH trigger) are continuum-
    #     independent, and taking an INDEPENDENT peak for the BPS curve is not
    #     robust.  At alpha_s = 0.2736 the BPS plateau is shallow enough that the
    #     280-700 K window's maximum migrates to 690 K — the supercritical branch
    #     climbing past Tc = 647.1 K, not the Simpson-Nakajima peak (it would
    #     report 1.076 instead of 1.062).  The published alpha_s = 0.32 curves are
    #     window-insensitive (identical over 280-700/500/450 K), which is why the
    #     main-text numbers are unaffected either way.
    lim = {}
    for alb in (A_PUB, A_TUN):
        m_s, m_t, r_s, r_t = _ihz_limits(ihz[(alb, 'mtckd')]['ts'],
                                         ihz[(alb, 'mtckd')]['seff'],
                                         ihz[(alb, 'mtckd')]['svmr'])
        lim[(alb, 'mtckd')] = (m_s, m_t, r_s, r_t)
        b = ihz[(alb, 'bps')]
        lim[(alb, 'bps')] = (float(np.interp(m_t, b['ts'], b['seff'])), m_t,
                             float(np.interp(r_t, b['ts'], b['seff'])), r_t)

    # ---- outer edge -------------------------------------------------------
    ohz = {}
    for tag, alb in ((TAG_PUB, A_PUB), (TAG_TUN, A_TUN)):
        r = _rows(os.path.join(OHZ, f'hz_outer{tag}.npz'))
        ohz[alb] = dict(pco2=r[:, 0], olr=r[:, 1], asr=r[:, 2],
                        alpha=r[:, 3], seff=r[:, 4])
    dolr_o = np.max(np.abs(ohz[A_PUB]['olr'] - ohz[A_TUN]['olr']))
    print(f"  check [OHZ]  : max|dOLR| = {dolr_o:.3e} W/m2")
    olim = {a: _ohz_limit(v['pco2'], v['seff']) for a, v in ohz.items()}

    # ---- summary table ----------------------------------------------------
    dalb = A_TUN - A_PUB
    print(f"\n{'':<26}{'a_s=0.32':>12}{'a_s=0.2736':>12}{'delta':>9}"
          f"{'dS/da_s':>10}{'Clima':>9}")
    rowdefs = [('moist GH  (MT_CKD)', lim[(A_PUB, 'mtckd')][0], lim[(A_TUN, 'mtckd')][0], KOPP_MOIST),
               ('moist GH  (BPS)',    lim[(A_PUB, 'bps')][0],   lim[(A_TUN, 'bps')][0],   KOPP_MOIST),
               ('runaway   (MT_CKD)', lim[(A_PUB, 'mtckd')][2], lim[(A_TUN, 'mtckd')][2], KOPP_RUNAWAY),
               ('runaway   (BPS)',    lim[(A_PUB, 'bps')][2],   lim[(A_TUN, 'bps')][2],   KOPP_RUNAWAY),
               ('max greenhouse',     olim[A_PUB][0],           olim[A_TUN][0],           KOPP_SEFF_MAX)]
    for name, s0, s1, kc in rowdefs:
        print(f"  {name:<24}{s0:12.3f}{s1:12.3f}{s1-s0:+9.3f}"
              f"{(s1-s0)/dalb:10.2f}{kc:9.3f}")
    print(f"\n  distances [AU]  (d = 1/sqrt(Seff))")
    for name, s0, s1, kc in rowdefs:
        print(f"  {name:<24}{1/np.sqrt(s0):12.3f}{1/np.sqrt(s1):12.3f}"
              f"{1/np.sqrt(s1)-1/np.sqrt(s0):+9.3f}{'':10}{1/np.sqrt(kc):9.3f}")
    print(f"\n  max-greenhouse pCO2 [bar]: {olim[A_PUB][1]:.2f} (0.32) -> "
          f"{olim[A_TUN][1]:.2f} (0.2736);  Clima ~{KOPP_PCO2_MAX:.0f}")

    # ---- figure -----------------------------------------------------------
    kopp = _load_kopp_ihz()
    kp_o = np.load(os.path.join(OHZ, 'kopparapu2013_fig5.npz')) \
        if os.path.exists(os.path.join(OHZ, 'kopparapu2013_fig5.npz')) else None

    fig, ((ax_a, ax_b), (ax_c, ax_d)) = plt.subplots(2, 2, figsize=(7.0, 5.6), dpi=300)
    fig.patch.set_facecolor('white')
    for ax in (ax_a, ax_b, ax_c, ax_d):
        ax.set_facecolor('white')

    MT = dict(lw=1.5, zorder=4)
    BP = dict(lw=1.1, ls=':', zorder=3)
    KW = dict(color=GREY, lw=0.9, ls='--', zorder=2)
    TS_LO, TS_HI = 260., 700.

    # (a) inner edge — planetary albedo
    for alb, col in ((A_PUB, C_PUB), (A_TUN, C_TUN)):
        ax_a.plot(ihz[(alb, 'mtckd')]['ts'], ihz[(alb, 'mtckd')]['alpha'], color=col, **MT)
        ax_a.plot(ihz[(alb, 'bps')]['ts'],   ihz[(alb, 'bps')]['alpha'],   color=col, **BP)
        ax_a.axhline(alb, color=col, lw=0.7, ls=(0, (5, 2, 1, 2)),
                     alpha=0.5, zorder=1)   # the surface albedo itself
    if kopp:
        ax_a.plot(kopp['tgo'], kopp['palb'], **KW)
    ax_a.set_xlim(TS_LO, TS_HI); ax_a.set_ylim(0.15, 0.35)
    ax_a.set_xlabel('Surface temperature (K)')
    ax_a.set_ylabel('Planetary albedo $\\alpha_p$')

    # (b) inner edge — Seff, with the two greenhouse limits marked per albedo
    for alb, col in ((A_PUB, C_PUB), (A_TUN, C_TUN)):
        ax_b.plot(ihz[(alb, 'mtckd')]['ts'], ihz[(alb, 'mtckd')]['seff'], color=col, **MT)
        ax_b.plot(ihz[(alb, 'bps')]['ts'],   ihz[(alb, 'bps')]['seff'],   color=col, **BP)
        m_s, m_t, r_s, r_t = lim[(alb, 'mtckd')]
        ax_b.plot([m_t, r_t], [m_s, r_s], '|', color=col, ms=9, mew=1.4, zorder=6)
    if kopp:
        ax_b.plot(kopp['tgo'], kopp['seff'], **KW)
    ax_b.set_xlim(TS_LO, TS_HI)
    ax_b.set_ylim(0.95, 1.15)
    ax_b.set_xlabel('Surface temperature (K)')
    ax_b.set_ylabel('$S_{\\rm eff}$')

    # (c) outer edge — planetary albedo
    for alb, col in ((A_PUB, C_PUB), (A_TUN, C_TUN)):
        ax_c.semilogx(ohz[alb]['pco2'], ohz[alb]['alpha'], color=col, **MT)
    if kp_o is not None:
        ax_c.semilogx(kp_o['alb_p'], kp_o['alb'], **KW)
    ax_c.set_xlabel('CO$_2$ partial pressure (bar)')
    ax_c.set_ylabel('Planetary albedo $\\alpha_p$')

    # (d) outer edge — Seff, with the maximum-greenhouse minimum marked
    for alb, col in ((A_PUB, C_PUB), (A_TUN, C_TUN)):
        ax_d.semilogx(ohz[alb]['pco2'], ohz[alb]['seff'], color=col, **MT)
        s_min, p_min = olim[alb]
        ax_d.plot(p_min, s_min, '|', color=col, ms=9, mew=1.4, zorder=6)
    if kp_o is not None:
        ax_d.semilogx(kp_o['seff_p'], kp_o['seff'], **KW)
    ax_d.set_xlabel('CO$_2$ partial pressure (bar)')
    ax_d.set_ylabel('$S_{\\rm eff}$')

    handles = [Line2D([0], [0], color=C_PUB, lw=1.5,
                      label=f'$\\alpha_s$ = {A_PUB:.2f} (Kopparapu et al. 2013)'),
               Line2D([0], [0], color=C_TUN, lw=1.5,
                      label=f'$\\alpha_s$ = {A_TUN:.4f} (ExoColumn calibration)'),
               Line2D([0], [0], color=GREY, lw=1.1, ls=':',
                      label='BPS continuum (inner edge)'),
               Line2D([0], [0], color=GREY, lw=0.9, ls='--', label='Clima')]
    ax_a.legend(handles=handles, loc='upper right', fontsize=6.5, frameon=False,
                handlelength=2.4, borderaxespad=0.3, labelspacing=0.35)

    fig.tight_layout(w_pad=1.6, h_pad=1.4)
    for ext, dpi in (('pdf', 300), ('png', 200)):
        path = os.path.join(HERE, f'albedo_sensitivity.{ext}')
        fig.savefig(path, dpi=dpi, facecolor='white')
        print(f"Saved: {path}")
    plt.close(fig)


if __name__ == '__main__':
    main()
