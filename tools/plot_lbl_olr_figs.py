#!/usr/bin/env python3
"""
plot_lbl_olr_figs.py — the two clear-sky OLR line-by-line benchmark figures,
one for each end of the habitable zone.  Each figure carries two columns: the
ExoColumn HZ-limit column on the left, and the matching Kopparapu et al. (2013)
radiative-transfer validation case on the right.

INNER HZ  -> reference/moist_runaway/lbl_olr_benchmark_ihz.{pdf,png}
  left  : moist column, Ts = 300 K (1 bar N2 + 330 ppm CO2 + H2O)
          — reference/moist_runaway/lbl_olr_benchmark_ts300.npz (+ Clima bands).
  right : the Kopparapu et al. (2013) Figure-2 configuration — dense H2O,
          Ts = 400 K, 200 K isothermal stratosphere, 4 bar N2 background
          (ps = 6.46 bar), Earth gravity.  Their published totals for this
          case are 285 W/m2 (Clima) and 297 W/m2 (SMART).

OUTER HZ  -> reference/max_greenhouse/lbl_olr_benchmark_ohz.{pdf,png}
  left  : maximum-greenhouse column, Ts = 273 K, pCO2 = 8.87 bar over 1 bar N2
          (the Seff minimum of reference/max_greenhouse).
  right : the Kopparapu et al. (2013) Figure-1 configuration — early Mars,
          Mars-mass planet (g = 3.73 m/s2), 2 bar 95% CO2 / 5% N2, Ts = 250 K,
          167 K isothermal stratosphere.  Kopparapu's own SMART line-by-line
          spectrum for this column is overlaid; their published totals are
          86 W/m2 (Clima) and 88.4 W/m2 (SMART).

Both Kopparapu-configuration columns are integrated by every model from the
SAME atmosphere: Clima's own profile for the dense-H2O case, and the profile
Kopparapu handed to SMART for early Mars.  Dense-CO2 columns use CO2 lines +
the Perrin & Hartmann (1989) sub-Lorentzian chi-factor + HITRAN-2024 CO2-CO2
CIA + trace H2O; moist columns use H2O + CO2 lines + the MT_CKD continuum.

Each column: spectral OLR density (LBL fine, LBL averaged onto the n68 bands,
ExoRT n68, and Clima in BOTH k-coefficient generations -- the 2013-era set used
by Kopparapu et al. 2013 and the post-2014 Wolf HITRAN-2016 set adopted by the
atmos repo in 2021) over a faint surface-T blackbody envelope, with the
model - LBL band residuals beneath.

Usage: python tools/plot_lbl_olr_figs.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
H, C, KB = 6.62607015e-34, 2.99792458e8, 1.380649e-23


def planck_nu(wn_cm, T):
    nu = wn_cm * 100.0
    B = 2 * H * C**2 * nu**3 / np.expm1(H * C * nu / (KB * T))
    return np.pi * B * 100.0                       # W/m2/cm-1


def band_integrate(wn, oc, edges):
    out = np.full(len(edges) - 1, np.nan)
    for i in range(len(out)):
        m = (wn >= edges[i]) & (wn < edges[i + 1])
        if m.sum() > 2:
            out[i] = np.trapezoid(oc[m], wn[m])
    return out


def panel(axes, npz, title, xlim, ylim, clima_txt=None, smart_txt=None,
          rlim=None, smooth_cm=6.0, legend_loc='upper right',
          rlegend_loc='lower right'):
    """One column: grey fine LBL + black LBL n68-band averages + red ExoRT n68
    + green/blue Clima (2013-era and Wolf-2016 k); residual = model - LBL
    beneath.  The dense-CO2 npz files carry both wing bounds (olr_nu_full =
    PH89 chi, olr_nu_full_nochi = pure-Lorentz); per the figure spec the
    reference LBL is the PH89-chi case (olr_nu_full), the same sub-Lorentzian
    convention as Kopparapu/Clima.  smart_txt, when given, overlays Kopparapu's
    own SMART line-by-line spectrum for the same column."""
    d = np.load(npz)
    wn = d['wn'].astype(float)
    edges, exo = d['band_edges'], d['band_exo']
    ts = float(d['ts'])
    w_b = np.diff(edges)
    a, b = axes

    oc = (d['olr_nu_full'] if 'olr_nu_full' in d.files else d['olr_nu_cont']).astype(float)
    lbl_b = band_integrate(wn, oc, edges)

    # surface-T blackbody envelope, then grey LBL, black LBL n68-band avg, red ExoRT
    wnbb = np.linspace(max(xlim[0], 1.0), xlim[1], 600)
    a.plot(wnbb, planck_nu(wnbb, ts), color='0.8', lw=1.0, zorder=1,
           label=f'{ts:.0f} K blackbody')
    k = max(1, int(round(smooth_cm / np.median(np.diff(wn)))) | 1)
    a.plot(wn[::5], np.convolve(oc, np.ones(k) / k, 'same')[::5], color='0.55',
           lw=0.6, zorder=2, label='LBL (line-by-line)')
    a.stairs(lbl_b / w_b, edges, color='k', lw=1.4, zorder=4,
             label='LBL, n68-band averages')
    a.stairs(exo / w_b, edges, color='C3', lw=1.2, zorder=3,
             label='ExoRT n68 (this work)')

    smart_b = None
    if smart_txt and os.path.isfile(smart_txt):
        sm = np.loadtxt(smart_txt)
        s_wn, s_f = sm[:, 0], sm[:, 1]
        ks = max(1, int(round(smooth_cm / np.median(np.diff(s_wn)))) | 1)
        a.plot(s_wn, np.convolve(s_f, np.ones(ks) / ks, 'same'), color='C4',
               lw=0.7, zorder=2, label='SMART LBL')
        smart_b = band_integrate(s_wn, s_f, edges)

    clima_b = clima_b16 = clima_e = None
    if clima_txt and os.path.isfile(clima_txt):
        clima = np.loadtxt(clima_txt)
        clima_e = np.append(clima[:, 0], clima[-1, 1])
        cw = np.diff(clima_e)
        # col 2 = 2013-era k (Kopparapu 2013 default);
        # col 3 = Wolf HITRAN-2016 k, adopted by the atmos repo in 2021.
        a.stairs(clima[:, 2] / cw, clima_e, color='C2', lw=1.0,
                 zorder=3, label='Clima (Kopparapu et al. 2013 $k$)')
        a.stairs(clima[:, 3] / cw, clima_e, color='C0', lw=1.0,
                 zorder=3, label='Clima (HITRAN-2016 $k$)')
        clima_b, clima_b16 = clima[:, 2], clima[:, 3]

    a.set_xlim(*xlim); a.set_ylim(*ylim)
    a.set_ylabel('F$_{IR}$ spectral density (W m$^{-2}$ / cm$^{-1}$)')
    a.legend(loc=legend_loc, fontsize=7.5, frameon=False)
    a.set_title(title, fontsize=9.5)
    a.set_facecolor('white')

    # residual panel: model - LBL (LBL = PH89-chi reference)
    ctr = 0.5 * (edges[:-1] + edges[1:])
    b.bar(ctr, (exo - lbl_b) / w_b, width=w_b * 0.92, color='C3', alpha=0.75,
          label='ExoRT − LBL')
    if smart_b is not None:
        b.step(np.append(edges[0], edges[1:]),
               np.append((smart_b - lbl_b) / w_b, np.nan),
               where='post', color='C4', lw=1.0, label='SMART − LBL')
    if clima_b is not None:
        lbl_c = band_integrate(wn, oc, clima_e)
        cw = np.diff(clima_e)
        for cb, col, lab in ((clima_b, 'C2', 'Clima 2013 $k$ − LBL'),
                             (clima_b16, 'C0', 'Clima 2016 $k$ − LBL')):
            b.step(np.append(clima_e[0], clima_e[1:]),
                   np.append((cb - lbl_c) / cw, np.nan),
                   where='post', color=col, lw=1.0, label=lab)
    b.axhline(0, color='k', lw=0.6)
    b.set_ylabel('model − LBL\n(W m$^{-2}$ / cm$^{-1}$)')
    b.legend(fontsize=7, frameon=False, loc=rlegend_loc, ncol=2)
    b.set_xlabel('Wavenumber (cm$^{-1}$)')
    b.set_xlim(*xlim); b.set_facecolor('white')
    if rlim:
        b.set_ylim(*rlim)

    inr = np.isfinite(lbl_b)
    ctot = (float(np.nansum(clima_b)), float(np.nansum(clima_b16))) \
        if clima_b is not None else (float('nan'), float('nan'))
    stot = float(np.nansum(smart_b)) if smart_b is not None else float('nan')
    return float(np.nansum(lbl_b[inr])), float(np.nansum(exo[inr])), ctot, stot


def report(tag, r):
    lbl, exo, cl, smart = r
    line = (f"{tag}: LBL={lbl:.1f}  ExoRT n68={exo:.1f}  diff={exo-lbl:+.1f}"
            f"  |  Clima 2013 k={cl[0]:.1f}  Clima HITRAN-2016 k={cl[1]:.1f}")
    if np.isfinite(smart):
        line += f"  |  SMART={smart:.1f}"
    print(line)


def make_figure(specs, outdir, stem):
    fig, ax = plt.subplots(2, 2, figsize=(12.0, 6.2), dpi=300,
                           height_ratios=[2.2, 1])
    fig.patch.set_facecolor('white')
    for j, s in enumerate(specs):
        report(s.pop('tag'), panel((ax[0, j], ax[1, j]), **s))
    fig.tight_layout()
    for ext, dpi in (('png', 200), ('pdf', 300)):
        out = os.path.join(ROOT, 'reference', outdir, f'{stem}.{ext}')
        fig.savefig(out, facecolor='white', dpi=dpi)
        print(f"wrote {out}")
    plt.close(fig)


def main():
    mr = os.path.join(ROOT, 'reference', 'moist_runaway')
    mg = os.path.join(ROOT, 'reference', 'max_greenhouse')

    # ------------------------------------------------------------ inner edge
    make_figure([
        dict(tag='IHZ  Ts=300 (1 bar N2 + 330 ppm CO2)',
             npz=os.path.join(mr, 'lbl_olr_benchmark_ts300.npz'),
             clima_txt=os.path.join(mr, 'clima_band_olr_ts300.txt'),
             title='Inner HZ: moist $T_s$ = 300 K '
                   '(1 bar N$_2$ + 330 ppm CO$_2$ + H$_2$O)',
             xlim=(10, 2000), ylim=(0, 0.42)),
        dict(tag='IHZ  Ts=400 (Kopparapu 2013 Fig. 2)',
             npz=os.path.join(mr, 'lbl_olr_benchmark_ts400.npz'),
             clima_txt=os.path.join(mr, 'clima_band_olr_ts400.txt'),
             title='Inner HZ: dense H$_2$O $T_s$ = 400 K '
                   '(4 bar N$_2$ + H$_2$O, $p_s$ = 6.46 bar)',
             xlim=(10, 2000), ylim=(0, 1.25)),
    ], 'moist_runaway', 'lbl_olr_benchmark_ihz')

    # ------------------------------------------------------------ outer edge
    make_figure([
        dict(tag='OHZ  Ts=273 (max greenhouse)',
             npz=os.path.join(mg, 'lbl_olr_co2_maxgh.npz'),
             clima_txt=os.path.join(mg, 'clima_band_olr_maxgh.txt'),
             title='Outer HZ: max-greenhouse $T_s$ = 273 K '
                   '(pCO$_2$ = 8.87 bar + 1 bar N$_2$)',
             xlim=(10, 1600), ylim=(0, 0.32),
             rlim=(-0.050, 0.030), rlegend_loc='upper right'),
        dict(tag='OHZ  Ts=250 (Kopparapu 2013 Fig. 1)',
             npz=os.path.join(mg, 'lbl_olr_co2_earlymars.npz'),
             clima_txt=os.path.join(mg, 'clima_band_olr_earlymars.txt'),
             smart_txt=os.path.join(mg, 'smart_earlymars_olr.txt'),
             title='Outer HZ: early Mars $T_s$ = 250 K '
                   '(2 bar 95% CO$_2$ / 5% N$_2$, Mars gravity)',
             xlim=(10, 1600), ylim=(0, 0.32),
             rlim=(-0.025, 0.075), rlegend_loc='upper right'),
    ], 'max_greenhouse', 'lbl_olr_benchmark_ohz')


if __name__ == '__main__':
    main()
