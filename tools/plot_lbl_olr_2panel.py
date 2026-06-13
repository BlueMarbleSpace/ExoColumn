#!/usr/bin/env python3
"""
plot_lbl_olr_2panel.py — side-by-side clear-sky OLR line-by-line benchmark figure:
  LEFT  : inner-HZ moist column, Ts = 300 K (1 bar N2 + 330 ppm CO2 + H2O)
          — reference/moist_runaway/lbl_olr_benchmark_ts300.npz (+ CLIMA bands).
          The Kopparapu et al. (2013) Figure-2 analogue (dense H2O).
  RIGHT : outer-HZ maximum-greenhouse column, Ts = 273 K, pCO2 = 8.87 bar over
          1 bar N2 (Seff-minimum of reference/max_greenhouse)
          — reference/max_greenhouse/lbl_olr_co2_maxgh.npz.
          The Kopparapu et al. (2013) Figure-1 analogue (dense CO2 + faint sun),
          with CO2 lines + Perrin-Hartmann 1989 sub-Lorentzian chi-factor +
          HITRAN-2024 CO2-CO2 CIA + trace H2O.

Each column: spectral OLR density (LBL fine, LBL averaged onto the n68 bands,
ExoRT n68) over a faint surface-T blackbody envelope, with the ExoRT - LBL band
residual beneath.

Usage: python tools/plot_lbl_olr_2panel.py
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


def panel(axes, npz, title, xlim, ylim, clima_txt=None, smooth_cm=6.0):
    """One column: grey fine LBL + black LBL n68-band averages + red ExoRT n68
    + green Clima (Kopparapu 2013); residual = model - LBL beneath.  The OHZ npz
    carries both wing bounds (olr_nu_full = PH89 chi, olr_nu_full_nochi =
    pure-Lorentz); per the figure spec the reference LBL is the PH89-chi case
    (olr_nu_full), the same sub-Lorentzian convention as Kopparapu/CLIMA."""
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

    clima_b = clima_e = None
    if clima_txt and os.path.isfile(clima_txt):
        clima = np.loadtxt(clima_txt)
        clima_e = np.append(clima[:, 0], clima[-1, 1])
        a.stairs(clima[:, 2] / np.diff(clima_e), clima_e, color='C2', lw=1.0,
                 ls='--', zorder=3, label='Clima (Kopparapu et al. 2013)')
        clima_b = clima[:, 2]

    a.set_xlim(*xlim); a.set_ylim(*ylim)
    a.set_ylabel('F$_{IR}$ spectral density (W m$^{-2}$ / cm$^{-1}$)')
    a.legend(loc='upper right', fontsize=7.5, frameon=False)
    a.set_title(title, fontsize=9.5)
    a.set_facecolor('white')

    # residual panel: model - LBL (LBL = PH89-chi reference)
    ctr = 0.5 * (edges[:-1] + edges[1:])
    b.bar(ctr, (exo - lbl_b) / w_b, width=w_b * 0.92, color='C3', alpha=0.75,
          label='ExoRT − LBL')
    if clima_b is not None:
        lbl_c = band_integrate(wn, oc, clima_e)
        b.step(np.append(clima_e[0], clima_e[1:]),
               np.append((clima_b - lbl_c) / np.diff(clima_e), np.nan),
               where='post', color='C2', lw=1.0, label='Clima − LBL')
    b.axhline(0, color='k', lw=0.6)
    b.set_ylabel('model − LBL\n(W m$^{-2}$ / cm$^{-1}$)')
    b.legend(fontsize=7, frameon=False, loc='lower right')
    b.set_xlabel('Wavenumber (cm$^{-1}$)')
    b.set_xlim(*xlim); b.set_facecolor('white')

    inr = np.isfinite(lbl_b)
    return float(np.nansum(lbl_b[inr])), float(np.nansum(exo[inr]))


def main():
    fig, ax = plt.subplots(2, 2, figsize=(12.0, 6.2), dpi=300,
                           height_ratios=[2.2, 1])
    fig.patch.set_facecolor('white')

    ihz = os.path.join(ROOT, 'reference', 'moist_runaway', 'lbl_olr_benchmark_ts300.npz')
    clima = os.path.join(ROOT, 'reference', 'moist_runaway', 'clima_band_olr_ts300.txt')
    ohz = os.path.join(ROOT, 'reference', 'max_greenhouse', 'lbl_olr_co2_maxgh.npz')
    clima_o = os.path.join(ROOT, 'reference', 'max_greenhouse', 'clima_band_olr_maxgh.txt')

    lbl_i, exo_i = panel((ax[0, 0], ax[1, 0]), ihz,
        'Inner HZ: moist $T_s$ = 300 K (1 bar N$_2$ + 330 ppm CO$_2$ + H$_2$O)',
        xlim=(10, 2000), ylim=(0, 0.42), clima_txt=clima)
    lbl_o, exo_o = panel((ax[0, 1], ax[1, 1]), ohz,
        'Outer HZ: max-greenhouse $T_s$ = 273 K (pCO$_2$ = 8.87 bar + 1 bar N$_2$)',
        xlim=(10, 1600), ylim=(0, 0.32), clima_txt=clima_o)

    print(f"IHZ (Ts=300): LBL={lbl_i:.1f}  ExoRT n68={exo_i:.1f}  diff={exo_i-lbl_i:+.1f}")
    print(f"OHZ (Ts=273): LBL(PH89)={lbl_o:.1f}  ExoRT n68={exo_o:.1f}  diff={exo_o-lbl_o:+.1f}")

    fig.tight_layout()
    for ext, dpi in (('png', 200), ('pdf', 300)):
        out = os.path.join(ROOT, 'reference', 'max_greenhouse',
                           f'lbl_olr_benchmark_2panel.{ext}')
        fig.savefig(out, facecolor='white', dpi=dpi)
        print(f"wrote {out}")
    plt.close(fig)


if __name__ == '__main__':
    main()
