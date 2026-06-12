#!/usr/bin/env python3
"""
plot_lbl_olr.py — render the clear-sky OLR line-by-line benchmark figure
(lbl_olr_benchmark_ts300.png) for the saturated Ts = 300 K IHZ column, from
the two co-located data files:

  * lbl_olr_benchmark_ts300.npz — LBL spectrum (RADIS/HITRAN H2O+CO2 lines +
    AER MT_CKD continuum, diffusivity-1.66 Schwarzschild) and the ExoRT n68
    band_lwup_toa for the same column.  Regenerate with the compute engine:
        run/exocol.exe on the hz_inner Ts=300 namelist (co2_vmr_total,
        ice cold trap, nonideal EOS), then
        python tools/lbl_olr_benchmark.py iofiles/exocol_out.nc
    (~25 min; HITRAN via RADIS, cached in ~/.radisdb; see the engine's
    docstring for physics details and caveats).
  * clima_band_olr_ts300.txt — CLIMA 55-interval TOA OLR on the SAME column,
    for both k-coefficient generations (2013-era = Kopparapu 2013, and the
    Wolf HITRAN2016 set adopted by the atmos repo in 2021).  Generated with
    the patched public CLIMA at /models/atmos (see the data file header and
    the README "LW offset" bullet).

Result (totals, 10-3000 cm-1): LBL 272.4, ExoRT n68 269.7, CLIMA-2021 267.5,
CLIMA-2013-era 251.8 (Kopparapu's tabulated value: 250.2) W/m2 — the
2013-era k-coefficient data, since superseded by CLIMA's own maintainers,
carries the F_IR offset between this work and Kopparapu et al. (2013).

Usage: python reference/moist_runaway/plot_lbl_olr.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    d = np.load(os.path.join(HERE, 'lbl_olr_benchmark_ts300.npz'))
    wn, oc = d['wn'].astype(float), d['olr_nu_cont'].astype(float)
    edges, exo = d['band_edges'], d['band_exo']
    ts = float(d['ts'])
    clima = np.loadtxt(os.path.join(HERE, 'clima_band_olr_ts300.txt'))

    # band-integrate the LBL onto the n68 bands
    lbl_b = np.full(len(exo), np.nan)
    for i in range(len(exo)):
        m = (wn >= edges[i]) & (wn < edges[i + 1])
        if m.sum() > 2:
            lbl_b[i] = np.trapezoid(oc[m], wn[m])
    w_b = np.diff(edges)

    fig, ax = plt.subplots(2, 1, figsize=(7.0, 6.2), dpi=300,
                           height_ratios=[2.2, 1], sharex=True)
    fig.patch.set_facecolor('white')
    a = ax[0]
    k = 501
    ker = np.ones(k) / k
    a.plot(wn[::10], np.convolve(oc, ker, 'same')[::10], color='0.55', lw=0.6,
           label='LBL: HITRAN lines + MT_CKD continuum (smoothed)')
    a.stairs(lbl_b / w_b, edges, color='k', lw=1.4, label='LBL, n68-band averages')
    a.stairs(exo / w_b, edges, color='C3', lw=1.2, label='ExoRT n68 (this work)')
    ce = np.append(clima[:, 0], clima[-1, 1])
    cw = np.diff(ce)
    a.stairs(clima[:, 3] / cw, ce, color='C0', lw=1.0, ls='--',
             label='CLIMA, Wolf HITRAN2016 k (2021)')
    a.stairs(clima[:, 2] / cw, ce, color='C2', lw=1.0, ls='--',
             label='CLIMA, 2013-era k (= Kopparapu 2013)')
    totals = ('totals (10–3000 cm$^{-1}$)\n'
              'LBL                272.4 W m$^{-2}$\n'
              'ExoRT n68          269.7 W m$^{-2}$\n'
              f'CLIMA 2021 k       {clima[:, 3].sum():.1f} W m$^{{-2}}$\n'
              f'CLIMA 2013-era k   {clima[:, 2].sum():.1f} W m$^{{-2}}$\n'
              'Kopparapu 2013 tab 250.2 W m$^{-2}$')
    a.set_ylabel('OLR spectral density (W m$^{-2}$ / cm$^{-1}$)')
    a.set_xlim(10, 2000)
    a.set_ylim(0, 0.42)
    a.legend(loc='upper left', fontsize=7.5, frameon=False)
    a.set_title(f'Clear-sky OLR, saturated $T_s$ = {ts:.0f} K column '
                '(1 bar N$_2$ + 330 ppm CO$_2$ + H$_2$O)', fontsize=10)
    a.text(0.985, 0.975, totals, transform=a.transAxes, ha='right', va='top',
           fontsize=8, family='monospace')
    a.set_facecolor('white')

    b = ax[1]
    ctr = 0.5 * (edges[:-1] + edges[1:])
    b.bar(ctr, (exo - lbl_b) / w_b, width=w_b * 0.92, color='C3', alpha=0.75,
          label='ExoRT − LBL')
    lbl_c = np.full(len(cw), np.nan)
    for i in range(len(cw)):
        m = (wn >= ce[i]) & (wn < ce[i + 1])
        if m.sum() > 2:
            lbl_c[i] = np.trapezoid(oc[m], wn[m])
    b.step(np.append(ce[0], ce[1:]), np.append((clima[:, 2] - lbl_c) / cw,
           np.nan), where='post', color='C2', lw=1.0,
           label='CLIMA 2013-era − LBL')
    b.legend(fontsize=7, frameon=False, loc='lower right')
    b.axhline(0, color='k', lw=0.6)
    b.set_ylabel('model − LBL\n(W m$^{-2}$ / cm$^{-1}$)')
    b.set_xlabel('Wavenumber (cm$^{-1}$)')
    b.set_facecolor('white')
    fig.tight_layout()
    out = os.path.join(HERE, 'lbl_olr_benchmark_ts300.png')
    fig.savefig(out, facecolor='white', dpi=200)
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == '__main__':
    main()
