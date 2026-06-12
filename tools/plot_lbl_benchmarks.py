#!/usr/bin/env python3
"""
plot_lbl_benchmarks.py — render the LBL benchmark figures (LW OLR and SW
albedo, from tools/lbl_{olr,sw}_benchmark_ts300.npz), spectral axes in
wavenumber.

Usage: python tools/plot_lbl_benchmarks.py   (renders whatever npz exists)
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))


def plot_lw():
    path = os.path.join(HERE, 'lbl_olr_benchmark_ts300.npz')
    if not os.path.isfile(path):
        print(f"skip LW ({path} missing)")
        return
    d = np.load(path)
    wn, oc = d['wn'], d['olr_nu_cont']
    edges, exo = d['band_edges'], d['band_exo']
    ts = float(d['ts'])

    # band-integrate the LBL
    lbl_b = np.full(len(exo), np.nan)
    for i in range(len(exo)):
        m = (wn >= edges[i]) & (wn < edges[i + 1])
        if m.sum() > 2:
            lbl_b[i] = np.trapezoid(oc[m], wn[m])
    w_b = np.diff(edges)

    # CLIMA band OLR on the matched column (both k-coefficient generations,
    # computed with the patched /models/atmos ir.f — see the data file header)
    clima_path = os.path.join(HERE, 'clima_band_olr_ts300.txt')
    clima = np.loadtxt(clima_path) if os.path.isfile(clima_path) else None

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
    totals = ('totals (10–3000 cm$^{-1}$)\n'
              'LBL                272.4 W m$^{-2}$\n'
              'ExoRT n68          269.7 W m$^{-2}$')
    if clima is not None:
        ce = np.append(clima[:, 0], clima[-1, 1])
        cw = np.diff(ce)
        a.stairs(clima[:, 3] / cw, ce, color='C0', lw=1.0, ls='--',
                 label='CLIMA, Wolf HITRAN2016 k (2021)')
        a.stairs(clima[:, 2] / cw, ce, color='C2', lw=1.0, ls='--',
                 label='CLIMA, 2013-era k (= Kopparapu 2013)')
        totals += (f'\nCLIMA 2021 k       {clima[:,3].sum():.1f} W m$^{{-2}}$'
                   f'\nCLIMA 2013-era k   {clima[:,2].sum():.1f} W m$^{{-2}}$'
                   '\nKopparapu 2013 tab 250.2 W m$^{-2}$')
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
    if clima is not None:
        # CLIMA(2013-era) − LBL on the common 55-interval grid
        lbl_c = np.full(len(cw), np.nan)
        for i in range(len(cw)):
            m = (wn >= ce[i]) & (wn < ce[i + 1])
            if m.sum() > 2:
                lbl_c[i] = np.trapezoid(oc[m], wn[m])
        cctr = 0.5 * (ce[:-1] + ce[1:])
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


def plot_sw():
    path = os.path.join(HERE, 'lbl_sw_benchmark_ts300.npz')
    if not os.path.isfile(path):
        print(f"skip SW ({path} missing)")
        return
    d = np.load(path)
    edges, alb_lbl, alb_exo = d['edges'], d['alb_lbl'], d['alb_exo']
    bswdn = d['bswdn']
    ts = float(d['ts'])
    w_b = np.diff(edges)

    m = np.isfinite(alb_lbl) & (bswdn > 1e-6)
    A_lbl = np.sum(alb_lbl[m] * bswdn[m]) / np.sum(bswdn[m])
    A_exo = np.sum(alb_exo[m] * bswdn[m]) / np.sum(bswdn[m])

    fig, ax = plt.subplots(2, 1, figsize=(7.0, 6.2), dpi=300,
                           height_ratios=[2.2, 1], sharex=True)
    fig.patch.set_facecolor('white')
    a = ax[0]
    a.stairs(np.where(np.isfinite(alb_lbl), alb_lbl, np.nan), edges, color='k',
             lw=1.4, label='LBL: HITRAN + MT_CKD, Toon-quadrature two-stream')
    a.stairs(alb_exo, edges, color='C0', lw=1.2, label='ExoRT n68 (this work)')
    # incident-flux distribution for context (right axis)
    a2 = a.twinx()
    a2.stairs(bswdn / w_b, edges, color='0.75', lw=0.7)
    a2.set_ylabel('incident flux (W m$^{-2}$ / cm$^{-1}$)', color='0.55', fontsize=8)
    a2.tick_params(axis='y', colors='0.55', labelsize=7)
    a2.set_ylim(0, None)
    a.set_xlim(2000, 42100)
    a.set_ylabel('band albedo')
    a.set_ylim(0, 0.8)
    a.legend(loc='upper left', fontsize=8, frameon=False)
    a.set_title(f'Planetary albedo by band, $T_s$ = {ts:.0f} K column '
                '(surface albedo 0.32, 6-node zenith avg)', fontsize=10)
    a.text(0.985, 0.42,
           f'totals (band_swdn weights)\n'
           f'LBL        {A_lbl:.4f}\n'
           f'ExoRT n68  {A_exo:.4f}\n'
           f'(ExoRT − LBL = {A_exo - A_lbl:+.4f})',
           transform=a.transAxes, ha='right', va='top', fontsize=8.5,
           family='monospace')
    a.set_facecolor('white')
    a.set_zorder(2)
    a.patch.set_visible(False)

    b = ax[1]
    ctr = 0.5 * (edges[:-1] + edges[1:])
    diff = np.where(np.isfinite(alb_lbl), alb_exo - alb_lbl, 0.0)
    b.bar(ctr, diff, width=w_b * 0.92, color='C0', alpha=0.75)
    b.axhline(0, color='k', lw=0.6)
    b.set_ylabel('ExoRT − LBL\n(band albedo)')
    b.set_xlabel('Wavenumber (cm$^{-1}$)')
    b.set_facecolor('white')
    fig.tight_layout()
    out = os.path.join(HERE, 'lbl_sw_benchmark_ts300.png')
    fig.savefig(out, facecolor='white', dpi=200)
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == '__main__':
    plot_lw()
    plot_sw()
