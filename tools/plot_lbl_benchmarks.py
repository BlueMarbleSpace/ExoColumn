#!/usr/bin/env python3
"""
plot_lbl_benchmarks.py — render the SW (planetary albedo) LBL benchmark
figure from tools/lbl_sw_benchmark_ts300.npz, spectral axis in wavenumber.
The LW (OLR) benchmark figure lives with the moist_runaway reference case:
reference/moist_runaway/plot_lbl_olr.py.

Usage: python tools/plot_lbl_benchmarks.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))


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
    plot_sw()
