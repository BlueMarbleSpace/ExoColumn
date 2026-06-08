#!/usr/bin/env python3
"""
olr_vs_kopparapu.py  —  direct OLR(Ts) comparison: ExoColumn vs Kopparapu (2013).

OLR is pure longwave: it depends only on the (prescribed) moist-pseudoadiabat
temperature profile + composition, NOT on the surface albedo / shortwave.  So
plotting OLR vs surface temperature isolates the *radiative-transfer* difference
between ExoColumn's ExoRT-n68 correlated-k scheme and Kopparapu (2013)'s updated
CLIMA (HITEMP 2010 + BPS H2O continuum), with no albedo / insolation-geometry
ambiguity.

Primary reference: Kopparapu et al. (2013), ApJ 765, 131, Fig. 3a (net outgoing
IR flux F_IR vs surface temperature; fully saturated "Earth" model, isothermal
200 K stratosphere, Kasting-1988 App-A moist pseudoadiabat, 1 Earth-ocean H2O
inventory, surface albedo 0.32 — the SAME construction as our inner-HZ sweep).
Text: F_IR levels out at 291 W/m^2 (HITEMP+BPS), rising again beyond ~2000 K;
this "closely matches Pierrehumbert (2010) Fig. 4.37."  Digitized from Fig. 3a.

Light secondary reference: Kasting (1988) Fig. 7a (plateau ~310 W/m^2) — the
pre-HITEMP value, shown to mark the lineage (Kasting 310 -> Kopparapu 291 as the
H2O opacity/continuum was updated).

Our OLR(Ts) is read from an hz_inner.py sweep log (OLR is albedo-independent and
PVER-independent to the sawtooth level).  Pass the log path as argv[1]
(default /tmp/hz_albedo032.log).
"""

import os
import re
import sys
import numpy as np
import matplotlib.pyplot as plt

OUT_DIR = os.path.dirname(os.path.abspath(__file__))
LOG = sys.argv[1] if len(sys.argv) > 1 else '/tmp/hz_albedo032.log'

SN_LIMIT = 282.0   # W/m^2  Pierrehumbert/Nakajima analytic Simpson-Nakajima limit

# --- Kopparapu et al. (2013) Fig. 3a F_IR(Ts), digitized [K, W/m^2] ----------
# HITEMP 2010 + BPS continuum; plateau 291 W/m^2 (stated in text), post-runaway
# rise beyond ~2000 K.
KOPP_TS  = np.array([200, 240, 270, 290, 320, 360, 400, 500, 700, 1000,
                     1400, 1800, 2000, 2100, 2200], dtype=float)
KOPP_OLR = np.array([100, 215, 272, 292, 300, 296, 292, 290, 290, 290,
                     290, 291, 296, 350, 440], dtype=float)

# --- Kasting (1988) Fig. 7a F_IR(Ts), digitized [K, W/m^2] (lineage context) --
KAST_TS  = np.array([200, 230, 260, 285, 310, 335, 370, 450, 600, 700,
                     1300, 1400, 1500, 1600, 1700, 1800], dtype=float)
KAST_OLR = np.array([108, 175, 245, 290, 315, 320, 312, 309, 308, 308,
                     309, 316, 345, 375, 400, 418], dtype=float)


def read_our_olr(path):
    pat = re.compile(r'Ts=\s*([\d.]+)\s*K\s+OLR=\s*([\d.]+)')
    ts, olr = [], []
    with open(path, encoding='utf-8') as f:
        for ln in f:
            m = pat.search(ln)
            if m:
                ts.append(float(m.group(1)))
                olr.append(float(m.group(2)))
    return np.array(ts), np.array(olr)


def main():
    ts, olr = read_our_olr(LOG)
    if ts.size == 0:
        raise SystemExit(f"No OLR data parsed from {LOG}")
    pl = (ts >= 700) & (ts <= 1500)
    our_plateau = olr[pl].mean()
    print(f"Read {ts.size} points from {LOG}  (Ts {ts.min():.0f}-{ts.max():.0f} K)")
    print(f"ExoColumn OLR plateau (700-1500 K) = {our_plateau:.1f} W/m^2")
    print(f"Kopparapu (2013)  plateau          = 291 W/m^2  (Δ = {our_plateau-291:+.1f})")
    print(f"Kasting   (1988)  plateau          = 310 W/m^2")

    fig, ax = plt.subplots(figsize=(7.0, 4.5), dpi=300)
    fig.patch.set_facecolor('white'); ax.set_facecolor('white')

    ax.plot(KAST_TS, KAST_OLR, color='0.6', lw=1.0, ls='--', marker='^', ms=3,
            mfc='white', label='Kasting (1988) Fig. 7a [pre-HITEMP]')
    ax.plot(KOPP_TS, KOPP_OLR, color='k', lw=1.4, ls='-', marker='o', ms=4,
            mfc='white', label='Kopparapu (2013) Fig. 3a [HITEMP+BPS]')
    ax.plot(ts, olr, color='C3', lw=2.0, label='ExoColumn (ExoRT-n68)')
    ax.axhline(SN_LIMIT, color='C0', lw=0.8, ls=':',
               label=f'Simpson–Nakajima limit ({SN_LIMIT:.0f} W m$^{{-2}}$)')

    # annotate the plateaus
    ax.annotate(f'{our_plateau:.0f}', xy=(1100, our_plateau), xytext=(0, -13),
                textcoords='offset points', color='C3', ha='center',
                fontsize=8, fontweight='bold')
    ax.annotate('291', xy=(1100, 291), xytext=(0, 7), textcoords='offset points',
                color='k', ha='center', fontsize=8, fontweight='bold')

    ax.set_xlim(200, 2200)
    ax.set_ylim(80, 460)
    ax.set_xlabel('Surface temperature $T_s$ (K)')
    ax.set_ylabel('Outgoing longwave radiation (W m$^{-2}$)')
    ax.legend(fontsize=7.5, framealpha=0.9, loc='upper left')

    fig.tight_layout()
    for ext, dpi in [('pdf', 300), ('png', 150)]:
        p = os.path.join(OUT_DIR, f'olr_vs_kopparapu.{ext}')
        fig.savefig(p, dpi=dpi, facecolor='white')
        print(f"Saved: {p}")
    plt.close(fig)


if __name__ == '__main__':
    main()
