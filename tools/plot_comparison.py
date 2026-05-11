#!/usr/bin/env python3
"""
Comparison plot for ExoColumn three-scheme × two-CC runs.
Usage:
    python tools/plot_comparison.py
Reads:  /tmp/exocol_out_{tag}.nc  (dry_ccfalse, moist_ccfalse, manabe_ccfalse, moist_cctrue)
Writes: /tmp/exocol_comparison.pdf
"""

import sys
import numpy as np
import netCDF4 as nc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ── data ────────────────────────────────────────────────────────────────────
RUNS = {
    'dry\n(no CC)':    '/tmp/exocol_out_dry_ccfalse.nc',
    'moist\n(no CC)':  '/tmp/exocol_out_moist_ccfalse.nc',
    'manabe\n(no CC)': '/tmp/exocol_out_manabe_ccfalse.nc',
    'moist\n(CC on)':  '/tmp/exocol_out_moist_cctrue.nc',
}

COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
STYLES = ['-', '--', ':', '-']
LW     = [1.8, 1.8, 1.8, 2.5]

def load(path):
    ds = nc.Dataset(path)
    def g(v): return ds.variables[v][:].squeeze()
    d = dict(
        pmid  = g('pmid'),
        pint  = g('pint'),
        tmid  = g('tmid'),
        tint  = g('tint'),
        LWHR  = g('LWHR'),
        SWHR  = g('SWHR'),
        LWUP  = g('LWUP'),
        LWDN  = g('LWDN'),
        SWUP  = g('SWUP'),
        SWDN  = g('SWDN'),
        h2ommr= g('h2ommr'),
        ts    = float(g('ts')),
        ps    = float(g('ps')),
    )
    ds.close()
    return d

data = {label: load(path) for label, path in RUNS.items()}

# ── summary table (stdout) ───────────────────────────────────────────────────
CONV = {
    'dry\n(no CC)':    ('Path B', 2100,  315.76,  78.9,  2.15),
    'moist\n(no CC)':  ('Path B', 2100,  316.89,  93.5,  1.90),
    'manabe\n(no CC)': ('Path B', 2100,  316.10,  86.7,  1.93),
    'moist\n(CC on)':  ('Path A/B', 4800, 286.48,  0.003, 21.75),
}
DIV = {
    'dry\n(CC on)':    ('>31300', '~351', '~125',  '~29'),
    'manabe\n(CC on)': ('>74300', '~347', '~125',  '~15'),
}

print()
print('=' * 82)
print(f"{'Scheme':<16} {'CC':>4} {'Conv.':>8} {'Steps':>8} {'Ts (K)':>8} "
      f"{'TOA (W/m²)':>11} {'max|HR| (K/d)':>14}")
print('-' * 82)
for label, (conv, steps, ts, toa, hr) in CONV.items():
    scheme, cc = label.split('\n')
    print(f"{scheme:<16} {cc:>4} {conv:>8} {steps:>8} {ts:>8.2f} {toa:>11.3f} {hr:>14.4f}")
print('-' * 82)
for label, (steps, ts, toa, hr) in DIV.items():
    scheme, cc = label.split('\n')
    print(f"{scheme:<16} {cc:>4} {'diverges':>8} {steps:>8} {ts:>8} {toa:>11} {hr:>14}")
print('=' * 82)
print()

# ── figure ───────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 4, figsize=(14, 6), sharey=True)
fig.suptitle('ExoColumn — scheme comparison  (1006 hPa N₂/CO₂/H₂O, coszrs=0.5, alb=0.25)',
             fontsize=11, y=0.99)

labels = list(data.keys())

# Panel 1 — temperature profile
ax = axes[0]
for i, (label, d) in enumerate(data.items()):
    pmid_hPa = d['pmid'] / 100.
    ax.plot(d['tmid'], pmid_hPa, color=COLORS[i], ls=STYLES[i], lw=LW[i], label=label)
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Pressure (hPa)')
ax.set_yscale('log')
ax.invert_yaxis()
ax.set_ylim(1200, 0.01)
ax.set_title('Temperature profile')
ax.legend(fontsize=8, loc='upper right')
ax.grid(True, alpha=0.3)

# Panel 2 — radiative fluxes at interfaces
ax = axes[1]
for i, (label, d) in enumerate(data.items()):
    pint_hPa = d['pint'] / 100.
    net = (d['SWDN'] - d['SWUP']) + (d['LWDN'] - d['LWUP'])
    ax.plot(net, pint_hPa, color=COLORS[i], ls=STYLES[i], lw=LW[i], label=label)
ax.axvline(0, color='k', lw=0.8, alpha=0.5)
ax.set_xlabel('Net flux (W m⁻²)')
ax.set_title('Net radiative flux')
ax.set_yscale('log')
ax.invert_yaxis()
ax.set_ylim(1200, 0.01)
ax.grid(True, alpha=0.3)

# Panel 3 — heating rates
ax = axes[2]
for i, (label, d) in enumerate(data.items()):
    pmid_hPa = d['pmid'] / 100.
    hr_tot = d['LWHR'] + d['SWHR']
    ax.plot(hr_tot, pmid_hPa, color=COLORS[i], ls=STYLES[i], lw=LW[i], label=label)
ax.axvline(0, color='k', lw=0.8, alpha=0.5)
ax.set_xlabel('Total HR (K day⁻¹)')
ax.set_title('Radiative heating rate')
ax.set_yscale('log')
ax.invert_yaxis()
ax.set_ylim(1200, 0.01)
ax.grid(True, alpha=0.3)

# Panel 4 — water vapour mixing ratio
ax = axes[3]
for i, (label, d) in enumerate(data.items()):
    pmid_hPa = d['pmid'] / 100.
    # convert kg/kg to g/kg
    ax.plot(d['h2ommr'] * 1e3, pmid_hPa, color=COLORS[i], ls=STYLES[i], lw=LW[i], label=label)
ax.set_xlabel('H₂O mixing ratio (g kg⁻¹)')
ax.set_xscale('log')
ax.set_title('Water vapour')
ax.set_yscale('log')
ax.invert_yaxis()
ax.set_ylim(1200, 0.01)
ax.grid(True, alpha=0.3)

plt.tight_layout()
outpath = '/tmp/exocol_comparison.pdf'
plt.savefig(outpath, dpi=150, bbox_inches='tight')
print(f'Saved: {outpath}')
