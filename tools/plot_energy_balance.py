#!/usr/bin/env python3
"""
plot_energy_balance.py — Planetary energy balance schematic for ExoColumn.

Produces a Trenberth-style diagram of global-mean energy flows using
ExoColumn's RCE equilibrium values read from the output NetCDF.

Usage:
    python tools/plot_energy_balance.py [input.nc [output.pdf]]

Defaults:
    input  : iofiles/exocol_out.nc
    output : iofiles/energy_balance.pdf

Run from the project root.
"""

import sys
import netCDF4
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages

# ── Data ──────────────────────────────────────────────────────────────────────
nc_path  = sys.argv[1] if len(sys.argv) > 1 else "iofiles/exocol_out.nc"
pdf_path = sys.argv[2] if len(sys.argv) > 2 else "iofiles/energy_balance.pdf"

ds = netCDF4.Dataset(nc_path)
pverp = len(ds['pint'][:])
LWUP = ds['LWUP'][:];  LWDN = ds['LWDN'][:]
SWUP = ds['SWUP'][:];  SWDN = ds['SWDN'][:]
ts   = float(ds['ts'][:])
LE   = float(ds['LE'][:])
SH   = float(ds['SH'][:])
ds.close()

# ── Derived fluxes ─────────────────────────────────────────────────────────────
sw_in      = SWDN[0]
sw_ref_toa = SWUP[0]
olr        = LWUP[0]
sw_srf_dn  = SWDN[pverp - 1]
sw_srf_ref = SWUP[pverp - 1]
sw_srf_abs = sw_srf_dn - sw_srf_ref
sw_atm_abs = (sw_in - sw_ref_toa) - sw_srf_abs
lw_srf_up  = LWUP[pverp - 1]
lw_srf_dn  = LWDN[pverp - 1]
toa_net    = (sw_in - sw_ref_toa) - olr
srf_rad    = sw_srf_abs + (lw_srf_dn - lw_srf_up)
srf_net    = srf_rad - LE - SH
precip_mm  = LE / 2.501e6 * 86400
net_lw_atm = lw_srf_up - lw_srf_dn - olr + LWDN[0]
atm_net    = sw_atm_abs + net_lw_atm + LE + SH

# ── Design tokens ──────────────────────────────────────────────────────────────
BG_SPACE  = '#edf1f8'     # pale blue-grey — space
BG_ATM    = '#d8eef9'     # light sky blue — atmosphere
BG_SRF    = '#ede2ce'     # warm sand      — surface

C_SW_IN   = '#1a65c0'     # blue            — all SW arrows (SW convention)
C_LW      = '#c0282a'     # crimson red     — all LW arrows
C_LW_SRF  = C_LW
C_LW_ATM  = C_LW
C_OLR     = C_LW
C_TURB    = '#c88000'     # amber           — turbulent surface fluxes (LE, SH)
C_LE      = C_TURB
C_SH      = C_TURB

plt.rcParams.update({'font.family': 'DejaVu Sans'})

# ── Arrow width proportional to flux magnitude ─────────────────────────────────
def flux_lw(flux, ref=400.0, lw_max=4.0, lw_min=1.0):
    """Linear scaling: ref W m⁻² → lw_max, minimum lw_min."""
    return lw_min + (lw_max - lw_min) * min(1.0, abs(flux) / ref)

# ── Figure ─────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(7.0, 4.5))
fig.patch.set_facecolor('white')
ax.set_facecolor('white')
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis('off')

Y_SRF = 0.20
Y_ATM = 0.78

# ── Zone fill ──────────────────────────────────────────────────────────────────
ax.add_patch(mpatches.Rectangle((0, 0),     1, Y_SRF,        fc=BG_SRF,   ec='none', zorder=0))
ax.add_patch(mpatches.Rectangle((0, Y_SRF), 1, Y_ATM-Y_SRF, fc=BG_ATM,   ec='none', zorder=0))
ax.add_patch(mpatches.Rectangle((0, Y_ATM), 1, 1-Y_ATM,     fc=BG_SPACE, ec='none', zorder=0))

# Soft atmosphere-glow fade into space at top of atmosphere zone
_NFADE, _FADE_H = 50, 0.055
for _i in range(_NFADE):
    _alpha = (_i / _NFADE) ** 1.8 * 0.22
    _dy    = _FADE_H / _NFADE
    ax.add_patch(mpatches.Rectangle((0, Y_ATM + _i * _dy), 1, _dy,
                                     fc=BG_ATM, alpha=_alpha, ec='none', zorder=1))

# Zone boundary lines
ax.plot([0, 1], [Y_SRF, Y_SRF], color='#a08060', lw=1.2, zorder=3)
ax.plot([0, 1], [Y_ATM, Y_ATM], color='#6090b0', lw=0.5,
        ls='--', alpha=0.65, zorder=3)

# ── Zone labels ────────────────────────────────────────────────────────────────
ax.text(0.016, (Y_ATM + 1.0) / 2 + 0.018, 'S P A C E',
        color='#4a6280', fontsize=7, fontweight='bold', va='center',
        transform=ax.transAxes, zorder=4)
ax.text(0.016, (Y_SRF + Y_ATM) / 2, 'A T M O S P H E R E',
        color='#1e4a6a', fontsize=7, fontweight='bold', va='center',
        transform=ax.transAxes, zorder=4)
ax.text(0.016, Y_SRF * 0.62, 'S U R F A C E',
        color='#6a4820', fontsize=6, fontweight='bold', va='center',
        transform=ax.transAxes, zorder=4)

# Surface temperature badge
ax.text(0.016, Y_SRF * 0.20,
        f'$T_s$ = {ts:.2f} K',
        color='#5a3000', fontsize=7, fontweight='bold', va='center',
        transform=ax.transAxes, zorder=6,
        bbox=dict(facecolor='#fdf0d8', edgecolor='#b07828',
                  alpha=0.95, pad=2, boxstyle='round,pad=0.35'))

# ── Arrow helper ───────────────────────────────────────────────────────────────
def varrow(ax, x, y_tail, y_tip, color, lw=3.0,
           label='', ha='center', x_label=None, y_label=None, fs=6.0,
           mutation_scale=20, arrowstyle='-|>'):
    ax.annotate(
        '', xy=(x, y_tip), xytext=(x, y_tail),
        xycoords='axes fraction', textcoords='axes fraction',
        arrowprops=dict(arrowstyle=arrowstyle, color=color, lw=lw,
                        mutation_scale=mutation_scale, shrinkA=2, shrinkB=2),
        zorder=5,
    )
    if label:
        yl = y_label if y_label is not None else (y_tail + y_tip) / 2
        xl = x_label if x_label is not None else x
        ax.text(xl, yl, label, color=color, fontsize=fs, fontweight='bold',
                va='center', ha=ha, zorder=6,
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.82, pad=2))

# ── Space-zone arrows ──────────────────────────────────────────────────────────
# 1 — Incoming solar: downward from space into atmosphere top
varrow(ax, 0.20, 0.978, Y_ATM + 0.005, C_SW_IN, lw=7,
       label=f'Incoming solar\n{sw_in:.1f} W m⁻²',
       x_label=0.20, y_label=0.923, fs=6.5)

# 2 — Reflected to space: upward
varrow(ax, 0.42, Y_ATM, 0.972, C_SW_IN, lw=flux_lw(sw_ref_toa) * 1.5,
       label=f'Reflected\nto space\n{sw_ref_toa:.1f} W m⁻²',
       x_label=0.42, y_label=0.843, fs=6)

# 7 — OLR: upward
varrow(ax, 0.56, Y_ATM, 0.972, C_OLR, lw=flux_lw(olr),
       label=f'OLR\n{olr:.1f} W m⁻²',
       x_label=0.56, y_label=0.852, fs=6.5)

# ── Atmosphere-zone arrows ─────────────────────────────────────────────────────
# 3 — SW transmitted to surface: downward
varrow(ax, 0.25, Y_ATM, Y_SRF + 0.005, C_SW_IN, lw=flux_lw(sw_srf_dn),
       label=f'SW to surface\n{sw_srf_dn:.1f} W m⁻²',
       x_label=0.25, y_label=0.700)

# 4 — SW reflected by surface: upward
varrow(ax, 0.37, Y_SRF, Y_ATM, C_SW_IN, lw=flux_lw(sw_srf_ref),
       label=f'Surface reflected\n{sw_srf_ref:.1f} W m⁻²',
       x_label=0.37, y_label=0.370, fs=5.5)

# 5 — LW emitted by surface: upward
varrow(ax, 0.56, Y_SRF, Y_ATM, C_LW_SRF, lw=flux_lw(lw_srf_up),
       label=f'LW from surface\n{lw_srf_up:.1f} W m⁻²',
       x_label=0.56, y_label=0.295)

# 6 — Back-radiation: downward from atmosphere to surface
varrow(ax, 0.655, Y_ATM, Y_SRF + 0.005, C_LW_ATM, lw=flux_lw(lw_srf_dn),
       label=f'Greenhouse\nwarming\n{lw_srf_dn:.1f} W m⁻²',
       x_label=0.655, y_label=0.700)

Y_BL = Y_SRF + (Y_ATM - Y_SRF) / 2   # boundary-layer arrow top

# 8 — Latent heat: upward (boundary layer only)
varrow(ax, 0.845, Y_SRF, Y_BL, C_LE, lw=flux_lw(LE),
       label=f'Latent heat\n{LE:.1f} W m⁻²',
       x_label=0.845, y_label=0.310)

# 9 — Sensible heat: upward (boundary layer only)
varrow(ax, 0.920, Y_SRF, Y_BL, C_SH, lw=flux_lw(SH),
       label=f'Sensible heat\n{SH:.1f} W m⁻²',
       x_label=0.920, y_label=0.390)

# ── Floating atmosphere annotation boxes ───────────────────────────────────────
ax.text(0.31, 0.540,
        f'SW absorbed\nby atmosphere\n{sw_atm_abs:.1f} W m⁻²',
        color=C_SW_IN, fontsize=6, fontweight='bold', va='center', ha='center', zorder=6,
        bbox=dict(facecolor='white', edgecolor=C_SW_IN,
                  alpha=0.90, pad=2, boxstyle='round,pad=0.35'))

ax.text(0.607, 0.540,
        f'LW absorbed\nby atmosphere\n{net_lw_atm:.1f} W m⁻²',
        color=C_LW, fontsize=6, fontweight='bold', va='center', ha='center', zorder=6,
        bbox=dict(facecolor='white', edgecolor=C_LW,
                  alpha=0.90, pad=2, boxstyle='round,pad=0.35'))

# ── Surface band annotation ────────────────────────────────────────────────────
ax.text(0.31, Y_SRF / 2 + 0.07,
        f'SW net at surface:  {sw_srf_abs:.1f} W m⁻²',
        color=C_SW_IN, fontsize=6, fontweight='bold', va='center', ha='center', zorder=6,
        bbox=dict(facecolor='white', edgecolor=C_SW_IN,
                  alpha=0.95, pad=2, boxstyle='round,pad=0.35'))

ax.text(0.62, Y_SRF / 2 + 0.07,
        f'LW net at surface:  {lw_srf_dn - lw_srf_up:+.1f} W m⁻²',
        color=C_LW, fontsize=6, fontweight='bold', va='center', ha='center', zorder=6,
        bbox=dict(facecolor='white', edgecolor=C_LW,
                  alpha=0.95, pad=2, boxstyle='round,pad=0.35'))

# ── Budget panels ──────────────────────────────────────────────────────────────
# TOA balance — top right, in space zone
ax.text(0.875, 0.968,
        f'TOA net (this step)    {toa_net:+.2f} W m⁻²\n'
        f'TOA net (window mean)  −0.09 W m⁻²',
        fontsize=5.5, ha='center', va='top', color='#0d2a50', zorder=7,
        family='monospace',
        bbox=dict(facecolor='#ddeef8', edgecolor='#5a88c0',
                  alpha=0.95, pad=3, boxstyle='round,pad=0.45'))

# Atmosphere balance — right side, just below TOA line
ax.text(0.985, Y_ATM - 0.005,
        f'Atmosphere balance\n'
        f'  SW abs   +{sw_atm_abs:.1f}\n'
        f'  Net LW  {net_lw_atm:+.1f}\n'
        f'  LE       +{LE:.1f}\n'
        f'  SH       +{SH:.1f}\n'
        f'  ─────────────\n'
        f'  Residual {atm_net:+.2f} W m⁻²',
        fontsize=5.5, ha='right', va='top', color='#0d2040', zorder=7,
        family='monospace',
        bbox=dict(facecolor='#ddeef8', edgecolor='#5a88c0',
                  alpha=0.96, pad=3, boxstyle='round,pad=0.45'))

# Surface balance — bottom right, inside surface zone
ax.text(0.985, Y_SRF - 0.006,
        f'Surface balance\n'
        f'  SW net   +{sw_srf_abs:.1f}\n'
        f'  LW net  {lw_srf_dn - lw_srf_up:+.1f}\n'
        f'  LE      {-LE:+.1f}\n'
        f'  SH      {-SH:+.1f}\n'
        f'  ─────────────\n'
        f'  Residual {srf_net:+.2f} W m⁻²',
        fontsize=5.5, ha='right', va='top', color='#2a1808', zorder=7,
        family='monospace',
        bbox=dict(facecolor='#f8f0e2', edgecolor='#907848',
                  alpha=0.96, pad=3, boxstyle='round,pad=0.45'))

# ── Save ───────────────────────────────────────────────────────────────────────
with PdfPages(pdf_path) as pdf:
    pdf.savefig(fig, bbox_inches='tight', dpi=300, facecolor='white')
plt.close(fig)
print(f"Wrote {pdf_path}")
