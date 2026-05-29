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

# --------------------------------------------------------------------------
# Read data
# --------------------------------------------------------------------------
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

# ---- Derived fluxes --------------------------------------------------------
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

# --------------------------------------------------------------------------
# Figure
# --------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(14, 9))
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis('off')

Y_SRF = 0.21
Y_ATM = 0.78

# Backgrounds
ax.add_patch(mpatches.Rectangle((0, 0),     1, Y_SRF,        color='#4e3820', zorder=0))
ax.add_patch(mpatches.Rectangle((0, Y_SRF), 1, Y_ATM-Y_SRF, color='#b8dff5', zorder=0))
ax.add_patch(mpatches.Rectangle((0, Y_ATM), 1, 1-Y_ATM,     color='#06101e', zorder=0))

# Band boundary lines
ax.plot([0,1], [Y_SRF, Y_SRF], color='#2e1c0a', lw=2.0, transform=ax.transAxes, zorder=2)
ax.plot([0,1], [Y_ATM, Y_ATM], color='#1a3a5c', lw=1.0, ls='--',
        alpha=0.6, transform=ax.transAxes, zorder=2)

# Band labels — left margin, clear of all arrows/boxes
ax.text(0.013, (Y_ATM + 1.0) / 2, 'SPACE',
        color='white', fontsize=13, fontweight='bold', va='center')
ax.text(0.013, (Y_SRF + Y_ATM) / 2, 'ATMOSPHERE',
        color='#0d2e4a', fontsize=11, fontweight='bold', va='center')
ax.text(0.013, Y_SRF * 0.68, 'SURFACE',
        color='white', fontsize=10, fontweight='bold', va='center')
ax.text(0.013, Y_SRF * 0.28, f'Ts = {ts:.2f} K',
        color='#f0c080', fontsize=9.5, fontweight='bold', va='center')

# --------------------------------------------------------------------------
# Arrow helper — arrowhead always at y_tip (xy); tail at y_tail (xytext).
# For downward arrows: y_tip < y_tail.  '->' puts head at xy = y_tip.
# --------------------------------------------------------------------------
def varrow(ax, x, y_tail, y_tip, color, lw=3.0,
           label='', ha='left', x_label=None, y_label=None, fs=9.0):
    ax.annotate(
        '', xy=(x, y_tip), xytext=(x, y_tail),
        xycoords='axes fraction', textcoords='axes fraction',
        arrowprops=dict(arrowstyle='->', color=color, lw=lw,
                        mutation_scale=18, shrinkA=2, shrinkB=2),
        zorder=5,
    )
    if label:
        yl = y_label if y_label is not None else (y_tail + y_tip) / 2
        xl = x_label if x_label is not None else x
        ax.text(xl, yl, label, color=color, fontsize=fs, fontweight='bold',
                va='center', ha=ha, zorder=6,
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.78, pad=1.5))

# --------------------------------------------------------------------------
# Colours
# --------------------------------------------------------------------------
C_SW_IN  = '#c87808'   # incoming / transmitted SW
C_SW_REF = '#d4a010'   # reflected SW
C_SW_ABS = '#a04000'   # absorbed SW label
C_LW_UP  = '#bb1a1a'   # LW up from surface
C_LW_DN  = '#d86040'   # back-radiation
C_OLR    = '#880606'   # OLR
C_LE     = '#1a58cc'   # latent heat
C_SH     = '#cc5808'   # sensible heat

# --------------------------------------------------------------------------
# Space-zone labels: centred above each arrow, clear of each other.
# Three space arrows: incoming (x=0.15), reflected (x=0.24), OLR (x=0.76).
# Put incoming label at y=0.91, reflected at y=0.84, OLR at y=0.88.
# --------------------------------------------------------------------------

# 1 — Incoming solar: downward, space → atmosphere top
varrow(ax, 0.15, 0.975, Y_ATM + 0.005, C_SW_IN, lw=6,
       label=f'Incoming solar\n{sw_in:.1f} W/m²',
       ha='center', x_label=0.15, y_label=0.912, fs=10)

# 2 — Reflected to space: upward
varrow(ax, 0.24, Y_ATM, 0.970, C_SW_REF, lw=3,
       label=f'Reflected\nto space\n{sw_ref_toa:.1f} W/m²',
       ha='center', x_label=0.24, y_label=0.840, fs=9)

# 7 — OLR: upward, atmosphere top → space  (placed here to keep space labels together)
varrow(ax, 0.76, Y_ATM, 0.970, C_OLR, lw=4.5,
       label=f'OLR\n{olr:.1f} W/m²',
       ha='center', x_label=0.76, y_label=0.880, fs=10)

# --------------------------------------------------------------------------
# Atmosphere-zone labels: three height tiers to avoid overlap.
#   Top tier    y ≈ 0.70  (just below TOA)
#   Middle tier y ≈ 0.50
#   Bottom tier y ≈ 0.29  (just above surface)
# --------------------------------------------------------------------------

# 3 — SW to surface: downward, atmosphere top → surface
#     label at top tier, to the LEFT of arrow
varrow(ax, 0.33, Y_ATM, Y_SRF + 0.005, C_SW_IN, lw=4,
       label=f'SW to surface\n{sw_srf_dn:.1f} W/m²',
       ha='right', x_label=0.315, y_label=0.700, fs=9)

# 4 — SW reflected by surface: upward
#     label at middle tier, to the RIGHT
varrow(ax, 0.42, Y_SRF, Y_ATM, C_SW_REF, lw=2.5,
       label=f'Surface reflected\n{sw_srf_ref:.1f} W/m²  (α={sw_srf_ref/sw_srf_dn:.3f})',
       ha='left', x_label=0.435, y_label=0.500, fs=8.5)

# 5 — LW emitted by surface: upward
#     label at bottom tier, to the LEFT
varrow(ax, 0.56, Y_SRF, Y_ATM, C_LW_UP, lw=4.5,
       label=f'LW emitted\nby surface\n{lw_srf_up:.1f} W/m²  (≈ σT⁴)',
       ha='right', x_label=0.545, y_label=0.295, fs=9)

# 6 — Back-radiation: downward, atmosphere top → surface
#     label at top tier, to the RIGHT
varrow(ax, 0.655, Y_ATM, Y_SRF + 0.005, C_LW_DN, lw=4,
       label=f'Back-radiation\n(greenhouse)\n{lw_srf_dn:.1f} W/m²',
       ha='left', x_label=0.670, y_label=0.700, fs=9)

# 8 — Latent heat: upward
#     label at bottom tier, to the LEFT
varrow(ax, 0.845, Y_SRF, Y_ATM, C_LE, lw=3,
       label=f'Latent heat\n{LE:.1f} W/m²\n({precip_mm:.2f} mm/day)',
       ha='right', x_label=0.828, y_label=0.295, fs=9)

# 9 — Sensible heat: upward
#     label at middle tier, to the RIGHT (short label to stay within bounds)
varrow(ax, 0.920, Y_SRF, Y_ATM, C_SH, lw=2.5,
       label=f'Sensible heat\n{SH:.1f} W/m²',
       ha='left', x_label=0.932, y_label=0.500, fs=9)

# --------------------------------------------------------------------------
# Floating annotation boxes inside the atmosphere band
# --------------------------------------------------------------------------
# SW absorbed by atmosphere
ax.text(0.255, 0.610,
        f'SW abs. by atmosphere\n{sw_atm_abs:.1f} W/m²',
        color=C_SW_ABS, fontsize=9, fontweight='bold', va='center', ha='center', zorder=6,
        bbox=dict(facecolor='white', edgecolor=C_SW_ABS, alpha=0.90,
                  pad=4, boxstyle='round,pad=0.3'))

# Net LW cooling of atmosphere
ax.text(0.720, 0.540,
        f'Net LW atm cooling\n{net_lw_atm:.1f} W/m²',
        color='#8a1818', fontsize=8.5, fontweight='bold', va='center', ha='center', zorder=6,
        bbox=dict(facecolor='white', edgecolor='#cc3030', alpha=0.88,
                  pad=3, boxstyle='round,pad=0.3'))

# --------------------------------------------------------------------------
# Surface band annotation
# --------------------------------------------------------------------------
ax.text(0.50, Y_SRF / 2,
        f'SW absorbed at surface:  {sw_srf_abs:.1f} W/m²',
        color='#e8a020', fontsize=9, fontweight='bold', va='center', ha='center', zorder=6,
        bbox=dict(facecolor='#4e3820', edgecolor='#e8a020', alpha=0.90,
                  pad=4, boxstyle='round,pad=0.3'))

# --------------------------------------------------------------------------
# Budget boxes
# --------------------------------------------------------------------------
# TOA (space zone, top-right)
ax.text(0.875, 0.965,
        f'TOA net (this step):   {toa_net:+.2f} W/m²\n'
        f'TOA net (window mean): −0.09 W/m²',
        fontsize=8.5, ha='center', va='top', color='white', zorder=7,
        bbox=dict(facecolor='#1a3a60', edgecolor='#6090cc',
                  alpha=0.92, pad=5, boxstyle='round,pad=0.4'))

# Surface balance (surface band, right side)
ax.text(0.985, Y_SRF - 0.005,
        f'Surface balance\n'
        f'  SW net    +{sw_srf_abs:.1f}\n'
        f'  LW net   {lw_srf_dn - lw_srf_up:+.1f}\n'
        f'  LE       {-LE:+.1f}\n'
        f'  SH       {-SH:+.1f}\n'
        f'  ─────────────\n'
        f'  Residual  {srf_net:+.2f} W/m²',
        fontsize=8.5, ha='right', va='top', color='white', zorder=7,
        family='monospace',
        bbox=dict(facecolor='#2a1a0a', edgecolor='#8a6030',
                  alpha=0.92, pad=5, boxstyle='round,pad=0.4'))

# --------------------------------------------------------------------------
# Title and caption
# --------------------------------------------------------------------------
ax.set_title(
    'ExoColumn  ·  Planetary Energy Balance\n'
    f'ZM convection (τ = 3600 s)  ·  α = 0.266  ·  Earth-like composition  ·  '
    f'$T_s$ = {ts:.2f} K',
    fontsize=12, fontweight='bold', pad=10,
)
ax.text(0.5, 0.003,
        'ExoColumn + ExoRT n68equiv  |  pver = 70  |  coszrs = 0.50  |  '
        'ps = 1000 hPa  |  Units: W m⁻²',
        ha='center', va='bottom', fontsize=8, color='#888888',
        transform=ax.transAxes)

# --------------------------------------------------------------------------
# Save
# --------------------------------------------------------------------------
with PdfPages(pdf_path) as pdf:
    pdf.savefig(fig, bbox_inches='tight', dpi=150)
plt.close(fig)
print(f"Wrote {pdf_path}")
