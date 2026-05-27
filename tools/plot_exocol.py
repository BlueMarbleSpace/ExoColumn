#!/usr/bin/env python3
"""
plot_exocol.py  —  Visualise ExoColumn RCE output.

Usage:
    python tools/plot_exocol.py [input.nc [output.pdf]]

Defaults:
    input  : iofiles/exocol_out.nc
    output : iofiles/exocol_out.pdf

Run from the project root.
"""

import sys
import numpy as np
import netCDF4
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages

# ---------------------------------------------------------------------------
# US Standard Atmosphere 1976 reference profile
# ---------------------------------------------------------------------------
def _us_std_atm_T(p_hpa):
    """Return USSA-1976 temperature (K) interpolated onto pressure levels (hPa)."""
    h_b = np.array([0,      11,     20,     32,     47,     51,     71    ])
    L_b = np.array([-6.5,   0.0,    1.0,    2.8,    0.0,   -2.8,   -2.0  ])
    T_b = np.array([288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65])
    P_b = np.array([1013.25, 226.32, 54.749, 8.6802, 1.1091, 0.66939, 0.039564])

    gMR = 9.80665 * 0.0289644 / 8.314462

    z = np.linspace(0, 86, 20000)
    T = np.empty_like(z)
    P = np.empty_like(z)
    for i, zi in enumerate(z):
        k = min(int(np.searchsorted(h_b, zi, side='right')) - 1, len(h_b) - 2)
        dz = zi - h_b[k]
        Ti = T_b[k] + L_b[k] * dz
        T[i] = Ti
        if L_b[k] == 0.0:
            P[i] = P_b[k] * np.exp(-gMR * dz * 1e3 / T_b[k])
        else:
            P[i] = P_b[k] * (T_b[k] / Ti) ** (gMR / (L_b[k] * 1e-3))

    return np.interp(np.log(p_hpa), np.log(P[::-1]), T[::-1])

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
nc_path  = sys.argv[1] if len(sys.argv) > 1 else "iofiles/exocol_out.nc"
pdf_path = sys.argv[2] if len(sys.argv) > 2 else "iofiles/exocol_out.pdf"

# ---------------------------------------------------------------------------
# Read data
# ---------------------------------------------------------------------------
ds = netCDF4.Dataset(nc_path)

pmid = ds["pmid"][:] / 100.0   # Pa → hPa
pint = ds["pint"][:] / 100.0   # Pa → hPa

tmid = ds["tmid"][:]
tint = ds["tint"][:]

LWUP = ds["LWUP"][:]
LWDN = ds["LWDN"][:]
SWUP = ds["SWUP"][:]
SWDN = ds["SWDN"][:]

LWHR   = ds["LWHR"][:]
SWHR   = ds["SWHR"][:]
totHR  = LWHR + SWHR

h2o_gkg = ds["h2ommr"][:] * 1e3

ts = float(ds["ts"][0])
ps = float(ds["ps"][0]) / 100.0

LE_diag = float(ds["LE"][0])
SH_diag = float(ds["SH"][0])

ds.close()

# ---------------------------------------------------------------------------
# Derived budget quantities (at interface levels)
# ---------------------------------------------------------------------------
SW_net  = SWDN - SWUP          # net downward SW at each interface [W/m²]
LW_net  = LWDN - LWUP          # net downward LW at each interface [W/m²]
Rad_net = SW_net + LW_net      # total radiative net (the "sum" of rad terms)

Rad_toa = float(Rad_net[0])    # TOA: ASR − OLR  (≈ 0 at equilibrium)
Rad_srf = float(Rad_net[-1])   # Surface: SW_net + LW_net before turbulent fluxes
residual = Rad_srf - LE_diag - SH_diag   # full surface budget residual (≈ 0 at equilibrium)

# ---------------------------------------------------------------------------
# Figure layout: 3 rows × 2 cols; bottom panel spans both columns
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(11, 12))
fig.suptitle(
    f"ExoColumn RCE equilibrium  |  Ts = {ts:.2f} K,  ps = {ps:.1f} hPa",
    fontsize=13, fontweight="bold"
)

gs = gridspec.GridSpec(3, 2, figure=fig, hspace=0.40, wspace=0.30,
                       top=0.94, bottom=0.06, left=0.08, right=0.97)

ax_T   = fig.add_subplot(gs[0, 0])
ax_fl  = fig.add_subplot(gs[0, 1])
ax_hr  = fig.add_subplot(gs[1, 0])
ax_q   = fig.add_subplot(gs[1, 1])
ax_bgt = fig.add_subplot(gs[2, :])   # energy budget — full width

def setup_yaxis(ax, pmin=None, pmax=None):
    ax.set_yscale("log")
    ax.set_ylim(pmax or pmid.max() * 1.05, pmin or pmid.min() * 0.95)
    ax.set_ylabel("Pressure (hPa)")

# ---- Panel 1: Temperature profile -----------------------------------------
ax = ax_T
T_ref = _us_std_atm_T(pmid)
ax.plot(T_ref, pmid, color="gray",     lw=1.2, ls=":", label="USSA-1976")
ax.plot(tmid,  pmid, color="firebrick",lw=1.5, label="T$_{mid}$")
ax.plot(tint,  pint, color="salmon",   lw=0.8, ls="--", label="T$_{int}$")
ax.axhline(ps, color="gray", lw=0.6, ls=":")
ax.set_xlabel("Temperature (K)")
ax.set_title("Temperature profile")
ax.legend(fontsize=9)
setup_yaxis(ax)

# ---- Panel 2: Radiative fluxes --------------------------------------------
ax = ax_fl
ax.plot(LWUP, pint, color="tomato",          lw=1.5, label="LW↑")
ax.plot(LWDN, pint, color="salmon",          lw=1.5, ls="--", label="LW↓")
ax.plot(SWDN, pint, color="steelblue",       lw=1.5, label="SW↓")
ax.plot(SWUP, pint, color="cornflowerblue",  lw=1.5, ls="--", label="SW↑")
ax.set_xlabel("Flux (W m$^{-2}$)")
ax.set_title("Radiative fluxes")
ax.legend(fontsize=9)
setup_yaxis(ax)

# ---- Panel 3: Heating rates -----------------------------------------------
ax = ax_hr
ax.axvline(0, color="gray", lw=0.6, ls=":")
ax.plot(LWHR,  pmid, color="tomato",    lw=1.5, label="LW")
ax.plot(SWHR,  pmid, color="steelblue", lw=1.5, label="SW")
ax.plot(totHR, pmid, color="black",     lw=1.8, label="Total")
ax.set_xlabel("Heating rate (K day$^{-1}$)")
ax.set_title("Radiative heating rates")
ax.legend(fontsize=9)
setup_yaxis(ax)

# ---- Panel 4: Water vapour ------------------------------------------------
ax = ax_q
ax.plot(h2o_gkg, pmid, color="royalblue", lw=1.5)
ax.set_xlabel("Specific humidity (g kg$^{-1}$)")
ax.set_title("Water vapour")
setup_yaxis(ax)

# ---- Panel 5: Column energy budget ----------------------------------------
# At each interface, the net downward radiative flux is:
#   SW_net = SWDN − SWUP   (positive downward)
#   LW_net = LWDN − LWUP   (negative = net upward in most of column)
#   Rad_net = SW_net + LW_net  (the sum of all radiative terms)
#
# At the surface, turbulent fluxes close the budget:
#   Rad_srf − LE − SH  =  surface residual  (≈ 0 at equilibrium)
#
# TOA balance requires Rad_toa ≈ 0.
ax = ax_bgt
ax.axvline(0, color="gray", lw=0.6, ls=":", zorder=0)
ax.axhline(pint[-1], color="gray", lw=0.5, ls="--", alpha=0.4, zorder=0)

# Three main profiles
ax.plot(SW_net,  pint, color="steelblue", lw=1.5,
        label="SW net  (SWDN − SWUP)")
ax.plot(LW_net,  pint, color="tomato",   lw=1.5,
        label="LW net  (LWDN − LWUP)")
ax.plot(Rad_net, pint, color="black",    lw=2.2,
        label="Rad net  = SW + LW  ← sum of rad terms")

# Surface budget decomposition:
#   start at Rad_srf, subtract LE, subtract SH, arrive at residual
p_srf = pint[-1]
x0 = Rad_srf                    # radiative input to surface
x1 = x0 - LE_diag              # after latent heat flux leaves
x2 = x1 - SH_diag              # after sensible heat flux leaves = residual

# Connecting line at surface level showing the budget steps
ax.plot([x0, x1, x2], [p_srf, p_srf, p_srf],
        color="gray", lw=1.2, ls="-", zorder=4)

# Markers at each budget step
ax.scatter([x0], [p_srf], color="black",      s=55, zorder=6,
           label=f"Rad$_{{srf}}$ = {x0:+.1f} W m$^{{-2}}$")
ax.scatter([x1], [p_srf], color="mediumorchid", marker="s", s=55, zorder=6,
           label=f"− LE ({LE_diag:.0f} W m$^{{-2}}$) → {x1:+.1f}")
ax.scatter([x2], [p_srf], color="darkorange",  marker="^", s=70, zorder=6,
           label=f"− SH ({SH_diag:.0f} W m$^{{-2}}$) → residual = {x2:+.1f}")

# TOA annotation
ax.scatter([Rad_toa], [pint[0]], color="navy", marker="*", s=100, zorder=6,
           label=f"TOA imbalance = {Rad_toa:+.2f} W m$^{{-2}}$")

ax.set_xlabel("Net downward flux (W m$^{-2}$)")
ax.set_title(
    "Column energy budget — net fluxes at each interface  "
    f"[TOA = {Rad_toa:+.2f}  |  Surface residual = {residual:+.1f} W m$^{{-2}}$]"
)
ax.legend(fontsize=8.5, ncol=2, loc="upper left")
setup_yaxis(ax)

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
with PdfPages(pdf_path) as pdf:
    pdf.savefig(fig)
plt.close(fig)

print(f"Wrote {pdf_path}")
