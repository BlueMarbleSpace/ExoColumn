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
from matplotlib.backends.backend_pdf import PdfPages

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
nc_path  = sys.argv[1] if len(sys.argv) > 1 else "iofiles/exocol_out.nc"
pdf_path = sys.argv[2] if len(sys.argv) > 2 else "iofiles/exocol_out.pdf"

# ---------------------------------------------------------------------------
# Read data
# ---------------------------------------------------------------------------
ds = netCDF4.Dataset(nc_path)

pmid = ds["pmid"][:] / 100.0   # Pa → hPa, shape (pver,)
pint = ds["pint"][:] / 100.0   # Pa → hPa, shape (pverp,)

tmid = ds["tmid"][:]            # K
tint = ds["tint"][:]            # K

LWUP = ds["LWUP"][:]            # W/m²
LWDN = ds["LWDN"][:]
SWUP = ds["SWUP"][:]
SWDN = ds["SWDN"][:]

LWHR = ds["LWHR"][:]            # K/day
SWHR = ds["SWHR"][:]
totHR = LWHR + SWHR

h2o_gkg = ds["h2ommr"][:] * 1e3   # kg/kg → g/kg

ts = float(ds["ts"][0])
ps = float(ds["ps"][0]) / 100.0    # hPa

ds.close()

# Index 0 = TOA, index -1 = surface (highest pressure).
# For vertical profiles plotted against log-p: surface at bottom → invert y-axis.

# ---------------------------------------------------------------------------
# Figure layout: 2 × 2
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))
fig.suptitle(
    f"ExoColumn RCE equilibrium  |  Ts = {ts:.2f} K,  ps = {ps:.1f} hPa",
    fontsize=13, fontweight="bold"
)

def setup_yaxis(ax, pmin=None, pmax=None):
    """Log-pressure y-axis, surface at bottom."""
    ax.set_yscale("log")
    ax.set_ylim(pmax or pmid.max() * 1.05, pmin or pmid.min() * 0.95)
    ax.set_ylabel("Pressure (hPa)")

# ---- Panel 1: Temperature profile -----------------------------------------
ax = axes[0, 0]
ax.plot(tmid, pmid, color="firebrick", lw=1.5, label="T$_{mid}$")
ax.plot(tint, pint, color="salmon",    lw=0.8, ls="--", label="T$_{int}$")
ax.axhline(ps, color="gray", lw=0.6, ls=":")
ax.set_xlabel("Temperature (K)")
ax.set_title("Temperature profile")
ax.legend(fontsize=9)
setup_yaxis(ax)

# ---- Panel 2: Radiative fluxes --------------------------------------------
ax = axes[0, 1]
ax.plot(LWUP, pint, color="tomato",      lw=1.5, label="LW↑")
ax.plot(LWDN, pint, color="salmon",      lw=1.5, ls="--", label="LW↓")
ax.plot(SWDN, pint, color="steelblue",   lw=1.5, label="SW↓")
ax.plot(SWUP, pint, color="cornflowerblue", lw=1.5, ls="--", label="SW↑")
ax.set_xlabel("Flux (W m$^{-2}$)")
ax.set_title("Radiative fluxes")
ax.legend(fontsize=9)
setup_yaxis(ax)

# ---- Panel 3: Heating rates -----------------------------------------------
ax = axes[1, 0]
ax.axvline(0, color="gray", lw=0.6, ls=":")
ax.plot(LWHR,  pmid, color="tomato",    lw=1.5, label="LW")
ax.plot(SWHR,  pmid, color="steelblue", lw=1.5, label="SW")
ax.plot(totHR, pmid, color="black",     lw=1.8, label="Total")
ax.set_xlabel("Heating rate (K day$^{-1}$)")
ax.set_title("Radiative heating rates")
ax.legend(fontsize=9)
setup_yaxis(ax)

# ---- Panel 4: Water vapour ------------------------------------------------
ax = axes[1, 1]
ax.plot(h2o_gkg, pmid, color="royalblue", lw=1.5)
ax.set_xlabel("Specific humidity (g kg$^{-1}$)")
ax.set_title("Water vapour")
setup_yaxis(ax)

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
fig.tight_layout(rect=[0, 0, 1, 0.95])
with PdfPages(pdf_path) as pdf:
    pdf.savefig(fig)
plt.close(fig)

print(f"Wrote {pdf_path}")
