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
import os
import numpy as np
import netCDF4
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
matplotlib.rcParams.update({
    "font.size":       13,
    "axes.titlesize":  14,
    "axes.labelsize":  13,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
})
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D

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
nc_path      = sys.argv[1] if len(sys.argv) > 1 else "iofiles/exocol_out.nc"
pdf_path     = sys.argv[2] if len(sys.argv) > 2 else "iofiles/exocol_out.pdf"
konrad_path  = sys.argv[3] if len(sys.argv) > 3 else "iofiles/konrad_rce_ref.npz"

# ---------------------------------------------------------------------------
# Read data
# ---------------------------------------------------------------------------
ds = netCDF4.Dataset(nc_path)

pmid = ds["pmid"][:] / 100.0   # Pa → hPa
pint = ds["pint"][:] / 100.0   # Pa → hPa

tmid = ds["tmid"][:]

LWUP = ds["LWUP"][:]
LWDN = ds["LWDN"][:]
SWUP = ds["SWUP"][:]
SWDN = ds["SWDN"][:]

LWHR   = ds["LWHR"][:]
SWHR   = ds["SWHR"][:]

cond_heating = ds["cond_heating"][:]   # K/day, per layer

h2ommr  = ds["h2ommr"][:]
co2mmr  = ds["co2mmr"][:]
ch4mmr  = ds["ch4mmr"][:]
o2mmr   = ds["o2mmr"][:]
o3mmr   = ds["o3mmr"][:]
n2mmr   = ds["n2mmr"][:]
# Argon is radiatively inert and not written to the output file;
# compute from the known VMR (0.00934) and molecular weight (39.948 g/mol).
armmr   = np.full_like(n2mmr, 0.00934 * 39.948 / 28.97)

ts = float(ds["ts"][0])
ps = float(ds["ps"][0]) / 100.0

LE_diag = float(ds["LE"][0])
SH_diag = float(ds["SH"][0])
cp_col  = float(ds["cp"][0])           # J/kg/K
pint_pa = ds["pint"][:]                # Pa (needed for pdel)

ds.close()

# ---------------------------------------------------------------------------
# SH heating rate profile: applied only to the bottom model layer
# SH [W/m²] * g [m/s²] / (cp [J/kg/K] * Δp [Pa]) * 86400 → K/day
# ---------------------------------------------------------------------------
_g = 9.80665
SH_HR = np.zeros_like(LWHR)
pdel_bot = float(pint_pa[-1] - pint_pa[-2])   # Pa, bottom layer thickness
SH_HR[-1] = SH_diag * _g / (cp_col * pdel_bot) * 86400.0

# ---------------------------------------------------------------------------
# Figure layout: 2 rows × 2 cols
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(11, 12))
fig.suptitle(
    f"ExoColumn Earth calibration",
    fontsize=15, fontweight="bold"
)

gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.28, wspace=0.30,
                       top=0.91, bottom=0.06, left=0.08, right=0.97)

ax_T   = fig.add_subplot(gs[0, 0])
ax_fl  = fig.add_subplot(gs[0, 1])
ax_hr  = fig.add_subplot(gs[1, 0])
ax_q   = fig.add_subplot(gs[1, 1])

def setup_yaxis(ax, pmin=None, pmax=None):
    ax.set_yscale("log")
    ax.set_ylim(pmax or ps * 1.02, pmin or pmid.min() * 0.95)
    ax.set_ylabel("Pressure (hPa)")

# ---- Panel 1: Temperature profile -----------------------------------------
ax = ax_T
T_ref = _us_std_atm_T(pmid)
ax.plot(T_ref, pmid, color="dimgray",   lw=1.5, ls=":")
if os.path.isfile(konrad_path):
    _kd = np.load(konrad_path)
    _kT = np.append(_kd["T_K"], float(_kd["Ts_K"]))
    _kp = np.append(_kd["plev_hpa"], float(_kd["phlev_hpa"].max()))
    ax.plot(_kT, _kp, color="gray", lw=1.0, ls="-.")
ax.plot(np.append(tmid, ts), np.append(pmid, ps), color="firebrick", lw=1.5)
ax.axhline(ps, color="gray", lw=0.6, ls=":")

def _iT(T_arr, p_arr, p_target):
    return float(np.interp(np.log(p_target), np.log(p_arr[::-1]), T_arr[::-1]))

# Direct labels: (T_array, p_array, label, color, p_hPa, dx_K, ha)
_T_lbls = [(T_ref, pmid, "US Std Atm, 1976", "dimgray", 120.0, -51, "left"),
           (tmid,  pmid, "ExoColumn", "firebrick", 20.0, -30, "left")]
if os.path.isfile(konrad_path):
    _T_lbls.append((_kT, _kp, "konrad", "gray", 400.0, -30, "left"))
for T_arr, p_arr, lbl, color, p_lbl, dx, ha in _T_lbls:
    ax.text(_iT(T_arr, p_arr, p_lbl) + dx, p_lbl, lbl,
            color=color, ha=ha, va="center", fontsize=12)

ax.set_xlabel("Temperature (K)")
ax.set_title("Temperature profile")
setup_yaxis(ax)

# ---- Panel 2: Radiative fluxes --------------------------------------------
ax = ax_fl
ax.plot(LWUP, pint, color="lightcoral",  lw=1.5)
ax.plot(LWDN, pint, color="crimson",     lw=1.5)
ax.plot(SWDN, pint, color="steelblue",   lw=1.5)
ax.plot(SWUP, pint, color="cornflowerblue", lw=1.5)

def _iflux(flux, p_arr, p_target):
    return float(np.interp(np.log(p_target), np.log(p_arr[::-1]), flux[::-1]))

# (flux, label, color, p_hPa, dx W/m², ha)
# ExoColumn labels and konrad labels staggered in pressure so they don't collide.
_flx_lbls = [
    (LWUP, "LW↑",   "lightcoral",     30.0, -50, "left"),
    (LWDN, "LW↓",   "crimson",        30.0, +55, "right"),
    (SWDN, "SW↓",   "steelblue",      30.0, +50, "right"),
    (SWUP, "SW↑",   "cornflowerblue", 30.0, +10, "left"),
]
for flux, lbl, color, p_lbl, dx, ha in _flx_lbls:
    x = _iflux(flux, pint, p_lbl)
    ax.text(x + dx, p_lbl, lbl, color=color, ha=ha, va="center", fontsize=12)

_style_handles = [
    Line2D([0], [0], color="gray", lw=1.5, ls="-",  label="ExoColumn"),
    Line2D([0], [0], color="gray", lw=1.0, ls="-.", label="konrad"),
]
if os.path.isfile(konrad_path):
    _kf = np.load(konrad_path)
    if "lw_flxu" in _kf.files and np.any(np.isfinite(_kf["lw_flxu"])):
        _kfp = _kf["phlev_hpa"]
        ax.plot(_kf["lw_flxu"], _kfp, color="lightcoral",    lw=1.0, ls="-.")
        ax.plot(_kf["lw_flxd"], _kfp, color="crimson",       lw=1.0, ls="-.")
        ax.plot(_kf["sw_flxd"], _kfp, color="steelblue",     lw=1.0, ls="-.")
        ax.plot(_kf["sw_flxu"], _kfp, color="cornflowerblue", lw=1.0, ls="-.")
    ax.legend(handles=_style_handles, loc="center",
              bbox_to_anchor=(0.44, 0.92), bbox_transform=ax_fl.transAxes)

ax.set_xlabel("Flux (W m$^{-2}$)")
ax.set_title("Radiative fluxes")
setup_yaxis(ax)

# ---- Panel 3: Heating rates -----------------------------------------------
ax = ax_hr
ax.axvline(0, color="gray", lw=0.6, ls=":")
grandHR = LWHR + SWHR + cond_heating + SH_HR
ax.plot(LWHR,         pmid, color="tomato",       lw=1.5)
ax.plot(SWHR,         pmid, color="steelblue",    lw=1.5)
ax.plot(cond_heating, pmid, color="mediumorchid", lw=1.5)
ax.plot(SH_HR,        pmid, color="darkorange",   lw=1.5)
ax.plot(grandHR,      pmid, color="black",        lw=2.2)

def _ihr(data, p_target):
    return float(np.interp(np.log(p_target), np.log(pmid[::-1]), data[::-1]))

# (data, label, color, p_hPa, dx K/day, ha)
_hr_lbls = [
    (LWHR,         "Longwave\nradiation",         "tomato",       13.0,         +1.0, "left"),
    (SWHR,         "Shortwave\nradiation",         "steelblue",   13.0,         -5.0, "left"),
    (cond_heating, "Latent heating",  "mediumorchid", 200.0,        +1.3, "left"),
    (SH_HR,        "Sensible heating", "darkorange",   800, +1.3, "left"),
    (grandHR,      "Total",    "black",        5.0,        +0.3, "left"),
]
for data, lbl, color, p_lbl, dx, ha in _hr_lbls:
    x = _ihr(data, p_lbl)
    ax.text(x + dx, p_lbl, lbl, color=color, ha=ha, va="center", fontsize=12)

ax.set_xlabel("Heating rate (K day$^{-1}$)")
ax.set_title("Heating rates")
setup_yaxis(ax)

# ---- Panel 4: Atmospheric composition ----------------------------------------
ax = ax_q
_ppm = 1e6
# (mmr_array, label, color, label_p_hPa, side)
# p levels chosen so no two labels share the same y; offsets kept small so
# the label sits tight against the curve rather than floating mid-panel.
_gases = [
    # (mmr, label, color, p_hPa, ha,      mult)
    (n2mmr,  "N$_2$",   "dimgray",        100.0, "left",  1.5),
    (o2mmr,  "O$_2$",   "steelblue",       15.0, "right", 0.7),
    (armmr,  "Ar",      "mediumpurple",   100.0, "left",  1.5),
    (h2ommr, "H$_2$O",  "royalblue",      100.0, "left",  1.5),
    (h2ommr, "H$_2$O\n(konrad)",  "royalblue",      400.0, "left",  0.6),
    (co2mmr, "CO$_2$",  "tomato",          15.0, "left",  1.5),
    (o3mmr,  "O$_3$",   "mediumseagreen",  15.0, "left",  4.0),  # nudged right past peak
    (ch4mmr, "CH$_4$",  "darkorange",     300.0, "right", 0.7),
]
for mmr, lbl, color, p_lbl, ha, mult in _gases:
    if mmr.max() > 0:
        ax.plot(mmr * _ppm, pmid, color=color, lw=1.5)
        valid = mmr > 0
        x_lbl = np.exp(np.interp(np.log(p_lbl),
                                  np.log(pmid[valid][::-1]),
                                  np.log(mmr[valid][::-1] * _ppm)))
        ax.text(x_lbl * mult, p_lbl, lbl, color=color,
                ha=ha, va="center", fontsize=12)

if os.path.isfile(konrad_path):
    _kd = np.load(konrad_path)
    if "h2o_mmr" in _kd.files and np.any(np.isfinite(_kd["h2o_mmr"])):
        ax.plot(_kd["h2o_mmr"] * _ppm, _kd["plev_hpa"], color="royalblue",
                lw=1.0, ls="-.")
    if "o3_mmr" in _kd.files and np.any(np.isfinite(_kd["o3_mmr"])):
        ax.plot(_kd["o3_mmr"] * _ppm, _kd["plev_hpa"], color="mediumseagreen",
                lw=1.0, ls="-.")

ax.set_xscale("log")
# Clamp: konrad O3 underflows to ~1e-32 ppm near the surface; clamping keeps
# the other constituents visible.
ax.set_xlim(1e-3, 5e6)
ax.set_xlabel("Mass mixing ratio (ppm)")
ax.set_title("Atmospheric composition")
setup_yaxis(ax)


# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
with PdfPages(pdf_path) as pdf:
    pdf.savefig(fig)
plt.close(fig)

print(f"Wrote {pdf_path}")
