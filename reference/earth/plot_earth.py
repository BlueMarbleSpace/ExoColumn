#!/usr/bin/env python3
"""
plot_earth.py — Publication figure for the ExoColumn Earth reference case.

Produces a 4-panel PDF:
  Panel 1: Temperature profile  (ExoColumn, konrad, CLIMA, US Std Atm 1976)
  Panel 2: Radiative fluxes     (ExoColumn, konrad, CLIMA)
  Panel 3: Heating rates        (ExoColumn only)
  Panel 4: Atmospheric composition (ExoColumn + konrad O3/H2O + CLIMA O3)

Run from any directory:
    python /path/to/plot_earth.py

All input paths are resolved relative to this script's own directory so the
figure can be reproduced from the reference snapshot without any arguments.
"""

import os
import numpy as np
import netCDF4
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D

# ---------------------------------------------------------------------------
# Paths — all relative to this script's directory
# ---------------------------------------------------------------------------
HERE      = os.path.dirname(os.path.abspath(__file__))
NC_FILE   = os.path.join(HERE, "exocol_out.nc")
KONRAD    = os.path.join(HERE, "konrad_rce_ref.npz")
CLIMA     = os.path.join(HERE, "clima_allout_presentEARTH.tab")
OUT_PDF   = os.path.join(HERE, "earth_rce.pdf")
OUT_PNG   = os.path.join(HERE, "earth_rce.png")

# ---------------------------------------------------------------------------
# Figure style
# ---------------------------------------------------------------------------
matplotlib.rcParams.update({
    "font.family":     "sans-serif",
    "font.size":        13,
    "axes.titlesize":   14,
    "axes.labelsize":   13,
    "xtick.labelsize":  12,
    "ytick.labelsize":  12,
    "legend.fontsize":  12,
    "figure.facecolor": "white",
    "axes.facecolor":   "white",
    "savefig.facecolor":"white",
})

FIGSIZE = (11, 12)    # inches
DPI     = 300

# Model line styles
LS_EXOCOL = dict(lw=1.5, ls="-")
LS_KONRAD = dict(lw=1.0, ls="-.")
LS_CLIMA  = dict(lw=1.2, ls="--")

# ---------------------------------------------------------------------------
# US Standard Atmosphere 1976
# ---------------------------------------------------------------------------
def _us_std_atm_T(p_hpa):
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
# CLIMA allout parser (Kasting/Kopparapu)
# ---------------------------------------------------------------------------
def read_clima(path):
    """Parse the last (converged) snapshot from a CLIMA allout file.

    Profile columns (midpoints, 101 levels):
        J  P[atm]  ALT[km]  T[K]  CONVEC  DT  TOLD  FH2O[VMR]  FSAVE  FO3[VMR]  TCOOL  THEAT

    Flux columns (flux levels, 101 levels):
        J  PF[atm]  ALT[km]  FTOTAL  FTIR  FDNIR  FUPIR  FTSOL  FDNSOL  FUPSOL  DIVF
        (all fluxes in erg cm⁻² s⁻¹; divide by 1000 → W m⁻²)

    Returns dict of arrays on the CLIMA grid:
        p_hpa   : midpoint pressure (hPa)
        T       : temperature (K)
        h2o_mmr : H2O mass mixing ratio (kg/kg)
        o3_mmr  : O3 mass mixing ratio (kg/kg)
        pf_hpa  : flux-level pressure (hPa)
        lw_flxd : LW downward flux (W/m²)
        lw_flxu : LW upward flux (W/m²)
        sw_flxd : SW downward flux (W/m²)
        sw_flxu : SW upward flux (W/m²)
    """
    ATM_TO_HPA = 1013.25
    ERG_TO_WM2 = 1e-3      # erg cm⁻² s⁻¹  →  W m⁻²
    MW_H2O, MW_O3, MW_AIR = 18.015, 47.998, 28.97

    with open(path) as f:
        lines = f.readlines()

    # Use the last converged snapshot
    last_time_idx = max(i for i, l in enumerate(lines) if "TIME=" in l)

    def find_header(keyword, start):
        for i in range(start, len(lines)):
            if keyword in lines[i]:
                return i
        raise ValueError(f"Header '{keyword}' not found after line {start}")

    def parse_block(hdr_idx, ncols):
        rows = []
        i = hdr_idx + 1
        while i < len(lines):
            s = lines[i].strip()
            if not s or not s[0].isdigit():
                break
            vals = s.split()
            if len(vals) >= ncols:
                rows.append([float(v) for v in vals[:ncols]])
            i += 1
        return np.array(rows)

    prof_hdr = find_header("  J     P ", last_time_idx)
    flux_hdr = find_header("  J    PF ", prof_hdr + 1)

    prof = parse_block(prof_hdr, 10)
    flux = parse_block(flux_hdr, 10)

    return dict(
        p_hpa   = prof[:, 1] * ATM_TO_HPA,
        T       = prof[:, 3],
        h2o_mmr = prof[:, 7] * (MW_H2O / MW_AIR),
        o3_mmr  = prof[:, 9] * (MW_O3  / MW_AIR),
        pf_hpa  = flux[:, 1] * ATM_TO_HPA,
        lw_flxd = flux[:, 5] * ERG_TO_WM2,   # FDNIR
        lw_flxu = flux[:, 6] * ERG_TO_WM2,   # FUPIR
        sw_flxd = flux[:, 8] * ERG_TO_WM2,   # FDNSOL
        sw_flxu = flux[:, 9] * ERG_TO_WM2,   # FUPSOL
    )

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
ds = netCDF4.Dataset(NC_FILE)

pmid = ds["pmid"][:] / 100.0
pint = ds["pint"][:] / 100.0
tmid = ds["tmid"][:]

LWUP = ds["LWUP"][:]
LWDN = ds["LWDN"][:]
SWUP = ds["SWUP"][:]
SWDN = ds["SWDN"][:]
LWHR = ds["LWHR"][:]
SWHR = ds["SWHR"][:]
cond_heating = ds["cond_heating"][:]

h2ommr = ds["h2ommr"][:]
co2mmr = ds["co2mmr"][:]
ch4mmr = ds["ch4mmr"][:]
o2mmr  = ds["o2mmr"][:]
o3mmr  = ds["o3mmr"][:]
n2mmr  = ds["n2mmr"][:]
armmr  = np.full_like(n2mmr, 0.00934 * 39.948 / 28.97)

ts      = float(ds["ts"][0])
ps      = float(ds["ps"][0]) / 100.0
LE_diag = float(ds["LE"][0])
SH_diag = float(ds["SH"][0])
cp_col  = float(ds["cp"][0])
pint_pa = ds["pint"][:]
ds.close()

# SH heating rate — applied to bottom layer only
_g = 9.80665
SH_HR = np.zeros_like(LWHR)
SH_HR[-1] = SH_diag * _g / (cp_col * float(pint_pa[-1] - pint_pa[-2])) * 86400.0

# Konrad reference
_kd = np.load(KONRAD) if os.path.isfile(KONRAD) else None
if _kd is not None:
    _kT  = np.concatenate([[float(_kd["Ts_K"])], _kd["T_K"]])
    _kp  = np.concatenate([[float(_kd["phlev_hpa"].max())], _kd["plev_hpa"]])

# CLIMA reference
_cl = read_clima(CLIMA) if os.path.isfile(CLIMA) else None

# ---------------------------------------------------------------------------
# Figure layout
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=FIGSIZE)

gs = gridspec.GridSpec(2, 2, figure=fig,
                       hspace=0.28, wspace=0.30,
                       top=0.96, bottom=0.06, left=0.08, right=0.97)
ax_T  = fig.add_subplot(gs[0, 0])
ax_fl = fig.add_subplot(gs[0, 1])
ax_hr = fig.add_subplot(gs[1, 0])
ax_q  = fig.add_subplot(gs[1, 1])

def setup_yaxis(ax, pmin=None, pmax=None):
    ax.set_yscale("log")
    ax.set_ylim(pmax or ps * 1.02, pmin or pmid.min() * 0.95)
    ax.set_ylabel("Pressure (hPa)")

def _iT(T_arr, p_arr, p_target):
    return float(np.interp(np.log(p_target), np.log(p_arr[::-1]), T_arr[::-1]))

def _iflux(flux, p_arr, p_target):
    return float(np.interp(np.log(p_target), np.log(p_arr[::-1]), flux[::-1]))

def _ihr(data, p_target):
    return float(np.interp(np.log(p_target), np.log(pmid[::-1]), data[::-1]))

# ---------------------------------------------------------------------------
# Panel 1: Temperature profile
# ---------------------------------------------------------------------------
ax = ax_T

T_ref = _us_std_atm_T(pmid)
ax.plot(T_ref, pmid, color="dimgray", lw=1.5, ls=":")

if _kd is not None:
    ax.plot(_kT, _kp, color="royalblue", **LS_KONRAD)

if _cl is not None:
    ax.plot(_cl["T"], _cl["p_hpa"], color="seagreen", **LS_CLIMA)

ax.plot(np.append(tmid, ts), np.append(pmid, ps), color="firebrick", **LS_EXOCOL)
ax.axhline(ps, color="gray", lw=0.6, ls=":")

_T_lbls = [
    (T_ref, pmid,        "US Std Atm, 1976", "dimgray",   120.0, -51, "left"),
    (tmid,  pmid,        "ExoColumn",        "firebrick",  10.0, -22, "left"),
]
if _kd is not None:
    _T_lbls.append((_kT, _kp, "konrad", "royalblue", 400.0, -30, "left"))
if _cl is not None:
    _T_lbls.append((_cl["T"], _cl["p_hpa"], "CLIMA", "seagreen", 20.0, -15, "left"))

for T_arr, p_arr, lbl, color, p_lbl, dx, ha in _T_lbls:
    ax.text(_iT(T_arr, p_arr, p_lbl) + dx, p_lbl, lbl,
            color=color, ha=ha, va="center", fontsize=12)

ax.set_xlabel("Temperature (K)")
ax.set_title("Temperature profile")
setup_yaxis(ax)

# ---------------------------------------------------------------------------
# Panel 2: Radiative fluxes
# ---------------------------------------------------------------------------
ax = ax_fl

ax.plot(LWUP, pint, color="lightcoral",      **LS_EXOCOL)
ax.plot(LWDN, pint, color="crimson",         **LS_EXOCOL)
ax.plot(SWDN, pint, color="steelblue",       **LS_EXOCOL)
ax.plot(SWUP, pint, color="cornflowerblue",  **LS_EXOCOL)

if _kd is not None and "lw_flxu" in _kd.files and np.any(np.isfinite(_kd["lw_flxu"])):
    _kfp = _kd["phlev_hpa"]
    ax.plot(_kd["lw_flxu"], _kfp, color="lightcoral",     **LS_KONRAD)
    ax.plot(_kd["lw_flxd"], _kfp, color="crimson",        **LS_KONRAD)
    ax.plot(_kd["sw_flxd"], _kfp, color="steelblue",      **LS_KONRAD)
    ax.plot(_kd["sw_flxu"], _kfp, color="cornflowerblue", **LS_KONRAD)

if _cl is not None:
    ax.plot(_cl["lw_flxu"], _cl["pf_hpa"], color="lightcoral",     **LS_CLIMA)
    ax.plot(_cl["lw_flxd"], _cl["pf_hpa"], color="crimson",        **LS_CLIMA)
    ax.plot(_cl["sw_flxd"], _cl["pf_hpa"], color="steelblue",      **LS_CLIMA)
    ax.plot(_cl["sw_flxu"], _cl["pf_hpa"], color="cornflowerblue", **LS_CLIMA)

# Quantity labels (positioned on ExoColumn curves)
_flx_lbls = [
    (LWUP, pint, "LW↑", "lightcoral",     30.0, -50, "left"),
    (LWDN, pint, "LW↓", "crimson",        30.0, +55, "right"),
    (SWDN, pint, "SW↓", "steelblue",      30.0, +50, "right"),
    (SWUP, pint, "SW↑", "cornflowerblue", 30.0, +10, "left"),
]
for flux, p_arr, lbl, color, p_lbl, dx, ha in _flx_lbls:
    ax.text(_iflux(flux, p_arr, p_lbl) + dx, p_lbl, lbl,
            color=color, ha=ha, va="center", fontsize=12)

# Model style legend
_legend_handles = [
    Line2D([0], [0], color="gray", **LS_EXOCOL, label="ExoColumn"),
]
if _kd is not None:
    _legend_handles.append(Line2D([0], [0], color="gray", **LS_KONRAD, label="konrad"))
if _cl is not None:
    _legend_handles.append(Line2D([0], [0], color="gray", **LS_CLIMA, label="CLIMA"))
ax.legend(handles=_legend_handles, loc="center",
          bbox_to_anchor=(0.44, 0.90), bbox_transform=ax_fl.transAxes)

ax.set_xlabel("Flux (W m$^{-2}$)")
ax.set_title("Radiative fluxes")
setup_yaxis(ax)

# ---------------------------------------------------------------------------
# Panel 3: Heating rates
# ---------------------------------------------------------------------------
ax = ax_hr

ax.axvline(0, color="gray", lw=0.6, ls=":")
grandHR = LWHR + SWHR + cond_heating + SH_HR
ax.plot(LWHR,         pmid, color="tomato",       lw=1.5)
ax.plot(SWHR,         pmid, color="steelblue",    lw=1.5)
ax.plot(cond_heating, pmid, color="mediumorchid", lw=1.5)
ax.plot(SH_HR,        pmid, color="darkorange",   lw=1.5)
ax.plot(grandHR,      pmid, color="black",        lw=2.2)

_hr_lbls = [
    (LWHR,         "Longwave\nradiation",  "tomato",       13.0,  +1.0, "left"),
    (SWHR,         "Shortwave\nradiation", "steelblue",    13.0,  -5.0, "left"),
    (cond_heating, "Latent heating",       "mediumorchid", 200.0, +1.3, "left"),
    (SH_HR,        "Sensible heating",     "darkorange",   800.0, +1.3, "left"),
    (grandHR,      "Total",               "black",          5.0, +0.3, "left"),
]
for data, lbl, color, p_lbl, dx, ha in _hr_lbls:
    ax.text(_ihr(data, p_lbl) + dx, p_lbl, lbl,
            color=color, ha=ha, va="center", fontsize=12)

ax.set_xlabel("Heating rate (K day$^{-1}$)")
ax.set_title("Heating rates")
setup_yaxis(ax)

# ---------------------------------------------------------------------------
# Panel 4: Atmospheric composition
# ---------------------------------------------------------------------------
ax = ax_q
_ppm = 1e6

_gases = [
    (n2mmr,  "N$_2$",  "dimgray",        100.0, "left",  1.5),
    (o2mmr,  "O$_2$",  "steelblue",       15.0, "right", 0.7),
    (armmr,  "Ar",     "mediumpurple",   100.0, "left",  1.5),
    (h2ommr, "H$_2$O", "royalblue",      100.0, "left",  1.5),
    (co2mmr, "CO$_2$", "tomato",          15.0, "left",  1.5),
    (o3mmr,  "O$_3$",  "mediumseagreen",  15.0, "left",  4.0),
    (ch4mmr, "CH$_4$", "darkorange",      15.0, "right", 0.7),
]
for mmr, lbl, color, p_lbl, ha, mult in _gases:
    if mmr.max() > 0:
        ax.plot(mmr * _ppm, pmid, color=color, **LS_EXOCOL)
        valid = mmr > 0
        x_lbl = np.exp(np.interp(np.log(p_lbl),
                                  np.log(pmid[valid][::-1]),
                                  np.log(mmr[valid][::-1] * _ppm)))
        ax.text(x_lbl * mult, p_lbl, lbl, color=color,
                ha=ha, va="center", fontsize=12)

if _kd is not None:
    if "h2o_mmr" in _kd.files and np.any(np.isfinite(_kd["h2o_mmr"])):
        ax.plot(_kd["h2o_mmr"] * _ppm, _kd["plev_hpa"],
                color="royalblue", **LS_KONRAD)
    if "o3_mmr" in _kd.files and np.any(np.isfinite(_kd["o3_mmr"])):
        ax.plot(_kd["o3_mmr"] * _ppm, _kd["plev_hpa"],
                color="mediumseagreen", **LS_KONRAD)

if _cl is not None:
    ax.plot(_cl["h2o_mmr"] * _ppm, _cl["p_hpa"], color="royalblue", **LS_CLIMA)
    ax.plot(_cl["o3_mmr"]  * _ppm, _cl["p_hpa"], color="mediumseagreen", **LS_CLIMA)
    # Label for CLIMA H2O — positioned using CLIMA's own profile at 700 hPa
    _p_lbl = 700.0
    _valid  = _cl["h2o_mmr"] > 0
    _x_lbl  = np.exp(np.interp(np.log(_p_lbl),
                                np.log(_cl["p_hpa"][_valid][::-1]),
                                np.log(_cl["h2o_mmr"][_valid][::-1] * _ppm)))

ax.set_xscale("log")
ax.set_xlim(1e-3, 5e6)
ax.set_xlabel("Mass mixing ratio (ppm)")
ax.set_title("Atmospheric composition")
setup_yaxis(ax)

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
with PdfPages(OUT_PDF) as pdf:
    pdf.savefig(fig, dpi=DPI)
fig.savefig(OUT_PNG, dpi=DPI)
plt.close(fig)
print(f"Wrote {OUT_PDF}")
print(f"Wrote {OUT_PNG}")
