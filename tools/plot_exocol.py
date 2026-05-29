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
# Saturated moist adiabat
# ---------------------------------------------------------------------------
def _moist_adiabat_T(ts_K, ps_hpa, p_hpa):
    """
    Saturated moist adiabat from (ts_K, ps_hpa) onto pressure levels p_hpa.
    Clausius-Clapeyron with phase-aware latent heat (Lv above 273.16 K, Ls below).
    Integrated in pressure coordinates via RK4. Returns NaN above ps.
    """
    Rd = 287.058; Rv = 461.5; cp = 1005.0
    Lv = 2.501e6; Ls = 2.834e6; eps = Rd / Rv

    def _dtdp(p_pa, T):
        L  = Lv if T >= 273.16 else Ls
        es = 611.2 * np.exp(L / Rv * (1.0/273.16 - 1.0/T))
        qs = eps * es / max(p_pa - (1.0 - eps) * es, 1.0)
        return (Rd * T / p_pa * (1.0 + L*qs/(Rd*T))) / (cp + L**2*qs/(Rv*T**2))

    n = 4000
    p0 = ps_hpa * 100.0
    p1 = float(np.min(p_hpa)) * 100.0
    p_path = np.linspace(p0, p1, n)
    T_path = np.empty(n); T_path[0] = ts_K
    for i in range(1, n):
        dp = p_path[i] - p_path[i-1]
        pi, Ti = p_path[i-1], T_path[i-1]
        k1 = _dtdp(pi,        Ti)
        k2 = _dtdp(pi+dp/2,   Ti+k1*dp/2)
        k3 = _dtdp(pi+dp/2,   Ti+k2*dp/2)
        k4 = _dtdp(pi+dp,     Ti+k3*dp)
        T_path[i] = Ti + dp*(k1 + 2*k2 + 2*k3 + k4)/6.0

    out = np.full(len(p_hpa), np.nan)
    for i, p in enumerate(p_hpa):
        if p * 100.0 <= p0:
            out[i] = np.interp(p*100.0, p_path[::-1], T_path[::-1])
    return out

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
tint = ds["tint"][:]

LWUP = ds["LWUP"][:]
LWDN = ds["LWDN"][:]
SWUP = ds["SWUP"][:]
SWDN = ds["SWDN"][:]

LWHR   = ds["LWHR"][:]
SWHR   = ds["SWHR"][:]
totHR  = LWHR + SWHR

cond_heating = ds["cond_heating"][:]   # K/day, per layer

h2ommr  = ds["h2ommr"][:]
co2mmr  = ds["co2mmr"][:]
ch4mmr  = ds["ch4mmr"][:]
o2mmr   = ds["o2mmr"][:]
o3mmr   = ds["o3mmr"][:]
n2mmr   = ds["n2mmr"][:]

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
SH_HR = np.zeros_like(totHR)
pdel_bot = float(pint_pa[-1] - pint_pa[-2])   # Pa, bottom layer thickness
SH_HR[-1] = SH_diag * _g / (cp_col * pdel_bot) * 86400.0

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
T_ref   = _us_std_atm_T(pmid)
T_moist = _moist_adiabat_T(ts, ps, pmid)
ax.plot(T_ref,   pmid, color="gray",          lw=1.2, ls=":",  label="USSA-1976")
ax.plot(T_moist, pmid, color="mediumseagreen", lw=1.2, ls="--", label="Moist adiabat")
if os.path.isfile(konrad_path):
    _kd = np.load(konrad_path)
    ax.plot(_kd["T_K"], _kd["plev_hpa"], color="steelblue", lw=1.2, ls="-.",
            label="konrad/RRTMG")
ax.plot(tmid,    pmid, color="firebrick",      lw=1.5, label="T$_{mid}$")
ax.plot(tint,  pint, color="salmon",   lw=0.8, ls="--", label="T$_{int}$")
ax.axhline(ps, color="gray", lw=0.6, ls=":")
ax.set_xlabel("Temperature (K)")
ax.set_xlim(left=175)
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
grandHR = LWHR + SWHR + cond_heating + SH_HR
ax.plot(LWHR,         pmid, color="tomato",       lw=1.5, label="LW rad")
ax.plot(SWHR,         pmid, color="steelblue",    lw=1.5, label="SW rad")
ax.plot(totHR,        pmid, color="black",        lw=1.8, ls="--", label="Rad total")
ax.plot(cond_heating, pmid, color="mediumorchid", lw=1.5, label="Cond (latent)")
ax.plot(SH_HR,        pmid, color="darkorange",   lw=1.5, label="Sensible (sfc)")
ax.plot(grandHR,      pmid, color="black",        lw=2.2, label="Grand total")
ax.set_xlabel("Heating rate (K day$^{-1}$)")
ax.set_title("Heating rates")
ax.legend(fontsize=9)
setup_yaxis(ax)

# ---- Panel 4: Atmospheric composition ----------------------------------------
ax = ax_q
_gases = [
    (n2mmr,  "N$_2$",                 "dimgray"),
    (o2mmr,  "O$_2$",                 "steelblue"),
    (h2ommr, "H$_2$O (moist-air)",    "royalblue"),
    (co2mmr, "CO$_2$",                "tomato"),
    (o3mmr,  "O$_3$",                 "mediumseagreen"),
    (ch4mmr, "CH$_4$",                "darkorange"),
]
for mmr, label, color in _gases:
    if mmr.max() > 0:
        ax.plot(mmr, pmid, color=color, lw=1.5, label=label)
ax.set_xscale("log")
ax.set_xlabel("Mass mixing ratio (kg kg$^{-1}$)")
ax.set_title("Atmospheric composition")
ax.legend(fontsize=8.5, loc="lower left")
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
