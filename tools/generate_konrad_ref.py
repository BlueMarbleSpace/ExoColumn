#!/usr/bin/env python3
"""
generate_konrad_ref.py  —  Run konrad to RCE equilibrium and save reference profile.

Usage (from project root, konrad conda env):
    conda activate konrad
    python tools/generate_konrad_ref.py [Ts_K [co2_ppm [output.npz]]]

Defaults:
    Ts_K      : 288.0   (match ExoColumn Earth-like equilibrium)
    co2_ppm   : 400.0   (modern Earth)
    output    : iofiles/konrad_rce_ref.npz

The output .npz contains:
    plev_hpa  : pressure levels (hPa), shape (N,)
    T_K       : equilibrium temperature (K), shape (N,)
    z_km      : geopotential height (km), shape (N,)
    rh        : relative humidity (fraction), shape (N,)
    Ts_K      : surface temperature used (scalar)
    co2_ppm   : CO2 concentration used (scalar)
    ch4_ppm   : CH4 concentration used (scalar)

Comparison fidelity vs ExoColumn
--------------------------------
This reference is configured to be as apples-to-apples with ExoColumn's
canonical cold-start RCE as konrad's framework allows:

  * Lapse rate   : MoistLapseRate (a true moist adiabat integrated by konrad's
                   ODE solver), NOT a prescribed 6.5 K/km Manabe slope.  This is
                   the physically correct convective profile and matches
                   ExoColumn's moist/zm convective-adjustment target in the
                   saturated limit.  (The numpy/scipy incompatibility that
                   forced FixedLapseRate in earlier konrad envs is absent in the
                   pinned env below.)
  * Composition  : CH4 = 1.8 ppm; N2O, CO, CFCs, CCl4 zeroed so konrad's set of
                   radiatively-active gases matches ExoColumn's
                   {H2O, CO2, CH4, O3, O2}.  O3 keeps konrad's default
                   mid-latitude climatology (comparable to ExoColumn's
                   o3_profile='earth'); O2 keeps its default ~0.21 vmr.
  * Surface      : FixedTemperature(Ts_K), matching ExoColumn's prognostic Ts
                   at equilibrium (~287.96 K for the canonical run).
  * Humidity     : FixedRH, vertically uniform 0.75.  konrad's standard
                   validation closure (Kluft et al. 2019).  Note ExoColumn's
                   prognostic moisture is NOT uniform — it decreases aloft and
                   dries out near the cold point — so some residual tropospheric
                   temperature difference vs ExoColumn is expected and physical,
                   not a configuration artifact.

Working environment (Python 3.11):
    konrad 1.0.2, numpy 1.26.4, scipy 1.13.1
    (numpy 2.x breaks HardAdjustment.stabilize; scipy >=1.14 pulls numpy >=2.1.)
"""

import sys
import warnings
import numpy as np
import konrad

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
Ts_K    = float(sys.argv[1]) if len(sys.argv) > 1 else 288.0
co2_ppm = float(sys.argv[2]) if len(sys.argv) > 2 else 400.0
outpath = sys.argv[3] if len(sys.argv) > 3 else "iofiles/konrad_rce_ref.npz"

ch4_ppm = 1.8   # match ExoColumn cold-start default (1.8e-6 vmr)

print(f"konrad version : {konrad.__version__}")
print(f"Surface Ts     : {Ts_K} K")
print(f"CO2            : {co2_ppm} ppm")
print(f"CH4            : {ch4_ppm} ppm")
print(f"Output         : {outpath}")

# ---------------------------------------------------------------------------
# Pressure grid: 128 half-levels from 1000 hPa -> 1 hPa  (quadratic spacing
# packs more resolution in the upper troposphere / lower stratosphere)
# ---------------------------------------------------------------------------
_, phlev = konrad.utils.get_pressure_grids(1000e2, 1e2, 128)

# ---------------------------------------------------------------------------
# Atmosphere: initialise with standard profiles, then set composition to match
# ExoColumn's radiatively-active gas set {H2O, CO2, CH4, O3, O2}.
# O3 keeps konrad's default mid-latitude climatology; O2 keeps its default.
# ---------------------------------------------------------------------------
atmosphere = konrad.atmosphere.Atmosphere(phlev)
atmosphere["CO2"][:] = co2_ppm * 1e-6   # volume mixing ratio
atmosphere["CH4"][:] = ch4_ppm * 1e-6   # volume mixing ratio

# Zero the gases ExoColumn does not model, so the greenhouse forcing comes from
# the same set of absorbers in both models.
for _gas in ("N2O", "CO", "CFC11", "CFC12", "CFC22", "CCl4"):
    if _gas in atmosphere.data_vars:
        atmosphere[_gas][:] = 0.0

# ---------------------------------------------------------------------------
# RCE: fixed surface temperature, hard convective adjustment to a true MOIST
#       ADIABAT (MoistLapseRate), fixed relative humidity.
# ---------------------------------------------------------------------------
rce = konrad.RCE(
    atmosphere,
    surface=konrad.surface.FixedTemperature(temperature=Ts_K),
    convection=konrad.convection.HardAdjustment(),
    lapserate=konrad.lapserate.MoistLapseRate(),
    humidity=konrad.humidity.FixedRH(
        rh_func=konrad.humidity.VerticallyUniform(0.75)
    ),
    timestep="12h",
    max_duration="500d",   # moist-adiabat RCE; ample margin for convergence
)

print("Running konrad RCE loop (MoistLapseRate) ...")
rce.run()
print("Converged (or hit max_duration).")

# ---------------------------------------------------------------------------
# Extract final state
# ---------------------------------------------------------------------------
T_eq     = np.asarray(rce.atmosphere["T"][-1, :])      # K, shape (N,)
plev_hpa = np.asarray(rce.atmosphere["plev"][:]) / 100.0   # Pa -> hPa
z_km     = np.asarray(rce.atmosphere["z"][-1, :]) / 1000.0  # m -> km
# Relative humidity diagnosed from the equilibrium H2O vmr (robust to konrad
# version differences in how the humidity model exposes its state).
try:
    h2o_vmr = np.asarray(rce.atmosphere["H2O"][-1, :])
    es = konrad.physics.e_eq_water_mk(T_eq)        # saturation vapour pressure [Pa]
    rh = h2o_vmr * np.asarray(rce.atmosphere["plev"][:]) / es
except Exception as _e:
    print(f"(RH diagnostic unavailable: {_e!r})")
    rh = np.full_like(T_eq, np.nan)

ts_final = float(rce.surface["temperature"][-1])
icp      = int(np.argmin(T_eq))
print(f"Ts final          : {ts_final:.2f} K")
print(f"T @ surface layer : {T_eq[0]:.2f} K   (p = {plev_hpa[0]:.1f} hPa)")
print(f"T @ TOA           : {T_eq[-1]:.2f} K   (p = {plev_hpa[-1]:.2f} hPa)")
print(f"Cold point        : {T_eq[icp]:.2f} K @ {plev_hpa[icp]:.1f} hPa ({z_km[icp]:.1f} km)")

np.savez(
    outpath,
    plev_hpa=plev_hpa, T_K=T_eq, z_km=z_km, rh=rh,
    Ts_K=Ts_K, co2_ppm=co2_ppm, ch4_ppm=ch4_ppm,
)
print(f"Saved -> {outpath}")
