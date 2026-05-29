#!/usr/bin/env python3
"""
generate_konrad_rad.py — konrad RCE reference WITH radiation diagnostics.

Same configuration as generate_konrad_ref.py (Ts=288, CO2=400 ppm, CH4=1.8 ppm,
MoistLapseRate, HardAdjustment, FixedRH 0.75, n68-equivalent absorber set) but
saves the full converged profile AND the RRTMG longwave heating rates / fluxes
and the gas volume mixing ratios, on konrad's native grid.  Used for the offline
ExoRT-vs-RRTMG radiation comparison (feed this exact profile through ExoRT and
compare LW cooling at matched humidity).

Output: iofiles/konrad_rad_ref.npz
"""
import sys, warnings
import numpy as np
import konrad
warnings.filterwarnings("ignore")

Ts_K    = float(sys.argv[1]) if len(sys.argv) > 1 else 288.0
co2_ppm = 400.0
ch4_ppm = 1.8
outpath = "iofiles/konrad_rad_ref.npz"

_, phlev = konrad.utils.get_pressure_grids(1000e2, 1e2, 128)
atm = konrad.atmosphere.Atmosphere(phlev)
atm["CO2"][:] = co2_ppm * 1e-6
atm["CH4"][:] = ch4_ppm * 1e-6
for g in ("N2O", "CO", "CFC11", "CFC12", "CFC22", "CCl4"):
    if g in atm.data_vars:
        atm[g][:] = 0.0

rce = konrad.RCE(
    atm,
    surface=konrad.surface.FixedTemperature(temperature=Ts_K),
    convection=konrad.convection.HardAdjustment(),
    lapserate=konrad.lapserate.MoistLapseRate(),
    humidity=konrad.humidity.FixedRH(rh_func=konrad.humidity.VerticallyUniform(0.75)),
    timestep="12h", max_duration="500d",
)
print("Running konrad RCE (radiation reference) ...", flush=True)
rce.run()

a = rce.atmosphere
r = rce.radiation
g1 = lambda x: np.asarray(x[-1, :]) if np.asarray(x).ndim == 2 else np.asarray(x)

np.savez(
    outpath,
    plev_pa=np.asarray(a["plev"][:]),
    phlev_pa=np.asarray(a["phlev"][:]),
    T_K=g1(a["T"]),
    z_m=g1(a["z"]),
    h2o_vmr=g1(a["H2O"]),
    o3_vmr=g1(a["O3"]),
    o2_vmr=g1(a["O2"]),
    co2_vmr=g1(a["CO2"]),
    ch4_vmr=g1(a["CH4"]),
    lw_htngrt=g1(r["lw_htngrt"]),
    sw_htngrt=g1(r["sw_htngrt"]),
    lw_flxu=g1(r["lw_flxu"]),
    lw_flxd=g1(r["lw_flxd"]),
    Ts_K=float(rce.surface["temperature"][-1]),
)
icp = int(np.argmin(g1(a["T"])))
print(f"Ts={float(rce.surface['temperature'][-1]):.2f}  "
      f"cold point {g1(a['T'])[icp]:.2f}K @ {np.asarray(a['plev'][:])[icp]/100:.1f}hPa", flush=True)
print(f"Saved -> {outpath}", flush=True)
