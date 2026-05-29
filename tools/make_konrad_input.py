#!/usr/bin/env python3
"""
make_konrad_input.py — build an ExoColumn input (RTprofile_in.nc) from konrad's
converged profile, on konrad's NATIVE grid, for the ExoRT-vs-RRTMG radiation test.

Reads iofiles/konrad_rad_ref.npz (from generate_konrad_rad.py) and writes an
ExoColumn-format input NetCDF with konrad's T, H2O, O3, O2, CO2, CH4 and N2
(as the remainder).  ExoColumn must be built with PVER = (number of konrad full
levels) so the dimension check passes; the script prints that number.
"""
import numpy as np
from netCDF4 import Dataset

MW = dict(N2=28.013, O2=31.999, CO2=44.010, CH4=16.043, O3=47.998, H2O=18.016)
CP = dict(N2=1039.0, O2=918.0, CO2=846.0, CH4=2226.0, O3=820.0)
GRAV = 9.80616

k = np.load("iofiles/konrad_rad_ref.npz")
pmid = k["plev_pa"].astype(float)          # (N,)  full levels
pint = k["phlev_pa"].astype(float)         # (N+1,) half levels
T    = k["T_K"].astype(float)
z    = k["z_m"].astype(float)
h2o  = k["h2o_vmr"].astype(float)
o3   = k["o3_vmr"].astype(float)
o2   = k["o2_vmr"].astype(float)
co2  = k["co2_vmr"].astype(float)
ch4  = k["ch4_vmr"].astype(float)
ts   = float(k["Ts_K"])
N    = pmid.size
print(f"konrad grid: {N} full levels, {pint.size} interfaces; build ExoColumn with PVER={N}")

# konrad indexing is surface->TOA (plev[0] = 1000 hPa); ExoColumn is TOA->surface.
if pmid[0] > pmid[-1]:
    pmid, T, z, h2o, o3, o2, co2, ch4 = [a[::-1] for a in (pmid,T,z,h2o,o3,o2,co2,ch4)]
    pint = pint[::-1]

pdel = np.diff(pint)                        # (N,)  positive, TOA->srf
ps   = float(pint[-1])

# Dry composition: N2 = remainder (per level).  mwdry from a representative dry
# mix (use column mean of the dry vmrs; O3 varies but is tiny).
n2 = 1.0 - (o2 + co2 + ch4 + o3)
def col_mean(x): return float(np.mean(x))
vmr_dry = dict(N2=col_mean(n2), O2=col_mean(o2), CO2=col_mean(co2),
               CH4=col_mean(ch4), O3=col_mean(o3))
mwdry = sum(vmr_dry[s]*MW[s] for s in vmr_dry)
cpdry = sum(vmr_dry[s]*MW[s]*CP[s] for s in vmr_dry) / mwdry
eps   = MW["H2O"] / mwdry
print(f"mwdry={mwdry:.3f}  cpdry={cpdry:.1f}  eps={eps:.4f}")

# Mass mixing ratios referenced to dry air (ExoColumn convention: w = (vmr*Mi)/mwdry).
def mmr(vmr_arr, sp): return vmr_arr * MW[sp] / mwdry
co2mmr = mmr(co2, "CO2"); ch4mmr = mmr(ch4, "CH4"); o2mmr = mmr(o2, "O2")
o3mmr  = mmr(o3, "O3");   n2mmr  = mmr(n2, "N2");   h2mmr = np.zeros(N)
c2h6mmr = np.zeros(N)
# Water: w = (vmr/(1-vmr))*eps  (mass water / mass dry air), matching qsat=eps*es/(p-es).
h2ommr = (h2o / (1.0 - h2o)) * eps

# Interface temperatures: log-p interpolation of T; pin surface to ts, extrapolate TOA.
lpm = np.log(pmid)
tint = np.empty(N+1)
for i in range(1, N):
    tint[i] = T[i-1] + (T[i]-T[i-1])*(np.log(pint[i])-lpm[i-1])/(lpm[i]-lpm[i-1])
tint[0]  = T[0] + (T[1]-T[0])*(np.log(pint[0])-lpm[0])/(lpm[1]-lpm[0])
tint[-1] = ts

# Interface heights: log-p interpolation of konrad z; surface = 0.
zint = np.interp(np.log(pint), lpm, z)   # lpm increasing TOA->srf
zint[-1] = 0.0

out = "iofiles/konrad_as_exocol_in.nc"
with Dataset(out, "w", clobber=True, format="NETCDF3_CLASSIC") as ds:
    ds.createDimension("pver", N)
    ds.createDimension("pverp", N+1)
    ds.createDimension("one", 1)
    def s(name, val):
        v = ds.createVariable(name, "f8", ("one",)); v[:] = val
    def a(name, arr, dim):
        v = ds.createVariable(name, "f8", (dim,)); v[:] = arr
    s("ts", ts); s("ps", ps); s("coszrs", 0.5); s("mw", mwdry); s("cp", cpdry)
    for nm in ("asdir","asdif","aldir","aldif"): s(nm, 0.2873)
    a("tmid", T, "pver"); a("pmid", pmid, "pver"); a("pdel", pdel, "pver")
    a("h2ommr", h2ommr, "pver")
    a("co2mmr", co2mmr, "pver"); a("ch4mmr", ch4mmr, "pver"); a("c2h6mmr", c2h6mmr, "pver")
    a("o2mmr", o2mmr, "pver"); a("o3mmr", o3mmr, "pver")
    a("n2mmr", n2mmr, "pver"); a("h2mmr", h2mmr, "pver")
    a("tint", tint, "pverp"); a("pint", pint, "pverp"); a("zint", zint, "pverp")
print(f"Wrote {out}")
