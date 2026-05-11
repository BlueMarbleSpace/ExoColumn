#!/usr/bin/env python3
"""Generate an Earth-like RTprofile_in.nc-format input for ExoColumn.

Replaces the IDL makeColumn.pro for Earth configurations.  Two presets
('preindustrial', 'present_day') set CO2 / CH4; everything else is held
to a clean Earth atmosphere with realistic O2, mwdry, cp_dry, surface
albedo, and a Manabe-Wetherald RH(p) initialisation.  Ozone is left at
zero (per the current development plan); add it later via --o3 once a
profile is selected.

Usage:
    python tools/make_earth_column.py --preset preindustrial \
        --out iofiles/exocol_in_preindustrial.nc
    python tools/make_earth_column.py --preset present_day \
        --out iofiles/exocol_in_present_day.nc
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

# ---------------------------------------------------------------------------
# Physical constants (SI)
# ---------------------------------------------------------------------------
GRAV = 9.80616          # m s-2 (Earth, matches ExoRT)
R_UNIV = 8.31446        # J mol-1 K-1
EPS_WV = 18.016 / 28.97 # mw_h2o / mw_air, used only for qsat-from-esat

# Per-species molar masses (g/mol) and isobaric specific heats (J/kg/K)
MW = dict(N2=28.013, O2=31.999, Ar=39.948,
          CO2=44.010, CH4=16.043, O3=47.998, H2=2.016, H2O=18.016)
CP = dict(N2=1039.0, O2=918.0, Ar=520.3,
          CO2=846.0, CH4=2226.0, H2=14320.0)

# ---------------------------------------------------------------------------
# Presets
# ---------------------------------------------------------------------------
PRESETS = {
    # 'control_n2_co2' reproduces the previously-validated baseline that
    # converged at Ts=286.48 K, TOA=0.003 W/m^2 (gap-analysis memo).  Use it
    # as the anchor when walking forward to add gases / change albedo: if the
    # control no longer converges, the regression is in our setup pipeline,
    # not the science.
    'control_n2_co2': dict(
        co2_vmr=376.0e-6,
        ch4_vmr=0.0,
        o2_vmr=0.0,
        ar_vmr=0.0,
        albedo=0.25,
        rh_surf=0.97,             # TS273K profile is near-saturated near surface
        ts=273.0,
    ),
    'preindustrial': dict(
        co2_vmr=280.0e-6,
        ch4_vmr=722.0e-9,
        o2_vmr=0.20946,
        ar_vmr=0.00934,
        albedo=0.30,
        rh_surf=0.77,
        ts=250.0,                 # cold-start; equilibrium is approached from below
    ),
    'present_day': dict(
        co2_vmr=420.0e-6,
        ch4_vmr=1900.0e-9,
        o2_vmr=0.20946,
        ar_vmr=0.00934,
        albedo=0.30,
        rh_surf=0.77,
        ts=250.0,                 # cold-start; equilibrium is approached from below
    ),
}
# NOTE: Ts here is the INITIAL surface temperature only — the model's slab
# ocean evolves it.  Starting hot (e.g. 287 K) with full Manabe-Wetherald
# RH seeds enough TCWV to tip the column into the runaway-greenhouse branch
# of the bistable RCE solution; starting cold (~250 K) approaches the
# climatological equilibrium from below, which is the previously-validated
# behaviour.


def esat_cc(T):
    """Saturation vapour pressure (Pa) over liquid water, Clausius-Clapeyron.

    Same formula as exocol_rce_loop.F90 so initialisation matches the
    relaxation target inside the model."""
    Lv = 2.5e6
    Rv = 461.5
    es0 = 611.2
    T0 = 273.15
    return es0 * np.exp(Lv / Rv * (1.0 / T0 - 1.0 / T))


def manabe_wetherald_rh(p, ps, rh_surf=0.77, p_min_frac=0.02):
    """RH profile from Manabe & Wetherald 1967.

    RH(p) = rh_surf * (p/ps - p_min_frac) / (1 - p_min_frac), clamped to >= 0.
    Yields rh_surf at the surface, zero above ~p_min_frac*ps.
    """
    h = (p / ps - p_min_frac) / (1.0 - p_min_frac)
    return np.clip(rh_surf * h, 0.0, rh_surf)


def build_grid(pver=300, p_top=10.0, ps=101325.0):
    """Log-spaced pressure grid; pmid = geometric mean of bracketing pints."""
    pverp = pver + 1
    pint = np.logspace(np.log10(p_top), np.log10(ps), pverp)
    pmid = np.sqrt(pint[:-1] * pint[1:])
    pdel = pint[1:] - pint[:-1]
    return pint, pmid, pdel


def initial_T(pmid, pint, ps, ts, R_dry, t_strat=200.0, gamma=6.5e-3):
    """6.5 K/km tropospheric lapse rate joined to an isothermal stratosphere.

    Matches the TS273K reference profile used by the previously-validated
    moist+CC run, and matches the Manabe-Wetherald 1967 convective scheme.
    Using a dry adiabat (kappa=R/cp) for init was ~10 K too cold in the
    upper troposphere, which the 5-day virtual timestep could not damp out
    before the column tipped past Komabayashi-Ingersoll.

    z(p) is taken as a constant-scale-height approximation -H ln(p/ps) with
    H = R_dry * 250 K / g (good enough for an init profile that the model
    will refine in step 1).
    """
    H = R_dry * 250.0 / GRAV
    z_mid = -H * np.log(pmid / ps)
    z_int = -H * np.log(pint / ps)
    tmid = np.maximum(t_strat, ts - gamma * z_mid)
    tint = np.maximum(t_strat, ts - gamma * z_int)
    tint[-1] = ts
    return tmid, tint


def hydrostatic_zint(pint, tmid, R_dry):
    """Geometric height at interfaces from the hypsometric equation.

    Ground (k=pverp-1) is z=0.  Builds upward.
    """
    pverp = pint.size
    zint = np.zeros(pverp)
    for k in range(pverp - 2, -1, -1):
        # layer k spans pint[k] (top) to pint[k+1] (bottom); use tmid[k]
        dz = R_dry * tmid[k] / GRAV * np.log(pint[k + 1] / pint[k])
        zint[k] = zint[k + 1] + dz
    return zint


def build_column(preset_name, pver=300, p_top=10.0, ps=101325.0,
                 albedo=None, coszrs=0.5, rh_surf=None, ts_init=None,
                 co2_ppmv=None, ch4_ppbv=None, no_o2=False, no_ch4=False):
    if preset_name not in PRESETS:
        raise ValueError(f'unknown preset {preset_name!r}; '
                         f'choose from {list(PRESETS)}')
    cfg = PRESETS[preset_name]
    co2_vmr = cfg['co2_vmr'] if co2_ppmv is None else co2_ppmv * 1e-6
    ch4_vmr = cfg['ch4_vmr'] if ch4_ppbv is None else ch4_ppbv * 1e-9
    if no_ch4:
        ch4_vmr = 0.0
    ts = cfg['ts'] if ts_init is None else ts_init
    if albedo is None:
        albedo = cfg['albedo']
    if rh_surf is None:
        rh_surf = cfg['rh_surf']

    # Earth dry composition (mole fractions before water).  Argon does not
    # appear as an absorber in the n68equiv tables, so we lump its mass into
    # n2mmr but keep the correct mwdry.
    o2_vmr = 0.0 if no_o2 else cfg['o2_vmr']
    ar_vmr = cfg['ar_vmr']
    o3_vmr = 0.0
    h2_vmr = 0.0
    n2_vmr = 1.0 - o2_vmr - ar_vmr - co2_vmr - ch4_vmr - o3_vmr - h2_vmr

    # Dry molecular weight and dry specific heat (mass-weighted)
    vmr = dict(N2=n2_vmr, O2=o2_vmr, Ar=ar_vmr,
               CO2=co2_vmr, CH4=ch4_vmr, H2=h2_vmr)
    mwdry = sum(vmr[s] * MW[s] for s in vmr)
    cp_dry = sum(vmr[s] * MW[s] * CP[s] for s in vmr) / mwdry
    R_dry = R_UNIV / (mwdry * 1.0e-3)              # J/kg/K
    kappa = R_dry / cp_dry

    # Mass mixing ratios (referenced to dry mass).  Argon gets folded into
    # n2mmr because ExoRT only carries N2 as the radiatively-inert background.
    def to_mmr(species):
        return vmr[species] * MW[species] / mwdry

    n2mmr_val = (vmr['N2'] * MW['N2'] + vmr['Ar'] * MW['Ar']) / mwdry
    o2mmr_val = to_mmr('O2')
    co2mmr_val = to_mmr('CO2')
    ch4mmr_val = to_mmr('CH4')
    h2mmr_val = to_mmr('H2')

    # Vertical grid + initial T
    pint, pmid, pdel = build_grid(pver=pver, p_top=p_top, ps=ps)
    tmid, tint = initial_T(pmid, pint, ps, ts, R_dry)

    # Initial water vapour: Manabe-Wetherald RH * qsat(T,p), with esat
    # capped at 0.99*p (matching the bug fix in exocol_rce_loop.F90).
    rh = manabe_wetherald_rh(pmid, ps, rh_surf=rh_surf)
    es = np.minimum(esat_cc(tmid), 0.99 * pmid)
    qsat = EPS_WV * es / (pmid - es)
    h2ommr = rh * qsat
    h2ommr = np.clip(h2ommr, 0.0, None)

    zint = hydrostatic_zint(pint, tmid, R_dry)

    # Tile mass mixing ratios constant in p
    def constant(val):
        return np.full(pver, val, dtype=np.float64)

    fields = dict(
        ts=ts, ps=ps, coszrs=coszrs,
        mw=mwdry, cp=cp_dry,
        asdir=albedo, asdif=albedo, aldir=albedo, aldif=albedo,
        tmid=tmid, tint=tint,
        pmid=pmid, pint=pint, pdel=pdel, zint=zint,
        h2ommr=h2ommr,
        co2mmr=constant(co2mmr_val),
        ch4mmr=constant(ch4mmr_val),
        o2mmr=constant(o2mmr_val),
        o3mmr=constant(0.0),
        n2mmr=constant(n2mmr_val),
        h2mmr=constant(h2mmr_val),
    )
    meta = dict(preset=preset_name, mwdry=mwdry, cp_dry=cp_dry, R_dry=R_dry,
                kappa=kappa, ts=ts, ps=ps, albedo=albedo,
                co2_ppmv=co2_vmr * 1e6, ch4_ppbv=ch4_vmr * 1e9,
                o2_pct=o2_vmr * 100,
                tcwv=float(np.sum(h2ommr * pdel) / GRAV))   # kg/m2
    return fields, meta


def write_nc(path, fields):
    pver = fields['pmid'].size
    pverp = fields['pint'].size
    assert pverp == pver + 1, 'grid mismatch'

    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Dataset(path, 'w', clobber=True, format='NETCDF3_CLASSIC') as ds:
        ds.createDimension('pver', pver)
        ds.createDimension('pverp', pverp)
        ds.createDimension('one', 1)

        scalars = dict(
            ts=('Surface temperature', 'K'),
            ps=('Surface pressure', 'Pa'),
            coszrs=('cosine of zenith angle', 'none'),
            mw=('molecular weight of dry air', 'AMU'),
            cp=('specific heat of dry air', 'J/kg K'),
            asdir=('albedo shortwave direct', 'albedo'),
            asdif=('albedo shortwave diffuse', 'albedo'),
            aldir=('albedo longwave direct', 'albedo'),
            aldif=('albedo longwave diffuse', 'albedo'),
        )
        for name, (title, units) in scalars.items():
            v = ds.createVariable(name, 'f4', ('one',))
            v.title = title
            v.units = units
            v[:] = np.float32(fields[name])

        midvars = dict(
            tmid=('Mid-layer temperatures', 'K'),
            pmid=('Mid-layer pressures', 'Pa'),
            pdel=('Pressures difference', 'Pa'),
            h2ommr=('H2O MMR', 'mmr'),
            co2mmr=('CO2 MMR', 'mmr'),
            ch4mmr=('CH4 MMR', 'mmr'),
            o2mmr=('O2 MMR', 'mmr'),
            o3mmr=('O3 MMR', 'mmr'),
            n2mmr=('N2 MMR', 'mmr'),
            h2mmr=('H2 MMR', 'mmr'),
        )
        for name, (title, units) in midvars.items():
            v = ds.createVariable(name, 'f4', ('pver',))
            v.title = title
            v.units = units
            v[:] = fields[name].astype(np.float32)

        intvars = dict(
            tint=('Interface temperatures', 'K'),
            pint=('Interface pressures', 'Pa'),
            zint=('Interface heights', 'm'),
        )
        for name, (title, units) in intvars.items():
            v = ds.createVariable(name, 'f4', ('pverp',))
            v.title = title
            v.units = units
            v[:] = fields[name].astype(np.float32)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--preset', choices=list(PRESETS), required=True)
    p.add_argument('--out', default='iofiles/exocol_in.nc',
                   help='output NetCDF path (default: iofiles/exocol_in.nc)')
    p.add_argument('--pver', type=int, default=300,
                   help='number of mid layers (must match ExoRT exo_pver)')
    p.add_argument('--ps', type=float, default=101325.0,
                   help='surface pressure [Pa] (default 101325)')
    p.add_argument('--p-top', type=float, default=10.0,
                   help='top interface pressure [Pa] (default 10)')
    p.add_argument('--albedo', type=float, default=None,
                   help='surface albedo override (default: preset value)')
    p.add_argument('--coszrs', type=float, default=0.5,
                   help='cosine of zenith angle (default 0.5)')
    p.add_argument('--rh-surf', type=float, default=None,
                   help='surface RH override for Manabe-Wetherald init (default: preset value)')
    p.add_argument('--ts-init', type=float, default=None,
                   help='override initial surface T (default: preset value). '
                        'Approach equilibrium from below to avoid the runaway-greenhouse branch.')
    p.add_argument('--co2-ppmv', type=float, default=None,
                   help='override CO2 in ppmv (default: preset value)')
    p.add_argument('--ch4-ppbv', type=float, default=None,
                   help='override CH4 in ppbv (default: preset value)')
    p.add_argument('--no-o2', action='store_true',
                   help='zero out O2 (ablation diagnostic)')
    p.add_argument('--no-ch4', action='store_true',
                   help='zero out CH4 (ablation diagnostic)')
    args = p.parse_args(argv)

    fields, meta = build_column(
        preset_name=args.preset, pver=args.pver, p_top=args.p_top,
        ps=args.ps, albedo=args.albedo, coszrs=args.coszrs,
        rh_surf=args.rh_surf, ts_init=args.ts_init,
        co2_ppmv=args.co2_ppmv, ch4_ppbv=args.ch4_ppbv,
        no_o2=args.no_o2, no_ch4=args.no_ch4,
    )
    write_nc(args.out, fields)

    print(f"wrote {args.out}")
    print(f"  preset      : {meta['preset']}")
    print(f"  Ts (init)   : {meta['ts']:.2f} K")
    print(f"  Ps          : {meta['ps']/100:.2f} hPa")
    print(f"  CO2         : {meta['co2_ppmv']:.1f} ppmv")
    print(f"  CH4         : {meta['ch4_ppbv']:.1f} ppbv")
    print(f"  O2          : {meta['o2_pct']:.2f} %")
    print(f"  mwdry       : {meta['mwdry']:.4f} g/mol")
    print(f"  cp_dry      : {meta['cp_dry']:.2f} J/kg/K")
    print(f"  R_dry       : {meta['R_dry']:.2f} J/kg/K")
    print(f"  kappa (R/cp): {meta['kappa']:.4f}")
    print(f"  albedo      : {meta['albedo']:.2f}")
    print(f"  TCWV (init) : {meta['tcwv']:.2f} kg/m^2")


if __name__ == '__main__':
    sys.exit(main())
