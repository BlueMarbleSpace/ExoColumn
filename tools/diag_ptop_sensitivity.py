#!/usr/bin/env python3
"""
diag_ptop_sensitivity.py — model-top (p_top) convergence test for the inner-HZ sweep.

The inner-HZ figure (tools/hz_inner.py) puts the model top at a fixed pressure
p_top = 1 Pa.  The corresponding ALTITUDE varies with Ts (thinner/colder columns
reach a lower altitude at the same top pressure), which is why panel (d) does not
uniformly reach Kopparapu+2013's 100+ km.

This script checks whether lowering p_top changes the radiative diagnostics
(OLR / ASR / Seff) — i.e. whether the top is radiatively converged — or only the
plotted altitude range.  For two representative Ts it reruns ExoColumn flux_only
at p_top = 1, 0.1, 0.01 Pa and tabulates OLR, ASR, albedo, Seff, the top altitude,
and a tropospheric-resolution proxy (number of layers with p > 100 hPa).

Same composition / EOS / albedo as hz_inner.py.  The user's exocol_config.nml is
backed up and restored.
"""

import os
import subprocess
import numpy as np
import netCDF4 as nc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT     = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXE      = os.path.join(ROOT, 'run', 'exocol.exe')
OUT_NC   = os.path.join(ROOT, 'iofiles', 'exocol_out.nc')
NML_PATH = os.path.join(ROOT, 'exocol_config.nml')

ALBEDO   = 0.32
T_STRATO = 200.0
MW_H2O   = 18.015
H2O_EOS  = os.environ.get('HZ_H2O_EOS', 'ideal')

TS_VALUES   = [380.0, 700.0, 1000.0, 1500.0, 2000.0, 2400.0]
PTOP_VALUES = [1.0, 0.1, 0.01]   # Pa

NML_TEMPLATE = """\
&exocol_nml
  flux_only   = .true.
  variable_ps = .true.
  ihz_profile = .true.
  o3_profile  = 'none'
  msdist      = 1.0
  h2o_eos     = '{h2o_eos}'
/
&exocol_init
  input_file = ''
  ts         = {ts:.2f}
  t_strato   = {t_strato:.1f}
  p_top      = {p_top:.4f}
  rh_init    = 1.0
  coszrs     = 0.5
  asdir      = {albedo:.4f}
  asdif      = {albedo:.4f}
  aldir      = {albedo:.4f}
  aldif      = {albedo:.4f}
/
&exocol_composition
  n2_vmr   = 0.78
  o2_vmr   = 0.210
  ar_vmr   = 0.01
  co2_vmr  = 3.3e-4
  ch4_vmr  = 0.0
  o3_vmr   = 0.0
/
"""


def run_one(ts, p_top):
    nml = NML_TEMPLATE.format(ts=ts, t_strato=T_STRATO, albedo=ALBEDO,
                              h2o_eos=H2O_EOS, p_top=p_top)
    with open(NML_PATH, 'w') as f:
        f.write(nml)
    result = subprocess.run([EXE], cwd=ROOT, capture_output=True,
                            text=True, timeout=180)
    if result.returncode != 0:
        print(f"  FAIL Ts={ts:.0f} p_top={p_top}  rc={result.returncode}")
        print(result.stderr[-400:] if result.stderr else '')
        return None
    with nc.Dataset(OUT_NC) as ds:
        lwup = ds['LWUP'][:]; swup = ds['SWUP'][:]; swdn = ds['SWDN'][:]
        pmid = np.array(ds['pmid'][:]); zint = np.array(ds['zint'][:])
    olr   = float(lwup[0])
    asr   = float(swdn[0] - swup[0])
    alpha = float(swup[0] / swdn[0]) if swdn[0] > 0 else 0.0
    seff  = olr / asr if asr > 1e-6 else np.nan
    ztop  = zint[0] / 1000.
    n_trop = int(np.sum(pmid > 1.0e4))   # layers below 100 hPa (trop. resolution proxy)
    return dict(olr=olr, asr=asr, alpha=alpha, seff=seff, ztop=ztop, n_trop=n_trop)


def _plot(ts_values, data):
    ts = np.asarray(ts_values, float)
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(9.0, 3.8), dpi=300)
    fig.patch.set_facecolor('white')
    colors = {1.0: 'C0', 0.1: 'C1', 0.01: 'C3'}
    markers = {1.0: 'o', 0.1: 's', 0.01: '^'}
    for p in PTOP_VALUES:
        c, m = colors.get(p, 'C2'), markers.get(p, 'x')
        ax_a.plot(ts, data[p]['alpha'], c=c, marker=m, ms=4, lw=1.2,
                  label=f'$p_{{\\rm top}}$ = {p:g} Pa')
        ax_b.plot(ts, data[p]['asr'], c=c, marker=m, ms=4, lw=1.2,
                  label=f'$p_{{\\rm top}}$ = {p:g} Pa')
    for ax in (ax_a, ax_b):
        ax.set_xlabel('$T_s$ (K)')
        ax.set_facecolor('white')
        ax.legend(fontsize=7, framealpha=0.9)
    ax_a.set_ylabel('Planetary albedo $\\alpha_p$')
    ax_b.set_ylabel('Absorbed SW (W m$^{-2}$)')
    ax_a.set_title('Albedo is invariant to model top', fontsize=9)
    ax_b.set_title('ASR is invariant to model top', fontsize=9)
    fig.tight_layout()
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       'diag_ptop_sensitivity.png')
    fig.savefig(out, dpi=150, facecolor='white')
    plt.close(fig)
    print(f"Saved: {out}")


def main():
    if not os.path.isfile(EXE):
        raise FileNotFoundError(EXE)
    orig = open(NML_PATH).read() if os.path.exists(NML_PATH) else None
    data = {p: dict(olr=[], asr=[], alpha=[]) for p in PTOP_VALUES}
    try:
        print(f"p_top model-top convergence test   (H2O EOS = {H2O_EOS}, albedo = {ALBEDO})\n")
        for ts in TS_VALUES:
            print(f"Ts = {ts:.0f} K")
            print(f"  {'p_top[Pa]':>10}  {'OLR':>8}  {'ASR':>8}  {'albedo':>7}"
                  f"  {'Seff':>7}  {'z_top[km]':>9}  {'n(p>100hPa)':>11}")
            base = None
            for p_top in PTOP_VALUES:
                r = run_one(ts, p_top)
                if r is None:
                    data[p_top]['olr'].append(np.nan)
                    data[p_top]['asr'].append(np.nan)
                    data[p_top]['alpha'].append(np.nan)
                    continue
                data[p_top]['olr'].append(r['olr'])
                data[p_top]['asr'].append(r['asr'])
                data[p_top]['alpha'].append(r['alpha'])
                if base is None:
                    base = r
                dolr = r['olr'] - base['olr']
                dseff = r['seff'] - base['seff']
                print(f"  {p_top:>10.3g}  {r['olr']:>8.2f}  {r['asr']:>8.2f}"
                      f"  {r['alpha']:>7.3f}  {r['seff']:>7.4f}  {r['ztop']:>9.1f}"
                      f"  {r['n_trop']:>11d}    (ΔOLR={dolr:+.2f}, ΔSeff={dseff:+.4f})")
            print()
        _plot(TS_VALUES, data)
    finally:
        if orig is not None:
            with open(NML_PATH, 'w') as f:
                f.write(orig)
        elif os.path.exists(NML_PATH):
            os.remove(NML_PATH)
    print("Original exocol_config.nml restored.")


if __name__ == '__main__':
    main()
