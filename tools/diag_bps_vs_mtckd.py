#!/usr/bin/env python3
"""A/B comparison: BPS vs MT_CKD H2O continuum on the inner-HZ sweep.

Reuses the moist-/runaway-greenhouse inner-HZ configuration of
reference/moist_runaway/hz_inner.py (flux_only, non-ideal H2O EOS, albedo 0.32,
6-point hemispheric solar-zenith quadrature, variable_ps, p_top=0.002 Pa) but
sweeps a COARSE surface-temperature grid for BOTH continuum models and reports
OLR, absorbed SW (ASR), planetary albedo, and Seff = OLR/ASR side by side.

The motivating question (see memory project_hz_albedo_offset / project_hz_roadmap):
our IHZ albedo sits ~0.02 above Kopparapu et al. (2013) and the runaway Seff is
~1.11 vs his 1.06, the residual attributed largely to near-IR H2O shortwave
absorption.  Kopparapu's CLIMA used the BPS continuum; ExoColumn defaults to
MT_CKD 3.3.  This script measures whether switching to BPS moves albedo/Seff
toward Kopparapu.

Requires an ExoColumn binary built at PVER>=200 (same as hz_inner.py).

Usage:
    python3 tools/diag_bps_vs_mtckd.py
Env overrides:
    AB_TS_MIN, AB_TS_MAX, AB_TS_STEP   coarse Ts grid (default 280 440 20)
"""

import os
import sys
import subprocess
import numpy as np
import netCDF4 as nc

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
EXE = os.path.join(ROOT, 'run', 'exocol.exe')
NML_PATH = os.path.join(ROOT, 'exocol_config.nml')
OUT_NC = os.path.join(ROOT, 'iofiles', 'exocol_out.nc')
KOPP_DAT = os.path.join(ROOT, 'reference', 'moist_runaway', 'waterloss_IHZ_present.dat')
OUT_PNG = os.path.join(HERE, 'diag_bps_vs_mtckd.png')

MW_H2O = 18.015
ALBEDO = 0.32          # Kopparapu surface albedo (matches hz_inner.py)
T_STRATO = 200.0

TS_MIN = float(os.environ.get('AB_TS_MIN', 280.0))
TS_MAX = float(os.environ.get('AB_TS_MAX', 440.0))
TS_STEP = float(os.environ.get('AB_TS_STEP', 20.0))

# Identical to hz_inner.py's NML_TEMPLATE, plus the h2o_continuum knob.
NML_TEMPLATE = """\
&exocol_nml
  flux_only      = .true.
  variable_ps    = .true.
  ihz_profile    = .true.
  o3_profile     = 'none'
  msdist         = 1.0
  h2o_eos        = 'nonideal'
  h2o_continuum  = '{continuum}'
  sw_zenith_quad = .true.
  sw_nquad       = 6
/
&exocol_init
  input_file = ''
  ts         = {ts:.2f}
  t_strato   = {t_strato:.1f}
  p_top      = 0.002
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


def run_one(ts, continuum):
    """Run ExoColumn flux_only at (ts, continuum); return diagnostics dict or None."""
    nml = NML_TEMPLATE.format(ts=ts, t_strato=T_STRATO, albedo=ALBEDO, continuum=continuum)
    with open(NML_PATH, 'w') as f:
        f.write(nml)
    result = subprocess.run([EXE], cwd=ROOT, capture_output=True, text=True, timeout=180)
    if result.returncode != 0:
        print(f"  FAIL Ts={ts:.0f} {continuum}  rc={result.returncode}")
        print((result.stderr or result.stdout)[-300:])
        return None
    with nc.Dataset(OUT_NC) as ds:
        lwup = ds['LWUP'][:]
        swup = ds['SWUP'][:]
        swdn = ds['SWDN'][:]
    olr = float(lwup[0])
    asr = float(swdn[0] - swup[0])
    alpha = float(swup[0] / swdn[0]) if swdn[0] > 0 else np.nan
    seff = olr / asr if asr > 1e-6 else np.nan
    return dict(olr=olr, asr=asr, alpha=alpha, seff=seff)


def load_kopparapu():
    """Return (Ts, Seff, albedo, OLR, ASR) arrays from Kopparapu's .dat, or None."""
    if not os.path.exists(KOPP_DAT):
        return None
    rows = []
    with open(KOPP_DAT) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            try:
                rows.append([float(x) for x in parts[:6]])
            except ValueError:
                continue
    if not rows:
        return None
    a = np.array(rows)
    # columns: TGO SEFF PALB FH2O FTIR[OLR] FTSO[absorbed SW]
    return dict(ts=a[:, 0], seff=a[:, 1], alpha=a[:, 2], olr=a[:, 4], asr=a[:, 5])


def main():
    if not os.path.exists(EXE):
        print(f"ERROR: {EXE} not found — build first (make PVER=200).")
        return 1
    ts_grid = np.arange(TS_MIN, TS_MAX + 0.5 * TS_STEP, TS_STEP)
    orig = open(NML_PATH).read() if os.path.exists(NML_PATH) else None
    results = {'mtckd': {}, 'bps': {}}
    try:
        for ts in ts_grid:
            for cont in ('mtckd', 'bps'):
                r = run_one(ts, cont)
                if r:
                    results[cont][ts] = r
                    print(f"  ran Ts={ts:6.1f}  {cont:5s}  OLR={r['olr']:7.2f}  "
                          f"alpha={r['alpha']:.4f}  Seff={r['seff']:.4f}")
    finally:
        if orig is not None:
            with open(NML_PATH, 'w') as f:
                f.write(orig)
        elif os.path.exists(NML_PATH):
            os.remove(NML_PATH)

    kopp = load_kopparapu()

    def kval(field, ts):
        if kopp is None:
            return np.nan
        return float(np.interp(ts, kopp['ts'], kopp[field]))

    print("\n" + "=" * 92)
    print(" BPS vs MT_CKD H2O continuum  —  inner-HZ flux_only sweep "
          "(nonideal EOS, albedo 0.32, 6-pt zenith)")
    print("=" * 92)
    print(f"{'Ts':>6} | {'OLR_mt':>7} {'OLR_bps':>7} {'OLR_K':>6} | "
          f"{'alb_mt':>6} {'alb_bps':>7} {'alb_K':>6} | "
          f"{'Sef_mt':>6} {'Sef_bps':>7} {'Sef_K':>6} | {'dAlb':>7} {'dSeff':>7}")
    print("-" * 92)
    for ts in ts_grid:
        m = results['mtckd'].get(ts)
        b = results['bps'].get(ts)
        if not m or not b:
            continue
        dalb = b['alpha'] - m['alpha']
        dseff = b['seff'] - m['seff']
        print(f"{ts:6.0f} | {m['olr']:7.2f} {b['olr']:7.2f} {kval('olr', ts):6.1f} | "
              f"{m['alpha']:6.4f} {b['alpha']:7.4f} {kval('alpha', ts):6.4f} | "
              f"{m['seff']:6.4f} {b['seff']:7.4f} {kval('seff', ts):6.4f} | "
              f"{dalb:+7.4f} {dseff:+7.4f}")
    print("=" * 92)
    print("dAlb/dSeff = bps - mtckd.  Negative dAlb means BPS lowers albedo (toward Kopparapu).")

    # ---- figure ----
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        tsm = np.array([t for t in ts_grid if t in results['mtckd'] and t in results['bps']])
        olr_m = np.array([results['mtckd'][t]['olr'] for t in tsm])
        olr_b = np.array([results['bps'][t]['olr'] for t in tsm])
        alb_m = np.array([results['mtckd'][t]['alpha'] for t in tsm])
        alb_b = np.array([results['bps'][t]['alpha'] for t in tsm])
        sef_m = np.array([results['mtckd'][t]['seff'] for t in tsm])
        sef_b = np.array([results['bps'][t]['seff'] for t in tsm])
        fig, ax = plt.subplots(1, 3, figsize=(11.0, 3.6), dpi=300)
        fig.patch.set_facecolor('white')
        for a in ax:
            a.set_facecolor('white')
            a.set_xlabel('surface temperature  $T_s$  [K]')
        ax[0].plot(tsm, olr_m, '-o', color='C3', ms=3, label='MT_CKD (ExoColumn)')
        ax[0].plot(tsm, olr_b, '-s', color='C1', ms=3, label='BPS')
        ax[1].plot(tsm, alb_m, '-o', color='C3', ms=3, label='MT_CKD')
        ax[1].plot(tsm, alb_b, '-s', color='C1', ms=3, label='BPS')
        ax[2].plot(tsm, sef_m, '-o', color='C3', ms=3, label='MT_CKD')
        ax[2].plot(tsm, sef_b, '-s', color='C1', ms=3, label='BPS')
        if kopp is not None:
            sel = (kopp['ts'] >= tsm.min() - 5) & (kopp['ts'] <= tsm.max() + 5)
            ax[0].plot(kopp['ts'][sel], kopp['olr'][sel], '--', color='k', lw=1, label='Kopparapu 2013')
            ax[1].plot(kopp['ts'][sel], kopp['alpha'][sel], '--', color='k', lw=1, label='Kopparapu')
            ax[2].plot(kopp['ts'][sel], kopp['seff'][sel], '--', color='k', lw=1, label='Kopparapu')
        ax[0].set_ylabel('OLR  [W m$^{-2}$]')
        ax[1].set_ylabel('planetary albedo')
        ax[2].set_ylabel('$S_{\\rm eff}=$ OLR/ASR')
        ax[0].legend(fontsize=7, frameon=False)
        fig.tight_layout()
        fig.savefig(OUT_PNG, facecolor='white')
        print(f"\nwrote {OUT_PNG}")
    except Exception as exc:  # plotting is best-effort
        print(f"(plot skipped: {exc})")
    return 0


if __name__ == '__main__':
    sys.exit(main())
