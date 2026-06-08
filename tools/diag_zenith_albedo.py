#!/usr/bin/env python3
"""
diag_zenith_albedo.py — does the single-zenith-angle approximation explain the
ExoColumn-vs-Kopparapu(2013) planetary-albedo offset (panel b of hz_inner)?

hz_inner.py runs the SW radiation at a SINGLE solar zenith angle coszrs = 0.5
(60 deg).  Kopparapu et al. (2013) instead average the solar flux over SIX
zenith angles (11, 25.3, 39.6, 54, 68.4, 82.8 deg) by Gaussian quadrature — i.e.
a proper hemispheric (Bond) average.  Because the plane albedo alpha(mu0) of a
Rayleigh-scattering atmosphere RISES toward grazing incidence (small mu0), a
single 60 deg point overweights grazing rays and biases the planetary albedo HIGH.

This script reconstructs alpha(mu0) by rerunning ExoColumn flux_only over a grid
of mu0 for several Ts, then forms the hemispheric Bond albedo

    A_Bond = 2 * integral_0^1 alpha(mu0) * mu0 dmu0

(the flux-weighted mean over the illuminated disk).  It prints, for each Ts:
  - alpha(0.5)      : ExoColumn's current single-60deg value (= hz_inner panel b)
  - A_Bond          : the hemispheric average (Kopparapu-equivalent)
  - Kopparapu Fig 3b reference value
so we can see how much of the offset the zenith treatment accounts for.

Same composition / EOS / surface albedo (0.32) as hz_inner.py.  The user's
exocol_config.nml is backed up and restored.
"""

import os
import subprocess
import numpy as np
import netCDF4 as nc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# numpy 2.x renamed trapz -> trapezoid; support both.
_trapz = getattr(np, 'trapezoid', getattr(np, 'trapz', None))

ROOT     = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXE      = os.path.join(ROOT, 'run', 'exocol.exe')
OUT_NC   = os.path.join(ROOT, 'iofiles', 'exocol_out.nc')
NML_PATH = os.path.join(ROOT, 'exocol_config.nml')

ALBEDO   = 0.32
T_STRATO = 200.0
H2O_EOS  = os.environ.get('HZ_H2O_EOS', 'ideal')

# Ts values: 400 K = Kopparapu albedo MINIMUM; >=1000 = plateau.
TS_VALUES = [300.0, 400.0, 500.0, 700.0, 1000.0, 1500.0, 2000.0]
# mu0 grid spanning [0,1] for the Bond-albedo integral (denser near the ends).
MU_VALUES = [0.02, 0.05, 0.125, 0.25, 0.368, 0.5, 0.588, 0.70, 0.771, 0.90, 0.99]
# Kopparapu Fig 3b reference (digitized): min ~0.16 @400K, plateau ~0.193.
KOPP_REF  = {300.0: 0.255, 400.0: 0.160, 500.0: 0.172, 700.0: 0.185,
             1000.0: 0.190, 1500.0: 0.193, 2000.0: 0.193}

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
  p_top      = 1.0
  rh_init    = 1.0
  coszrs     = {coszrs:.5f}
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


def alpha_at(ts, mu0):
    nml = NML_TEMPLATE.format(ts=ts, t_strato=T_STRATO, albedo=ALBEDO,
                              h2o_eos=H2O_EOS, coszrs=mu0)
    with open(NML_PATH, 'w') as f:
        f.write(nml)
    r = subprocess.run([EXE], cwd=ROOT, capture_output=True, text=True, timeout=180)
    if r.returncode != 0:
        print(f"  FAIL Ts={ts:.0f} mu0={mu0}  rc={r.returncode}")
        print((r.stderr or '')[-300:])
        return np.nan
    with nc.Dataset(OUT_NC) as ds:
        swup = float(ds['SWUP'][:][0]); swdn = float(ds['SWDN'][:][0])
    return swup / swdn if swdn > 0 else np.nan


def _plot(ts, a_half, a_bond, kopp):
    fig, ax = plt.subplots(figsize=(7.0, 4.5), dpi=300)
    fig.patch.set_facecolor('white'); ax.set_facecolor('white')
    ax.plot(ts, a_half, 'o-', color='C3', lw=1.6, ms=5,
            label='ExoColumn, single 60° (current panel b)')
    ax.plot(ts, a_bond, 's-', color='C1', lw=1.6, ms=5,
            label='ExoColumn, hemispheric Bond avg')
    ax.plot(ts, kopp, '^--', color='k', lw=1.4, ms=5,
            label='Kopparapu+2013 Fig 3b')
    # shade the two contributions at the plateau as a visual guide
    ax.annotate('', xy=(2000, a_bond[-1]), xytext=(2000, a_half[-1]),
                arrowprops=dict(arrowstyle='<->', color='C0', lw=1.2))
    ax.text(2010, 0.5*(a_half[-1]+a_bond[-1]), 'zenith\nangle', color='C0',
            fontsize=7, va='center')
    ax.annotate('', xy=(2000, kopp[-1]), xytext=(2000, a_bond[-1]),
                arrowprops=dict(arrowstyle='<->', color='C2', lw=1.2))
    ax.text(2010, 0.5*(a_bond[-1]+kopp[-1]), 'near-IR\nH₂O abs.', color='C2',
            fontsize=7, va='center')
    ax.set_xlabel('$T_s$ (K)'); ax.set_ylabel('Planetary albedo')
    ax.set_xlim(250, 2200); ax.set_ylim(0.14, 0.30)
    ax.legend(fontsize=8, framealpha=0.9)
    fig.tight_layout()
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       'diag_zenith_albedo.png')
    fig.savefig(out, dpi=150, facecolor='white')
    plt.close(fig)
    print(f"Saved: {out}")


def main():
    if not os.path.isfile(EXE):
        raise FileNotFoundError(EXE)
    orig = open(NML_PATH).read() if os.path.exists(NML_PATH) else None
    mu = np.array(MU_VALUES)
    try:
        print(f"Zenith-angle vs planetary albedo   (surface albedo {ALBEDO}, EOS {H2O_EOS})")
        print("A_Bond = 2*integral alpha(mu0) mu0 dmu0  (hemispheric / Kopparapu-equivalent)\n")
        print(f"  {'Ts[K]':>6}  {'alpha(0.5)':>10}  {'A_Bond':>8}  {'Kopparapu':>9}  {'a(0.5)-A_Bond':>13}")
        curves = {}
        a_half_l, a_bond_l, kopp_l = [], [], []
        for ts in TS_VALUES:
            alpha = np.array([alpha_at(ts, m) for m in mu])
            curves[ts] = alpha
            # single-angle value at mu0=0.5 (exact grid point)
            a_half = float(alpha[np.argmin(np.abs(mu - 0.5))])
            # hemispheric Bond albedo: 2 * integral_0^1 alpha(mu) mu dmu
            integ = _trapz(alpha * mu, mu)
            norm  = _trapz(mu, mu)        # = 0.5 over [0,1] if endpoints 0..1
            a_bond = integ / norm
            kref = KOPP_REF.get(ts)
            a_half_l.append(a_half); a_bond_l.append(a_bond); kopp_l.append(kref)
            kstr = f"{kref:.3f}" if kref is not None else "   -"
            print(f"  {ts:>6.0f}  {a_half:>10.3f}  {a_bond:>8.3f}  {kstr:>9}"
                  f"  {a_half - a_bond:>+13.3f}")
        print("\n  alpha(mu0) curves (rows=Ts, cols=mu0):")
        print("   mu0 ->  " + "  ".join(f"{m:5.3f}" for m in mu))
        for ts in TS_VALUES:
            print(f"   {ts:5.0f}K  " + "  ".join(f"{a:5.3f}" for a in curves[ts]))
        _plot(np.array(TS_VALUES), np.array(a_half_l), np.array(a_bond_l),
              np.array([k if k is not None else np.nan for k in kopp_l]))
    finally:
        if orig is not None:
            with open(NML_PATH, 'w') as f:
                f.write(orig)
        elif os.path.exists(NML_PATH):
            os.remove(NML_PATH)
    print("\nOriginal exocol_config.nml restored.")


if __name__ == '__main__':
    main()
