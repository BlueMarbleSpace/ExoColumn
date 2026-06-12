#!/usr/bin/env python3
"""
lbl_sw_benchmark.py — line-by-line clear-sky SHORTWAVE benchmark for an
ExoColumn column: the SW twin of lbl_olr_benchmark.py, built to referee the
IHZ planetary-albedo offset vs Kopparapu et al. (2013) (~+0.012 with the BPS
continuum at the limits; the 220 K dry-end albedo agreement of ±0.003 already
pins the gap on near-IR H2O absorption — see reference/moist_runaway/README).

Design — everything except the gas opacity replicates ExoRT, so differences
vs the model's band_swup_toa isolate the H2O/CO2 absorption data:
  * Lines: RADIS/HITRAN H2O + CO2 (2000-20000 cm-1, i.e. 0.5-5 um, where all
    significant H2O/CO2 SW absorption lives), plus the MT_CKD continuum (same
    AER port as the LW benchmark).
  * Rayleigh: ExoRT's own formulas verbatim (Vardavas & Carver/Allen for N2 &
    CO2, Bucholtz-index for H2O; rayleigh_data.F90 constants), evaluated at
    the BAND MIDPOINT and held constant within each band, exactly as
    calc_gasopd does — so the Rayleigh discretization bucket is zero.
  * Two-stream: Toon et al. (1989) QUADRATURE coefficients (mu1 = 1/sqrt(3)),
    the choice ExoRT makes for the solar stream (exo_init_ref.F90: U1Isol =
    sqrt(3)); homogeneous-layer analytic R/T + adding method with a Lambertian
    surface (asdir = asdif = 0.32).  Conservative-scattering guard:
    omega <= 1 - 1e-6.
  * Zenith treatment: the same 6-point Gauss-Legendre quadrature in mu over
    [0,1] as &exocol_nml::sw_zenith_quad with sw_nquad=6, flux-weighted.
  * Solar spectrum: ExoRT's own per-band incident fluxes (band_swdn_toa from
    the input file) are used as band weights for BOTH models, with a 5778-K
    blackbody SHAPE inside each band for the LBL — so the stellar input is
    identical at band resolution by construction and only the within-band
    shape (a second-order effect over a few-hundred-cm-1 band) differs.

Output: per-n68-band albedo, LBL vs ExoRT, and totals weighted by ExoRT's
band_swdn_toa.  Bands above 20000 cm-1 are treated as gas-free
(Rayleigh + surface only) — H2O/CO2 absorption there is negligible.

Usage:
  python tools/lbl_sw_benchmark.py [exocol_out.nc] [--nlay 50] [--wstep 0.02]
                                   [--no-gas]   # solver validation: Rayleigh only
"""
import os, sys, argparse
import numpy as np
import netCDF4 as nc

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from lbl_olr_benchmark import load_column, MTCKD, H, C, KB

SQRT3 = np.sqrt(3.0)
ALB_SFC = 0.32
OMEGA_MAX = 1.0 - 1.0e-6
TSUN = 5778.0

# 6-point Gauss-Legendre on mu in [0,1] (zenith_quad_nodes convention)
GL_X = np.array([-0.9324695142031521, -0.6612093864662645, -0.2386191860831969,
                  0.2386191860831969,  0.6612093864662645,  0.9324695142031521])
GL_W = np.array([0.1713244923791704, 0.3607615730481386, 0.4679139345726910,
                 0.4679139345726910, 0.3607615730481386, 0.1713244923791704])
MU_NODES = (GL_X + 1.0) / 2.0
W_NODES = GL_W / 2.0

# ExoRT rayleigh_data.F90 constants
DEL_N2, DEL_CO2, DEL_H2O = 0.0305, 0.0805, 0.17
RAYLA_N2, RAYLB_N2 = 29.06, 7.70
RAYLA_CO2, RAYLB_CO2 = 43.90, 6.40


def sigma_rayleigh_band(wn_mid):
    """ExoRT calc_gasopd Rayleigh cross sections [cm2/molec] at band-midpoint
    wavenumber wn_mid [cm-1].  Returns (sigma_N2, sigma_CO2, sigma_H2O)."""
    wl = 1.0e4 / wn_mid                      # microns
    depolN2 = (6 + 3 * DEL_N2) / (6 - 7 * DEL_N2)
    depolCO2 = (6 + 3 * DEL_CO2) / (6 - 7 * DEL_CO2)
    depolH2O = (6 + 3 * DEL_H2O) / (6 - 7 * DEL_H2O)
    allenN2 = (1.0e-5 * RAYLA_N2 * (1.0 + 1.0e-3 * RAYLB_N2 / wl**2))**2
    allenCO2 = (1.0e-5 * RAYLA_CO2 * (1.0 + 1.0e-3 * RAYLB_CO2 / wl**2))**2
    sigN2 = 4.577e-21 / wl**4 * depolN2 * allenN2
    sigCO2 = 4.577e-21 / wl**4 * depolCO2 * allenCO2
    ns = (5791817. / (238.0185 - (1. / wl)**2)
          + 167909. / (57.362 - (1. / wl)**2)) / 1.0e8     # Bucholtz (1995)
    r = 0.85 * ns
    sigH2O = 4.577e-21 * depolH2O * r**2 / wl**4
    return sigN2, sigCO2, sigH2O


def planck_shape(wn, T=TSUN):
    """Un-normalized blackbody flux density vs wavenumber (for within-band
    weighting of the stellar spectrum)."""
    nu = wn * 100.0
    return nu**3 / np.expm1(H * C * nu / (KB * T))


# ----------------------------------------------------------------------------
# Toon89-quadrature two-stream: homogeneous-layer R/T + adding
# ----------------------------------------------------------------------------
def layer_rt(tau, omega, mu0):
    """Reflection/transmission of one homogeneous layer (g = 0, Rayleigh) for
    diffuse and collimated (mu0) incidence, Toon89 quadrature coefficients.
    Returns (Rdif, Tdif, Rcol, Tcol_dif, tdir); all vectorized over tau/omega."""
    w = np.minimum(omega, OMEGA_MAX)
    g1 = SQRT3 * (2.0 - w) / 2.0
    g2 = SQRT3 * w / 2.0
    g3 = 0.5                                  # (1 - sqrt(3) g mu0)/2 with g=0
    g4 = 0.5
    lam = np.sqrt(np.maximum(g1**2 - g2**2, 1e-30))
    Gam = g2 / (g1 + lam)
    E = np.exp(-np.minimum(lam * tau, 200.0))
    den = 1.0 - (Gam * E)**2
    Rdif = Gam * (1.0 - E**2) / den
    Tdif = (1.0 - Gam**2) * E / den

    # particular (direct-beam) solution; unit incident DIRECT FLUX mu0*F0*pi=1
    # Toon89 eq 23/24 with pi*F0 = 1/mu0:
    mu0i = 1.0 / mu0
    denom = lam**2 - mu0i**2
    # resonance guard
    denom = np.where(np.abs(denom) < 1e-12, np.sign(denom + 1e-30) * 1e-12, denom)
    piF0 = mu0i                                # so that mu0*pi*F0 = 1
    Cp0 = w * piF0 * ((g1 - mu0i) * g3 + g4 * g2) / denom
    Cm0 = w * piF0 * ((g1 + mu0i) * g4 + g2 * g3) / denom
    tdir = np.exp(-np.minimum(tau * mu0i, 200.0))
    CpB = Cp0 * tdir
    CmB = Cm0 * tdir

    # homogeneous coefficients from BCs F-(0) = 0, F+(tau_b) = 0
    # F+ = A' e^{-lam(tau-t)} + B Gam e^{-lam t} + Cp(t)
    # F- = A' Gam e^{-lam(tau-t)} + B e^{-lam t} + Cm(t)
    det = 1.0 - (Gam * E)**2
    Ap = (-CpB + Gam * E * Cm0) / det
    Bh = (-Cm0 + Gam * E * CpB) / det
    Rcol = Ap * E + Bh * Gam + Cp0            # F+(0)
    Tcol = Ap * Gam + Bh * E + CmB            # F-(tau_b), diffuse part
    return Rdif, Tdif, np.maximum(Rcol, 0.0), np.maximum(Tcol, 0.0), tdir


def adding_albedo(tau_layers, omega_layers, mu0):
    """TOA albedo for collimated unit flux at mu0 over a Lambertian surface,
    layers TOA->surface, via the adding method.  Vectorized over wavenumber."""
    nlay = tau_layers.shape[0]
    # cumulative slab (initialized to the TOP layer)
    Rd_a, Td_a, Rc_a, Tc_a, t_a = layer_rt(tau_layers[0], omega_layers[0], mu0)
    Rd_b_a = Rd_a            # from-below diffuse reflectance (homog: same)
    for k in range(1, nlay):
        Rd, Td, Rc, Tc, td = layer_rt(tau_layers[k], omega_layers[k], mu0)
        D = 1.0 / (1.0 - Rd_b_a * Rd)
        # collimated: beam t_a hits layer k; diffuse Tc_a hits from above
        u0 = t_a * Rc + Tc_a * Rd
        u = u0 * D
        d = Tc_a + Rd_b_a * u
        Rc_ab = Rc_a + Td_a * u
        Tc_ab = t_a * Tc + Td * d
        t_ab = t_a * td
        # diffuse from above / from below
        Rd_ab = Rd_a + Td_a**2 * Rd * D
        Td_ab = Td_a * Td * D
        Rd_b_ab = Rd + Td**2 * Rd_b_a * D
        Rc_a, Tc_a, t_a = Rc_ab, Tc_ab, t_ab
        Rd_a, Td_a, Rd_b_a = Rd_ab, Td_ab, Rd_b_ab
    # surface: Lambertian, reflects direct and diffuse with ALB_SFC
    D = 1.0 / (1.0 - Rd_b_a * ALB_SFC)
    u = (t_a * ALB_SFC + Tc_a * ALB_SFC) * D
    return Rc_a + Td_a * u


# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('ncfile', nargs='?',
                    default=os.path.join(os.path.dirname(HERE), 'iofiles', 'exocol_out.nc'))
    ap.add_argument('--nlay', type=int, default=50)
    ap.add_argument('--wstep', type=float, default=0.02)
    ap.add_argument('--gasmax', type=float, default=20000.0)
    ap.add_argument('--gasmin', type=float, default=2000.0)
    ap.add_argument('--no-gas', action='store_true',
                    help='Rayleigh-only (two-stream solver validation)')
    ap.add_argument('--out', default=os.path.join(HERE, 'lbl_sw_benchmark'))
    args = ap.parse_args()

    lay = load_column(args.ncfile, args.nlay)
    band = lay['band']
    with nc.Dataset(args.ncfile) as ds:
        bswup = np.array(ds['band_swup_toa'][:])
        bswdn = np.array(ds['band_swdn_toa'][:])
    edges = band['edges']
    nb = len(bswup)
    N_n2 = lay['N_dry'] * 0.99967           # pure-N2+CO2 background of the sweep

    print(f"column: Ts={lay['ts']:.1f} K, {len(lay['T'])} layers")

    # --- gas optical depths on the fine grid (2000-20000 cm-1) ---
    if not args.no_gas:
        from lbl_olr_benchmark import layer_tau
        wn, tau_h2o, tau_co2, tau_cont = layer_tau(lay, args.gasmin, args.gasmax,
                                                   args.wstep)
        tau_gas = (tau_h2o + tau_co2 + tau_cont).astype(np.float32)
        del tau_h2o, tau_co2, tau_cont
    else:
        wn = np.arange(args.gasmin, args.gasmax, 10.0)
        tau_gas = np.zeros((len(lay['T']), len(wn)), np.float32)

    shape = planck_shape(wn)

    print("\nband-by-band two-stream (6-node GL zenith average):")
    print(f"{'band cm-1':>15} {'alb_LBL':>8} {'alb_ExoRT':>9} {'dalb':>7} {'Sdn_w':>8}")
    alb_lbl = np.full(nb, np.nan)
    alb_exo = bswup / np.maximum(bswdn, 1e-30)
    for i in range(nb):
        n1, n2 = edges[i], edges[i+1]
        wmid = 0.5 * (n1 + n2)
        sigN2, sigCO2, sigH2O = sigma_rayleigh_band(wmid)
        tau_ray_lay = (sigN2 * N_n2 + sigCO2 * lay['N_co2']
                       + sigH2O * lay['N_h2o'])           # per layer, const in band
        if n2 <= args.gasmin:
            # thermal-IR band: negligible solar flux; excluded from totals
            continue
        elif n1 >= args.gasmax:
            # gas-free band (vis/UV): Rayleigh + surface only; flat in band
            tau_tot = tau_ray_lay[:, None].astype(np.float64)
            omega = np.ones_like(tau_tot) * 1.0           # pure scattering
            R = np.zeros(1)
            for mu0, wq in zip(MU_NODES, W_NODES):
                R = R + wq * mu0 * adding_albedo(tau_tot, omega, mu0)
            alb_lbl[i] = float(R[0] / np.sum(W_NODES * MU_NODES))
        else:
            m = (wn >= n1) & (wn < n2)
            if m.sum() < 3:
                continue
            tg = tau_gas[:, m].astype(np.float64)
            tau_tot = tg + tau_ray_lay[:, None]
            omega = tau_ray_lay[:, None] / np.maximum(tau_tot, 1e-30)
            wgt = shape[m]
            Rnum = np.zeros(m.sum())
            for mu0, wq in zip(MU_NODES, W_NODES):
                Rnum = Rnum + wq * mu0 * adding_albedo(tau_tot, omega, mu0)
            Rb = Rnum / np.sum(W_NODES * MU_NODES)         # flux-weighted albedo(nu)
            alb_lbl[i] = float(np.sum(Rb * wgt) / np.sum(wgt))
        if bswdn[i] > 1e-3:
            print(f"{n1:7.0f}-{n2:7.0f} {alb_lbl[i]:8.4f} {alb_exo[i]:9.4f} "
                  f"{alb_exo[i]-alb_lbl[i]:+7.4f} {bswdn[i]:8.2f}")

    # totals weighted by ExoRT's incident band fluxes (identical input spectrum)
    mfin = np.isfinite(alb_lbl) & (bswdn > 1e-6)
    A_lbl = np.sum(alb_lbl[mfin] * bswdn[mfin]) / np.sum(bswdn[mfin])
    A_exo = np.sum(alb_exo[mfin] * bswdn[mfin]) / np.sum(bswdn[mfin])
    asr_lbl = np.sum((1 - alb_lbl[mfin]) * bswdn[mfin])
    asr_exo = np.sum((1 - alb_exo[mfin]) * bswdn[mfin])
    print(f"\n=== totals (ExoRT band_swdn weights, surface albedo {ALB_SFC}) ===")
    print(f"  planetary albedo:  LBL = {A_lbl:.4f}   ExoRT n68 = {A_exo:.4f}"
          f"   (ExoRT - LBL = {A_exo-A_lbl:+.4f})")
    print(f"  absorbed SW:       LBL = {asr_lbl:.2f}  ExoRT n68 = {asr_exo:.2f} W/m2")

    np.savez(args.out + '.npz', edges=edges, alb_lbl=alb_lbl, alb_exo=alb_exo,
             bswdn=bswdn, bswup=bswup, ts=lay['ts'])
    print(f"saved {args.out}.npz")


if __name__ == '__main__':
    main()
