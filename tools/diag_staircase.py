#!/usr/bin/env python3
"""
diag_staircase.py — reproduce ExoColumn's cold-start column construction in pure
Python (NO ExoRT radiation call) to test whether the hz_inner staircase is a
vertical-grid discretization artifact.

For each Ts in the hz_inner sweep it rebuilds exactly what exocol_coldstart.F90
builds (variable_ps grid + Kasting IHZ dry+moist adiabat + cold-trap water) and
records the *integer-valued* quantities that snap to model layers:
  - k_top_conv : topmost tropospheric layer index (tmid > t_strato)
  - k_cond     : IHZ condensation level index
  - p_tropo    : pint(k_top_conv)  (the snapped tropopause pressure)
  - strat VMR  : q_cold_trap expressed as H2O volume mixing ratio

If the staircase in OLR/albedo/Seff is a grid artifact, these integer indices
will step at the same Ts values as the plotted sawtooth.
"""
import numpy as np
import matplotlib.pyplot as plt

# ---- constants (mirror shr_const_mod + exocol_convadj) ----
RGAS   = 6.02214e26 * 1.38065e-23   # J/K/kmole  = 8314.x
MWWV   = 18.016
RWV    = RGAS / MWWV
TKFRZ  = 273.16
LATVAP = 2.501e6
LATSUB = 3.337e5 + LATVAP
ES0    = 611.2
G      = 9.80616
NSUB   = 20

# composition (hz_inner): N2 .78 / O2 .21 / Ar .01 / CO2 3.3e-4
MW = dict(CO2=44.010, CH4=16.043, C2H6=30.069, H2=2.016,
          N2=28.013, O3=47.998, O2=31.999, AR=39.948)
CP = dict(CO2=844.0, N2=1039.0, O2=919.0, AR=520.3)
vmr = dict(CO2=3.3e-4, N2=0.78, O2=0.210, AR=0.01)
mwdry = sum(vmr[g] * MW[g] for g in vmr)
RD    = RGAS / mwdry
EPS   = MWWV / mwdry

T_STRATO = 200.0
P_TOP    = 1.0
PVER     = 70
PVERP    = PVER + 1


def esat_cc(T):
    Tu = min(max(T, 50.0), 5000.0)
    L = LATVAP if Tu >= TKFRZ else LATSUB
    return ES0 * np.exp((L / RWV) * (1.0 / TKFRZ - 1.0 / Tu))


def Lvap(T):
    return LATVAP if T >= TKFRZ else LATSUB


def malr(T, p):
    es = esat_cc(T)
    L = Lvap(T)
    ws = EPS * es / (p - es) if p > es else 1e-10
    cp = mwdry_cp
    gm = (G / cp) * (1.0 + L * ws / (RD * T)) / \
         (1.0 + L**2 * ws / (cp * RWV * T**2))
    return max(gm, 1e-3)


# mass-weighted cp (argon included)
mmr = {g: vmr[g] * MW[g] / mwdry for g in vmr}
mwdry_cp = sum(mmr[g] * CP[g] for g in vmr)


def build_column(ts):
    e_sfc = esat_cc(ts)
    ps = 1.0e5 + e_sfc                       # variable_ps
    pint = np.exp(np.log(P_TOP) +
                  np.arange(PVERP) / (PVERP - 1) * np.log(ps / P_TOP))
    pmid = 0.5 * (pint[:-1] + pint[1:])

    T_at_int = np.empty(PVERP)
    T_at_int[-1] = ts

    # --- IHZ dry+moist adiabat (mirror exocol_coldstart) ---
    e_s = min(esat_cc(ts), 0.99 * ps)
    f_vmr = e_s / ps
    ws_sfc = EPS * e_s / max(ps - e_s, 1.0)
    k_cond = 0
    if ts <= T_STRATO:
        T_at_int[:-1] = T_STRATO
    elif ws_sfc > 1.0:
        Mw_mix = f_vmr * MWWV + (1 - f_vmr) * mwdry
        w_h2o = f_vmr * MWWV / Mw_mix
        R_mix = w_h2o * RWV + (1 - w_h2o) * RD
        cp_mix = w_h2o * 1870.0 + (1 - w_h2o) * mwdry_cp
        kappa = R_mix / cp_mix
        k_cond = 1
        in_dry = True
        for k in range(PVER - 1, -1, -1):   # k index 0..pver-1 -> interfaces
            kk = k  # 0-based interface index
            if not in_dry:
                if T_at_int[kk + 1] <= T_STRATO:
                    T_at_int[kk] = T_STRATO
                    continue
                T_lev = T_at_int[kk + 1]; p_lev = pint[kk + 1]
                dl = np.log(pint[kk] / pint[kk + 1]) / NSUB
                for _ in range(NSUB):
                    Gm = malr(T_lev, p_lev)
                    es_ = min(esat_cc(T_lev), 0.99 * p_lev)
                    r_ = EPS * es_ / (p_lev - es_)
                    Tv = T_lev * (1 + r_ / EPS) / (1 + r_)
                    T_lev += Gm * RD * Tv / G * dl
                    p_lev *= np.exp(dl)
                    if T_lev <= T_STRATO:
                        T_lev = T_STRATO; break
                T_at_int[kk] = T_lev
            else:
                T_at_int[kk] = ts * (pint[kk] / ps) ** kappa
                if T_at_int[kk] <= T_STRATO:
                    T_at_int[kk] = T_STRATO
                es_ = min(esat_cc(T_at_int[kk]), 0.99 * pint[kk])
                r_ = EPS * es_ / max(pint[kk] - es_, 1.0)
                if r_ <= 1.0:
                    k_cond = kk
                    in_dry = False
    else:
        for k in range(PVER - 1, -1, -1):
            kk = k
            if T_at_int[kk + 1] <= T_STRATO:
                T_at_int[kk] = T_STRATO; continue
            T_lev = T_at_int[kk + 1]; p_lev = pint[kk + 1]
            dl = np.log(pint[kk] / pint[kk + 1]) / NSUB
            for _ in range(NSUB):
                Gm = malr(T_lev, p_lev)
                es_ = min(esat_cc(T_lev), 0.99 * p_lev)
                r_ = EPS * es_ / (p_lev - es_)
                Tv = T_lev * (1 + r_ / EPS) / (1 + r_)
                T_lev += Gm * RD * Tv / G * dl
                p_lev *= np.exp(dl)
                if T_lev <= T_STRATO:
                    T_lev = T_STRATO; break
            T_at_int[kk] = T_lev

    tmid = 0.5 * (T_at_int[:-1] + T_at_int[1:])
    pdel = pint[1:] - pint[:-1]

    # k_top_conv: first (from TOA) layer with tmid > t_strato  (1-based in F90)
    k_top = 0
    for k in range(PVER):
        if tmid[k] > T_STRATO:
            k_top = k + 1  # 1-based
            break
    if k_top > 0:
        es_tropo = esat_cc(T_STRATO)
        p_tropo = pint[k_top - 1]        # pint(k_top_conv), 1-based->0-based
        q_ct = EPS * es_tropo / (p_tropo - es_tropo)
        strat_vmr = q_ct * mwdry / MWWV
    else:
        p_tropo = np.nan; strat_vmr = 0.0; q_ct = 0.0

    # ---- reconstruct the actual h2ommr profile exactly as coldstart does ----
    h2o = np.empty(PVER)
    for k in range(PVER):
        if tmid[k] > T_STRATO:
            if ws_sfc > 1.0 and (k + 1) >= k_cond and k_cond > 0:
                # dry layer: constant mass fraction w_h2o
                Mw_mix = f_vmr * MWWV + (1 - f_vmr) * mwdry
                h2o[k] = f_vmr * MWWV / Mw_mix
            else:
                es_k = min(esat_cc(tmid[k]), 0.99 * pmid[k])
                h2o[k] = 1.0 * EPS * es_k / (pmid[k] - es_k)  # rh_init=1
        else:
            h2o[k] = 1.0 * q_ct   # rh_init=1

    # column water path [kg/m2] = sum q*pdel/g  (q≈mmr here)
    tcwv = np.sum(h2o * pdel) / G
    # water path ABOVE 100 hPa (the part that controls OLR most strongly)
    upper = pmid < 1.0e4
    wp_upper = np.sum(h2o[upper] * pdel[upper]) / G

    return dict(ps=ps, k_top=k_top, k_cond=k_cond, p_tropo=p_tropo,
                strat_vmr=strat_vmr, n_warm=int(np.sum(tmid > T_STRATO)),
                tcwv=tcwv, wp_upper=wp_upper)


def continuous_coldpoint(ts):
    """True cold-point pressure where the *continuous* adiabat = t_strato,
    integrated on a fine log-p grid (independent of the 70-layer interfaces).
    This is what the cold trap SHOULD use — a smooth function of Ts."""
    e_sfc = esat_cc(ts)
    ps = 1.0e5 + e_sfc
    if ts <= T_STRATO:
        return np.nan
    N = 4000
    p = np.exp(np.linspace(np.log(ps), np.log(P_TOP), N))
    # IHZ dry base if supercritical
    e_s = min(esat_cc(ts), 0.99 * ps)
    f_vmr = e_s / ps
    ws_sfc = EPS * e_s / max(ps - e_s, 1.0)
    T = ts
    if ws_sfc > 1.0:
        Mw_mix = f_vmr * MWWV + (1 - f_vmr) * mwdry
        w_h2o = f_vmr * MWWV / Mw_mix
        R_mix = w_h2o * RWV + (1 - w_h2o) * RD
        cp_mix = w_h2o * 1870.0 + (1 - w_h2o) * mwdry_cp
        kappa = R_mix / cp_mix
        dry = True
    else:
        dry = False
    p_prev, T_prev = ps, ts
    for i in range(1, N):
        if dry:
            T = ts * (p[i] / ps) ** kappa
            es_ = min(esat_cc(T), 0.99 * p[i])
            r_ = EPS * es_ / max(p[i] - es_, 1.0)
            if r_ <= 1.0:
                dry = False
        else:
            dl = np.log(p[i] / p[i - 1])
            Gm = malr(T, p[i - 1])
            es_ = min(esat_cc(T), 0.99 * p[i - 1])
            r_ = EPS * es_ / (p[i - 1] - es_)
            Tv = T * (1 + r_ / EPS) / (1 + r_)
            T = T + Gm * RD * Tv / G * dl
        if T <= T_STRATO:
            # linear interp in log-p between (p_prev,T_prev) and (p[i],T)
            f = (T_prev - T_STRATO) / (T_prev - T)
            return np.exp(np.log(p_prev) + f * (np.log(p[i]) - np.log(p_prev)))
        p_prev, T_prev = p[i], T
    return P_TOP


def main():
    TS = np.arange(200, 2205, 5, dtype=float)
    keys = ['k_top', 'k_cond', 'p_tropo', 'strat_vmr', 'n_warm',
            'tcwv', 'wp_upper']
    D = {k: [] for k in keys}
    for ts in TS:
        c = build_column(ts)
        for k in keys:
            D[k].append(c[k])
    for k in keys:
        D[k] = np.array(D[k], dtype=float)

    fig, ax = plt.subplots(2, 2, figsize=(9, 6.5), dpi=120)
    ax[0, 0].plot(TS, D['p_tropo'] / 100.0, color='C2')
    ax[0, 0].set_ylabel('p_tropo = pint(k_top_conv) [hPa]')
    ax[0, 0].set_title('Snapped tropopause pressure (the discontinuity)')

    ax[0, 1].plot(TS, D['k_cond'], color='C1', label='k_cond (dry/moist switch)')
    ax[0, 1].plot(TS, D['k_top'], color='C0', label='k_top_conv')
    ax[0, 1].set_ylabel('snapped layer index')
    ax[0, 1].set_title('Snapped layer indices vs Ts')
    ax[0, 1].legend(fontsize=8)

    ax[1, 0].semilogy(TS, D['tcwv'], color='C3')
    ax[1, 0].set_ylabel('total column water [kg m$^{-2}$]')
    ax[1, 0].set_title('Total water column (OLR proxy) vs Ts')

    # ---- the fix: cold trap from the CONTINUOUS cold-point pressure ----
    es_tropo = esat_cc(T_STRATO)
    q_snap = np.where(np.isfinite(D['p_tropo']),
                      EPS * es_tropo / (D['p_tropo'] - es_tropo), 0.0)
    p_cp = np.array([continuous_coldpoint(ts) for ts in TS])
    q_fix = np.where(np.isfinite(p_cp),
                     EPS * es_tropo / (p_cp - es_tropo), 0.0)
    ax[1, 1].semilogy(TS, q_snap, color='C4', lw=1.2,
                      label='current: pint(k_top_conv) (snapped)')
    ax[1, 1].semilogy(TS, q_fix, color='k', lw=1.2,
                      label='fix: interpolated cold-point p')
    ax[1, 1].set_ylabel('stratospheric H2O mmr (cold trap)')
    ax[1, 1].set_title('Cold-trap water: snapped vs continuous')
    ax[1, 1].legend(fontsize=8)
    ax[1, 1].set_xlim(700, 2200)

    for a in ax.flat:
        a.set_xlabel('Ts (K)')
        a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig('tools/diag_staircase.png', dpi=120)
    print('saved tools/diag_staircase.png')

    # report step count and locations in the 700-2200 K range
    seg = (TS >= 700)
    tsseg = TS[seg]
    for name in ['k_top', 'k_cond']:
        jumps = np.where(np.diff(D[name][seg]) != 0)[0]
        print(f'\n{name}: range {D[name].min():.0f}..{D[name].max():.0f}; '
              f'{len(jumps)} jumps in 700-2200 K', end='')
        if len(jumps) > 1:
            print(f'  (mean spacing {np.diff(tsseg[jumps]).mean():.0f} K)', end='')
    # fractional jump in upper-atm water at each k_top step
    dwp = np.abs(np.diff(D['wp_upper'])) / (D['wp_upper'][:-1] + 1e-30)
    big = dwp[(TS[:-1] >= 700)]
    print(f'\n\nupper-atm water column: typical step-to-step fractional jump '
          f'(700-2200 K) = {np.median(big)*100:.1f}%, max = {big.max()*100:.0f}%')


if __name__ == '__main__':
    main()
