#!/usr/bin/env python3
"""
check_steam_adiabat.py — validate the Kasting Appendix-A adiabat evaluator.

Two checks for the SATURATED branch (A4/A5):
  (1) transcription: an independent Python reimplementation of the same
      A4/A5 assembly (using the `iapws` package for water sat properties)
      must reproduce the Fortran dlnT/dlnP.
  (2) physics: Kasting's ideal-gas analogs A7/A8 must agree with the
      non-ideal A4/A5 in the near-ideal regime (low T) and diverge in the
      steam-dominated regime (high T / near critical) — confirming both the
      implementation and that non-ideality matters.
And for the UNSATURATED branch (A11/A12): transcription check.

Usage:  /tmp/.../test_steam_adiabat | python tools/check_steam_adiabat.py
"""
import sys
from math import log
from iapws import IAPWS95

R   = 8.314462618        # J/mol/K
MV  = 0.018015268        # kg/mol  (water)
MN  = 0.02897            # kg/mol  (mwdry 28.97)
RV  = R/MV               # 461.5
CC  = 4218.0             # liquid-water cp [J/kg/K] (condensate)

TC  = 647.096

def sat_props(T):
    g = IAPWS95(T=T, x=1.0); l = IAPWS95(T=T, x=0.0)
    return g.P*1e6, g.rho, g.s*1e3, l.s*1e3          # Psat, rhov, sv, sc

def sat_dln(T):
    dT  = max(1e-2, 1e-3*T)
    Thi = min(T+dT, TC-1e-3); Tlo = max(T-dT, 50.0)
    Ph, rvh, svh, sch = sat_props(Thi)
    Pl, rvl, svl, scl = sat_props(Tlo)
    den = log(Thi)-log(Tlo)
    return ((log(Ph)-log(Pl))/den, (log(rvh)-log(rvl))/den,
            (svh-svl)/den, (sch-scl)/den)

def latent(T):
    Ps, rhov, sv, sc = sat_props(T)
    return T*(sv-sc)                                  # L = T (sv - sc) [J/kg]

# ---- (1) independent A4/A5 reimplementation ----
def dlnTdlnP_sat_py(T, P, Rd, cpdry):
    Psat, rhov, sv, sc = sat_props(T)
    dlnPsat, dlnrhov, dsv, dsc = sat_dln(T)
    Pv = Psat; Pn = max(P-Pv, 1e-6*P)
    rhon = Pn/(Rd*T); alfav = rhov/rhon; Cvn = cpdry-Rd
    num = Rd*dlnrhov - Cvn - alfav*dsv
    den = alfav*(sv-sc) + Rd
    dlnalfav = num/den
    dlnPn = 1.0 + dlnrhov - dlnalfav
    dlnP = (Pv/P)*dlnPsat + (Pn/P)*dlnPn
    return 1.0/dlnP

# ---- (2) Kasting ideal-gas analogs A7/A8 ----
def dlnTdlnP_sat_ideal(T, P, Rd, cpdry):
    Psat = IAPWS95(T=T, x=1.0).P*1e6
    Pv = Psat; Pn = max(P-Pv, 1e-6*P)
    alfav = (Pv*MV)/(Pn*MN)                           # ideal mass ratio rho_v/rho_n
    L  = latent(T)
    dL = (latent(min(T+0.5,TC-1e-3)) - latent(max(T-0.5,50)))/ \
         (min(T+0.5,TC-1e-3) - max(T-0.5,50))          # dL/dT
    gamma = cpdry + alfav*(CC - L/T + dL)
    Rn = R/MN
    dadlnT = ((MV/MN)*L/T - gamma)/(L/T + Rn/alfav)    # A7  (d alpha_v/d ln T)
    dlnalfav = dadlnT/alfav
    dlnPdlnT = L/(RV*T) - dlnalfav/(1.0 + alfav*MN/MV) # A8
    return 1.0/dlnPdlnT

# ---- unsaturated A11/A12 reimplementation ----
def dlnTdlnP_dry_py(T, P, alfav_fixed, Rd, cpdry):
    # solve rho_v(Pv,T) = alfav_fixed*(P-Pv)/(Rd T)
    Pcap = 0.999*P
    if T < TC:
        Pcap = min(Pcap, 0.999*IAPWS95(T=T, x=1.0).P*1e6)
    lo, hi = 1e-3, Pcap
    for _ in range(100):
        mid = 0.5*(lo+hi)
        st = IAPWS95(P=mid/1e6, T=T)
        g = st.rho - alfav_fixed*(P-mid)/(Rd*T)
        if g > 0: hi = mid
        else:     lo = mid
        if hi-lo < 1e-7*P: break
    Pv = 0.5*(lo+hi); Pn = P-Pv
    st = IAPWS95(P=Pv/1e6, T=T)
    dlnPv_s = st.cp*1e3*st.rho/(Pv*st.alfav)
    dlnPn_s = cpdry/Rd
    dlnP = (Pv/P)*dlnPv_s + (Pn/P)*dlnPn_s
    return 1.0/dlnP

def main():
    Rd, cpdry = 287.0, 1004.0
    worst = 0.0; nfail = 0
    print(f"{'state':>22} | {'Fortran':>11} {'PyA4/A5':>11} {'relerr':>9} | "
          f"{'idealA7/A8':>11} {'nonid/ideal-1':>13}")
    for line in (open(sys.argv[1]) if len(sys.argv)>1 else sys.stdin):
        t = line.split()
        if not t: continue
        if t[0] == 'SATAD':
            T, P, dlF, lapseF = map(float, t[1:5])
            dlPy = dlnTdlnP_sat_py(T, P, Rd, cpdry)
            dlId = dlnTdlnP_sat_ideal(T, P, Rd, cpdry)
            e = abs(dlF-dlPy)/abs(dlPy)
            worst = max(worst, e)
            dev = dlF/dlId - 1.0
            print(f"sat T={T:6.1f} P={P:8.2e} | {dlF:11.5f} {dlPy:11.5f} {e:9.1e} | "
                  f"{dlId:11.5f} {dev:+13.2%}")
            if e > 1e-5: nfail += 1
        elif t[0] == 'DRYAD':
            T, P, alfa, dlF, lapseF = map(float, t[1:6])
            dlPy = dlnTdlnP_dry_py(T, P, alfa, Rd, cpdry)
            e = abs(dlF-dlPy)/abs(dlPy)
            worst = max(worst, e)
            print(f"dry T={T:6.1f} P={P:8.2e} | {dlF:11.5f} {dlPy:11.5f} {e:9.1e} | "
                  f"{'(alpha_v='+format(alfa,'.3f')+')':>26}")
            if e > 1e-4: nfail += 1
    print("-"*80)
    print(f"worst Fortran-vs-Python relerr: {worst:.2e}")
    print("Physics: nonideal/ideal deviation should be ~0 at low T and grow "
          "to tens of % near the 647 K critical point.")
    if nfail:
        print(f"FAILED: {nfail} transcription mismatch(es)"); sys.exit(1)
    print("TRANSCRIPTION OK")

if __name__ == '__main__':
    main()
