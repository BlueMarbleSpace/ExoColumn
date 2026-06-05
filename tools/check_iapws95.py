#!/usr/bin/env python3
"""
check_iapws95.py — validate the native Fortran IAPWS-95 port.

Reads the stdout of test/test_iapws95 (piped in or via a file argument),
recomputes each state with the reference `iapws` Python package, and reports
the relative error per property.  Exit status 0 iff every property is within
tolerance.

Usage:
    /tmp/test_iapws95 | python tools/check_iapws95.py
    python tools/check_iapws95.py fort_out.txt
"""
import sys
import numpy as np
from iapws import IAPWS95

# IAPWS-95 specific gas constant differs slightly between the canonical release
# (R = 0.46151805 kJ/kgK) and the iapws package (R_molar/M).  That ~4e-5 offset
# bounds the achievable agreement on R-scaled properties.
RTOL_PRIMARY = 3e-4    # P, s, h, u, cv, cp, w, alfav, Z, densities
RTOL_DERIV   = 1e-3    # latent heat (difference of two large enthalpies)

def relerr(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    scale = np.maximum(np.abs(b), 1e-30)
    return np.abs(a - b) / scale

def sc(x, fac=1.0):
    """Scale a reference property, tolerating None (two-phase / undefined)."""
    return None if x is None else x*fac

def main():
    src = open(sys.argv[1]) if len(sys.argv) > 1 else sys.stdin
    worst = 0.0
    nfail = 0
    for line in src:
        tok = line.split()
        if not tok:
            continue
        tag = tok[0]
        if tag == 'RHOT':
            rho, T, P, s, h, u, cv, cp, w, alfav, Z = map(float, tok[1:12])
            ref = IAPWS95(rho=rho, T=T)
            cmp = [
                ('P',    P,     sc(ref.P, 1e6)),
                ('s',    s,     sc(ref.s, 1e3)),
                ('h',    h,     sc(ref.h, 1e3)),
                ('u',    u,     sc(ref.u, 1e3)),
                ('cv',   cv,    sc(ref.cv, 1e3)),
                ('cp',   cp,    sc(ref.cp, 1e3)),
                ('w',    w,     sc(ref.w)),
                ('alfav',alfav, sc(ref.alfav)),
                ('Z',    Z,     sc(ref.Z)),
            ]
            label = f"RHOT rho={rho:9.4g} T={T:7.2f}"
            tol = RTOL_PRIMARY
        elif tag == 'SAT':
            ok = tok[1]
            T, Psat, rv, rl, sv, sl, L = map(float, tok[2:9])
            g = IAPWS95(T=T, x=1.0); l = IAPWS95(T=T, x=0.0)
            cmp = [
                ('Psat', Psat, g.P*1e6),
                ('rhov', rv,   g.rho),
                ('rhol', rl,   l.rho),
                ('sv',   sv,   g.s*1e3),
                ('sl',   sl,   l.s*1e3),
                ('L',    L,    (g.h-l.h)*1e3),
            ]
            label = f"SAT  T={T:7.2f} ok={ok}"
            tol = RTOL_DERIV
        elif tag == 'PT':
            ok = tok[1]
            P, T, rho, Pback = map(float, tok[2:6])
            ref = IAPWS95(P=P/1e6, T=T)
            cmp = [('rho', rho, ref.rho), ('P(roundtrip)', Pback, P)]
            label = f"PT   P={P:9.3g} T={T:7.2f} ok={ok}"
            tol = RTOL_PRIMARY
        else:
            continue

        errs = []
        bad = []
        for name, got, exp in cmp:
            if exp is None:
                print(f"         (skip {name}: reference is None — two-phase boundary)")
                continue
            e = float(relerr(got, exp))
            errs.append(e)
            if e > tol:
                bad.append(f"{name}: got {got:.7g} exp {exp:.7g} (relerr {e:.2e})")
        mx = max(errs) if errs else 0.0
        worst = max(worst, mx)
        flag = 'OK ' if not bad else 'FAIL'
        print(f"[{flag}] {label}  max relerr {mx:.2e}")
        for b in bad:
            print(f"         {b}")
        if bad:
            nfail += 1

    print("-"*60)
    print(f"worst relative error over all states: {worst:.3e}")
    if nfail:
        print(f"FAILED: {nfail} state(s) out of tolerance")
        sys.exit(1)
    print("ALL STATES WITHIN TOLERANCE")

if __name__ == '__main__':
    main()
