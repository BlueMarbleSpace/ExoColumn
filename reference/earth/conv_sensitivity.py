#!/usr/bin/env python3
"""
conv_sensitivity.py  —  Sensitivity of the ExoColumn Earth equilibrium to the
choice of convective-adjustment scheme.

Holds the modern-Earth cold-start configuration (reference/earth/exocol_config.nml)
completely fixed — same insolation, surface albedo, composition, surface-flux
scheme, moisture treatment, and vertical grid (PVER=70, matching Fig. 1) — and
varies ONLY &exocol_nml::conv_scheme.  This isolates how much the equilibrium
temperature profile and surface temperature depend on the convective closure,
demonstrating the robustness of the production result.

Schemes compared:
    sbm    : simplified Betts-Miller (production default)
    zm     : Zhang-McFarlane soft adjustment (ExoCAM-consistent)
    moist  : RH-weighted local moist-adiabatic adjustment
    manabe : fixed 6.5 K/km lapse rate (Manabe & Wetherald 1967)

Run (sweep + plot):   python reference/earth/conv_sensitivity.py
Re-plot from cache:   CONV_REPLOT=1 python reference/earth/conv_sensitivity.py

Requires the PVER=70 binary (run/exocol.exe) and the Intel runtime on PATH
(source /opt/intel/oneapi/setvars.sh).  The sweep temporarily overwrites the
project-root exocol_config.nml and restores it on exit.
"""

import os
import re
import shutil
import subprocess
import numpy as np
import netCDF4
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
BASE_NML = os.path.join(HERE, "exocol_config.nml")      # fixed Earth reference config
ROOT_NML = os.path.join(ROOT, "exocol_config.nml")
EXE = os.path.join(ROOT, "run", "exocol.exe")
OUT_NC = os.path.join(ROOT, "iofiles", "exocol_out.nc")
SWEEPDIR = os.path.join(HERE, "conv_sensitivity")
CACHE = os.path.join(HERE, "conv_sensitivity.npz")
FIG_PDF = os.path.join(HERE, "conv_sensitivity.pdf")
FIG_PNG = os.path.join(HERE, "conv_sensitivity.png")

# (scheme key, display label, colour).  sbm drawn last/heaviest as the reference.
SCHEMES = [
    ("manabe", "Manabe 6.5 K/km", "#2ca02c"),
    ("moist",  "Moist adiabat",   "#1f77b4"),
    ("zm",     "Zhang-McFarlane", "#ff7f0e"),
    ("sbm",    "Betts-Miller (default)", "#d62728"),
]


def run_scheme(scheme):
    """Write the root config with conv_scheme=scheme, run the model, return the
    output NetCDF path (copied into SWEEPDIR) and the reported surface T."""
    with open(BASE_NML) as f:
        txt = f.read()
    txt = re.sub(r"conv_scheme\s*=\s*'[^']*'", f"conv_scheme     = '{scheme}'", txt, count=1)
    with open(ROOT_NML, "w") as f:
        f.write(txt)
    print(f"  running conv_scheme='{scheme}' ...", flush=True)
    res = subprocess.run([EXE], cwd=ROOT, capture_output=True, text=True)
    if res.returncode != 0:
        print(res.stdout[-2000:])
        print(res.stderr[-2000:])
        raise RuntimeError(f"exocol.exe failed for scheme {scheme}")
    dst = os.path.join(SWEEPDIR, f"{scheme}.nc")
    shutil.copyfile(OUT_NC, dst)
    # Pull the final Ts and step count from stdout if printed (diagnostic only).
    ts = None
    for line in res.stdout.splitlines()[::-1]:
        m = re.search(r"[Tt]s\s*=\s*([0-9.]+)", line)
        if m:
            ts = float(m.group(1))
            break
    return dst, ts


def load_profile(nc):
    d = netCDF4.Dataset(nc)
    pmid = d.variables["pmid"][:].astype(float)      # Pa
    tmid = d.variables["tmid"][:].astype(float)      # K
    ts = float(d.variables["ts"][:]) if "ts" in d.variables else float("nan")
    q = d.variables["h2ommr"][:].astype(float) if "h2ommr" in d.variables else None
    d.close()
    return pmid, tmid, ts, q


def sweep():
    os.makedirs(SWEEPDIR, exist_ok=True)
    backup = ROOT_NML + ".convsens_bak"
    shutil.copyfile(ROOT_NML, backup)
    data = {}
    try:
        for key, _, _ in SCHEMES:
            nc, ts_stdout = run_scheme(key)
            pmid, tmid, ts, q = load_profile(nc)
            data[key] = dict(pmid=pmid, tmid=tmid, ts=ts, q=q)
            print(f"    -> Ts = {ts:.3f} K", flush=True)
    finally:
        shutil.copyfile(backup, ROOT_NML)
        os.remove(backup)
    # Cache (npz can't hold dict-of-dict directly; flatten).
    flat = {}
    for key in data:
        for fld in ("pmid", "tmid", "q"):
            if data[key][fld] is not None:
                flat[f"{key}__{fld}"] = data[key][fld]
        flat[f"{key}__ts"] = np.array(data[key]["ts"])
    np.savez(CACHE, **flat)
    return data


def load_cache():
    z = np.load(CACHE)
    data = {}
    for key, _, _ in SCHEMES:
        data[key] = dict(
            pmid=z[f"{key}__pmid"], tmid=z[f"{key}__tmid"],
            ts=float(z[f"{key}__ts"]),
            q=z[f"{key}__q"] if f"{key}__q" in z.files else None,
        )
    return data


def plot(data):
    plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.8,
                         "figure.facecolor": "white", "savefig.facecolor": "white"})
    fig, (ax, axd) = plt.subplots(1, 2, figsize=(7.0, 4.5),
                                  gridspec_kw=dict(width_ratios=[2, 1], wspace=0.05))

    ref = data["sbm"]  # difference reference
    for key, label, color in SCHEMES:
        d = data[key]
        p_hpa = d["pmid"] / 100.0
        lw = 1.8 if key == "sbm" else 1.2
        ax.plot(d["tmid"], p_hpa, "-", color=color, lw=lw,
                label=f"{label}  ($T_s$={d['ts']:.1f} K)")
        # Temperature difference from the SBM reference, interpolated to its grid.
        dT = np.interp(np.log(ref["pmid"]), np.log(d["pmid"]), d["tmid"]) - ref["tmid"]
        axd.plot(dT, ref["pmid"] / 100.0, "-", color=color, lw=lw)

    ax.set_yscale("log")
    ax.invert_yaxis()
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel("Pressure (hPa)")
    ax.legend(frameon=False, fontsize=7.5, loc="upper right")
    ax.set_ylim(1000, p_hpa.min())

    axd.set_yscale("log")
    axd.invert_yaxis()
    axd.set_ylim(1000, p_hpa.min())
    axd.axvline(0, color="0.7", lw=0.7)
    axd.set_xlabel(r"$T - T_{\rm SBM}$ (K)")
    axd.tick_params(labelleft=False)

    fig.savefig(FIG_PDF, bbox_inches="tight")
    fig.savefig(FIG_PNG, dpi=200, bbox_inches="tight")
    print(f"wrote {FIG_PDF}\nwrote {FIG_PNG}")
    # Summary table.
    print("\n scheme           Ts [K]   dTs vs SBM")
    for key, label, _ in SCHEMES:
        print(f"  {label:<24s} {data[key]['ts']:7.3f}  {data[key]['ts']-ref['ts']:+7.3f}")


if __name__ == "__main__":
    if os.environ.get("CONV_REPLOT") and os.path.exists(CACHE):
        data = load_cache()
    else:
        data = sweep()
    plot(data)
