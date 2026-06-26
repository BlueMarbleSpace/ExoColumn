#!/usr/bin/env bash
#
# sweep_all_masses.sh — build ExoColumn at each Kopparapu+2014 planetary mass
# (0.1, 1, 5 M_Earth -> EXO_G via the mass-radius relation), run the inner
# (runaway) + outer (maximum-greenhouse) multi-stellar HZ sweep for each, then
# restore the default Earth-gravity binary and draw the Figure-3 analogue.
#
# Surface gravity is a COMPILE-TIME constant in ExoRT (exo_g -> SHR_CONST_G), so
# every mass needs its own binary.  Run from anywhere; needs the Intel OneAPI
# runtime for the build:
#
#     source /opt/intel/oneapi/setvars.sh
#     bash reference/planet_mass/sweep_all_masses.sh
#
# Thin wrapper around `python hz_mass.py all` (n68 M->G sweep) followed by
# `addf` (the n84 F 7200 K endpoint), kept separate so the build environment
# (setvars.sh) is obviously the caller's responsibility.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"
python3 reference/planet_mass/hz_mass.py all     # n68 cool/solar ladder ×3 masses
python3 reference/planet_mass/hz_mass.py addf    # n84 F 7200 K endpoint ×3 masses
