#!/usr/bin/env python3
"""Map the valid vs extrapolated regime of the hz_inner sweep to explain the
departures from Kopparapu+2013, independent of the (cosmetic) staircase.

Overlays on the OLR/ASR/albedo curves:
  - S-N runaway OLR limit (282 W/m²): OLR should asymptote to this, never exceed/dip.
  - ps(Ts) = 1 bar N2 + esat(Ts): variable_ps grows WITHOUT BOUND (unlimited water).
  - Earth ocean inventory (~270 bar of H2O): above this, variable_ps is unphysical
    (a real ocean fully evaporates; ps then saturates and Ts rises at fixed water mass).
  - Ts where the surface adiabat pushes layer T past the ExoRT k-table ceiling (500 K):
    above this OLR/SW opacity is frozen/extrapolated.
"""
import numpy as np
import matplotlib.pyplot as plt
from diag_staircase import esat_cc

d = np.loadtxt('/tmp/hz_data.txt')
TS, OLR, ASR, ALB, SEFF = d.T

ps_bar = np.array([(1.0e5 + esat_cc(t)) / 1.0e5 for t in TS])   # total ps [bar]
OCEAN_BAR = 270.0     # ~Earth ocean as surface pressure of steam
SN = 282.0

# Ts where esat alone = ocean inventory (variable_ps becomes unphysical)
ts_ocean = TS[np.argmin(np.abs(ps_bar - 1.0 - OCEAN_BAR))]
# Ts where ExoRT T-ceiling (500 K) is first exceeded at the surface
ts_ceiling = TS[np.argmin(np.abs(TS - 500.0))]

fig, ax = plt.subplots(3, 1, figsize=(8.5, 9), dpi=120, sharex=True)

ax[0].plot(TS, OLR, 'C3', label='OLR')
ax[0].plot(TS, ASR, 'C0', label='ASR')
ax[0].axhline(SN, color='k', ls='--', lw=0.8, label=f'S-N limit ({SN:.0f})')
ax[0].set_ylabel('flux (W/m²)'); ax[0].legend(fontsize=8)
ax[0].set_title('OLR should rise then PLATEAU at S-N (282); ours peaks→dips→rises')

ax[1].plot(TS, ALB, 'C2')
ax[1].set_ylabel('planetary albedo')
ax[1].set_title('Kopparapu (G-star): albedo LOW & ~flat/decreasing; ours dips to 0.17 then rises to 0.42')

ax[2].semilogy(TS, ps_bar, 'C4')
ax[2].axhline(OCEAN_BAR, color='k', ls=':', lw=1.0, label=f'~Earth ocean ({OCEAN_BAR:.0f} bar)')
ax[2].set_ylabel('surface pressure (bar)'); ax[2].set_xlabel('Ts (K)')
ax[2].legend(fontsize=8)
ax[2].set_title('variable_ps grows WITHOUT bound (unlimited water)')

for a in ax:
    a.axvspan(ts_ceiling, TS[-1], color='red', alpha=0.06)
    a.axvline(ts_ceiling, color='red', lw=0.8, alpha=0.6)
    a.axvline(ts_ocean, color='purple', lw=0.8, alpha=0.7)
ax[0].text(ts_ceiling+10, ax[0].get_ylim()[0]+15,
           'k-table T-ceiling\n(>500 K extrapolated)', color='red', fontsize=7.5)
ax[0].text(ts_ocean+10, 110, 'ocean fully\nevaporated', color='purple', fontsize=7.5)

fig.tight_layout(); fig.savefig('tools/diag_kopparapu_regime.png', dpi=120)
print('saved tools/diag_kopparapu_regime.png')

print(f"\nps reaches Earth-ocean inventory ({OCEAN_BAR} bar) at Ts ~ {ts_ocean:.0f} K")
print(f"k-table T-ceiling (500 K) exceeded at the surface for Ts > {ts_ceiling:.0f} K")
print(f"\nOLR peak = {OLR.max():.1f} W/m² at Ts={TS[np.argmax(OLR)]:.0f} K  "
      f"(S-N limit {SN}; overshoot {OLR.max()-SN:+.1f})")
print(f"OLR at Ts=400 K = {OLR[np.argmin(np.abs(TS-400))]:.1f};  "
      f"at 1000 K = {OLR[np.argmin(np.abs(TS-1000))]:.1f};  "
      f"at 2200 K = {OLR[-1]:.1f}  (should be ~flat at {SN})")
print(f"albedo: min={ALB.min():.3f} at Ts={TS[np.argmin(ALB)]:.0f} K; "
      f"max={ALB.max():.3f} at Ts={TS[np.argmax(ALB)]:.0f} K")
print(f"\nValid regime (below ceiling AND ocean-reasonable): Ts ≲ "
      f"{min(ts_ceiling, ts_ocean):.0f} K  — this brackets the actual inner-edge "
      f"definitions (moist GH ~340 K, runaway OLR limit).")
