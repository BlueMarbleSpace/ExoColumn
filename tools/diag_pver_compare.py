#!/usr/bin/env python3
"""Compare the OLR staircase at PVER=70 vs PVER=140 over 1000-1600 K.
If the staircase is a vertical-resolution artifact (tropopause snapping),
PVER=140 should show ~2x as many steps with ~half the amplitude."""
import numpy as np
import matplotlib.pyplot as plt

# PVER=70: full sweep (Ts OLR ASR alb Seff), slice 1000-1600
d70 = np.loadtxt('/tmp/hz_data.txt')
m = (d70[:, 0] >= 1000) & (d70[:, 0] <= 1600)
ts70, olr70, asr70 = d70[m, 0], d70[m, 1], d70[m, 2]

# PVER=140 subset (Ts OLR ASR albedo)
d140 = np.loadtxt('/tmp/pver140.csv')
ts140, olr140, asr140 = d140[:, 0], d140[:, 1], d140[:, 2]

def steps(ts, y, thr=0.6):
    dj = np.abs(np.diff(y))
    return ts[1:][dj > thr], dj[dj > thr]

j70, a70 = steps(ts70, olr70)
j140, a140 = steps(ts140, olr140)

fig, ax = plt.subplots(2, 1, figsize=(9, 7), dpi=120, sharex=True)
ax[0].plot(ts70, olr70, 'C3-o', ms=2.5, lw=1, label=f'PVER=70  ({len(j70)} OLR steps)')
ax[0].plot(ts140, olr140, 'C0-o', ms=2.5, lw=1, label=f'PVER=140 ({len(j140)} OLR steps)')
ax[0].set_ylabel('OLR (W/m²)'); ax[0].legend()
ax[0].set_title('OLR staircase: PVER=70 vs PVER=140 (1000-1600 K)')

ax[1].plot(ts70, asr70, 'C3-o', ms=2.5, lw=1, label='PVER=70')
ax[1].plot(ts140, asr140, 'C0-o', ms=2.5, lw=1, label='PVER=140')
ax[1].set_ylabel('ASR (W/m²)'); ax[1].set_xlabel('Ts (K)'); ax[1].legend()
fig.tight_layout(); fig.savefig('tools/diag_pver_compare.png', dpi=120)
print('saved tools/diag_pver_compare.png')

print(f"\nPVER=70 : {len(j70):2d} OLR steps, mean |jump| = {a70.mean():.2f} W/m², "
      f"peak-peak OLR = {olr70.max()-olr70.min():.2f}")
print(f"PVER=140: {len(j140):2d} OLR steps, mean |jump| = {a140.mean():.2f} W/m², "
      f"peak-peak OLR = {olr140.max()-olr140.min():.2f}")
# std around a smooth (savgol-like) trend = staircase roughness
def roughness(ts, y):
    p = np.polyfit(ts, y, 3); return np.std(y - np.polyval(p, ts))
print(f"\nstaircase roughness (std vs cubic trend):")
print(f"  PVER=70  = {roughness(ts70, olr70):.3f} W/m²")
print(f"  PVER=140 = {roughness(ts140, olr140):.3f} W/m²")
