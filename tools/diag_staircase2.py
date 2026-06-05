#!/usr/bin/env python3
"""Characterize the hz_inner staircase from the actual sweep output and
correlate the OLR/ASR jumps with candidate discrete quantities:
  (1) k_top_conv / cold-trap (from diag_staircase.build_column)
  (2) # of model layers whose pmid exceeds the k-table max (10 bar)
  (3) # of model layers with pmid above each k-table pressure node region
"""
import numpy as np
import matplotlib.pyplot as plt
from diag_staircase import build_column, esat_cc, P_TOP, PVER, PVERP

d = np.loadtxt('/tmp/hz_data.txt')
TS, OLR, ASR, ALB, SEFF = d.T

# candidate discrete quantities recomputed on the same grid as the model
n_over_10bar = []     # layers with pmid > 1e4 mb = 1e6 Pa  (k-table max)
k_top = []
for ts in TS:
    e = esat_cc(ts); ps = 1.0e5 + e
    pint = np.exp(np.log(P_TOP) + np.arange(PVERP)/(PVERP-1)*np.log(ps/P_TOP))
    pmid = 0.5*(pint[:-1]+pint[1:])
    n_over_10bar.append(int(np.sum(pmid > 1.0e6)))
    k_top.append(build_column(ts)['k_top'])
n_over_10bar = np.array(n_over_10bar)
k_top = np.array(k_top)

# locate OLR jumps (discrete drops) in the staircase regime
seg = TS >= 700
dOLR = np.diff(OLR)
dASR = np.diff(ASR)
# a "tooth" = a step where |dOLR| spikes well above the local ramp
thr = 0.8
olr_jumps = TS[1:][(np.abs(dOLR) > thr) & (TS[1:] >= 700)]

fig, ax = plt.subplots(3, 1, figsize=(9, 9), dpi=120, sharex=True)
ax[0].plot(TS, OLR, 'C3', label='OLR')
ax[0].plot(TS, ASR, 'C0', label='ASR')
for x in olr_jumps:
    ax[0].axvline(x, color='gray', lw=0.4, alpha=0.5)
ax[0].set_ylabel('flux (W/m²)'); ax[0].legend(); ax[0].set_xlim(700, 2200)
ax[0].set_title('OLR/ASR with detected jumps (vertical lines)')

ax[1].step(TS, n_over_10bar, 'C1', where='mid')
for x in olr_jumps:
    ax[1].axvline(x, color='gray', lw=0.4, alpha=0.5)
ax[1].set_ylabel('# layers pmid>10 bar')
ax[1].set_title('Layers above k-table max (10 bar) vs Ts — does each +1 line up with a jump?')

ax[2].step(TS, k_top, 'C2', where='mid')
for x in olr_jumps:
    ax[2].axvline(x, color='gray', lw=0.4, alpha=0.5)
ax[2].set_ylabel('k_top_conv'); ax[2].set_xlabel('Ts (K)')
ax[2].set_title('Cold-trap tropopause layer index vs Ts')

fig.tight_layout()
fig.savefig('tools/diag_staircase2.png', dpi=120)
print('saved tools/diag_staircase2.png')

# quantify correlation: do OLR jumps coincide with n_over_10bar steps?
n10_steps = TS[1:][(np.diff(n_over_10bar) != 0) & (TS[1:] >= 700)]
ktop_steps = TS[1:][(np.diff(k_top) != 0) & (TS[1:] >= 700)]
print(f'\n# OLR jumps (|dOLR|>{thr}) in 700-2200 K : {len(olr_jumps)}')
print(f'# n_over_10bar steps in 700-2200 K       : {len(n10_steps)}')
print(f'# k_top_conv steps in 700-2200 K         : {len(ktop_steps)}')

def coincidence(a, b, tol=7.5):
    if len(a) == 0: return np.nan
    return np.mean([np.min(np.abs(b - x)) <= tol for x in a]) if len(b) else 0.0
print(f'\nfraction of OLR jumps within 7.5 K of a n_over_10bar step: '
      f'{coincidence(olr_jumps, n10_steps):.2f}')
print(f'fraction of OLR jumps within 7.5 K of a k_top_conv step : '
      f'{coincidence(olr_jumps, ktop_steps):.2f}')
print(f'\nn_over_10bar range over 700-2200 K: '
      f'{n_over_10bar[seg].min()}..{n_over_10bar[seg].max()}')
