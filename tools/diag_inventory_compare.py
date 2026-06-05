#!/usr/bin/env python3
"""Before/after: old unbounded variable_ps + ws>1 switch vs new finite-inventory
Kasting switch.  /tmp/hz_data.txt = OLD, /tmp/hz_data2.txt = NEW (Ts OLR ASR alb Seff)."""
import numpy as np
import matplotlib.pyplot as plt

old = np.loadtxt('/tmp/hz_data.txt'); new = np.loadtxt('/tmp/hz_data2.txt')
to, oo, ao, lo, so = old.T
tn, on, an, ln, sn = new.T
SN = 282.0; KOPP = 1.0140

fig, ax = plt.subplots(2, 2, figsize=(10, 7.5), dpi=120)
(a, b), (c, d) = ax

a.plot(to, oo, 'C7', lw=1, alpha=0.7, label='OLD OLR (unbounded)')
a.plot(tn, on, 'C3', lw=1.3, label='NEW OLR (inventory)')
a.axhline(SN, color='k', ls='--', lw=0.8, label='S-N limit 282')
a.set_ylabel('OLR (W/m²)'); a.legend(fontsize=7); a.set_title('OLR: now plateaus at the S-N limit')

b.plot(to, ao, 'C7', lw=1, alpha=0.7, label='OLD ASR')
b.plot(tn, an, 'C0', lw=1.3, label='NEW ASR')
b.set_ylabel('ASR (W/m²)'); b.legend(fontsize=7); b.set_title('ASR')

c.plot(to, lo, 'C7', lw=1, alpha=0.7, label='OLD α')
c.plot(tn, ln, 'C2', lw=1.3, label='NEW α')
c.set_ylabel('albedo'); c.set_xlabel('Ts (K)'); c.legend(fontsize=7)
c.set_title('Albedo: no more runaway-Rayleigh climb to 0.43')

d.plot(to, so, 'C7', lw=1, alpha=0.7, label='OLD Seff')
d.plot(tn, sn, 'C1', lw=1.3, label='NEW Seff')
d.axhline(KOPP, color='C3', ls='--', lw=0.8, label='Kopparapu runaway 1.014')
d.set_ylabel('Seff'); d.set_xlabel('Ts (K)'); d.legend(fontsize=7)
d.set_title('Seff: 1.3 → ~1.05-1.14 (near Kopparapu)')
d.set_ylim(0.9, 1.4); d.set_xlim(300, 2200)

fig.tight_layout(); fig.savefig('tools/diag_inventory_compare.png', dpi=120)
print('saved tools/diag_inventory_compare.png')

def at(t, ts, y): return y[np.argmin(np.abs(ts - t))]
print("\n            OLD ->  NEW")
for T in (340, 400, 600, 1000, 2200):
    print(f"Ts={T:4d}  OLR {at(T,to,oo):6.1f}->{at(T,tn,on):6.1f}   "
          f"alb {at(T,to,lo):.3f}->{at(T,tn,ln):.3f}   "
          f"Seff {at(T,to,so):.3f}->{at(T,tn,sn):.3f}")
# staircase roughness (std vs cubic) over 700-2200
def rough(ts, y):
    m = ts >= 700; p = np.polyfit(ts[m], y[m], 4); return np.std(y[m]-np.polyval(p, ts[m]))
print(f"\nOLR staircase roughness (700-2200): OLD {rough(to,oo):.2f} -> NEW {rough(tn,on):.2f} W/m²")
print(f"Seff peak: OLD {so.max():.3f}  NEW {sn.max():.3f}  (Kopparapu 1.014)")
print(f"OLR plateau (1000-2200) mean: NEW {on[(tn>=1000)].mean():.1f} W/m² (S-N 282)")
