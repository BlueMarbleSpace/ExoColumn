#!/usr/bin/env python3
"""
digitize_kopparapu_fig5.py — extract the three Kopparapu et al. (2013) Figure 5
curves (outer HZ: F_IR & F_SOL, planetary albedo, Seff vs CO2 partial pressure)
from the paper PDF by pixel digitization, for direct overlay in the ExoColumn
outer-HZ sweep figure (reference/max_greenhouse/hz_outer.py).

Method: render page 8 at 300 dpi (pdftoppm), locate each MATLAB axes frame
(full-height/-width black lines), calibrate using the frame limits — which in
these MATLAB panels coincide with the labelled end ticks:
  (a) x: 10^0..10^2 bar (log), y: 50..250 W/m2;  F_IR red solid, F_SOL blue dashed
  (b) x: 10^0..10^2 bar (log), y: 0.2..0.7;      albedo blue solid
  (c) x: 10^0..10^2 bar (log), y: 0.32..0.46;    Seff blue solid
then collect pure-blue / pure-red curve pixels column by column.

Validation (run this script to see it): the digitized Seff minimum is ~0.338 at
~7-8 bar, matching the caption (1.70 AU -> Seff 0.346) / text (0.325 at ~8 bar)
to figure-reading accuracy, and F_IR at 1 bar is ~112 W/m2.

Output: reference/max_greenhouse/kopparapu2013_fig5.npz with arrays
  fir_p, fir / fsol_p, fsol / alb_p, alb / seff_p, seff
(_p = pCO2 in bar).

Requires: pdftoppm, PIL, the paper at ~/Documents/library/kopparapu_etal2013.pdf
"""

import os
import subprocess
import numpy as np
from PIL import Image

PDF  = os.path.expanduser('~/Documents/library/kopparapu_etal2013.pdf')
HERE = os.path.dirname(os.path.abspath(__file__))
OUT  = os.path.abspath(os.path.join(HERE, '..', 'reference', 'max_greenhouse',
                                    'kopparapu2013_fig5.npz'))
PAGE_PNG = '/tmp/kopp_fig5_page8.png'


def render_page():
    subprocess.run(['pdftoppm', '-png', '-r', '300', '-f', '8', '-l', '8',
                    PDF, '/tmp/kopp_fig5_page'], check=True)
    os.replace('/tmp/kopp_fig5_page-08.png', PAGE_PNG)


def find_frame(a, xr, yr):
    """Locate the axes box (full-height cols, full-width rows of black) inside
    the page-pixel search window xr=(x0,x1), yr=(y0,y1)."""
    blk = (a[:, :, 0] < 60) & (a[:, :, 1] < 60) & (a[:, :, 2] < 60)
    sub = blk[yr[0]:yr[1], xr[0]:xr[1]]
    rowcnt = sub.sum(axis=1)
    rows = np.where(rowcnt > 0.85 * rowcnt.max())[0]
    if len(rows) < 2:
        raise RuntimeError('frame rows not found')
    rtop, rbot = rows.min(), rows.max()
    band = sub[rtop:rbot + 1, :]
    cols = np.where(band.sum(axis=0) > 0.9 * (rbot - rtop))[0]
    if len(cols) < 2:
        raise RuntimeError('frame cols not found')
    cleft, cright = cols.min(), cols.max()
    return (xr[0] + cleft, xr[0] + cright, yr[0] + rtop, yr[0] + rbot)


def digitize(a, frame, color, ylim, exclude=None):
    """Column-by-column mean row of `color` pixels inside the frame.
    color: 'blue' or 'red'.  x is log-calibrated 10^0..10^2 bar.
    exclude: optional (x_bar_min, x_bar_max, y_val_min, y_val_max) region to
    drop (e.g. in-panel label text in the curve colour)."""
    xl, xr_, yt, yb = frame
    if color == 'blue':
        mask = (a[:, :, 2] > 120) & (a[:, :, 0] < 110) & (a[:, :, 1] < 110)
    else:
        mask = (a[:, :, 0] > 150) & (a[:, :, 2] < 110) & (a[:, :, 1] < 110)
    xs, vs = [], []
    for x in range(xl + 2, xr_ - 1):
        rr = np.where(mask[yt + 2:yb - 1, x])[0]
        if len(rr) == 0:
            continue
        r = rr.mean() + yt + 2
        pbar = 10.0 ** (2.0 * (x - xl) / (xr_ - xl))
        val  = ylim[1] + (ylim[0] - ylim[1]) * (r - yt) / (yb - yt)
        if exclude and exclude[0] <= pbar <= exclude[1] \
                   and exclude[2] <= val <= exclude[3]:
            continue
        xs.append(pbar)
        vs.append(val)
    return np.array(xs), np.array(vs)


def main():
    render_page()
    a = np.asarray(Image.open(PAGE_PNG).convert('RGB')).astype(int)
    h, w, _ = a.shape

    fr_a = find_frame(a, (int(0.05 * w), int(0.50 * w)), (int(0.04 * h), int(0.29 * h)))
    fr_b = find_frame(a, (int(0.50 * w), int(0.99 * w)), (int(0.04 * h), int(0.29 * h)))
    fr_c = find_frame(a, (int(0.20 * w), int(0.80 * w)), (int(0.27 * h), int(0.53 * h)))

    # (a) F_IR (red solid) + F_SOL (blue dashed); drop the black-text labels'
    # surroundings automatically (labels/arrows are black, not curve-coloured).
    fir_p,  fir  = digitize(a, fr_a, 'red',  (50.0, 250.0))
    fsol_p, fsol = digitize(a, fr_a, 'blue', (50.0, 250.0))
    alb_p,  alb  = digitize(a, fr_b, 'blue', (0.2, 0.7))
    seff_p, seff = digitize(a, fr_c, 'blue', (0.32, 0.46))

    np.savez(OUT, fir_p=fir_p, fir=fir, fsol_p=fsol_p, fsol=fsol,
             alb_p=alb_p, alb=alb, seff_p=seff_p, seff=seff)
    print(f'Saved: {OUT}')

    imin = np.argmin(seff)
    print(f'  check: F_IR(1 bar) = {fir[0]:.1f} W/m2 (expect ~112)')
    print(f'  check: albedo(1 bar) = {alb[0]:.3f} (expect ~0.277)')
    print(f'  check: Seff min = {seff[imin]:.4f} at pCO2 = {seff_p[imin]:.2f} bar '
          f'(expect ~0.338 at ~7-8 bar)')


if __name__ == '__main__':
    main()
