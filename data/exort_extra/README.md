# Supplementary ExoRT stellar spectra

The three BT-Settl spectral energy distributions in this directory are used by
the multi-stellar habitable-zone calculations (`reference/habitablezone/`,
Figures 6 and 7 of the paper) but are **not yet part of the public
[ExoRT](https://github.com/storyofthewolf/ExoRT) distribution**:

| File | Host star |
|------|-----------|
| `btsettl_T2000_g4.5_m0.0_n68.nc` | M dwarf, $T_{\rm eff}$ = 2000 K |
| `btsettl_T2200_g4.5_m0.0_n68.nc` | M dwarf, $T_{\rm eff}$ = 2200 K |
| `btsettl_T2400_g4.5_m0.0_n68.nc` | M dwarf, $T_{\rm eff}$ = 2400 K |

They extend ExoRT's host-star ladder below its previous 2600 K floor. Each is a
BT-Settl model atmosphere (log g = 4.5, [Fe/H] = 0) binned onto the 68-band
n68equiv grid with ExoRT's `makeStellarSpectrum_fromSED.py`, in the same format
as ExoRT's own `data/solar/*.nc` files (`wav_low`, `wav_high`, `S0`,
`solarflux`).

They are shipped here so that this release is self-contained and the coolest
three M-dwarf points in Figures 6 and 7 can be reproduced. **Copy them into your
ExoRT tree before running the multi-stellar sweeps:**

```bash
cp data/exort_extra/*.nc "$EXORT_ROOT"/data/solar/     # default EXORT_ROOT=/models/ExoRT
```

This is the only step in ExoColumn that writes into ExoRT. It adds new files and
modifies none, so the read-only policy for ExoRT source is preserved. Once these
SEDs are merged upstream into ExoRT, this directory can be removed.

All other host-star spectra used in the paper (2600–7200 K, plus the n84 F-star
SED) are already distributed with ExoRT and require no action.
