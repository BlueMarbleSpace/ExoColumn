#!/usr/bin/env python3
"""Generate ExoColumn-owned local copies of three ExoRT n68equiv source files.

ExoRT is read-only (CLAUDE.md guiding principle #1).  Rather than editing
/models/ExoRT, the build copies a few of its source files into ExoColumn's
src/ directory, applying small ExoColumn-owned patches, and compiles the
copies.  ExoRT's source tree is never touched; the generated files are
gitignored and regenerated on every build, so they auto-track ExoRT updates.

This is the same local-copy pattern the Makefile already uses for
exoplanet_mod.F90 (PVER substitution).  The single-line seds it used for
calc_opd_mod.F90 (k-table clamp) have been folded in here because the BPS
continuum toggle requires multi-line block insertions that are awkward and
fragile to express as inline Makefile seds.

Files generated (from <spec_dir> = ExoRT/source/src.n68equiv):

  calc_opd_mod.F90          k-table pressure CLAMP (unchanged behaviour) +
                            runtime MT_CKD <-> BPS H2O-continuum branch
  radgrid.F90               public `use_bps_continuum` runtime flag, the BPS
                            continuum filename, and the `solar_file_override`
                            runtime stellar-spectrum selector
  initialize_rad_mod_1D.F90 conditional read of bps_h20_continuum_n68.nc, plus
                            honouring solar_file_override in initialize_solar

Each patch is keyed to an explicit anchor string and asserts the anchor was
found; any miss prints `ERROR: ...` and exits non-zero so the Makefile fails
loudly if a future ExoRT update moves the anchor (mirrors the Makefile's
existing `$(error ...)` guard).

Usage:  patch_exort.py <spec_dir> <out_dir>
"""

import sys
import os


class PatchError(Exception):
    pass


def _require(condition, message):
    if not condition:
        raise PatchError(message)


def replace_once(text, old, new, where):
    """Replace exactly one occurrence of `old`; error if 0 or >1 found."""
    n = text.count(old)
    _require(n == 1, f"{where}: expected exactly 1 occurrence of anchor, found {n}: {old!r}")
    return text.replace(old, new)


def insert_after_line(text, anchor, block, where, start=0):
    """Insert `block` immediately after the first line containing `anchor`
    (searching from character offset `start`).  Returns (new_text, end_offset)."""
    idx = text.find(anchor, start)
    _require(idx != -1, f"{where}: anchor not found: {anchor!r}")
    line_end = text.find("\n", idx)
    _require(line_end != -1, f"{where}: anchor line has no newline: {anchor!r}")
    insert_at = line_end + 1
    return text[:insert_at] + block + text[insert_at:], insert_at + len(block)


def insert_before_line(text, anchor, block, where):
    """Insert `block` immediately before the start of the line containing `anchor`."""
    idx = text.find(anchor)
    _require(idx != -1, f"{where}: anchor not found: {anchor!r}")
    line_start = text.rfind("\n", 0, idx) + 1  # 0 if not found -> start of file
    return text[:line_start] + block + text[line_start:]


# --------------------------------------------------------------------------
# calc_opd_mod.F90
# --------------------------------------------------------------------------

# (1) k-table pressure clamp — identical effect to the historical Makefile sed.
CLAMP_OLD = "pressure = log10(pmid(ik))       ! log pressure"
CLAMP_NEW = ("pressure = min(log10(pmid(ik)), log10pgrid(kc_npress))  "
             "! log pressure; CLAMPED at k-table top by ExoColumn build")

# (2) scratch scalars for the BPS continuum, added after the amagats block.
CALC_DECL_ANCHOR = "    real(r8) :: amaN2, amaH2, amaCO2, amaCH4, amaH2O, amaFRGN, amaO2"
CALC_DECL_BLOCK = (
    "    ! ExoColumn: BPS H2O-continuum scratch (used when use_bps_continuum=.true.)\n"
    "    real(r8) :: arg1, arg2, radfield, tau_h2os_bps, tau_h2of_bps\n"
)

# (3a) open the runtime branch just before the MT_CKD continuum section.
CALC_MTCKD_ANCHOR = "      !===== MT_CKD =====!"
CALC_IF_OPEN = (
    "      ! ExoColumn: select the H2O self+foreign continuum model at runtime\n"
    "      ! (&exocol_nml::h2o_continuum).  Default .false. => MT_CKD, bit-identical\n"
    "      ! to stock n68equiv.  .true. => BPS layer-average continuum.\n"
    "      if (.not. use_bps_continuum) then\n"
)

# (3b) the MT_CKD branch's last statement (unique) — used to locate the enddo
#      that closes its band loop, after which we splice in the BPS else-branch.
CALC_MTCKD_LASTLINE = ("tau_gas(itc,ik) = tau_gas(itc,ik) + "
                       "(ans_h2os(ig,iw)*amaH2O + ans_h2of(ig,iw)*amaFRGN) * u_h2o")
CALC_BANDLOOP_CLOSE = "enddo    ! close band loop"

CALC_ELSE_BLOCK = """      else
        !===== BPS H2O continuum (ExoColumn; layer-average self+foreign) =====!
        ! Ported verbatim from ExoRT src.n68h2o/calc_opd_mod.F90 (Wolf, BPS
        ! continuum), with h2ovap_press -> ppH2O (the n68equiv H2O partial
        ! pressure, numerically identical).  Applied as a per-band layer
        ! average, added to every g-point in the band.
        itc = 0
        do iw = iwbeg, iwend     ! loop over bands
          arg1 = 1.4388_r8*wavenum_mid(iw)/(2._r8*tmid(ik))
          arg2 = 1.4388_r8*wavenum_mid(iw)/(2._r8*296.0_r8)
          radfield = (exp(2._r8*arg1)-1._r8)/(exp(2._r8*arg1)+1._r8) &
                   / (exp(2._r8*arg2)-1._r8)/(exp(2._r8*arg2)+1._r8)
          tau_h2os_bps = (296.0_r8/tmid(ik)) &
                       * (ppH2O/1013.25_r8) &
                       * (pathlength(ik)*ppH2O)/(SHR_CONST_BOLTZ*tmid(ik)) &
                       * (self(iw)*exp(TempCoeff(iw)*(296.0_r8-tmid(ik))) + radfield*base_self(iw)) &
                       / (100._r8**2)    ! cm^2 -> m^2
          tau_h2of_bps = (296.0_r8/tmid(ik)) &
                       * ((pmid(ik)-ppH2O)/1013.25_r8) &
                       * (pathlength(ik)*ppH2O)/(SHR_CONST_BOLTZ*tmid(ik)) &
                       * radfield * (foreign(iw) + base_foreign(iw)) &
                       / (100._r8**2)    ! cm^2 -> m^2
          do ig = 1, ngauss_pts(iw)
            itc = itc + 1
            tau_gas(itc,ik) = tau_gas(itc,ik) + tau_h2os_bps + tau_h2of_bps
          enddo
        enddo    ! close band loop
      endif
"""


def patch_calc_opd(text):
    text = replace_once(text, CLAMP_OLD, CLAMP_NEW, "calc_opd:clamp")
    text = insert_after_line(text, CALC_DECL_ANCHOR, CALC_DECL_BLOCK, "calc_opd:decls")[0]
    text = insert_before_line(text, CALC_MTCKD_ANCHOR, CALC_IF_OPEN, "calc_opd:if-open")

    # Locate the MT_CKD branch's closing `enddo ! close band loop`: it is the
    # first such close *after* the unique MT_CKD final tau_gas statement.
    last_idx = text.find(CALC_MTCKD_LASTLINE)
    _require(last_idx != -1, f"calc_opd:else: MT_CKD last-line anchor not found: {CALC_MTCKD_LASTLINE!r}")
    text, _ = insert_after_line(text, CALC_BANDLOOP_CLOSE, CALC_ELSE_BLOCK,
                                "calc_opd:else", start=last_idx)
    return text


# --------------------------------------------------------------------------
# radgrid.F90
# --------------------------------------------------------------------------

RADGRID_ANCHOR = "  !!!=== end bps definitions ===="
RADGRID_BLOCK = (
    "\n"
    "  ! ---- ExoColumn: runtime MT_CKD <-> BPS H2O-continuum toggle -------------\n"
    "  ! Set by exocol_driver from &exocol_nml::h2o_continuum BEFORE\n"
    "  ! initialize_kcoeff.  .false. (default) => MT_CKD 3.3 8-gpt path,\n"
    "  ! bit-identical to stock n68equiv.  .true. => BPS layer-average continuum\n"
    "  ! (bps_h20_continuum_n68.nc), for direct comparison with Kopparapu (2013).\n"
    "  logical :: use_bps_continuum = .false.\n"
    "  character(len=256), parameter :: kh2oself_bps_file = 'bps_h20_continuum_n68.nc'\n"
    "\n"
    "  ! ---- ExoColumn: runtime stellar-spectrum override ----------------------\n"
    "  ! Set by exocol_driver from &exocol_nml::solar_file BEFORE initialize_solar.\n"
    "  ! Empty (default) => use the compile-time exoplanet_mod::solar_file\n"
    "  ! (G2V_SUN_n68.nc), bit-identical to stock behaviour.  Non-empty => load\n"
    "  ! that n68-band file from data/solar/ instead (HZ multi-stellar sweeps).\n"
    "  character(len=256) :: solar_file_override = ''\n"
)


def patch_radgrid(text):
    return insert_after_line(text, RADGRID_ANCHOR, RADGRID_BLOCK, "radgrid:toggle")[0]


# --------------------------------------------------------------------------
# initialize_rad_mod_1D.F90
# --------------------------------------------------------------------------

INIT_ANCHOR = "      !! /mtckd"
INIT_BLOCK = """
      !! ExoColumn: optional BPS H2O continuum (toggled by &exocol_nml::h2o_continuum).
      !! Read only when active; uses dirct (kabs) + kh2oself_bps_file (radgrid).
      if (use_bps_continuum) then
        filename = trim(exort_rootdir)//trim(dirct)//trim(kh2oself_bps_file)
        call getfil(filename, locfn, 0)
        call wrap_open(locfn, 0, ncid)
        call wrap_inq_varid(ncid, 'self', keff_id)
        call wrap_get_var_realx(ncid, keff_id, self)
        call wrap_inq_varid(ncid, 'foreign', keff_id)
        call wrap_get_var_realx(ncid, keff_id, foreign)
        call wrap_inq_varid(ncid, 'base_self', keff_id)
        call wrap_get_var_realx(ncid, keff_id, base_self)
        call wrap_inq_varid(ncid, 'base_foreign', keff_id)
        call wrap_get_var_realx(ncid, keff_id, base_foreign)
        call wrap_inq_varid(ncid, 'TempCoeff', keff_id)
        call wrap_get_var_realx(ncid, keff_id, TempCoeff)
        write (6, '(2x, a)') 'ExoColumn: BPS H2O continuum active (bps_h20_continuum_n68.nc)'
      endif
"""


# initialize_solar builds the solar filename from exoplanet_mod::solar_file.
# Swap that single line so a non-empty radgrid::solar_file_override (set at
# runtime from &exocol_nml::solar_file) selects a different host-star spectrum.
SOLAR_FNAME_OLD = (
    "      filename = trim(exort_rootdir)//trim(dirsol)//trim(solar_file)"
)
SOLAR_FNAME_NEW = (
    "      ! ExoColumn: honour runtime &exocol_nml::solar_file override (radgrid).\n"
    "      ! Empty override => compile-time exoplanet_mod::solar_file (bit-identical).\n"
    "      if (len_trim(solar_file_override) > 0) then\n"
    "        filename = trim(exort_rootdir)//trim(dirsol)//trim(solar_file_override)\n"
    "      else\n"
    "        filename = trim(exort_rootdir)//trim(dirsol)//trim(solar_file)\n"
    "      end if"
)


def patch_initialize(text):
    text = replace_once(text, SOLAR_FNAME_OLD, SOLAR_FNAME_NEW,
                        "initialize:solar-override")
    return insert_after_line(text, INIT_ANCHOR, INIT_BLOCK, "initialize:bps-read")[0]


# --------------------------------------------------------------------------

JOBS = [
    ("calc_opd_mod.F90", patch_calc_opd),
    ("radgrid.F90", patch_radgrid),
    ("initialize_rad_mod_1D.F90", patch_initialize),
]


def main(argv):
    if len(argv) != 3:
        print("ERROR: usage: patch_exort.py <spec_dir> <out_dir>", file=sys.stderr)
        return 2
    spec_dir, out_dir = argv[1], argv[2]
    try:
        for fname, patch_fn in JOBS:
            src = os.path.join(spec_dir, fname)
            dst = os.path.join(out_dir, fname)
            _require(os.path.isfile(src), f"{fname}: ExoRT source not found: {src}")
            with open(src, "r") as fh:
                text = fh.read()
            patched = patch_fn(text)
            with open(dst, "w") as fh:
                fh.write(patched)
    except PatchError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    # Positive success sentinel: the Makefile treats its ABSENCE as failure, so
    # a missing python3 (empty log) or an exception are both caught.
    print("PATCH_OK: generated calc_opd_mod.F90, radgrid.F90, initialize_rad_mod_1D.F90 (clamp + BPS toggle)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
