#------------------------------------------------------------------------------
# ExoColumn Makefile
# Links ExoColumn src/ against ExoRT source/src.main + source/src.n68equiv
# and NetCDF-Fortran.
#
# Local settings (compiler path, NETCDF_ROOT, EXORT_ROOT) should be placed in
# config.mk at the project root — copy config.mk.example to get started.
# config.mk is included automatically if it exists; it is not committed to git.
#
# Usage:
#   make          build (reads ../config.mk if present)
#   make clean
#
# The executable is placed in ../run/exocol.exe.
# Before running, set working directory to the project root so that ExoRT can
# find its data files via the relative paths in sys_rootdir.F90 / exoplanet_mod.
#------------------------------------------------------------------------------

# Pull in local machine settings (silently skipped if file does not exist).
-include ../config.mk

# --- NetCDF -------------------------------------------------------------------
# If NETCDF_ROOT is set (in config.mk or on the command line) use it directly.
# Otherwise fall back to nf-config (works for gfortran builds where the system
# NetCDF matches the compiler).
ifdef NETCDF_ROOT
  LIB_NETCDF := -L$(NETCDF_ROOT)/lib -Wl,-Bstatic -lnetcdff -Wl,-Bdynamic -lnetcdf
  INC_NETCDF := -I$(NETCDF_ROOT)/include
else
  LIB_NETCDF := $(shell nf-config --flibs)
  INC_NETCDF := $(shell nf-config --fflags)
endif

# --- Compiler -----------------------------------------------------------------
# Default: Intel ifx at its well-known OneAPI path.
# Override via USER_FC in config.mk or on the command line.
IFX := /opt/intel/oneapi/compiler/2025.3/bin/ifx

ifeq ($(USER_FC),)
  FC := $(IFX)
else
  FC := $(USER_FC)
endif

# Compiler flags + the .mod output-directory flag (vendor-specific: Intel uses
# `-module <dir>`, gfortran uses `-J<dir>`; both also add <dir> to the module
# search path).  BUILD_DIR is expanded lazily, so its later definition is fine.
ifeq ($(FC),$(IFX))
  FFLAGS  = -O2 -traceback -fpe0
  MODFLAG = -module $(BUILD_DIR)
else ifeq ($(FC),gfortran)
  # -fallow-argument-mismatch / -fallow-invalid-boz: tolerate the legacy
  # CESM-lineage interface looseness in the ExoRT source (scalar/array netcdf
  # wrappers, BOZ bit patterns in infnan.F90) that modern gfortran rejects but
  # ifx accepts.  ExoRT is read-only, so these are handled via flags.
  FFLAGS  = -O2 -fbacktrace -ffpe-trap=invalid,zero,overflow -cpp \
            -fallow-argument-mismatch -fallow-invalid-boz -fno-range-check
  MODFLAG = -J$(BUILD_DIR)
else
  # Fallback: assume ifort-compatible flags
  FFLAGS  = -O2 -traceback -fpe0
  MODFLAG = -module $(BUILD_DIR)
endif

# --- ExoRT source directories -------------------------------------------------
# Override EXORT_ROOT in config.mk if ExoRT lives elsewhere on your machine.
EXORT_ROOT ?= /models/ExoRT
MAIN_DIR    = $(EXORT_ROOT)/source/src.main
MISC_DIR    = $(EXORT_ROOT)/source/src.misc
SPEC_DIR    = $(EXORT_ROOT)/source/src.n68equiv
SRC_DIR     = $(EXORT_ROOT)/source

# ExoColumn source directory — absolute path so the compiler can be invoked
# from a different working directory (see build rule below).
EXCOL_DIR   = $(abspath ../src)
BUILD_DIR   = $(abspath .)

# --- ExoColumn-owned model settings ------------------------------------------
# ExoRT exposes a few compile-time `parameter` constants (vertical levels,
# gravity, stellar spectrum, etc.).  Rather than edit /models/ExoRT, we let
# ExoColumn own these values: at build time we read ExoRT's exoplanet_mod.F90,
# substitute the values we control, and write a local copy into src/.  The
# build then compiles our local copy instead of ExoRT's, leaving ExoRT
# completely untouched.  Override in config.mk or on the command line, e.g.:
#
#   make PVER=120
#
PVER ?= 70

# Generate sys_rootdir.F90 with the local EXORT_ROOT path, replacing ExoRT's
# hardcoded machine-specific path.  Runs on every make invocation (fast).
# The generated file is machine-specific and is not committed to version control.
$(shell printf 'module sys_rootdir\n  implicit none\n  public\n  character(len=256), parameter :: exort_rootdir = "%s/"\nend module sys_rootdir\n' \
  "$(EXORT_ROOT)" > $(EXCOL_DIR)/sys_rootdir.F90)

# Generate src/exoplanet_mod.F90 from ExoRT's source with our chosen exo_pver
# value.  Matches the (single) uncommented `exo_pver = N` line and rewrites it;
# commented variants are left alone.  ExoRT's file is preserved.
$(shell sed -E 's|^([[:space:]]*integer[[:space:]]*,[[:space:]]*parameter[[:space:]]*::[[:space:]]*exo_pver[[:space:]]*=[[:space:]]*)[0-9]+.*|\1$(PVER)  ! set by ExoColumn build (PVER=$(PVER))|' \
  $(SRC_DIR)/exoplanet_mod.F90 > $(EXCOL_DIR)/exoplanet_mod.F90)

# Generate local copies of three n68equiv source files via patch_exort.py:
#   calc_opd_mod.F90          k-coefficient pressure interpolation CLAMPED at the
#                             table top (10 bar; above its grid ExoRT linearly
#                             EXTRAPOLATES k in log10(p), unphysical for >10 bar
#                             columns and can produce negative k — the clamp
#                             freezes pressure broadening at its 10-bar value for
#                             deeper layers, while partial pressures/amagats keep
#                             the true pmid; <10 bar columns are BIT-IDENTICAL)
#                             PLUS the runtime MT_CKD<->BPS H2O-continuum branch.
#   radgrid.F90               public `use_bps_continuum` flag + BPS filename.
#   initialize_rad_mod_1D.F90 conditional read of bps_h20_continuum_n68.nc.
# Same local-copy pattern as exoplanet_mod above; ExoRT's source is preserved.
# The patches are anchored to explicit lines and the script exits non-zero
# (printing ERROR) if any anchor disappears in a future ExoRT update.
PATCH_EXORT_LOG := $(shell python3 $(BUILD_DIR)/patch_exort.py $(SPEC_DIR) $(EXCOL_DIR) 2>&1)
# Success requires the positive PATCH_OK sentinel; its absence (anchor miss,
# python exception, or missing python3 → empty log) is a hard failure.
ifeq ($(findstring PATCH_OK,$(PATCH_EXORT_LOG)),)
$(error patch_exort.py failed (is python3 available?): $(PATCH_EXORT_LOG))
endif

# ---- ExoRT objects (same as n68equiv.exe target minus main.F90) ----
# Ordering matters: each file must follow the modules it USEs.
EXORT_SRCS = \
  $(MISC_DIR)/shr_kind_mod.F90        \
  $(EXCOL_DIR)/exoplanet_mod.F90      \
  $(EXCOL_DIR)/sys_rootdir.F90        \
  $(MISC_DIR)/shr_const_mod.F90       \
  $(MISC_DIR)/shr_sys_mod.F90         \
  $(MISC_DIR)/spmd_utils.F90          \
  $(MISC_DIR)/time_manager.F90        \
  $(MISC_DIR)/wrap_nf.F90             \
  $(MISC_DIR)/infnan.F90              \
  $(MISC_DIR)/ppgrid.F90              \
  $(MISC_DIR)/pmgrid.F90              \
  $(MISC_DIR)/ioFileMod.F90           \
  $(EXCOL_DIR)/radgrid.F90            \
  $(MAIN_DIR)/input.F90               \
  $(MISC_DIR)/physconst.F90           \
  $(MAIN_DIR)/output.F90              \
  $(SPEC_DIR)/cloud.F90               \
  $(SPEC_DIR)/kabs.F90                \
  $(MAIN_DIR)/mcica_random_numbers.F90 \
  $(MAIN_DIR)/mcica.F90               \
  $(SPEC_DIR)/rad_interp_mod.F90      \
  $(MAIN_DIR)/rayleigh_data.F90       \
  $(EXCOL_DIR)/calc_opd_mod.F90       \
  $(EXCOL_DIR)/initialize_rad_mod_1D.F90 \
  $(MAIN_DIR)/planck_mod.F90          \
  $(SPEC_DIR)/model_specific.F90      \
  $(MAIN_DIR)/exo_init_ref.F90        \
  $(MAIN_DIR)/exo_radiation_mod.F90

# ---- ExoColumn objects ----
EXCOL_SRCS = \
  $(EXCOL_DIR)/exocol_iapws95.F90     \
  $(EXCOL_DIR)/exocol_steam.F90       \
  $(EXCOL_DIR)/exocol_co2.F90         \
  $(EXCOL_DIR)/exocol_mod.F90         \
  $(EXCOL_DIR)/exocol_ozone.F90       \
  $(EXCOL_DIR)/exocol_config.F90      \
  $(EXCOL_DIR)/exocol_radiation.F90   \
  $(EXCOL_DIR)/exocol_io.F90          \
  $(EXCOL_DIR)/exocol_convadj.F90     \
  $(EXCOL_DIR)/exocol_surface.F90     \
  $(EXCOL_DIR)/exocol_coldstart.F90   \
  $(EXCOL_DIR)/exocol_rce_loop.F90    \
  $(EXCOL_DIR)/exocol_driver.F90

ALL_SRCS = $(EXORT_SRCS) $(EXCOL_SRCS)

# Build target
# Compile from $(SRC_DIR) so that ExoRT's input.F90 can resolve its relative
# include 'netcdf.inc' path.  $(MODFLAG) redirects generated .mod files back to
# build/.  All source paths are absolute so the CWD change is transparent.
exocol.exe: $(ALL_SRCS)
	cd $(SRC_DIR) && $(FC) $(FFLAGS) $(INC_NETCDF) -I$(MISC_DIR) \
	  $(MODFLAG) -o $(BUILD_DIR)/$@ $(ALL_SRCS) $(LIB_NETCDF)
	@mkdir -p ../run
	mv -f $(BUILD_DIR)/$@ ../run/
	@echo "Built ../run/exocol.exe"

# Phony targets
.PHONY: clean run

clean:
	$(RM) -f *.o *.mod
	$(RM) -f ../run/exocol.exe

# Convenience: run from the project root so ExoRT data paths resolve
run:
	cd .. && run/exocol.exe
