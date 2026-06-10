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

# Compiler flags
ifeq ($(FC),$(IFX))
  FFLAGS = -O2 -traceback -fpe0
else ifeq ($(FC),gfortran)
  FFLAGS = -O2 -fbacktrace -ffpe-trap=invalid,zero,overflow -cpp
else
  # Fallback: assume ifort-compatible flags
  FFLAGS = -O2 -traceback -fpe0
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

# Generate src/calc_opd_mod.F90 from ExoRT's n68equiv source with the
# k-coefficient pressure interpolation CLAMPED at the table top (10 bar).
# Above its pressure grid ExoRT linearly EXTRAPOLATES k in log10(p)
# (rad_interp_mod: "at the tops of grids, extrapolate"), which is unphysical
# for >10 bar columns and can even produce negative k.  The clamp freezes the
# pressure broadening at its 10-bar value for deeper layers; partial pressures
# and amagats for the CIA/continuum paths still use the true pmid.  Columns
# that stay below 10 bar are BIT-IDENTICAL (the min() is a no-op).  Same
# local-copy pattern as exoplanet_mod above; ExoRT's file is preserved.
$(shell sed 's|pressure = log10(pmid(ik))       ! log pressure|pressure = min(log10(pmid(ik)), log10pgrid(kc_npress))  ! log pressure; CLAMPED at k-table top by ExoColumn build|' \
  $(SPEC_DIR)/calc_opd_mod.F90 > $(EXCOL_DIR)/calc_opd_mod.F90)
# Fail loudly if the anchor line ever changes in a future ExoRT update (the
# sed above would otherwise silently produce an unpatched copy).
ifeq ($(shell grep -c 'CLAMPED at k-table top' $(EXCOL_DIR)/calc_opd_mod.F90),0)
$(error calc_opd_mod clamp patch failed: anchor line not found in $(SPEC_DIR)/calc_opd_mod.F90 — ExoRT source changed?)
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
  $(SPEC_DIR)/radgrid.F90             \
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
  $(SPEC_DIR)/initialize_rad_mod_1D.F90 \
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
# include 'netcdf.inc' path.  -module redirects generated .mod files back to
# build/.  All source paths are absolute so the CWD change is transparent.
exocol.exe: $(ALL_SRCS)
	cd $(SRC_DIR) && $(FC) $(FFLAGS) $(INC_NETCDF) -I$(MISC_DIR) \
	  -module $(BUILD_DIR) -o $(BUILD_DIR)/$@ $(ALL_SRCS) $(LIB_NETCDF)
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
