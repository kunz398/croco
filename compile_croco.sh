#!/bin/bash
# Dedicated CROCO compilation script for the Forecast_Niue1k config inside Docker.
# Called by run_forecast.sh.  All paths are container-absolute (volume-mounted).
# Config: MPI + gfortran, no AGRIF/XIOS/OASIS/OpenMP.
set -e

# BASE_DIR is the workspace root (host path or container path)
BASE_DIR="${BASE_DIR:-/home/croco/croco_pytools}"

OCEAN_SRC="${BASE_DIR}/croco-v2.1.0/OCEAN"
RUNDIR="${BASE_DIR}/croco-v2.1.0/CONFIGS/Forecast_Niue1k"
# Use /tmp so the scratch dir is always writable regardless of NTFS ownership
SCRDIR="/tmp/croco_compile"
ROOT_DIR="${BASE_DIR}/croco-v2.1.0"

echo "=== CROCO compile: RUNDIR=$RUNDIR ==="
cd "$RUNDIR"

# ---- Compiler and flag setup -------------------------------------------
FC=gfortran
# Use the apt-installed mpich mpif90 (installed via apt in the Docker image).
MPIF90=/usr/bin/mpif90
NETCDFINC=$(nf-config --includedir)
NETCDFLIB=$(nf-config --flibs)

CPP1="cpp -traditional -DLinux"
CFT1="$MPIF90"           # MPI is defined in cppdefs.h
FFLAGS1="-O2 -mcmodel=medium -fdefault-real-8 -fdefault-double-8 -std=legacy"
LDFLAGS1="$NETCDFLIB"
CPPFLAGS1="-I${NETCDFINC} -ICROCOFILES/AGRIF_INC"

NPROCS=$(nproc)

# ---- Clean and recreate scratch dir ------------------------------------
echo "Cleaning Compile dir..."
rm -rf "$SCRDIR"
mkdir -p "$SCRDIR/CROCOFILES/AGRIF_INC"
mkdir -p "$SCRDIR/CROCOFILES/AGRIF_MODELFILES"

# ---- Copy OCEAN sources ------------------------------------------------
echo "Copying OCEAN sources from $OCEAN_SRC ..."
cp "$OCEAN_SRC"/*.F   "$SCRDIR/"
cp "$OCEAN_SRC"/*.F90 "$SCRDIR/" 2>/dev/null || true
cp "$OCEAN_SRC"/*.h   "$SCRDIR/"
cp "$OCEAN_SRC"/*.py  "$SCRDIR/" 2>/dev/null || true
cp "$OCEAN_SRC"/Make* "$SCRDIR/"
cp "$OCEAN_SRC"/jobcomp "$SCRDIR/" 2>/dev/null || true
cp "$OCEAN_SRC"/amr.in  "$SCRDIR/" 2>/dev/null || true

# AGRIF
[ -d "$ROOT_DIR/AGRIF" ] && cp -r "$ROOT_DIR/AGRIF" "$SCRDIR/"

# PISCES — copy top-level files first, then flatten subdirectories
if [ -d "$ROOT_DIR/PISCES" ]; then
    cp -r "$ROOT_DIR/PISCES"/* "$SCRDIR/" 2>/dev/null || true
    # Flatten LIB/ into compile dir (oce_trc.F90, in_out_manager.F90, iom.F90, lib_mpp.F90, prtctl.F90, ...)
    [ -d "$ROOT_DIR/PISCES/LIB" ] && cp -r "$ROOT_DIR/PISCES/LIB"/* "$SCRDIR/" 2>/dev/null || true
    # Flatten SED/ into compile dir (sedadv.F90, sedbtb.F90, oce_sed.F90, par_sed.F90, ...)
    [ -d "$ROOT_DIR/PISCES/SED" ] && cp    "$ROOT_DIR/PISCES/SED"/* "$SCRDIR/" 2>/dev/null || true
fi

# MUSTANG
[ -d "$ROOT_DIR/MUSTANG" ] && cp -r "$ROOT_DIR/MUSTANG"/* "$SCRDIR/" 2>/dev/null || true

# OBSTRUCTION
[ -d "$ROOT_DIR/OBSTRUCTION" ] && cp -r "$ROOT_DIR/OBSTRUCTION"/* "$SCRDIR/" 2>/dev/null || true

# XIOS (send_xios_diags.F and any other Fortran stubs)
[ -d "$ROOT_DIR/XIOS" ] && cp "$ROOT_DIR/XIOS"/*.F "$SCRDIR/" 2>/dev/null || true

# ---- Overwrite with local config files ---------------------------------
echo "Applying local config overrides..."
for ext in F F90 h h90; do
    ls "$RUNDIR"/*.$ext > /dev/null 2>&1 && cp -f "$RUNDIR"/*.$ext "$SCRDIR/" || true
done
cp -f "$RUNDIR"/Make* "$SCRDIR/" 2>/dev/null || true

# ---- Build Makedefs from template --------------------------------------
cd "$SCRDIR"
sed \
    -e "s?\$(FFLAGS1)?${FFLAGS1}?g" \
    -e "s?\$(LDFLAGS1)?${LDFLAGS1}?g" \
    -e "s?\$(CPP1)?${CPP1}?g" \
    -e "s?\$(CFT1)?${CFT1}?g" \
    -e "s?\$(CPPFLAGS1)?${CPPFLAGS1}?g" \
    Makedefs.generic > Makedefs

# Suppress GNU make's built-in m2c pattern rule (%.o: %.mod).
# CROCO/PISCES drops circular .mod<-.o deps; without this override make tries
# to rebuild the .o from the .mod using m2c (a MATLAB tool) which isn't installed.
echo "" >> Makedefs
echo "%.o: %.mod ;" >> Makedefs

echo "--- Makedefs ---"
cat Makedefs
echo "----------------"

# ---- Compile -----------------------------------------------------------
echo "Running make depend..."
make depend

echo "Running make -j ${NPROCS}..."
make -j "${NPROCS}"

# ---- Install binary back to RUNDIR -------------------------------------
if [ -f croco ]; then
    cp croco "$RUNDIR/"
    chmod +x "$RUNDIR/croco"
    echo "=== croco binary installed to $RUNDIR/croco ==="
else
    echo "ERROR: croco binary not produced after make"
    exit 1
fi
[ -f partit ] && cp partit "$RUNDIR/"
[ -f ncjoin ] && cp ncjoin  "$RUNDIR/"
