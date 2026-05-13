#!/bin/bash
# Bash script to run CROCO forecast with proper directory setup

rundir=$(pwd)
folder=$(date +%d-%m-%Y)
yfolder=$(date -d "yesterday" +%d-%m-%Y)
targetDir="$rundir/$folder"

export PATH="/home/ekunals/anaconda3/bin:/home/ekunals/.local/bin:$PATH"

# ---------------------------------------------------------------------------
# Logging setup — tee all output (stdout + stderr) to the same log file that
# the Python forecast script uses, so the entire run is captured in one place.
# ---------------------------------------------------------------------------
LOG_DIR="$rundir/logs"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/${folder}_run_forecast.log"
exec > >(tee -a "$LOG_FILE") 2>&1

# Helpers for consistent, timestamped log lines
log()     { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
log_ok()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] OK    $*"; }
log_warn(){ echo "[$(date '+%Y-%m-%d %H:%M:%S')] WARN  $*"; }
log_err() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR $*"; }
section() {
    echo ""
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ============================================================"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ============================================================"
}

# ---------------------------------------------------------------------------
section "CROCO FORECAST RUN — $folder"
log "Log file: $LOG_FILE"

# Detect the uid:gid of the workspace directory owner (not of the calling user).
# This ensures container processes always match the filesystem owner, which is
# required on FUSE-mounted drives (NTFS/exFAT) where root cannot write to
# files owned by a different uid.
HOST_UID=$(stat -c '%u' "${rundir}")
HOST_GID=$(stat -c '%g' "${rundir}")
log "Running containers as uid:gid = ${HOST_UID}:${HOST_GID} (owner of workspace)"

# ---------------------------------------------------------------------------
# STEP 1: Compile CROCO ocean model using the  jobcomp script.
# jobcomp handles: source copying, PISCES/LIB/SED flattening, XIOS stubs,
# MUSTANG, OBSTRUCTION, Makedefs generation, make depend, make -j, and
# moves the croco binary to CONFIGS/Forecast_Niue1k/.
# ---------------------------------------------------------------------------
#comment out compile from here to -->>

# JOBCOMP_DIR="${rundir}/croco-v2.1.0/CONFIGS/Forecast_Niue1k"
# echo "Compiling CROCO ocean model (jobcomp)..."
# ( cd "$JOBCOMP_DIR" && NPROCS=$(nproc) bash jobcomp ) 2>&1 | tee /tmp/jobcomp.log
# if [ ${PIPESTATUS[0]} -ne 0 ]; then
#     echo "ERROR: CROCO compilation failed. See /tmp/jobcomp.log"
#     exit 1
# fi
# echo "CROCO compilation complete."
## <<-- end of compile section

# PREPRO Fortran/f2py extension is compiled inside the Docker container (Step 3 below)
# so that the .so is built for the container's Python ABI, not the host's.
# The croco binary produced by jobcomp lands in CONFIGS/Forecast_Niue1k/ (same as SRC_DIR below).

# ---------------------------------------------------------------------------
section "STEP 1: Update FORECAST_BASE"

# Install compiled model into FORECAST_BASE (backup existing files)
# croco binary is taken directly from the workspace root ($rundir/croco)
TARGET_BASE="$rundir/FORECAST_BASE"

log "Source binary : $rundir/croco"
log "Target base   : $TARGET_BASE"
mkdir -p "$TARGET_BASE"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
BACKUP_DIR="$TARGET_BASE/backup_$TIMESTAMP"
mkdir -p "$BACKUP_DIR"
log "Backup dir    : $BACKUP_DIR"
for f in croco; do
    if [ -f "$TARGET_BASE/$f" ]; then
        mv -v "$TARGET_BASE/$f" "$BACKUP_DIR/"
    fi
done
if [ -f "$rundir/croco" ]; then
    cp -pv "$rundir/croco" "$TARGET_BASE/"
    chmod +x "$TARGET_BASE/croco"
    log_ok "FORECAST_BASE croco binary updated."
else
    log_warn "$rundir/croco not found — keeping existing FORECAST_BASE binary."
fi

section "STEP 2: Set up forecast directory — $folder"
log "Target directory: $targetDir"

# Remove old directory if it exists
if [ -d "$targetDir" ]; then
    log "Removing existing directory..."
    rm -rf "$targetDir"
fi

# Create directory structure on host
log "Creating directory structure..."
mkdir -p "$targetDir"
mkdir -p "$targetDir/CROCO_FILES"
mkdir -p "$targetDir/DATA/Forecast"
mkdir -p "$targetDir/FORECAST"
mkdir -p "$targetDir/SCRATCH"

# Copy files from FORECAST_BASE
log "Copying files from FORECAST_BASE..."
forecastBase="$rundir/FORECAST_BASE"

if [ -f "$forecastBase/croco" ]; then
    cp "$forecastBase/croco" "$targetDir/"
fi
if [ -f "$forecastBase/croco_forecast.in" ]; then
    cp "$forecastBase/croco_forecast.in" "$targetDir/"
fi
if [ -d "$forecastBase/CROCO_FILES" ]; then
    cp -r "$forecastBase/CROCO_FILES"/* "$targetDir/CROCO_FILES/" 2>/dev/null || true
fi

log_ok "Directory setup complete. Structure:"
# Ensure directory is writable by container user (permissions fix for sudo usage)
chmod -R 777 "$targetDir"
ls -1 "$targetDir" | while read -r line; do log "  $line"; done

section "STEP 3: Configure parameters"
paramsFile="$rundir/cp_g_col_croco_pytools/Forecast_CROCO/croco_tools_params.py"
if [ -f "$paramsFile" ]; then
    newRunDir="/home/croco/croco_pytools/$folder/"
    sed -i "s|RUN_dir = '[^']*'|RUN_dir = '$newRunDir'|" "$paramsFile"
    log_ok "croco_tools_params.py: RUN_dir updated to $newRunDir"
fi

# Load Copernicus Marine credentials if available
credsFile="$rundir/cp_g_col_croco_pytools/Forecast_CROCO/config/.copernicusmarine-credentials"
if [ -f "$credsFile" ]; then
    log "Loading Copernicus Marine credentials..."
    decoded=$(base64 -d "$credsFile" 2>/dev/null)
    if [ -n "$decoded" ]; then
        export COPERNICUSMARINE_SERVICE_USERNAME=$(echo "$decoded" | grep username | cut -d'=' -f2)
        export COPERNICUSMARINE_SERVICE_PASSWORD=$(echo "$decoded" | grep password | cut -d'=' -f2)
        log_ok "Copernicus Marine credentials loaded."
    fi
fi

section "STEP 4: Docker container — preprocessing + MPI run"
log "Starting Docker container..."
# docker run --rm -i \
#     --user "${HOST_UID}:${HOST_GID}" \
#     --ipc=host \
#     --cpuset-cpus "0-17" \
#     -e HOME=/tmp \
#     -e COPERNICUSMARINE_SERVICE_USERNAME="$COPERNICUSMARINE_SERVICE_USERNAME" \
#     -e COPERNICUSMARINE_SERVICE_PASSWORD="$COPERNICUSMARINE_SERVICE_PASSWORD" \
#     -v "${rundir}:/home/croco/croco_pytools" \
#     --entrypoint "/bin/bash" \
#     croco-forecast \
#     -s <<'FORECAST_EOF'
docker run --rm -i \
    --user "${HOST_UID}:${HOST_GID}" \
    --ipc=host \
    --cpuset-cpus "0-35" \
    --cpu-shares 1024 \
    --ulimit stack=-1 \
    --ulimit memlock=-1 \
    -e HOME=/tmp \
    -e OMP_NUM_THREADS=2 \
    -e OMPI_MCA_btl=vader,self \
    -e COPERNICUSMARINE_SERVICE_USERNAME="$COPERNICUSMARINE_SERVICE_USERNAME" \
    -e COPERNICUSMARINE_SERVICE_PASSWORD="$COPERNICUSMARINE_SERVICE_PASSWORD" \
    -e FORECAST_FOLDER="$folder" \
    -e YESTERDAY_FOLDER="$yfolder" \
    -v "${rundir}:/home/croco/croco_pytools" \
    --entrypoint "/bin/bash" \
    croco-forecast \
    -s <<'FORECAST_EOF'
set -e
set -o pipefail
# Ensure the host-mounted logs directory exists inside the container
mkdir -p /home/croco/croco_pytools/logs

_ts() { date '+%Y-%m-%d %H:%M:%S'; }
_log()     { echo "[$(_ts)] $*"; }
_section() { echo ""; echo "[$(_ts)] ============================================================"; echo "[$(_ts)] $*"; echo "[$(_ts)] ============================================================"; }

# Step 3a: compile PREPRO toolsf.so for this container's Python version
# Commented out to skip recompilation on every run (saves significant time).
# Uncomment below if toolsf.so needs to be rebuilt 

# from here -->>>
#_section "Compiling PREPRO inside container"
#cd /home/croco/croco_pytools/cp_g_col_croco_pytools/prepro/Modules/tools_fort_routines
#make clean && make
# <<<-- till here

# Step 3b: run the forecast (preprocessing / boundary conditions / initial conditions)
_section "[Docker 3b] Python pre-processing — GFS + OGCM"
cd /home/croco/croco_pytools/cp_g_col_croco_pytools/Forecast_CROCO
python3 run_croco_forecast.py

# Step 3c: hot restart patching (before MPI run)
_section "[Docker 3c] Hot restart check"
FOLDER=${FORECAST_FOLDER:-$(date +%d-%m-%Y)}
YFOLDER=${YESTERDAY_FOLDER:-$(date -d "yesterday" +%d-%m-%Y)}
SCRATCH=/home/croco/croco_pytools/$FOLDER/SCRATCH
YRST=/home/croco/croco_pytools/$YFOLDER/SCRATCH/croco_rst_$YFOLDER.nc

if [ -f "$YRST" ]; then
    _log "OK    Hot restart: found $YRST — patching croco_forecast.in"
    # Replace ini file with yesterday's restart (relative path from SCRATCH)
    sed -i "s|croco_ini.nc|../../$YFOLDER/SCRATCH/croco_rst_$YFOLDER.nc|" $SCRATCH/croco_forecast.in
else
    _log "WARN  No yesterday restart found ($YRST) — cold start"
fi
# Always name today's restart with the date so tomorrow can pick it up
sed -i "s|croco_rst.nc|croco_rst_$FOLDER.nc|" $SCRATCH/croco_forecast.in

# Step 3d: run CROCO ocean model with MPI
_section "[Docker 3d] CROCO ocean model — MPI run"
cd /home/croco/croco_pytools/$FOLDER/SCRATCH
_log "Working dir : $(pwd)"
_log "MPI log     : /home/croco/croco_pytools/logs/${FOLDER}_mpi.log"
mpirun -np 18 \
    --allow-run-as-root \
    --map-by core \
    --bind-to core \
    --mca btl vader,self \
    --mca rtc ^hwloc \
    ./croco croco_forecast.in 2>&1 | tee /home/croco/croco_pytools/logs/${FOLDER}_mpi.log

# Step 3e: post-process CROCO output - convert sigma to z-levels and concatenate
_section "[Docker 3e] Sigma-to-z post-processing"
_log "Input  : /home/croco/croco_pytools/$FOLDER/SCRATCH"
_log "Log    : /home/croco/croco_pytools/logs/${FOLDER}_postprocess.log"
python3 /home/croco/croco_pytools/FromGaby/postprocess.py /home/croco/croco_pytools/$FOLDER/SCRATCH \
    2>&1 | tee /home/croco/croco_pytools/logs/${FOLDER}_postprocess.log
FORECAST_EOF

if [ $? -ne 0 ]; then
  log_err "Docker container failed (PREPRO compile or forecast error)."
  exit 1
fi

section "STEP 5: Transfer output to ocean portal"
DEST_DIR="/data/ocean_portal/datasets/model/country/spc/forecast/hourly/NIU_Currents"
SRC_FILE="$targetDir/SCRATCH/d1_temp_salt_uv_z_all.nc"
log "Source : $SRC_FILE"
log "Dest   : $DEST_DIR"
if [ -f "$SRC_FILE" ]; then
    mkdir -p "$DEST_DIR"
    cp -v "$SRC_FILE" "$DEST_DIR/"
    log_ok "Transfer complete."
else
    log_warn "$SRC_FILE not found — skipping transfer."
fi

section "CROCO FORECAST COMPLETE — $folder"
log "Full log: $LOG_FILE"
