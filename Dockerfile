# Dockerfile for CROCO Python Forecast System
# Image provides compilers + Python packages only.
# Project code and data are always supplied via volume mount at runtime.
# Nothing is copied into the image — builds stay fast and incremental.
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV SETUPTOOLS_USE_DISTUTILS=stdlib

# Install system dependencies (compilers, NetCDF, MPI, Python)
RUN apt-get update && apt-get install -y \
    python3 python3-pip python3-dev \
    ca-certificates curl \
    gfortran gcc make cmake \
    build-essential \
    libeigen3-dev \
    libgsl-dev \
    libnetcdf-dev netcdf-bin libnetcdff-dev \
    libopenblas-dev \
    liblapack-dev \
    libgeos-dev \
    libproj-dev \
    proj-data \
    libwxgtk3.0-gtk3-dev \
    libboost-all-dev \
    git \
    mpich libmpich-dev \
    libeccodes-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies (cached layer — only rebuilt when this list changes)
RUN pip3 install --upgrade pip setuptools wheel && \
    pip3 install numpy==1.23.5 matplotlib==3.6.2 scipy==1.9.3 cartopy==0.21.0 \
    traits==6.4.1 traitsui==7.4.2 netcdf4==1.6.2 dask==2022.11.1 'xarray>=2023.4.0' \
    geopandas==0.12.1 regionmask==0.9.0 pyinterp==2022.10.1 dill copernicusmarine pydap jinja2 \
    cfgrib eccodes requests s3fs boto3 fsspec h5py h5netcdf \
    xgcm cf_xarray intake pandas distributed pyamg xrft numba

# Install nco (provides ncrcat CLI tool used by postprocess.py)
RUN apt-get update && apt-get install -y nco && rm -rf /var/lib/apt/lists/*

# Fix 1: Repair the _sysconfigdata file corrupted by the previous bad patch run.
# Use sed (atomic in-place edit) — the Python list-comprehension approach
# truncated the file before reading it (open(f,'w') evaluated before open(f).read()).
RUN sed -i 's/-Wl,-Bsymbolic-functions//g' \
    /usr/lib/python3.10/_sysconfigdata__x86_64-linux-gnu.py

# Fix 2: Strip linker hardening flags from /usr/bin/nf-config itself.
# The flags (-Wl,-Bsymbolic-functions, -flto=auto, -Wl,-z,relro, etc.) are
# hardcoded in the nf-config shell script's flibs= variable. They get passed
# by PREPRO's Makefile NETCDFLIB=$(shell nf-config --flibs) to f2py, which
# passes them to gfortran as source-file arguments → "error: unknown file type".
RUN sed -i \
    's/ -Wl,-Bsymbolic-functions//g; s/ -flto=auto//g; s/ -ffat-lto-objects//g; s/ -Wl,-z,relro//g; s/ -Wl,-z,now//g' \
    /usr/bin/nf-config

# Working directory matches the volume mount point used at runtime
WORKDIR /home/croco/croco_pytools

ENV PYTHONPATH="/home/croco/croco_pytools:/home/croco/croco_pytools/cp_g_col_croco_pytools/Forecast_CROCO:/home/croco/croco_pytools/cp_g_col_croco_pytools/PREPRO:/home/croco/croco_pytools/cp_g_col_croco_pytools/PREPRO/Modules"

VOLUME ["/home/croco/croco_pytools"]

ENTRYPOINT ["/bin/bash", "-c"]
CMD ["bash /home/croco/croco_pytools/run_croco.sh"]
