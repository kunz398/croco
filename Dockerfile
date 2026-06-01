# Dockerfile for CROCO Niue Forecast System
# Installs all dependencies; CROCO is compiled at runtime on the target machine
# All forecast data lives outside the container, mounted at /data

FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV SETUPTOOLS_USE_DISTUTILS=stdlib

# System dependencies
RUN apt-get update && apt-get install -y \
    wget curl git make cmake \
    gfortran gcc g++ cpp \
    libnetcdf-dev netcdf-bin libnetcdff-dev \
    libopenblas-dev liblapack-dev \
    libopenmpi-dev openmpi-bin \
    libhdf5-dev \
    libeccodes-dev \
    libboost-all-dev \
    nco \
    bzip2 ca-certificates \
    rsync \
    && rm -rf /var/lib/apt/lists/*

# Create croco user and /data mount point
RUN useradd -ms /bin/bash croco && \
    mkdir -p /data && \
    chown croco:croco /data

USER croco
WORKDIR /home/croco

# Install Miniconda
RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
      -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /home/croco/miniconda3 && \
    rm /tmp/miniconda.sh

ENV PATH="/home/croco/miniconda3/bin:${PATH}"

# Create conda env with Python 3.9 and pyinterp
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r && \
    conda create -n croco_forecast python=3.9 -y && \
    conda install -n croco_forecast -c conda-forge -y pyinterp && \
    conda clean -afy

# Install remaining Python packages
RUN /home/croco/miniconda3/envs/croco_forecast/bin/pip install --no-cache-dir \
    numpy==1.24.3 \
    matplotlib==3.9.4 \
    scipy==1.13.1 \
    cartopy==0.22.0 \
    netcdf4==1.7.2 \
    dask==2024.8.0 \
    distributed==2024.8.0 \
    xarray==2024.7.0 \
    geopandas \
    regionmask \
    dill \
    mpi4py \
    pydap \
    jinja2 \
    pandas \
    boto3 \
    fsspec \
    zarr \
    h5netcdf \
    h5py \
    xgcm \
    intake \
    cf_xarray \
    pyamg \
    xrft \
    numba \
    pyarrow \
    bokeh \
    lxml \
    pyyaml \
    tqdm \
    cfgrib \
    eccodes \
    copernicusmarine

ENTRYPOINT ["/bin/bash", "-c"]
CMD ["source /home/croco/miniconda3/bin/activate croco_forecast && bash /data/run_croco_docker.sh"]
