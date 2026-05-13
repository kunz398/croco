#!/bin/bash
# Simplified version - assumes directory structure is already created by PowerShell script

rundir=/home/croco/croco_pytools/
folder=$(date +%d-%m-%Y)
yfolder=$(date -d "yesterday" +%d-%m-%Y)

echo "Working directory: $rundir$folder"
echo "Checking directory structure..."
ls -la "$rundir$folder" 2>&1 || echo "ERROR: Directory not found!"

# Verify required subdirectories exist
for dir in CROCO_FILES DATA FORECAST SCRATCH; do
    if [ ! -d "$rundir$folder/$dir" ]; then
        echo "WARNING: $dir directory missing, creating..."
        mkdir -p "$rundir$folder/$dir"
    fi
done





cd cp_g_col_croco_pytools/Forecast_CROCO/

echo $pwd
echo pwd

sed -i "s|RUN_dir = '.*'|RUN_dir = '$rundir$folder/'|" croco_tools_params.py

echo  sed -i "s|RUN_dir = '.*'|RUN_dir = '$rundir$folder/'|" croco_tools_params.py

cd $rundir$folder/SCRATCH
cp $rundir$yfolder/SCRATCH/croco_rst$yfolder.nc ./

echo sed -i "s|croco_ini.nc|croco_rst_$yfolder.nc|" croco_forecast.in
echo sed -i "s|croco_rst.nc|croco_rst_$folder.nc|" croco_forecast.in


#mpirun -np 18 croco croco_forecast.in > $(date +%d-%m-%Y).log



python3 run_croco_forecast.py
