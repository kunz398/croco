#bash

basedir=/DATA/CROCO/croco-v2.1.0/CONFIGS/
rundir=/media/claudiaa/Elements/eDRIVE/
#basedir=/
folder=$(date +%d-%m-%Y) 
yfolder=$(date -d "yesterday" +%d-%m-%Y)

edir+=$basedir"FORECAST_BASE"
echo "$edir"


cp -r $edir  $rundir$folder

source ~/miniconda3/bin/activate

conda activate croco_forecast

cd /DATA/CROCO/cp_g_col_croco_pytools/Forecast_CROCO/

sed -i "s|RUN_dir = '.*'|RUN_dir = '$rundir$folder/'|" croco_tools_params.py

python run_croco_forecast.py

conda deactivate
cd $rundir$folder/SCRATCH
cp $rundir$yfolder/SCRATCH/croco_rst$yfolder.nc ./

echo sed -i "s|croco_ini.nc|croco_rst_$yfolder.nc|" croco_forecast.in
echo sed -i "s|croco_rst.nc|croco_rst_$folder.nc|" croco_forecast.in


mpirun -np 18 croco croco_forecast.in > $(date +%d-%m-%Y).log 

