######################################################################
# Create and fill frc and bulk files with GFS data.
# for a forecast run
#
# The on-line reference to GFS is at
# http://nomad3.ncep.noaa.gov/
################################################################
#
# Common parameters
#
from datetime import datetime, timedelta
import numpy as np
from netCDF4 import Dataset
import pandas as pd
import Forecast_tools as ft
import Preprocessing_tools as ppt
from scipy.interpolate import RegularGridInterpolator
from croco_tools_params import *

# Helper function for interpolation
def interpolate_to_grid(lon_src, lat_src, var_src, lon_dst, lat_dst, method='cubic'):
    """
    Interpolate data from source grid to destination grid using RegularGridInterpolator
    """
    # Ensure var_src is a regular numpy array (not masked)
    if hasattr(var_src, 'mask'):
        var_src = np.ma.filled(var_src, fill_value=np.nan)
    var_src = np.asarray(var_src)
    
    # Create interpolator
    interp_func = RegularGridInterpolator((lat_src, lon_src), var_src, 
                                        method=method, bounds_error=False, fill_value=np.nan)
    
    # Create coordinate points for interpolation  
    points = np.column_stack([lat_dst.ravel(), lon_dst.ravel()])
    result = interp_func(points).reshape(lon_dst.shape)
    
    return result
#
it=2
keep_grib=True
#
frc_prefix=frc_prefix+'_GFS_'
blk_prefix=blk_prefix+'_GFS_'
#
################################################################
# end of user input  parameters
################################################################
#
# time (in matlab time)
#
rundate_str=datetime.today()
rundate=rundate_str-datetime(Yorig,1,1)
rundate=rundate.days
#
# GFS data name
#
gfs_name=FRCST_dir+'GFS_'+str(rundate)+'.nc'
#
#
if level==0:
  nc_suffix='.nc'
else:
  nc_suffix='.nc.'+str(level)
  grdname=grdname+'.'+str(level)
#
# Get the model grid
#
nc=Dataset(grdname, 'r',set_auto_maskandscale=False)
lon=nc['lon_rho'][:]
lat=nc['lat_rho'][:]
angle=nc['angle'][:]
h=nc['h'][:]
nc.close()
cosa=np.cos(angle)
sina=np.sin(angle)
#
# Extract data over the internet
#
#Download_data=0
if Download_data==1:
  #
  # Get model limits
  #
  lonmin=np.min(lon, axis=(0,1)) 
  lonmax=np.max(lon, axis=(0,1))
  latmin=np.min(lat, axis=(0,1))
  latmax=np.max(lat, axis=(0,1))
  #
  # Download data with DODS (the download matlab routine depends on the OGCM)
  # 
  print('Download data...')
  ft.download_GFS(rundate_str,lonmin,lonmax,latmin,latmax,FRCST_dir,Yorig,it,keep_grib=keep_grib)
#
#end
#
# Get the GFS grid 
# 
nc=Dataset(gfs_name)
lon1=nc['lon'][:]
lat1=nc['lat'][:]
time=nc['time'][:]
mask=nc['mask'][:]
tlen=len(time)
#
# bulk and forcing files
#
blkname=blk_prefix+str(rundate)+nc_suffix
print('Create a new bulk file: '+blkname)
ppt.create_bulk(blkname,grdname,CROCO_title,time,0)
nc_blk=Dataset(blkname,'r+')
frcname=frc_prefix+str(rundate)+nc_suffix
print('Create a new forcing file: '+frcname)
# Set proper time dimensions for forcing - we need momentum stress and surface flux time
ppt.create_forcing(frcname,grdname,CROCO_title,time,0,0,0,0,0,0,0,0,0,0,0)
nc_frc=Dataset(frcname,'a')  # Open in append mode

# Add wind variables to forcing file manually since they're not in create_forcing
nc_frc.createDimension('wind_time', len(time))
nw_wind_time = nc_frc.createVariable('wind_time','d',('wind_time',))
nw_wind_time.long_name = 'surface wind time'
nw_wind_time.units = 'days'
nw_wind_time[:] = time

nw_uwnd = nc_frc.createVariable('uwnd','d',('wind_time', 'eta_u', 'xi_u'))
nw_uwnd.long_name = 'u-wind component'
nw_uwnd.units = 'm s-1'

nw_vwnd = nc_frc.createVariable('vwnd','d',('wind_time', 'eta_v', 'xi_v'))
nw_vwnd.long_name = 'v-wind component'
nw_vwnd.units = 'm s-1'

nc_frc.close()
nc_frc=Dataset(frcname,'r+')  # Reopen for writing data

#
# Loop on time
#
missval=float('nan')
default=float('nan')

for l in range(tlen):
  print('time index: ',str(l+1),' of total: '+str(tlen))
  var=nc['tair'][l,:,:]
  
  # Ensure var is a proper numpy array and handle masked arrays
  if hasattr(var, 'mask'):
    var = np.ma.filled(var, fill_value=np.nan)
  var = np.asarray(var)
  
  if np.mean(np.isnan(var)!=1):
    var, interp_flag = ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    # Use helper function for interpolation
    nc_blk['tair'][l,:,:] = interpolate_to_grid(lon1, lat1, var, lon, lat, interp_method)
  else:
    var=nc['tair'][l-1,:,:] 
    # Handle masked arrays for previous timestep too
    if hasattr(var, 'mask'):
      var = np.ma.filled(var, fill_value=np.nan)
    var = np.asarray(var)
    var, interp_flag = ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    # Use helper function for interpolation
    nc_blk['tair'][l,:,:] = interpolate_to_grid(lon1, lat1, var, lon, lat, interp_method)
  #end

  var=nc['rhum'][l,:,:]
  if np.mean(np.isnan(var)!=1):
    var, interp_flag = ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    nc_blk['rhum'][l,:,:] = interpolate_to_grid(lon1, lat1, var, lon, lat, interp_method)
  else:
    var=nc['rhum'][l-1,:,:] 
    var, interp_flag = ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    nc_blk['rhum'][l,:,:] = interpolate_to_grid(lon1, lat1, var, lon, lat, interp_method)
  #end
  
  var=nc['prate'][l,:,:]
  if np.mean(np.isnan(var)!=1):
    var, interp_flag = ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    # Replaced by helper function
    nc_blk['prate'][l,:,:]= interpolate_to_grid(lon1, lat1, var, lon, lat, interp_method) 
  else:
    var=nc['prate'][l-1,:,:] 
    var, interp_flag = ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    # Replaced by helper function
    nc_blk['prate'][l,:,:]= interpolate_to_grid(lon1, lat1, var, lon, lat, interp_method)
  #end
  
  var=nc['wspd'][l,:,:]
  if np.mean(np.isnan(var)!=1):
    var, interp_flag = ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    # Replaced by helper function
    nc_blk['wspd'][l,:,:]= interpolate_to_grid(lon1, lat1, var, lon, lat, interp_method) 
  else:
    var=nc['wspd'][l-1,:,:] 
    var, interp_flag = ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    # Replaced by helper function
    nc_blk['wspd'][l,:,:]= interpolate_to_grid(lon1, lat1, var, lon, lat, interp_method)
  #end
  
  #Zonal wind speed
  var=nc['uwnd'][l,:,:]
  if np.mean(np.isnan(var)!=1):
    uwnd, interp_flag = ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    uwnd = interpolate_to_grid(lon1, lat1, uwnd, lon, lat, interp_method) 
  else:
    var=nc['uwnd'][l-1,:,:] 
    uwnd, interp_flag = ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    uwnd = interpolate_to_grid(lon1, lat1, uwnd, lon, lat, interp_method)
  #end
  
  #Meridian wind speed
  var=nc['vwnd'][l,:,:]
  if np.mean(np.isnan(var)!=1):
    vwnd, interp_flag = ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    # Removed old interp2d call
    vwnd = interpolate_to_grid(lon1, lat1, vwnd, lon, lat, interp_method) 
  else:
    var=nc['vwnd'][l-1,:,:] 
    vwnd, interp_flag = ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    # Removed old interp2d call
    vwnd = interpolate_to_grid(lon1, lat1, vwnd, lon, lat, interp_method)
  #end
  
  nc_frc['uwnd'][l,:,:]=ppt.rho2u_2d(uwnd*cosa+vwnd*sina)
  nc_frc['vwnd'][l,:,:]=ppt.rho2v_2d(vwnd*cosa-uwnd*sina)
  
  nc_blk['uwnd'][l,:,:]=ppt.rho2u_2d(uwnd*cosa+vwnd*sina)
  nc_blk['vwnd'][l,:,:]=ppt.rho2v_2d(vwnd*cosa-uwnd*sina)
  
  #Net longwave flux
  var=nc['radlw'][l,:,:]
  if np.mean(np.isnan(var)!=1):
    var, interp_flag = ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    # Replaced by helper function
    nc_blk['radlw'][l,:,:]= interpolate_to_grid(lon1, lat1, var, lon, lat, interp_method) 
  else:
    var=nc['radlw'][l-1,:,:] 
    var, interp_flag = ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    # Replaced by helper function
    nc_blk['radlw'][l,:,:]= interpolate_to_grid(lon1, lat1, var, lon, lat, interp_method)
  #end
  
  #Downward longwave flux
  var=nc['radlw_in'][l,:,:]
  if np.mean(np.isnan(var)!=1):
    var, interp_flag = ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    # Replaced by helper function
    nc_blk['radlw_in'][l,:,:]= interpolate_to_grid(lon1, lat1, var, lon, lat, interp_method)
  else:
    var=nc['radlw_in'][l-1,:,:] 
    var, interp_flag = ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    # Replaced by helper function
    nc_blk['radlw_in'][l,:,:]= interpolate_to_grid(lon1, lat1, var, lon, lat, interp_method)
  #end
   
  #Net solar short wave radiation
  var=nc['radsw'][l,:,:]
  if np.mean(np.isnan(var)!=1):
    var, interp_flag = ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    # Replaced by helper function
    nc_blk['radsw'][l,:,:]= interpolate_to_grid(lon1, lat1, var, lon, lat, interp_method)
  else:
    var=nc['radsw'][l-1,:,:] 
    var, interp_flag = ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    # Replaced by helper function
    nc_blk['radsw'][l,:,:]= interpolate_to_grid(lon1, lat1, var, lon, lat, interp_method)
  #end
  
  var=nc['tx'][l,:,:]
  if np.mean(np.isnan(var)!=1):
    tx, interp_flag = ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    # Removed old interp2d call
    tx = interpolate_to_grid(lon1, lat1, tx, lon, lat, interp_method)
  else:
    var=nc['tx'][l,:,:]
    tx, interp_flag = ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    # Removed old interp2d call
    tx = interpolate_to_grid(lon1, lat1, tx, lon, lat, interp_method)
  #end
  
  var=nc['ty'][l,:,:]
  if np.mean(np.isnan(var)!=1):
    ty, interp_flag = ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    # Removed old interp2d call
    ty = interpolate_to_grid(lon1, lat1, ty, lon, lat, interp_method)
  else:
    var=nc['ty'][l,:,:]
    ty, interp_flag = ppt.get_missing_val(lon1,lat1,var,missval,Roa,default)
    # Removed old interp2d call
    ty = interpolate_to_grid(lon1, lat1, ty, lon, lat, interp_method)
  #end
  
  nc_frc['sustr'][l,:,:]=ppt.rho2u_2d(tx*cosa+ty*sina)
  nc_frc['svstr'][l,:,:]=ppt.rho2v_2d(ty*cosa-tx*sina)
  
  nc_blk['sustr'][l,:,:]=ppt.rho2u_2d(tx*cosa+ty*sina)
  nc_blk['svstr'][l,:,:]=ppt.rho2v_2d(ty*cosa-tx*sina)
#end
# 
nc_frc.close()
nc_blk.close()
nc.close()