######################################################################
#
# Extract a subgrid from mercator to get a CROCO forcing
# Store that into monthly files.
# Take care of the Greenwitch Meridian.

############################# LYBRARIES ##############################
from datetime import datetime, timedelta
import numpy as np
from netCDF4 import Dataset 
import pandas as pd
import get_file_python_mercator as gfm
from pydap.client import open_url
import pprint
import subprocess
import Preprocessing_tools as ppt
import oforc_OGCM as ofgm
import os
import sys
from cp_g_col_croco_pytools.Forecast_CROCO.old__croco_tools_params import *
import boto3
from botocore import UNSIGNED
from botocore.config import Config
import xarray as xr
import tempfile

######################################################################
#
# Extract a subset from Marcator using python motu client (cls)
# Write it in a local file (keeping the classic SODA netcdf format)
######################################################################
def write_mercator_frcst(FRCST_dir,FRCST_prefix,raw_mercator_name,mercator_type,vars,time,Yorig):
    #
    print('    writing MERCATOR file')
    #
    # Get grid and time frame
    #
    nc = Dataset(raw_mercator_name, set_auto_maskandscale=True)
    if mercator_type==1:
        lon = nc.variables['longitude'][:]
        lat = nc.variables['latitude'][:]
        depth = nc.variables['depth'][:]
        time = nc.variables['time'][:]
        diff_dates = (datetime(1950,1,1) - datetime(Yorig,1,1)).days
        #time = time.astype('float32') / 24. + diff_dates.astype('float32')
        time = time/24. + diff_dates
    else:
        lon = nc.variables['longitude'][:]
        lat = nc.variables['latitude'][:]
        depth = nc.variables['depth'][:]
        time = nc.variables['time'][:]
        # Fix: Properly convert from Mercator time reference to CROCO time reference
        # Mercator uses 1950-01-01 as reference, CROCO uses Yorig (2005) as reference
        # Time is in hours, so convert to days and adjust reference
        diff_dates = (datetime(1950,1,1) - datetime(Yorig,1,1)).days
        time = time/24. + diff_dates
        
    #
    # Get SSH
    #
    print('    ...SSH')
    if 'zos' in nc.variables:
        ncc = nc['zos']
        ssh = ncc[:,:,:]
    else:
        ssh = np.zeros((len(nc.dimensions['time']), len(nc.dimensions['latitude']), len(nc.dimensions['longitude'])))
    #
    # Get U
    #
    print('    ...U')
    if 'uo' in nc.variables:
        ncc = nc['uo']
        u = ncc[:,:,:,:]
    else:
        u = np.zeros((len(nc.dimensions['time']), len(nc.dimensions['depth']), len(nc.dimensions['latitude']), len(nc.dimensions['longitude'])-1))
    #
    # Get V
    #
    print('    ...V')
    if 'vo' in nc.variables:
        ncc = nc['vo']
        v = ncc[:,:,:,:]
    else:
        v = np.zeros((len(nc.dimensions['time']), len(nc.dimensions['depth']), len(nc.dimensions['latitude'])-1, len(nc.dimensions['longitude'])))
    #
    # Get TEMP
    #
    print('    ...TEMP')
    if 'thetao' in nc.variables:
        ncc = nc['thetao']
        temp = ncc[:,:,:,:]
    else:
        temp = np.zeros((len(nc.dimensions['time']), len(nc.dimensions['depth']), len(nc.dimensions['latitude']), len(nc.dimensions['longitude'])))
    #
    # Get SALT
    #
    print('    ...SALT')
    if 'so' in nc.variables:
        ncc = nc['so']
        salt = ncc[:,:,:,:]
    else:
        salt = np.zeros((len(nc.dimensions['time']), len(nc.dimensions['depth']), len(nc.dimensions['latitude']), len(nc.dimensions['longitude'])))

    #
    # Create the Mercator file
    #
    rundate_str=datetime.today()
    rundate=rundate_str-datetime(Yorig,1,1)
    nc.close()


    ofgm.create_OGCM(FRCST_dir+FRCST_prefix+str(rundate.days)+'.cdf',lon,lat,lon,lat,lon,lat,depth,time,\
      (temp),(salt),(u),(v),(ssh),Yorig)

    #
    return


###############################################################################
def download_mercator(pathMotu,user,password,mercator_type,lh,lf,lonmin,lonmax,latmin,latmax,zmax,FRCST_dir,FRCST_prefix,raw_mercator_name,Yorig):
  # pathMotu is a deprecated parameter !
  raw_exist = os.path.exists(raw_mercator_name)
  if not raw_exist: 
    download_raw_data=1
    print(' ')
    print('Downloading Raw Mercator File')
  else:
    download_raw_data=0
    print(' ')
    print('Raw Mercator file already downloaded, converting to croco_tools')
    print(' ')

  
  convert_raw2crocotools=1 # convert -> crocotools format data
  #
  # Set variable names according to mercator type data
  #
  if mercator_type==1:
    vars = ' --variable zos --variable uo --variable vo --variable thetao --variable so '
    var_names = ['zos', 'uo', 'vo', 'thetao', 'so']
  else:
    vars = ' --variable zos --variable uo --variable vo --variable thetao --variable so '
    var_names = ['zos', 'uo', 'vo', 'thetao', 'so']

  #
  # Get dates
  #
  rundate_str=datetime.today()
  rundate=rundate_str-datetime(Yorig,1,1)
  rundate=rundate.days

  time1 = []
  for i in range(1,lh+1):
    time1.append(rundate_str-timedelta(days=lh+2-i))

  time2=rundate_str
  time3 = []
  for j in range(1,lf+2):
    time3.append(rundate_str+timedelta(days=j))

  time=[*time1,time2,*time3]
  tiempo_inicial = time1[0]
  tiempo_final = time3[-1]

  if(lonmin > 180):
    lonmin = lonmin - 360
    
  if(lonmax > 180):
    lonmax = lonmax - 360

  print(' ')
  print('Get data for '+str(rundate_str))
  print('Minimum Longitude: '+str(lonmin))
  print('Maximum Longitude: '+str(lonmax))
  print('Minimum Latitude:  '+str(latmin))
  print('Maximum Latitude:  '+str(latmax))
  print(' ')

  if download_raw_data==1:
    #
    # Get data 
    #
    gfm.get_mercator(mercator_type, vars,
    [str(lonmin-1), str(lonmax+1), str(latmin-1), str(latmax+1), str(0), str(zmax)],
    [str(tiempo_inicial), str(tiempo_final)],
    [user, password],
    raw_mercator_name)

  if convert_raw2crocotools==1: 
    #
    # Convert data format and write in a more CROCOTOOLS 
    # compatible input file 
    #
    print('Making output data directory '+FRCST_dir) # create directory
    frcst_ext = os.path.exists(FRCST_dir)
    if not frcst_ext:
      os.makedirs(FRCST_dir)
        #
    mercator_name=FRCST_dir+FRCST_prefix+str(rundate)+'.cdf'
    if mercator_name != None:
      print('Mercator file already exist => overwrite it')
    
    write_mercator_frcst(FRCST_dir,FRCST_prefix,raw_mercator_name,mercator_type,var_names,time,Yorig)

def interp_OGCM_frcst(OGCM_name,Roa,interp_method,lonU,latU,lonV,latV,lonT,latT,Z,tin,nc_clm,nc_bry,lon,lat,angle,h,pm,pn,rmask,tout,vtransform,obc):
  conserv=1 # same barotropic velocities as the OGCM
  #
  print(['  Horizontal interpolation: ',OGCM_name])
  #
  # CROCO grid angle
  #
  cosa=np.cos(angle)
  sina=np.sin(angle)
  #
  # Open the OGCM file
  #
  nc=Dataset(OGCM_name, 'r', set_auto_maskandscale=False)
  #
  # Interpole data on the OGCM Z grid and CROCO horizontal grid
  #
  #
  # Read and extrapole the 2D variables
  #
  zeta=ofgm.ext_data_OGCM(nc,lonT,latT,'ssh',tin,lon,lat,1,Roa,interp_method)
  u2d=ofgm.ext_data_OGCM(nc,lonU,latU,'ubar',tin,lon,lat,1,Roa,interp_method)
  v2d=ofgm.ext_data_OGCM(nc,lonV,latV,'vbar',tin,lon,lat,1,Roa,interp_method)
  ubar=ppt.rho2u_2d(u2d*cosa+v2d*sina)
  vbar=ppt.rho2v_2d(v2d*cosa-u2d*sina)
  #
  # Read and extrapole the 3D variables
  #
  NZ=len(Z)
  [M,L]=np.shape(lon)
  dz=np.gradient(Z)
  temp=np.zeros((NZ,M,L))
  salt=np.zeros((NZ,M,L))
  u=np.zeros((NZ,M,L-1))
  v=np.zeros((NZ,M-1,L))
  for k in range(NZ):
    if (k%10)==0:
      print('  Level ',str(k),' of ',str(NZ))
    #
    u2d=ofgm.ext_data_OGCM(nc,lonU,latU,'u',tin,lon,lat,k,Roa,interp_method)
    v2d=ofgm.ext_data_OGCM(nc,lonV,latV,'v',tin,lon,lat,k,Roa,interp_method)
    u[k]=ppt.rho2u_2d(u2d*cosa+v2d*sina)                                 
    v[k]=ppt.rho2v_2d(v2d*cosa-u2d*sina)
    temp[k]=ofgm.ext_data_OGCM(nc,lonT,latT,'temp',tin,lon,lat,k,Roa,interp_method)
    salt[k]=ofgm.ext_data_OGCM(nc,lonT,latT,'salt',tin,lon,lat,k,Roa,interp_method)
  #
  #
  # Close the OGCM file
  #
  nc.close()
  #
  # Get the CROCO vertical grid
  #
  print('  Vertical interpolations')
  if type(nc_clm)!=list:
    theta_s=nc_clm['theta_s'][:]
    theta_b=nc_clm['theta_b'][:]
    hc=nc_clm['hc'][:]
    N=nc_clm.dimensions['s_rho'].size
  #
  if type(nc_bry)!=list:
    theta_s=nc_bry['theta_s'][:]
    theta_b=nc_bry['theta_b'][:]
    hc=nc_bry['hc'][:]
    N=nc_bry.dimensions['s_rho'].size
  #
  zr=ppt.zlevs(h,zeta,theta_s,theta_b,hc,N,'r',vtransform)
  zu=ppt.rho2u_3d(zr)    
  zv=ppt.rho2v_3d(zr)
  zw=ppt.zlevs(h,zeta,theta_s,theta_b,hc,N,'w',vtransform)
  dzr=zw[1:,:,:]-zw[:-1,:,:]
  dzu=ppt.rho2u_3d(dzr)
  dzv=ppt.rho2v_3d(dzr)
  #
  # Add an extra bottom layer (-100000m) and an extra surface layer (+100m)
  # to prevent vertical extrapolations
  #
  pre_z = np.zeros(len(Z)+2)
  pre_z[0] = 100.
  pre_z[-1] = -100000.
  pre_z[1:-1]=Z
  Z=pre_z
  u=np.concatenate((np.reshape(u[0,:,:], (1,u.shape[1], u.shape[2])),u))
  u=np.concatenate((u,np.reshape(u[-1,:,:], (1,u.shape[1], u.shape[2]))))
  v=np.concatenate((np.reshape(v[0,:,:], (1,v.shape[1], v.shape[2])),v))
  v=np.concatenate((v,np.reshape(v[-1,:,:], (1,v.shape[1], v.shape[2]))))
  temp=np.concatenate((np.reshape(temp[0,:,:], (1, temp.shape[1], temp.shape[2])),temp))
  temp=np.concatenate((temp,np.reshape(temp[-1,:,:], (1,temp.shape[1], temp.shape[2]))))
  salt=np.concatenate((salt,np.reshape(salt[-1,:,:], (1,salt.shape[1], salt.shape[2]))))
  salt=np.concatenate((np.reshape(salt[0,:,:], (1,salt.shape[1], salt.shape[2])),salt))
  # 
  # Perform the vertical interpolations 
  #
  temp=ppt.ztosigma(np.flip(temp,0),zr,np.flipud(Z))
  salt=ppt.ztosigma(np.flip(salt,0),zr,np.flipud(Z))
  u=ppt.ztosigma(np.flip(u,0),zu,np.flipud(Z))
  v=ppt.ztosigma(np.flip(v,0),zv,np.flipud(Z))
  #
  # Correct the horizontal transport 
  # i.e. remove the interpolated tranport and add 
  #      the OGCM transport
  #
  if conserv==1:
    u=u-ppt.tridim((np.sum(u*dzu, axis=0)/np.sum(dzu, axis=0)),N)
    v=v-ppt.tridim((np.sum(v*dzv, axis=0)/np.sum(dzv, axis=0)),N)
    u=u+ppt.tridim(ubar,N)
    v=v+ppt.tridim(vbar,N)
  #end
  #
  # Barotropic velocities
  #
  ubar=(np.sum(u*dzu, axis=0)/np.sum(dzu, axis=0))
  vbar=(np.sum(v*dzv, axis=0)/np.sum(dzv, axis=0))
  #
  #  fill the files
  #
  if type(nc_clm)!=list: 
    nc_clm['zeta'][tout,:,:]=zeta
    if 'SSH' in list(nc_clm.variables):
      nc_clm['SSH'][tout,:,:]=zeta
    nc_clm['temp'][tout,:,:,:]=temp
    nc_clm['salt'][tout,:,:,:]=salt
    nc_clm['u'][tout,:,:,:]=u
    nc_clm['v'][tout,:,:,:]=v
    nc_clm['ubar'][tout,:,:]=ubar
    nc_clm['vbar'][tout,:,:]=vbar
  #end
  if type(nc_bry)!=list:
    for obcndx in range(0,4):
      if obcndx==0:
        nc_bry['zeta_south'][tout,:]=zeta[0,:]
        nc_bry['temp_south'][tout,:,:]=temp[:,0,:]
        nc_bry['salt_south'][tout,:,:]=salt[:,0,:]
        nc_bry['u_south'][tout,:,:]=u[:,0,:]
        nc_bry['v_south'][tout,:,:]=v[:,0,:]
        nc_bry['ubar_south'][tout,:]=ubar[0,:]
        nc_bry['vbar_south'][tout,:]=vbar[0,:]
      elif obcndx==1:
        nc_bry['zeta_east'][tout,:]=zeta[:,-1]
        nc_bry['temp_east'][tout,:,:]=temp[:,:,-1]
        nc_bry['salt_east'][tout,:,:]=salt[:,:,-1]
        nc_bry['u_east'][tout,:,:]=u[:,:,-1]
        nc_bry['v_east'][tout,:,:]=v[:,:,-1]
        nc_bry['ubar_east'][tout,:]=ubar[:,-1]
        nc_bry['vbar_east'][tout,:]=vbar[:,-1]
      elif obcndx==2:
        nc_bry['zeta_north'][tout,:]=zeta[-1,:]
        nc_bry['temp_north'][tout,:,:]=temp[:,-1,:]
        nc_bry['salt_north'][tout,:,:]=salt[:,-1,:]
        nc_bry['u_north'][tout,:,:]=u[:,-1,:]
        nc_bry['v_north'][tout,:,:]=v[:,-1,:]
        nc_bry['ubar_north'][tout,:]=ubar[-1,:]
        nc_bry['vbar_north'][tout,:]=vbar[-1,:]
      elif obcndx==3:
        nc_bry['zeta_west'][tout,:]=zeta[:,0]
        nc_bry['temp_west'][tout,:,:]=temp[:,:,0]
        nc_bry['salt_west'][tout,:,:]=salt[:,:,0]
        nc_bry['u_west'][tout,:,:]=u[:,:,0]
        nc_bry['v_west'][tout,:,:]=v[:,:,0]
        nc_bry['ubar_west'][tout,:]=ubar[:,0]
        nc_bry['vbar_west'][tout,:]=vbar[:,0]
      #end
    #end
  return nc_clm, nc_bry

######################################################################
#  Reproduce the old readattribute behavior when using loaddap library
#  But uses Built-in Support for OPeNDAP from Matlab >= 2012a 
#
#  Get the attribute of an OPENDAP dataset
######################################################################
def loaddap(url):
  x={}
  #
  try:
    ncid = open_url(url)
    
    nvar=len(ncid.keys())
    for ii in range(nvar):
      varname = list(ncid.keys())[ii]
      nattvar = len(ncid[varname].attributes)
      varname2=varname
      varname2=varname.replace('-','_2d')
      for jj in range(nattvar):
        attname = list(ncid[varname].attributes)[jj]
        attdict = ncid[varname].attributes
        
        if attname.find('missing_value')!=-1 or attname.find('FillValue')!=-1:
          attval = attdict[attname]
          #x[varname2].setncattr('missing_value',attval)
          #x[varname2].setncattr('ml__FillValue',attval)
          x['missing_value'] = attval
          x['ml__FillValue'] = attval
          #end
        if attname.find('scale_factor')!=-1:
          attval = attdict[attname]
          #x[varname2].setncattr('scale_factor',attval)
          x['scale_factor'] = attval
          #end
        if attname.find('add_offset')!=-1:
          attval = attdict[attname]
          #x[varname2].setncattr('add_offset',attval)
          x['add_offset'] = attval
          #end
        if attname.find('units')!=-1:
          attval = attdict[attname]
          #x[varname2].setncattr('units',attval)
          x['units'] = attval
          #end
        if attname.find('DODS_ML_Size')!=-1:
          attval = attdict[attname]
          #x[varname2].setncattr('DODS_ML_Size',attval)
          x['DODS_ML_Size'] = attval
  except:
    print('PB with loaddap')
    x=[]

  return x

def find_all(a_str, sub):
  start = 0
  while True:
    start = a_str.find(sub, start)
    if start == -1: return
    yield start
    start += len(sub)

######################################################################
#  Reproduce the old readdap behavior when using loaddap library
#  But uses Built-in Support for OPeNDAP from Matlab >= 2012a 
#  Retry (100 times) in case of network failure.
######################################################################
#
def readdap(url,varname,query):
  nmax=100
  data=[]
  ntry=0
  nargin = readdap.__code__.co_argcount
  #
  if nargin <2:
    print('not engough input argments')
  elif nargin <3 or len(query)==0:
    print('READDAP_New: Extract : '+ varname)
    while len(data)==0:
      if ntry>nmax:
        sys.exit('READDAP_New: repeated failures after '+str(nmax)+' queries')
      #end
      ntry=ntry+1
      ncid = open_url(url)
      if len(ncid)!=0:  
        raw_data = ncid[varname][:]
        # Handle PyDAP objects properly
        if hasattr(raw_data, 'array') and hasattr(raw_data.array, 'data'):
          # This is a GridType with array attribute  
          data = raw_data.array.data
        elif hasattr(raw_data, 'data'):
          # This is a BaseType or simple data
          data = raw_data.data
        elif hasattr(raw_data, 'mask'):
          # This is a masked array
          data = np.ma.filled(raw_data, np.nan)
        else:
          # Fallback
          data = np.array(raw_data)
      else:
        data=[]
        print('READDAP_New: did not work at '+str(ntry)+' try: lets try again.')
      #end
    #end
      
  else:
    print('READDAP_New: Extract : ', varname, query)
    ind1 = list(find_all(query, '['))
    ind2 = list(find_all(query, ']'))

    nb_dims=len(ind1)
    start2=[]
    count2=[]

    for ii in range(nb_dims):
      str_tmp=query[ind1[ii]+1:ind2[ii]]
      if str_tmp.find(':')==-1:
        start2.append(int(str_tmp))
      else:
        start2.append(int(str_tmp[:str_tmp.find(':')]))
        count2.append(int(str_tmp[str_tmp.find(':')+1:]))#-start2[ii]+1
      #end
    #end
    #start=np.flip(start2)
    #count=np.flip(count2)
      #[startcount]
      
  while len(data)==0:
    if ntry>nmax:
      sys.exit('READDAP_New: repeated failures after '+str(nmax)+' queries')
      #end
    ntry=ntry+1
    try:
      ncid = open_url(url)

      if len(start2)==1:
        if start2[0]==count2[0]:
          raw_data = ncid[varname][start2[0]]
        else:
          raw_data = ncid[varname][start2[0]:count2[0]]
      elif len(start2)==2:
        raw_data = ncid[varname].array[start2[0]:count2[0], start2[1]:count2[1]]
        #data=np.transpose(data,(1,0))
      elif len(start2)==3:
        if start2[0]==count2[0]:
          raw_data = ncid[varname].array[start2[0], start2[1]:count2[1], start2[2]:count2[2]]
        else:
          raw_data = ncid[varname].array[start2[0]:count2[0], start2[1]:count2[1], start2[2]:count2[2]]
        #data=np.transpose(data,(2, 1, 0))               
      elif len(start2)==4:
        raw_data = ncid[varname].array[start2[0]:count2[0], start2[1]:count2[1], start2[2]:count2[2], start2[3]:count2[3]]
        #data=np.transpose(data,(3, 2, 1, 0))
      
      # Handle PyDAP objects properly - BaseType from slicing or GridType from direct access
      if hasattr(raw_data, 'array') and hasattr(raw_data.array, 'data'):
        # This is a GridType with array attribute
        data = raw_data.array.data
      elif hasattr(raw_data, 'data'):
        # This is a BaseType from slicing or simple data
        data = raw_data.data
      elif hasattr(raw_data, 'mask'):
        # This is a masked array
        data = np.ma.filled(raw_data, np.nan)
      else:
        # Fallback - try to convert to array
        data = np.array(raw_data)
        
    except:
      data=[]
      print('READDAP_New: did not work at '+str(ntry)+' try: lets try again.')
  #
  return data

######################################################################
#
#  var=getdap(path,fname,vname,trange,krange,jrange,...
#             i1min,i1max,i2min,i2max,i3min,i3max)
#
#  Download a data subsets from a OPENDAP server.
#
#  Take care of the greenwitch meridian
#  (i.e. get 3 subgrids defined by i1min,i1max,i2min,i2max,i3min,i3max
#  and concatenate them).

######################################################################
def getdap(path,fname,vname,trange,krange,jrange,i1min,i1max,i2min,i2max,i3min,i3max):
#
  url=path+fname
  #
  var=[]
  #
  if i2min:
    irange='['+str(i2min)+':'+str(i2max)+']'
    var=readdap(url,vname,trange+krange+jrange+irange)
  #end
  #
  if i1min:
    irange='['+str(i1min)+':'+str(i1max)+']'
    var0=readdap(url,vname,trange+krange+jrange+irange)
    if var:
      var=np.concatenate((var0,var), axis=-1)
    else:
      var = var0
  #end
  #
  if i3min:
    irange='[',str(i3min)+':'+str(i3max)+']'
    var0=readdap(url,vname,trange+krange+jrange+irange)
    if var:
      var=np.concatenate((var,var0), axis=-1)
    else:
      var=var0
  #end
  #
  return var


######################################################################
#  Give the GFS url for a given date.
######################################################################
def get_GFS_fname(time,gfs_run_time,gfstype):
  #
  # set URL
  #
  url='https://nomads.ncep.noaa.gov'
  #url='https://nomads-cprk.ncep.noaa.gov'

  # set file types
  if gfstype==0:
    gfsname ='fnl'
    gfsname1='fnlflx'     # 1/2 GDAS data
  else:
    gfsname ='gfs'
    gfsname1='gfs_0p25'   # 1/4 deg res GFS data
    #gfsname1='gfs_0p50'  # 1/2 deg res GFS data
  #end
  #
  # Get the date
  #
  str_time = time.strftime('%Y%m%d')
  stry=str_time[:4]
  strm=str_time[4:6]
  strd=str_time[6:]
  #end
  if gfs_run_time < 10:
    strh='_0'
  else:
    strh='_'
  #end
  #
  # Get the grid
  #
  if gfstype==0:
    gfsdir =url+'/dods/'+gfsname+'/'+gfsname+stry+strm+strd+'/'
    fname=gfsdir+gfsname1+strh+str(gfs_run_time)+'z'
  else:
    gfsdir =url+'/dods/'+gfsname1+'/'+gfsname+stry+strm+strd+'/'
    fname=gfsdir+gfsname1+strh+str(gfs_run_time)+'z'
  #end
  return fname

         
######################################################################
# Get the indices for a GFS subgrid 
######################################################################
def get_GFS_grid(fname,lonmin,lonmax,latmin,latmax):
  dl=1
  lonmin=lonmin-dl
  lonmax=lonmax+dl
  latmin=latmin-dl
  latmax=latmax+dl
  #
  # Get the global grid
  #
  nc=open_url(fname)
  lon = nc['lon'][:].data
  lat = nc['lat'][:].data
  #
  # Get a subgrid
  #
  # 1 Longitude: take care of greenwitch
  #
  i1=np.where((lon-360>=lonmin) & (lon-360<=lonmax))
  i2=np.where((lon>=lonmin) & (lon<=lonmax))
  i3=np.where((lon+360>=lonmin) & (lon+360<=lonmax))
  #
  lon=np.concatenate((lon[i1]-360,lon[i2],lon[i3]+360), axis=0)
  #
  if i1[0].shape[0]!=0:
    i1min=np.min(i1)-1
    i1max=np.max(i1)
  else:
    i1min=[]
    i1max=[]
  #end

  if i2[0].shape[0]!=0:
    i2min=np.min(i2)-1
    i2max=np.max(i2)
  else:
    i2min=[]
    i2max=[]
  #end

  if i3[0].shape[0]!=0:
    i3min=np.min(i3)-1
    i3max=np.max(i3)
  else:
    i3min=[]
    i3max=[]
  #end
  #
  # 2 Latitude
  #
  j=np.where((lat>=latmin) & (lat<=latmax))
  lat=lat[j]
  jmin=np.min(j)-1
  jmax=np.max(j)
  jrange='['+str(jmin)+':'+str(jmax)+']'
  #
  return [i1min,i1max,i2min,i2max,i3min,i3max,jrange,lon,lat]


         
######################################################################
# Download one full subset of GFS for CROCO bulk for 1 time step
# Put them in the CROCO units
######################################################################
def get_GFS(fname,mask,tndx,jrange,i1min,i1max,i2min,i2max,i3min,i3max,missvalue):
  if type(tndx)==int:
    trange='['+str(tndx)+':'+str(tndx)+']'
  else:
    trange='['+str(min(tndx))+':'+str(max(tndx))+']'
  #
  # Get GFS variables for 1 time step
  #
  print(' ')
  print('====================================================')

  t=readdap(fname,'time',trange)
  print('TRANGE='+str(trange))
  print('GFS raw time='+str(t))
  #t=t+365 # put it in "matlab" time
  print('GFS: ', datetime(1,1,1)+t*timedelta(days=1)-timedelta(days=2))
  print('====================================================')

  #print('u...')
  u=mask*getdap('',fname,'ugrd10m',trange,'',jrange,i1min,i1max,i2min,i2max,i3min,i3max)
  u[abs(u)>=missvalue]=float('nan')

  #print('v...')
  v=mask*getdap('',fname,'vgrd10m',trange,'',jrange,i1min,i1max,i2min,i2max,i3min,i3max)
  v[abs(v)>=missvalue]=float('nan')

  #print('ty...')
  ty=mask*getdap('',fname,'vflxsfc',trange,'',jrange,i1min,i1max,i2min,i2max,i3min,i3max)
  ty[abs(ty)>=missvalue]=float('nan')

  #print('tx...')
  tx=mask*getdap('',fname,'uflxsfc',trange,'',jrange,i1min,i1max,i2min,i2max,i3min,i3max)
  tx[abs(tx)>=missvalue]=float('nan')


  #print('skt...')
  #skt=mask.*getdap('',fname,'tmpsfc',trange,'',jrange,...
  #                  i1min,i1max,i2min,i2max,i3min,i3max)
  #skt(abs(skt)>=missvalue)=float('nan')

  #print('tair...')
  tair=mask*getdap('',fname,'tmp2m',trange,'',jrange,i1min,i1max,i2min,i2max,i3min,i3max)
  tair[abs(tair)>=missvalue]=float('nan')

  #print('rhum...')
  rhum=mask*getdap('',fname,'rh2m',trange,'',jrange,i1min,i1max,i2min,i2max,i3min,i3max)
  rhum[abs(rhum)>=missvalue]=float('nan')

  #print('prate...')
  prate=mask*getdap('',fname,'pratesfc',trange,'',jrange,i1min,i1max,i2min,i2max,i3min,i3max)
  prate[abs(prate)>=missvalue]=float('nan')

  #print('down radlw...')
  dradlw=mask*getdap('',fname,'dlwrfsfc',trange,'',jrange,i1min,i1max,i2min,i2max,i3min,i3max)
  dradlw[abs(dradlw)>=missvalue]=float('nan')

  #print('up radlw')
  uradlw=mask*getdap('',fname,'ulwrfsfc',trange,'',jrange,i1min,i1max,i2min,i2max,i3min,i3max)
  uradlw[abs(uradlw)>=missvalue]=float('nan')

  #print('down radsw')
  dradsw=mask*getdap('',fname,'dswrfsfc',trange,'',jrange,i1min,i1max,i2min,i2max,i3min,i3max)
  dradsw[abs(dradsw)>=missvalue]=float('nan')

  #print('up radsw')
  uradsw=mask*getdap('',fname,'uswrfsfc',trange,'',jrange,i1min,i1max,i2min,i2max,i3min,i3max)
  uradsw[abs(uradsw)>=missvalue]=float('nan')
  #
  # Transform the variables
  #
  #
  # 1: Air temperature: Convert from Kelvin to Celsius
  #
  tair=tair-273.15
  #
  # 2: Relative humidity: Convert from # to fraction
  #
  rhum=rhum/100
  #
  # 3: Precipitation rate: Convert from [kg/m^2/s] to cm/day
  #
  prate=prate*0.1*(24*60*60.0)
  prate[abs(prate)<1.e-4]=0
  #
  # 4: Net shortwave flux: [W/m^2]
  #    CROCO convention: positive downward: same as GFS
  # ?? albedo ??
  #
  radsw=dradsw - uradsw
  radsw[radsw<1.e-10]=0
  #
  # 5: Net outgoing Longwave flux:  [W/m^2]
  #    CROCO convention: positive upward (opposite to nswrs)
  #    GFS convention: positive downward --> * (-1)
  #    input: downward longwave rad. and
  #    skin temperature.
  #
  #skt=skt-273.15 
  #radlw=-lwhf(skt,radlw) 
  radlw = uradlw - dradlw
  radlw_in=dradlw

  #
  # 6: Wind speed
  #
  wspd=(u**2+v**2)**0.5
  # 7:  Wind vectors
  uwnd=u # rho point
  vwnd=v # rho point
  # #
  # # 7: Compute the stress following large and pond
  # #
  # [Cd,uu]=cdnlp(wspd,10.)
  # rhoa=air_dens(tair,rhum*100)
  # tx=Cd.*rhoa.*u.*wspd
  # ty=Cd.*rhoa.*v.*wspd
  #
  return [t,tx,ty,tair,rhum,prate,wspd,uwnd,vwnd,radlw,radlw_in,radsw]



######################################################################
#  function write_GFS(fname,Yorig,lon,lat,mask,time,tx,ty,...
#                     tair,rhum,prate,wspd,radlw,radsw)
#  Write into a GFS file
######################################################################
def write_GFS(fname,Yorig,lon,lat,mask,time,tx,ty,tair,rhum,prate,wspd,uwnd,vwnd,radlw,radlw_in,radsw):
  #
  print(['Create ',fname])
  nc=Dataset(fname,'w')
  #
  nc.createDimension('lon',len(lon))
  nc.createDimension('lat',len(lat))
  #
  # nc('latu') = len(lat)
  # nc('latv') = len(lat)-1
  # nc('lonu') = len(lon)-1
  # nc('lonv') = len(lon)
  #
  nc.createDimension('time', len(time))
  #
  nc_lon = nc.createVariable('lon','f',('lon'))
  nc_lon.long_name = 'longitude of RHO-points'
  nc_lon.units = 'degree_east'
  # 
  nc_lat = nc.createVariable('lat','f',('lat'))
  nc_lat.long_name = 'latitude of RHO-points'
  nc_lat.units = 'degree_north'
  #
  nc_time = nc.createVariable('time','f',('time'))
  nc_time.long_name = 'Time'
  nc_time.units = 'days since 1-Jan-'+str(Yorig)+' 00:00:0.0'
  #
  nc_mask = nc.createVariable('mask','f',('lat','lon'))
  nc_tx = nc.createVariable('tx','f',('time','lat','lon'))
  nc_ty = nc.createVariable('ty','f',('time','lat','lon'))
  nc_tair = nc.createVariable('tair','f',('time','lat','lon'))
  nc_rhum = nc.createVariable('rhum','f',('time','lat','lon'))
  nc_prate = nc.createVariable('prate','f',('time','lat','lon'))
  nc_wspd = nc.createVariable('wspd','f',('time','lat','lon'))
  nc_radlw = nc.createVariable('radlw','f',('time','lat','lon'))
  nc_radsw = nc.createVariable('radsw','f',('time','lat','lon'))
  nc_radlw_in = nc.createVariable('radlw_in','f',('time','lat','lon'))
  nc_uwnd = nc.createVariable('uwnd','f',('time','lat','lon'))
  nc_vwnd = nc.createVariable('vwnd','f',('time','lat','lon'))
  #
  #endef(nc)
  #
  
  nc_lon[:]=lon.data
  nc_lat[:]=lat.data
  nc_time[:]=time.data
  nc_mask[:]=mask[0].data
  nc_tx[:]=tx.data
  nc_ty[:]=ty.data
  nc_tair[:]=tair.data
  nc_rhum[:]=rhum.data
  nc_prate[:]=prate.data
  nc_wspd[:]=wspd.data
  nc_uwnd[:]=uwnd.data
  nc_vwnd[:]=vwnd.data
  nc_radlw[:]=radlw.data
  nc_radlw_in[:]=radlw_in.data
  nc_radsw[:]=radsw.data
  #
  nc.close()


######################################################################
#  download_GFS_AWS(today,lonmin,lonmax,latmin,latmax,FRCST_dir,Yorig,it)
#  Extract a subgrid from GFS AWS S3 to get a CROCO forcing
#  Uses AWS Open Data instead of deprecated OpenDAP
#  Downloads GRIB2 files from S3 and processes with xarray/cfgrib
######################################################################
def download_GFS(today,lonmin,lonmax,latmin,latmax,FRCST_dir,Yorig,it,keep_grib=True):
  """
  Download GFS data from AWS S3 Open Data bucket
  Data is in GRIB2 format and accessed via boto3
  """
  #
  # Put the date in 'Yorig' time
  #
  rundate_str=datetime.today()
  rundate=rundate_str-datetime(Yorig,1,1)
  rundate=rundate.days

  # Target forcing start aligned with CROCO ini convention:
  # ini starts at (rundate - hdays) at 00:00 local day boundary.
  ini_start = datetime(today.year, today.month, today.day) - timedelta(days=hdays)
  #
  # GFS output name
  #
  gfs_name=FRCST_dir+'GFS_'+str(rundate)+'.nc'
  
  #
  # Setup AWS S3 client (no credentials needed for public data)
  # Note: verify=False may be needed in some institutional networks with SSL issues
  #
  import urllib3
  urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
  
  # Disable SSL verification for institutional networks if needed
  s3 = boto3.client('s3', config=Config(signature_version=UNSIGNED), verify=False)
  GFS_BUCKET_NAME = "noaa-gfs-bdp-pds"
  
  print(' ')
  print('Get GFS data from AWS S3 for ', str(today))
  print('Minimum Longitude: ',str(lonmin))
  print('Maximum Longitude: ',str(lonmax))
  print('Minimum Latitude:  ',str(latmin))
  print('Maximum Latitude:  ',str(latmax))
  print(' ')
  #
  # Create directory if needed
  #
  print('Making output data directory '+FRCST_dir)
  frcst_ext = os.path.exists(FRCST_dir)
  if not frcst_ext:
    os.makedirs(FRCST_dir)
  
  #
  # Find latest available GFS forecast
  # GFS forecasts are available 4 times per day: 00, 06, 12, 18 UTC
  #
  gfs_cycles = ['18', '12', '06', '00']
  forecast_date = today
  found_forecast = False
  
  # Try to find most recent forecast (up to 3 days back)
  # Only accept cycles that already have f000 available for download.
  def _object_exists(bucket, key):
    try:
      s3.head_object(Bucket=bucket, Key=key)
      return True
    except Exception:
      return False

  for days_back in range(6):
    check_date = forecast_date - timedelta(days=days_back)
    date_str = check_date.strftime('%Y%m%d')
    
    for cycle in gfs_cycles:
      folder = f"gfs.{date_str}/{cycle}/atmos/"
      cycle_dt = datetime(check_date.year, check_date.month, check_date.day, int(cycle))
      try:
        response = s3.list_objects_v2(Bucket=GFS_BUCKET_NAME, Prefix=folder, MaxKeys=1)
        if 'Contents' in response:
          # Do not select a cycle newer than ini_start, otherwise hindcast start drifts late.
          if cycle_dt > ini_start:
            continue
          f000_025 = f"{folder}gfs.t{cycle}z.pgrb2.0p25.f000"
          f000_050 = f"{folder}gfs.t{cycle}z.pgrb2.0p50.f000"
          if _object_exists(GFS_BUCKET_NAME, f000_025) or _object_exists(GFS_BUCKET_NAME, f000_050):
            found_forecast = True
            forecast_cycle = cycle
            forecast_date_str = date_str
            print(f'  Found GFS forecast with available f000: {date_str} cycle {cycle}Z (<= ini_start {ini_start})')
            break
          else:
            print(f'  Skipping {date_str} cycle {cycle}Z (folder exists but f000 not yet available)')
      except Exception as e:
        continue
    
    if found_forecast:
      break
  
  if not found_forecast:
    sys.exit('  Could not find any recent GFS forecast on AWS S3')
  
  #
  # Download GFS GRIB2 files for the forecast period
  # GFS provides forecasts at 3-hour intervals
  #
  temp_dir = tempfile.mkdtemp(prefix='gfs_tmp_', dir=FRCST_dir)
  print(f'  Using temporary directory: {temp_dir}')
  
  # Calculate number of forecast hours needed
  # For hindcast: use analysis (f000) and early forecast hours
  # For forecast: use forecast hours
  total_days = hdays + fdays + 1
  forecast_hours = list(range(0, total_days * 24, it * 3))  # it*3 hour intervals
  
  print(f'  Downloading {len(forecast_hours)} forecast hours: {forecast_hours[0]} to {forecast_hours[-1]}')
  
  # Download GRIB2 files
  # - pgrb2  : core atmospheric variables (10u/10v/2t/2r/prate/...)
  # - pgrb2b : additional fields including radiation in many GFS streams
  grib_files = []
  grib_files_b = []
  for fhour in forecast_hours:
    fhour_str = f"{fhour:03d}"
    object_key = f"gfs.{forecast_date_str}/{forecast_cycle}/atmos/gfs.t{forecast_cycle}z.pgrb2.0p25.f{fhour_str}"
    local_file = os.path.join(temp_dir, f"gfs_f{fhour_str}.grib2")
    
    try:
      print(f'    Downloading forecast hour {fhour_str}...')
      s3.download_file(GFS_BUCKET_NAME, object_key, local_file)
      grib_files.append(local_file)

      # Download companion pgrb2b file (radiation often stored here)
      object_key_b = f"gfs.{forecast_date_str}/{forecast_cycle}/atmos/gfs.t{forecast_cycle}z.pgrb2b.0p25.f{fhour_str}"
      local_file_b = os.path.join(temp_dir, f"gfsb_f{fhour_str}.grib2")
      try:
        s3.download_file(GFS_BUCKET_NAME, object_key_b, local_file_b)
        grib_files_b.append(local_file_b)
      except Exception:
        grib_files_b.append(None)
    except Exception as e:
      print(f'    Warning: Could not download f{fhour_str}: {e}')
      # Try alternative resolution if 0.25 degree not available
      try:
        object_key_alt = f"gfs.{forecast_date_str}/{forecast_cycle}/atmos/gfs.t{forecast_cycle}z.pgrb2.0p50.f{fhour_str}"
        s3.download_file(GFS_BUCKET_NAME, object_key_alt, local_file)
        grib_files.append(local_file)
        object_key_b_alt = f"gfs.{forecast_date_str}/{forecast_cycle}/atmos/gfs.t{forecast_cycle}z.pgrb2b.0p50.f{fhour_str}"
        local_file_b = os.path.join(temp_dir, f"gfsb_f{fhour_str}.grib2")
        try:
          s3.download_file(GFS_BUCKET_NAME, object_key_b_alt, local_file_b)
          grib_files_b.append(local_file_b)
        except Exception:
          grib_files_b.append(None)
        print(f'    Downloaded f{fhour_str} at 0.5 degree resolution')
      except Exception as e2:
        print(f'    Error: Could not download f{fhour_str} at either resolution')
  
  if len(grib_files) == 0:
    sys.exit('  Error: No GRIB2 files could be downloaded')
  
  n_b = len([f for f in grib_files_b if f is not None])
  print(f'  Successfully downloaded {len(grib_files)} pgrb2 files and {n_b} pgrb2b files')
  
  #
  # Process GRIB2 files and extract variables
  #
  print(' ... PROCESSING GRIB2 FILES ... ')
  
  # Lists to store data from all time steps
  all_data = []
  all_times = []
  lon = None
  lat = None
  report_keys = ['u10', 'v10', 't2m', 'r2', 'uflx', 'vflx', 'prate', 'dswrf', 'uswrf', 'dlwrf', 'ulwrf']
  extract_report = {
    key: {'ok': 0, 'miss': 0, 'src_main': 0, 'src_b': 0, 'step_avg': 0, 'step_instant': 0, 'step_none': 0}
    for key in report_keys
  }

  def _subset_roi(ds_in, lonmin_req, lonmax_req, latmin_req, latmax_req):
    lat_lo = min(latmin_req, latmax_req) - 1
    lat_hi = max(latmin_req, latmax_req) + 1
    lat_vals = ds_in['latitude'].values
    if lat_vals[0] > lat_vals[-1]:
      ds_lat = ds_in.sel(latitude=slice(lat_hi, lat_lo))
    else:
      ds_lat = ds_in.sel(latitude=slice(lat_lo, lat_hi))

    lon_vals = ds_lat['longitude'].values
    if np.nanmax(lon_vals) > 180:
      lon_lo = lonmin_req % 360
      lon_hi = lonmax_req % 360
    else:
      lon_lo = ((lonmin_req + 180) % 360) - 180
      lon_hi = ((lonmax_req + 180) % 360) - 180

    if lon_lo <= lon_hi:
      ds_lon = ds_lat.sel(longitude=slice(lon_lo - 1, lon_hi + 1))
    else:
      ds_w = ds_lat.sel(longitude=slice(lon_lo - 1, np.nanmax(lon_vals)))
      ds_e = ds_lat.sel(longitude=slice(np.nanmin(lon_vals), lon_hi + 1))
      ds_lon = xr.concat([ds_w, ds_e], dim='longitude')

    return ds_lon

  def _extract_grib_field(grib_path, short_name, type_of_level, level=None, step_types=None):
    if grib_path is None:
      raise ValueError(f'No GRIB path available for {short_name}')

    if step_types is None:
      step_types = [None]

    last_error = None
    for step_type in step_types:
      filter_keys = {'shortName': short_name, 'typeOfLevel': type_of_level}
      if level is not None:
        filter_keys['level'] = level
      if step_type is not None:
        filter_keys['stepType'] = step_type

      try:
        ds = xr.open_dataset(
          grib_path,
          engine='cfgrib',
          backend_kwargs={'filter_by_keys': filter_keys, 'indexpath': ''}
        )
        ds_sub = _subset_roi(ds, lonmin, lonmax, latmin, latmax)

        if ds_sub.sizes.get('latitude', 0) == 0 or ds_sub.sizes.get('longitude', 0) == 0:
          ds.close()
          raise ValueError(f'Empty subset for {short_name}')

        data_vars = list(ds_sub.data_vars)
        if len(data_vars) == 0:
          ds.close()
          raise ValueError(f'No data variables found for {short_name}')

        vname = data_vars[0]
        arr = np.asarray(ds_sub[vname].values).squeeze()

        if 'valid_time' in ds_sub.coords:
          tval = np.asarray(ds_sub['valid_time'].values).squeeze()
        elif 'time' in ds_sub.coords:
          tval = np.asarray(ds_sub['time'].values).squeeze()
        else:
          tval = None

        lon_out = np.asarray(ds_sub['longitude'].values)
        lat_out = np.asarray(ds_sub['latitude'].values)
        ds.close()

        return arr, tval, lon_out, lat_out, step_type
      except Exception as e:
        last_error = e
        continue

    raise ValueError(f'Could not read {short_name} from {os.path.basename(grib_path)}: {last_error}')

  def _align_to_target_grid(arr_in, lon_src, lat_src, lon_tgt, lat_tgt):
    arr = np.asarray(arr_in)
    if arr.shape == (len(lat_tgt), len(lon_tgt)):
      return arr

    # Sort source coordinates to satisfy interpolation requirements
    lon_src = np.asarray(lon_src)
    lat_src = np.asarray(lat_src)
    lon_order = np.argsort(lon_src)
    lat_order = np.argsort(lat_src)
    lon_src_sorted = lon_src[lon_order]
    lat_src_sorted = lat_src[lat_order]
    arr_sorted = arr[np.ix_(lat_order, lon_order)]

    da = xr.DataArray(
      arr_sorted,
      coords={'latitude': lat_src_sorted, 'longitude': lon_src_sorted},
      dims=('latitude', 'longitude')
    )

    da_interp = da.interp(
      latitude=np.asarray(lat_tgt),
      longitude=np.asarray(lon_tgt),
      method='linear'
    )
    return np.asarray(da_interp.values)
  
  for idx, grib_file in enumerate(grib_files):
    try:
      print(f'    Processing {os.path.basename(grib_file)}...')
      grib_file_b = grib_files_b[idx] if idx < len(grib_files_b) else None

      data_dict = {}
      time_val = None

      # 10 m winds
      try:
        arr, tval, lon_tmp, lat_tmp, step_used = _extract_grib_field(grib_file, '10u', 'heightAboveGround', 10)
        if lon is None:
          lon = np.where(lon_tmp > 180, lon_tmp - 360, lon_tmp)
          lat = lat_tmp
          data_dict['u10'] = arr
        else:
          lon_src = np.where(lon_tmp > 180, lon_tmp - 360, lon_tmp)
          data_dict['u10'] = _align_to_target_grid(arr, lon_src, lat_tmp, lon, lat)
        extract_report['u10']['ok'] += 1
        extract_report['u10']['src_main'] += 1
        if step_used == 'avg':
          extract_report['u10']['step_avg'] += 1
        elif step_used == 'instant':
          extract_report['u10']['step_instant'] += 1
        else:
          extract_report['u10']['step_none'] += 1
        if time_val is None:
          time_val = tval
      except Exception as e:
        extract_report['u10']['miss'] += 1
        print(f'      Warning: could not read 10u: {e}')

      try:
        arr, _, lon_tmp, lat_tmp, step_used = _extract_grib_field(grib_file, '10v', 'heightAboveGround', 10)
        if lon is None:
          lon = np.where(lon_tmp > 180, lon_tmp - 360, lon_tmp)
          lat = lat_tmp
          data_dict['v10'] = arr
        else:
          lon_src = np.where(lon_tmp > 180, lon_tmp - 360, lon_tmp)
          data_dict['v10'] = _align_to_target_grid(arr, lon_src, lat_tmp, lon, lat)
        extract_report['v10']['ok'] += 1
        extract_report['v10']['src_main'] += 1
        if step_used == 'avg':
          extract_report['v10']['step_avg'] += 1
        elif step_used == 'instant':
          extract_report['v10']['step_instant'] += 1
        else:
          extract_report['v10']['step_none'] += 1
      except Exception as e:
        extract_report['v10']['miss'] += 1
        print(f'      Warning: could not read 10v: {e}')

      try:
        arr, _, lon_tmp, lat_tmp, step_used = _extract_grib_field(grib_file, '2t', 'heightAboveGround', 2)
        if lon is None:
          lon = np.where(lon_tmp > 180, lon_tmp - 360, lon_tmp)
          lat = lat_tmp
          data_dict['t2m'] = arr
        else:
          lon_src = np.where(lon_tmp > 180, lon_tmp - 360, lon_tmp)
          data_dict['t2m'] = _align_to_target_grid(arr, lon_src, lat_tmp, lon, lat)
        extract_report['t2m']['ok'] += 1
        extract_report['t2m']['src_main'] += 1
        if step_used == 'avg':
          extract_report['t2m']['step_avg'] += 1
        elif step_used == 'instant':
          extract_report['t2m']['step_instant'] += 1
        else:
          extract_report['t2m']['step_none'] += 1
      except Exception as e:
        extract_report['t2m']['miss'] += 1
        print(f'      Warning: could not read 2t: {e}')

      try:
        arr, _, lon_tmp, lat_tmp, step_used = _extract_grib_field(grib_file, '2r', 'heightAboveGround', 2)
        if lon is None:
          lon = np.where(lon_tmp > 180, lon_tmp - 360, lon_tmp)
          lat = lat_tmp
          data_dict['r2'] = arr
        else:
          lon_src = np.where(lon_tmp > 180, lon_tmp - 360, lon_tmp)
          data_dict['r2'] = _align_to_target_grid(arr, lon_src, lat_tmp, lon, lat)
        extract_report['r2']['ok'] += 1
        extract_report['r2']['src_main'] += 1
        if step_used == 'avg':
          extract_report['r2']['step_avg'] += 1
        elif step_used == 'instant':
          extract_report['r2']['step_instant'] += 1
        else:
          extract_report['r2']['step_none'] += 1
      except Exception as e:
        extract_report['r2']['miss'] += 1
        print(f'      Warning: could not read 2r: {e}')

      # Surface momentum fluxes and other surface fields
      # Use alias lists because GFS short names vary by product/cycle.
      for short_names, key_name, step_types, src_kind in [
          (['uflx', 'utaua'], 'uflx', ['avg', 'instant', None], 'main'),
          (['vflx', 'vtaua'], 'vflx', ['avg', 'instant', None], 'main'),
          (['prate'], 'prate', ['avg', 'instant', None], 'main'),
          (['dswrf', 'sdswrf'], 'dswrf', ['avg', 'instant', None], 'main_then_b'),
          (['uswrf', 'suswrf'], 'uswrf', ['avg', 'instant', None], 'main_then_b'),
          (['dlwrf', 'sdlwrf'], 'dlwrf', ['avg', 'instant', None], 'main_then_b'),
          (['ulwrf', 'sulwrf'], 'ulwrf', ['avg', 'instant', None], 'main_then_b')]:
        if src_kind == 'main':
          source_candidates = [grib_file]
        elif src_kind == 'b_then_main':
          source_candidates = [f for f in [grib_file_b, grib_file] if f is not None]
        else:  # main_then_b
          source_candidates = [f for f in [grib_file, grib_file_b] if f is not None]
        got_var = False
        for src_file in source_candidates:
          for short_name in short_names:
            try:
              arr, _, lon_tmp, lat_tmp, step_used = _extract_grib_field(src_file, short_name, 'surface', step_types=step_types)
              if lon is None:
                lon = np.where(lon_tmp > 180, lon_tmp - 360, lon_tmp)
                lat = lat_tmp
                data_dict[key_name] = arr
              else:
                lon_src = np.where(lon_tmp > 180, lon_tmp - 360, lon_tmp)
                data_dict[key_name] = _align_to_target_grid(arr, lon_src, lat_tmp, lon, lat)
              extract_report[key_name]['ok'] += 1
              if src_file == grib_file_b:
                extract_report[key_name]['src_b'] += 1
              else:
                extract_report[key_name]['src_main'] += 1
              if step_used == 'avg':
                extract_report[key_name]['step_avg'] += 1
              elif step_used == 'instant':
                extract_report[key_name]['step_instant'] += 1
              else:
                extract_report[key_name]['step_none'] += 1
              got_var = True
              break
            except Exception:
              continue
          if got_var:
            break
        if not got_var:
          extract_report[key_name]['miss'] += 1

      if time_val is None:
        print('      Warning: no valid time found, skipping timestep')
        continue

      # Keep timestep only if critical fields exist
      if not all(k in data_dict for k in ['u10', 'v10', 't2m', 'r2']):
        print('      Warning: missing critical atmospheric variables, skipping timestep')
        continue

      all_data.append(data_dict)
      all_times.append(time_val)
      
    except Exception as e:
      print(f'    Error processing {grib_file}: {e}')
      continue
  
  print(f'  Processed {len(all_data)} time steps')
  print('  Extraction report [var: ok/miss | src(main,b) | step(avg,instant,none)]')
  for key in report_keys:
    rep = extract_report[key]
    print(f"    {key:6s}: {rep['ok']:2d}/{rep['miss']:2d} | src({rep['src_main']:2d},{rep['src_b']:2d}) | step({rep['step_avg']:2d},{rep['step_instant']:2d},{rep['step_none']:2d})")

  if len(all_data) == 0:
    sys.exit('  Error: no valid GFS timesteps after GRIB parsing/subsetting')

  if lon is None or lat is None:
    sys.exit('  Error: could not extract longitude/latitude grid from GRIB files')

  def _time_to_datetime(time_in):
    if time_in is None:
      raise ValueError('Missing time value for timestep')
    tarr = np.asarray(time_in)
    if tarr.size == 0:
      raise ValueError('Empty time value for timestep')
    tscalar = tarr.reshape(-1)[0]
    return pd.to_datetime(tscalar).to_pydatetime()

  # Align atmospheric forcing start time to CROCO ini day boundary.
  # ini_start = today 00:00 - hdays
  run_day_start = ini_start
  keep_indices = []
  for tidx, time_val in enumerate(all_times):
    try:
      if _time_to_datetime(time_val) >= run_day_start:
        keep_indices.append(tidx)
    except Exception:
      continue

  if len(keep_indices) > 0:
    if keep_indices[0] > 0:
      first_kept_time = _time_to_datetime(all_times[keep_indices[0]])
      print(f'  Trimming {keep_indices[0]} timestep(s) before {run_day_start} to match CROCO ini start')
      print(f'  First kept forcing time: {first_kept_time}')
    all_data = [all_data[k] for k in keep_indices]
    all_times = [all_times[k] for k in keep_indices]
  else:
    print('  Warning: no GFS timesteps on/after run-date midnight; keeping original time series')

  if len(all_times) > 0:
    first_final_time = _time_to_datetime(all_times[0])
    last_final_time = _time_to_datetime(all_times[-1])
    print(f'  Final forcing coverage: {first_final_time} -> {last_final_time} ({len(all_times)} records)')
  
  #
  # Convert extracted data to CROCO format and write output file
  #
  print(' ... CONVERTING TO CROCO FORMAT ... ')
  
  # Create 2D mesh grid for lon/lat
  lon_2d, lat_2d = np.meshgrid(lon, lat)
  
  # Initialize output arrays
  n_times = len(all_data)
  M, L = lon_2d.shape
  
  tx = np.full((n_times, M, L), np.nan)
  ty = np.full((n_times, M, L), np.nan)
  tair = np.full((n_times, M, L), np.nan)
  rhum = np.full((n_times, M, L), np.nan)
  prate = np.full((n_times, M, L), np.nan)
  wspd = np.full((n_times, M, L), np.nan)
  uwnd = np.full((n_times, M, L), np.nan)
  vwnd = np.full((n_times, M, L), np.nan)
  radlw = np.full((n_times, M, L), np.nan)
  radlw_in = np.full((n_times, M, L), np.nan)
  radsw = np.full((n_times, M, L), np.nan)
  gfstime = np.zeros(n_times)
  
  # Create land mask (assume all ocean for now, can be refined)
  mask = np.ones((1, M, L))

  def _fill_missing_timesteps(field_arr):
    nt = field_arr.shape[0]
    if nt == 0:
      return field_arr

    # Forward fill from previous valid timestep
    for tidx in range(1, nt):
      if np.all(~np.isfinite(field_arr[tidx])):
        field_arr[tidx] = field_arr[tidx-1]

    # Backward fill leading invalid timesteps
    first_valid = None
    for tidx in range(nt):
      if np.any(np.isfinite(field_arr[tidx])):
        first_valid = tidx
        break

    if first_valid is None:
      # No valid timestep exists for this field: keep original values (NaN)
      return field_arr

    for tidx in range(0, first_valid):
      field_arr[tidx] = field_arr[first_valid]
    return field_arr
  
  # Process each time step
  for i, (data, time_val) in enumerate(zip(all_data, all_times)):
    # Convert time to CROCO format (days since Yorig-01-01)
    time_dt = _time_to_datetime(time_val)
    gfstime[i] = (time_dt - datetime(Yorig, 1, 1)).total_seconds() / 86400.0
    
    # Extract and convert variables with proper fallbacks
    # Wind components
    if 'u10' in data and data['u10'] is not None:
      uwnd[i] = data['u10']
    elif i > 0:
      uwnd[i] = uwnd[i-1]  # Use previous value if missing
      
    if 'v10' in data and data['v10'] is not None:
      vwnd[i] = data['v10']
    elif i > 0:
      vwnd[i] = vwnd[i-1]
    
    # Temperature: Kelvin to Celsius
    if 't2m' in data and data['t2m'] is not None:
      tair[i] = data['t2m'] - 273.15
    elif i > 0:
      tair[i] = tair[i-1]
    
    # Relative humidity: Percent to fraction
    if 'r2' in data and data['r2'] is not None:
      rhum[i] = data['r2'] / 100.0
    elif i > 0:
      rhum[i] = rhum[i-1]
    
    # Calculate wind speed
    wspd[i] = np.sqrt(uwnd[i]**2 + vwnd[i]**2)
    
    # Wind stress - use direct values if available, otherwise calculate from winds
    if 'uflx' in data and data['uflx'] is not None:
      tx[i] = data['uflx']
    else:
      # Calculate using bulk formula: Cd ~ 1.3e-3 for moderate winds
      Cd = 1.3e-3
      rho_air = 1.225  # kg/m^3
      tx[i] = Cd * rho_air * uwnd[i] * wspd[i]
    
    if 'vflx' in data and data['vflx'] is not None:
      ty[i] = data['vflx']
    else:
      Cd = 1.3e-3
      rho_air = 1.225
      ty[i] = Cd * rho_air * vwnd[i] * wspd[i]
    
    # Precipitation rate: Convert from [kg/m^2/s] to cm/day if in correct units
    if 'prate' in data and data['prate'] is not None:
      prate[i] = data['prate'] * 0.1 * 86400.0  # kg/m^2/s to cm/day
    elif i > 0:
      prate[i] = prate[i-1]
    
    # Radiation fluxes
    # Net shortwave: downward - upward (CROCO convention: positive downward)
    dswrf_val = data.get('dswrf', 0) if 'dswrf' in data else 0
    uswrf_val = data.get('uswrf', 0) if 'uswrf' in data else 0
    if isinstance(dswrf_val, np.ndarray) and isinstance(uswrf_val, np.ndarray):
      radsw[i] = dswrf_val - uswrf_val
    elif i > 0:
      radsw[i] = radsw[i-1]
    
    # Net longwave: upward - downward (CROCO convention: positive upward)
    dlwrf_val = data.get('dlwrf', 0) if 'dlwrf' in data else 0
    ulwrf_val = data.get('ulwrf', 0) if 'ulwrf' in data else 0
    if isinstance(dlwrf_val, np.ndarray) and isinstance(ulwrf_val, np.ndarray):
      radlw[i] = ulwrf_val - dlwrf_val
      radlw_in[i] = dlwrf_val
    elif i > 0:
      radlw[i] = radlw[i-1]
      radlw_in[i] = radlw_in[i-1]

  if len(gfstime) > 0:
    first_written_time = datetime(Yorig, 1, 1) + timedelta(days=float(gfstime[0]))
    last_written_time = datetime(Yorig, 1, 1) + timedelta(days=float(gfstime[-1]))
    print(f'  Writing GFS file time range: {first_written_time} -> {last_written_time} ({len(gfstime)} records)')

  # Ensure all forcing arrays are finite to avoid downstream interpolation failures
  prate = _fill_missing_timesteps(prate)
  radsw = _fill_missing_timesteps(radsw)
  radlw = _fill_missing_timesteps(radlw)
  radlw_in = _fill_missing_timesteps(radlw_in)
  
  # Write output file
  mask[np.isnan(mask)] = 0
  write_GFS(gfs_name, Yorig, lon, lat, mask, gfstime, tx, ty, tair, rhum, prate, wspd, uwnd, vwnd, radlw, radlw_in, radsw)
  
  if keep_grib:
    print(' ... KEEPING DOWNLOADED GRIB FILES ... ')
    print(f'  GRIB files kept in: {temp_dir}')
  else:
    print(' ... CLEANING UP DOWNLOADED GRIB FILES ... ')
    cleanup_files = list(grib_files) + [f for f in grib_files_b if f is not None]
    for grib_file in cleanup_files:
      try:
        os.remove(grib_file)
      except Exception:
        pass
    try:
      os.rmdir(temp_dir)
    except Exception:
      pass
  
  print('Download GFS from AWS S3: done')
  return
