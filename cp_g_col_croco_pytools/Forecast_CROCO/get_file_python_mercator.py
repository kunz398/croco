from cp_g_col_croco_pytools.Forecast_CROCO.old__croco_tools_params import *
import copernicusmarine
from datetime import datetime, timedelta
import re
from netCDF4 import Dataset
import numpy as np
import os

def get_mercator(mercator_type, vars, geom, date, info, outname):
    """
    Generalized function to download Mercator data using the new Copernicus Marine API.
    Simplified signature: removed unnecessary parameters (pathmotu, url, sid, pid).
    
    Parameters:
    - mercator_type: 1 for Mercator 1/12°, 2 for UK Met Office 1/4°, 4 for independent variable datasets.
    - vars: String of variables, e.g., ' --variable thetao --variable so'.
    - geom: [min_lon, max_lon, min_lat, max_lat, min_depth, max_depth].
    - date: [start_date_str, end_date_str] in ISO format.
    - info: [username, password].
    - outname: Output filename.
    """
    
    # Authenticate without prompting
    username = info[0] if len(info) > 0 else None
    password = info[1] if len(info) > 1 else None
    if username and password:
        os.environ['COPERNICUSMARINE_USERNAME'] = username
        os.environ['COPERNICUSMARINE_PASSWORD'] = password
    else:
        print("Warning: No credentials provided. Ensure COPERNICUSMARINE_USERNAME and COPERNICUSMARINE_PASSWORD are set.")
    
    # Parse dates
    start_datetime = datetime.fromisoformat(date[0].replace('Z', '+00:00'))
    end_datetime = datetime.fromisoformat(date[1].replace('Z', '+00:00'))
    
    # Parse geometry
    minimum_longitude = float(geom[0])
    maximum_longitude = float(geom[1])
    minimum_latitude = float(geom[2])
    maximum_latitude = float(geom[3])
    minimum_depth = float(geom[4]) if len(geom) > 4 and geom[4] != '0' else None
    
    # Always use Mercator's full depth range
    maximum_depth = 5727.9  # Use Mercator's maximum depth to get all available levels
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(outname), exist_ok=True)
    
    if mercator_type == 4:
        # For type=4, handle separate datasets
        parsed_vars = re.findall(r'--variable\s+(\w+)', vars) if vars else []
        
        # Define groups
        groups = []
        if 'uo' in parsed_vars or 'vo' in parsed_vars:
            groups.append(('currents', "cmems_mod_glo_phy-cur_anfc_0.083deg_PT6H-i", ['uo', 'vo']))
        if 'so' in parsed_vars:
            groups.append(('salinity', "cmems_mod_glo_phy-so_anfc_0.083deg_PT6H-i", ['so']))
        if 'thetao' in parsed_vars:
            groups.append(('temperature', "cmems_mod_glo_phy-thetao_anfc_0.083deg_PT6H-i", ['thetao']))
        if 'zos' in parsed_vars or 'ssh' in parsed_vars or 'sea_surface_height' in parsed_vars:
            groups.append(('ssh', "cmems_mod_glo_phy_anfc_merged-sl_PT1H-i", ['sea_surface_height']))
        
        if not groups:
            raise ValueError(f"No matching dataset for variables {parsed_vars} in mercator_type=4.")
        
        if len(groups) == 1:
            # Single group, download directly
            dataset_id = groups[0][1]
            variables = groups[0][2]
            subset_kwargs = {
                "dataset_id": dataset_id,
                "variables": variables,
                "start_datetime": start_datetime,
                "end_datetime": end_datetime,
                "minimum_longitude": minimum_longitude,
                "maximum_longitude": maximum_longitude,
                "minimum_latitude": minimum_latitude,
                "maximum_latitude": maximum_latitude,
                "output_filename": outname,
            }
            if minimum_depth is not None:
                subset_kwargs["minimum_depth"] = minimum_depth
            if maximum_depth is not None:
                subset_kwargs["maximum_depth"] = maximum_depth
            print(f"Starting download for {dataset_id}...")
            copernicusmarine.subset(**subset_kwargs)
            print(f"Download completed: {outname}")
        else:
            # Multiple groups, download separately and merge
            temp_files = []
            for name, ds_id, vars_list in groups:
                temp_file = outname.replace('.nc', f'_{name}.nc')
                if os.path.exists(temp_file):
                    print(f"File {temp_file} already exists, skipping download")
                    temp_files.append(temp_file)
                    continue

                subset_kwargs = {
                    "dataset_id": ds_id,
                    "variables": vars_list,
                    "start_datetime": start_datetime,
                    "end_datetime": end_datetime,
                    "minimum_longitude": minimum_longitude,
                    "maximum_longitude": maximum_longitude,
                    "minimum_latitude": minimum_latitude,
                    "maximum_latitude": maximum_latitude,
                    "output_filename": temp_file,
                }
                if minimum_depth is not None:
                    subset_kwargs["minimum_depth"] = minimum_depth
                if maximum_depth is not None:
                    subset_kwargs["maximum_depth"] =  5727.9 #maximum_depth
                print(f"Starting download for {ds_id}...")
                copernicusmarine.subset(**subset_kwargs)
                temp_files.append(temp_file)
            
            # Merge datasets using netCDF4
            print("Merging datasets...")
            
            # Get reference time and dimensions from currents file (6-hourly)
            ref_file = next(f for f in temp_files if '_currents.nc' in f)
            with Dataset(ref_file) as nc:
                time = nc.variables['time'][:]
                depth = nc.variables['depth'][:]
                lat = nc.variables['latitude'][:]
                lon = nc.variables['longitude'][:]
            
            # Create output file
            print(f"Creating merged file: {outname}")
            with Dataset(outname, 'w') as nc_out:
                # Create dimensions
                nc_out.createDimension('time', len(time))
                nc_out.createDimension('depth', len(depth))
                nc_out.createDimension('latitude', len(lat))
                nc_out.createDimension('longitude', len(lon))
                
                # Create coordinate variables
                time_var = nc_out.createVariable('time', 'f8', ('time',))
                depth_var = nc_out.createVariable('depth', 'f8', ('depth',))
                lat_var = nc_out.createVariable('latitude', 'f8', ('latitude',))
                lon_var = nc_out.createVariable('longitude', 'f8', ('longitude',))
                
                # Copy coordinate data
                time_var[:] = time
                depth_var[:] = depth
                lat_var[:] = lat
                lon_var[:] = lon
                
                # Copy variables from each file
                for temp_file in temp_files:
                    print(f"Processing {temp_file}...")
                    with Dataset(temp_file) as nc_in:
                        # First check if this file has depths and validate them (except for SSH)
                        if 'depth' in nc_in.variables and 'ssh' not in temp_file:
                            file_depths = nc_in.variables['depth'][:]
                            if len(file_depths) != len(depth):
                                print(f"Warning: {temp_file} has {len(file_depths)} depth levels, expected {len(depth)}.")
                                print("This might indicate an issue with the downloaded data.")
                                
                        for var_name, var_in in nc_in.variables.items():
                            if var_name in ['time', 'depth', 'latitude', 'longitude']:
                                continue
                                
                            if var_name == 'sea_surface_height':
                                var_name = 'zos'  # Rename SSH to zos
                                
                                # Resample 1-hourly SSH to 6-hourly if needed
                                if len(var_in) > len(time):
                                    # Simple subsample every 6th point for now
                                    data = var_in[::6]
                                else:
                                    data = var_in[:]
                                
                                # Create SSH variable without depth dimension
                                dims = ('time', 'latitude', 'longitude')
                                var_out = nc_out.createVariable(var_name, var_in.dtype, dims)
                                var_out[:] = data
                            else:
                                # Handle potential dimension mismatches
                                src_shape = var_in.shape
                                if 'depth' in var_in.dimensions and src_shape[var_in.dimensions.index('depth')] != len(depth):
                                    print(f"Warning: {var_name} has different depth dimension. Interpolating to match reference depths...")
                                    # Create the variable with correct dimensions
                                    var_out = nc_out.createVariable(var_name, var_in.dtype, var_in.dimensions)
                                    
                                    # Get source depths
                                    src_depths = nc_in.variables['depth'][:]
                                    
                                    # Create interpolation for each timestep
                                    for t in range(len(time)):
                                        for lat_idx in range(len(lat)):
                                            for lon_idx in range(len(lon)):
                                                # Extract profile
                                                profile = var_in[t, :, lat_idx, lon_idx]
                                                # Interpolate to new depths
                                                interp_profile = np.interp(depth, src_depths, profile)
                                                var_out[t, :, lat_idx, lon_idx] = interp_profile
                                else:
                                    # Copy other variables directly
                                    var_out = nc_out.createVariable(var_name, var_in.dtype, var_in.dimensions)
                                    var_out[:] = var_in[:]
            
            print(f"Merged file created: {outname}")
            # Keep temp files for debugging
    
    else:
        # For type=1 and 2, use combined dataset
        if mercator_type == 1:
            dataset_id = "cmems_mod_glo_phy_anfc_0.083deg_P1D-m"
        elif mercator_type == 2:
            dataset_id = "cmems_mod_glo_phy_coupled_anfc_0.25deg_P1D-m"
        else:
            raise ValueError(f"Unsupported mercator_type: {mercator_type}")
        
        # Parse variables
        if vars:
            variables = re.findall(r'--variable\s+(\w+)', vars)
        else:
            variables = None
        
        subset_kwargs = {
            "dataset_id": dataset_id,
            "variables": variables,
            "start_datetime": start_datetime,
            "end_datetime": end_datetime,
            "minimum_longitude": minimum_longitude,
            "maximum_longitude": maximum_longitude,
            "minimum_latitude": minimum_latitude,
            "maximum_latitude": maximum_latitude,
            "output_filename": outname,
        }
        
        # Add depth if specified
        if minimum_depth is not None:
            subset_kwargs["minimum_depth"] = minimum_depth
        if maximum_depth is not None:
            subset_kwargs["maximum_depth"] = maximum_depth
        
        # Download the data
        print(f"Starting download for {dataset_id}...")
        copernicusmarine.subset(**subset_kwargs)
        print(f"Download completed: {outname}")

# Example usage (uncomment to test)
# get_mercator(
#     mercator_type=1,
#     vars=' --variable thetao --variable so --variable uo --variable vo',
#     geom=['-173', '-165', '-23', '-16', '0', '6000'],
#     date=['2023-01-01T00:00:00.000Z', '2023-01-02T00:00:00.000Z'],
#     info=['your_username', 'your_password'],
#     outname='test_mercator.nc'
# )
